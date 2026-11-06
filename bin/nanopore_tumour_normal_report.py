#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nanopore_tumour_normal_report.py
================================

Comprehensive tumour-normal reporting for the DETERMINE Oxford Nanopore
long-read somatic pipeline. Companion to ``nanopore_germline_report.py``.

Ingests the contents of a per-case ``Tissue/`` working directory and produces a
single integrated molecular-tumour-board (MTB) style report in three formats:

    <out>/<patient>.tn_report.html    interactive, self-contained
    <out>/<patient>.tn_report.pdf     print/archival (WeasyPrint)
    <out>/<patient>.tn_report.json    machine-readable, LIMS-friendly

Annotation source of truth:
    * Somatic SNV/indel, CNA, TMB, signatures -> PCGR  (pcgr/)
    * Germline cancer predisposition           -> CPSR  (cpsr/)
The script never *re-annotates*; it harvests, reconciles and presents.

Design goals
------------
* Defensive: any missing input degrades to "not available", never a crash.
* Reconciling: purity/ploidy is reported from ASCAT, wakhan and SAVANA side by
  side, with a consensus call and a divergence flag.
* Reproducible: every input file consumed is logged with size + mtime + sha256
  into the JSON provenance block.
* Extensible: every data source is a small ``parse_*`` function returning a
  plain dict; add a parser, add a section.

Author: built for A. Beggs (University of Birmingham / UHB), DETERMINE.
"""

from __future__ import annotations

import argparse
import base64
import datetime as _dt
import glob
import gzip
import hashlib
import io
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
import sys
import textwrap
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional, Tuple

# ----------------------------------------------------------------------------
# Soft dependencies. None of these are allowed to be fatal at import time.
# ----------------------------------------------------------------------------
try:
    import numpy as np
except Exception:  # pragma: no cover
    np = None

try:
    import pandas as pd
    _PD_ERR = None
except Exception as _exc:  # pragma: no cover
    pd = None
    _PD_ERR = _exc

try:
    import matplotlib
    matplotlib.use("Agg")          # headless: must be set before pyplot import
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator  # noqa: F401
    _HAVE_MPL = True
    _MPL_ERR = None
except Exception as _exc:  # pragma: no cover
    _HAVE_MPL = False
    _MPL_ERR = _exc

try:
    from jinja2 import Environment, BaseLoader, select_autoescape
    _HAVE_JINJA = True
except Exception:  # pragma: no cover
    _HAVE_JINJA = False

LOG = logging.getLogger("tn_report")

SCHEMA_VERSION = "1.0.0"
REPORT_GENERATOR = "nanopore_tumour_normal_report.py"


# ============================================================================
# Configuration
# ============================================================================
@dataclass
class Config:
    input_dir: str
    out_dir: str
    patient_id: str = "PATIENT05"
    display_id: str = ""                   # ID shown in the report title (defaults to patient_id)
    tumour_id: str = "T9"
    normal_id: str = "K6"
    genome_build: str = "grch38"          # PCGR build for somatic annotation
    sample_label: str = "PATIENT_05_T9"   # file stem used across the tree
    assay: str = "ONT WGS, tumour-normal"
    samtools: str = "samtools"
    mosdepth: str = "mosdepth"
    bcftools: str = "bcftools"
    skip_coverage: bool = False
    no_pdf: bool = False
    # Targeted promoter panel for methylation read-out. Coordinates are GRCh38.
    # Extend as needed; values that look like estimates are flagged in-report.
    promoter_panel: Dict[str, Tuple[str, int, int, bool]] = field(default_factory=lambda: {
        # gene: (chrom, start, end, verified?)
        "MLH1":  ("chr3", 36_993_226, 36_993_572, True),
        "MSH2":  ("chr2", 47_403_067, 47_403_437, True),
        "MGMT":  ("chr10", 129_466_683, 129_467_311, True),
        "CDKN2A":("chr9", 21_994_500, 21_995_300, False),
    })


# ============================================================================
# Small utilities
# ============================================================================
def _exists(path: Optional[str]) -> bool:
    return bool(path) and os.path.exists(path)


def _first(input_dir: str, *patterns: str) -> Optional[str]:
    """Return the first existing path matching any glob pattern (relative)."""
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(input_dir, pat)))
        for h in hits:
            if os.path.exists(h):
                return h
    return None


def _sha256(path: str, limit_bytes: int = 0) -> Optional[str]:
    try:
        h = hashlib.sha256()
        read = 0
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
                read += len(chunk)
                if limit_bytes and read >= limit_bytes:
                    break
        return h.hexdigest()
    except Exception:
        return None


def _open_text(path: str):
    """Open plain or gzip text transparently."""
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _read_table(path: str, sep: str = "\t", comment: Optional[str] = None,
                header: Any = "infer") -> Optional["pd.DataFrame"]:
    if pd is None or not _exists(path):
        return None
    try:
        return pd.read_csv(path, sep=sep, comment=comment, header=header,
                           low_memory=False, dtype=str)
    except Exception as exc:
        LOG.warning("Failed to read table %s: %s", path, exc)
        return None


def _to_num(x: Any) -> Optional[float]:
    try:
        if x is None or x == "" or str(x).upper() in {"NA", "NAN", ".", "NONE"}:
            return None
        return float(x)
    except Exception:
        return None


def _pick_col(df: "pd.DataFrame", *candidates: str) -> Optional[str]:
    """Case-insensitive fuzzy column picker; returns the first match found."""
    if df is None:
        return None
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        cl = cand.lower()
        if cl in lower:
            return lower[cl]
    # substring fallback
    for cand in candidates:
        cl = cand.lower()
        for c in df.columns:
            if cl in c.lower():
                return c
    return None


def _png_b64(fig) -> Optional[str]:
    if not _HAVE_MPL:
        return None
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
        plt.close(fig)
        return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")
    except Exception as exc:
        LOG.warning("Figure render failed: %s", exc)
        return None


def _img_file_b64(path: str) -> Optional[str]:
    """Embed an existing PNG file (e.g. ASCAT plots) as a data URI.

    When Pillow is available the image is decoded and re-encoded to a clean
    PNG. This both validates it (a corrupt or placeholder image is dropped
    rather than embedded) and guarantees the PDF backend receives a
    well-formed stream. Without Pillow we fall back to a lightweight
    signature/size sanity check and pass the raw bytes through.
    """
    if not _exists(path):
        return None
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except Exception:
        return None
    if not data:
        return None
    try:
        from PIL import Image  # optional dependency
        with Image.open(io.BytesIO(data)) as im:
            im.load()
            # drop 1x1 / placeholder tiles that carry no information
            if im.width < 8 or im.height < 8:
                return None
            buf = io.BytesIO()
            im.convert("RGBA").save(buf, format="PNG")
            data = buf.getvalue()
    except ImportError:
        # Pillow absent: basic magic-byte + size guard (real ASCAT PNGs pass)
        is_png = data[:8] == b"\x89PNG\r\n\x1a\n"
        is_jpg = data[:3] == b"\xff\xd8\xff"
        if not (is_png or is_jpg) or len(data) < 256:
            return None
        mime = "image/png" if is_png else "image/jpeg"
        return f"data:{mime};base64," + base64.b64encode(data).decode("ascii")
    except Exception:
        return None
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def _run(cmd: List[str], timeout: int = 1800) -> Tuple[int, str, str]:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr
    except Exception as exc:
        return 1, "", str(exc)


# ============================================================================
# VCF mini-reader (no pysam dependency)
# ============================================================================
def iter_vcf(path: str, max_records: int = 0):
    """Yield (chrom, pos, id, ref, alt, qual, filt, info_dict, fmt_keys, samples)."""
    if not _exists(path):
        return
    n = 0
    sample_names: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_names = cols[9:] if len(cols) > 9 else []
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            info_d: Dict[str, str] = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_d[k] = v
                elif kv:
                    info_d[kv] = "True"
            fmt_keys = parts[8].split(":") if len(parts) > 8 else []
            samples = {}
            for i, sname in enumerate(sample_names):
                idx = 9 + i
                if idx < len(parts):
                    vals = parts[idx].split(":")
                    samples[sname] = dict(zip(fmt_keys, vals))
            yield (chrom, int(pos) if pos.isdigit() else pos, vid, ref, alt,
                   qual, filt, info_d, fmt_keys, samples)
            n += 1
            if max_records and n >= max_records:
                return


def count_vcf_records(path: str, pass_only: bool = False) -> int:
    n = 0
    for rec in iter_vcf(path):
        if pass_only and rec[6] not in ("PASS", ".", ""):
            continue
        n += 1
    return n


# ============================================================================
# Parsers — each returns a plain dict, tolerant of missing inputs
# ============================================================================
def parse_provenance(cfg: Config) -> Dict[str, Any]:
    """Hash and stat every input file we care about for the audit trail."""
    files = []
    interesting = [
        "pcgr/*", "cpsr/*", "ClairS/*.vcf.gz", "ascat/*.txt", "ascatpurity.txt",
        "ascatploidy.txt", "savana/*.tsv", "savana/*.vcf", "severus_out/*",
        "wakhan/*.tsv", "wakhan/*.html", "nanoDx/*", "*.msisensor*",
        "*.maf", "*.bedMethyl", "*.bedmethyl", "*.modcall.vcf",
        "pvacseq/*/*", "SpecImmune/*",
    ]
    seen = set()
    for pat in interesting:
        for path in sorted(glob.glob(os.path.join(cfg.input_dir, pat))):
            if not os.path.isfile(path) or path in seen:
                continue
            seen.add(path)
            try:
                st = os.stat(path)
                size = st.st_size
                mtime = _dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")
            except Exception:
                size, mtime = None, None
            # Only hash files < 200 MB to keep this quick; otherwise prefix-hash.
            digest = _sha256(path, limit_bytes=0 if (size or 0) < 200 * 1 << 20 else 64 << 20)
            files.append({
                "path": os.path.relpath(path, cfg.input_dir),
                "bytes": size, "mtime": mtime, "sha256": digest,
            })
    return {"n_files": len(files), "files": files}


def parse_coverage(cfg: Config) -> Dict[str, Any]:
    """Coverage QC from existing mosdepth summaries (read regardless of
    --skip-coverage); live samtools coverage only as a fallback when no summary
    exists and computation is permitted."""
    out: Dict[str, Any] = {"available": False, "samples": {}, "method": None}

    # mosdepth summaries already produced by the pipeline (tumour + normal).
    # New layout puts them under mosdepth_tumour/ and mosdepth_germline/.
    mos_found = False
    for role, sid, subdir in (("tumour", cfg.tumour_id, "mosdepth_tumour"),
                              ("normal", cfg.normal_id, "mosdepth_germline")):
        msum = _first(cfg.input_dir, "%s/*mosdepth.summary.txt" % subdir,
                      "%s/*.mosdepth.summary.txt" % subdir,
                      "*%s*mosdepth.summary.txt" % sid)
        if msum and pd is not None:
            try:
                df = pd.read_csv(msum, sep="\t")
                total = df[df["chrom"] == "total"]
                if len(total):
                    r = total.iloc[0]
                    out["samples"][role] = {"mean_depth": _to_num(r.get("mean")),
                                            "min_depth": _to_num(r.get("min")),
                                            "max_depth": _to_num(r.get("max")),
                                            "source": os.path.basename(msum)}
                    mos_found = True
            except Exception as exc:
                LOG.warning("mosdepth summary parse failed (%s): %s", msum, exc)
    if mos_found:
        out["available"] = True
        out["method"] = "mosdepth"
        return out

    if cfg.skip_coverage:
        out["note"] = "no mosdepth summary found; live computation skipped (--skip-coverage)"
        return out

    bam_map = {
        "tumour": _first(cfg.input_dir, "bams/*%s*.sorted.bam" % cfg.tumour_id),
        "normal": _first(cfg.input_dir, "bams/*%s*.sorted.bam" % cfg.normal_id),
    }

    # samtools coverage fallback (mean depth across primary contigs)
    any_ok = False
    for role, bam in bam_map.items():
        if not _exists(bam):
            continue
        rc, so, se = _run([cfg.samtools, "coverage", bam], timeout=3600)
        if rc != 0:
            LOG.warning("samtools coverage failed for %s: %s", bam, se[:200])
            continue
        depths, widths = [], []
        for ln in so.splitlines():
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) >= 7 and re.match(r"^(chr)?([0-9]{1,2}|X|Y)$", f[0]):
                d = _to_num(f[6]); w = _to_num(f[2])
                if d is not None and w is not None:
                    depths.append(d * w); widths.append(w)
        mean_depth = (sum(depths) / sum(widths)) if widths else None
        out["samples"][role] = {"mean_depth": mean_depth, "source": os.path.basename(bam)}
        any_ok = True
    out["available"] = any_ok
    out["method"] = "samtools coverage" if any_ok else None
    return out


def parse_pcgr(cfg: Config) -> Dict[str, Any]:
    """Harvest PCGR somatic outputs (the annotation source of truth)."""
    pdir = os.path.join(cfg.input_dir, "pcgr")
    out: Dict[str, Any] = {"available": os.path.isdir(pdir), "html": None}
    if not out["available"]:
        return out
    out["html"] = _first(cfg.input_dir, "pcgr/*.html")

    # --- SNV/indel annotated table -----------------------------------------
    snv_tsv = _first(cfg.input_dir, "pcgr/*snv_indel_ann.tsv.gz", "pcgr/*.pass.tsv.gz")
    snvs: List[Dict[str, Any]] = []
    tier_counts: Dict[str, int] = {}
    if snv_tsv and pd is not None:
        df = _read_table(snv_tsv)
        if df is not None and len(df):
            c_gene = _pick_col(df, "SYMBOL", "GENE_SYMBOL", "GENE")
            c_cons = _pick_col(df, "CONSEQUENCE", "VARIANT_CLASS")
            c_hgvsp = _pick_col(df, "HGVSp_short", "HGVSp", "PROTEIN_CHANGE")
            c_hgvsc = _pick_col(df, "HGVSc", "CDS_CHANGE")
            c_tier = _pick_col(df, "TIER", "ACTIONABILITY_TIER", "CPSR_CLASSIFICATION")
            c_onc = _pick_col(df, "ONCOGENICITY", "CLINPATH_ONCOGENICITY")
            c_vaf = _pick_col(df, "VAF_TUMOR", "TVAF", "AF_TUMOR", "VAF")
            c_dp = _pick_col(df, "DP_TUMOR", "TDP", "DP")
            c_chrom = _pick_col(df, "CHROM", "Chromosome")
            c_pos = _pick_col(df, "POS", "Start_Position", "POSITION")
            c_ref = _pick_col(df, "REF", "Reference_Allele")
            c_alt = _pick_col(df, "ALT", "Tumor_Seq_Allele2", "ALLELE")
            c_evi = _pick_col(df, "CLINICAL_EVIDENCE", "BIOMARKER_EVIDENCE", "ACTIONABILITY")
            for _, r in df.iterrows():
                tier = str(r.get(c_tier, "")) if c_tier else ""
                tier_counts[tier] = tier_counts.get(tier, 0) + 1
                snvs.append({
                    "gene": r.get(c_gene) if c_gene else None,
                    "consequence": r.get(c_cons) if c_cons else None,
                    "hgvsp": r.get(c_hgvsp) if c_hgvsp else None,
                    "hgvsc": r.get(c_hgvsc) if c_hgvsc else None,
                    "tier": tier or None,
                    "oncogenicity": r.get(c_onc) if c_onc else None,
                    "vaf": _to_num(r.get(c_vaf)) if c_vaf else None,
                    "dp": _to_num(r.get(c_dp)) if c_dp else None,
                    "chrom": r.get(c_chrom) if c_chrom else None,
                    "pos": r.get(c_pos) if c_pos else None,
                    "ref": r.get(c_ref) if c_ref else None,
                    "alt": r.get(c_alt) if c_alt else None,
                    "evidence": r.get(c_evi) if c_evi else None,
                })
    out["snv_indel"] = {
        "source": os.path.basename(snv_tsv) if snv_tsv else None,
        "n_total": len(snvs),
        "tier_counts": tier_counts,
        "variants": snvs,
    }

    # --- TMB ----------------------------------------------------------------
    tmb_tsv = _first(cfg.input_dir, "pcgr/*tmb.tsv")
    tmb: Dict[str, Any] = {"available": bool(tmb_tsv)}
    if tmb_tsv and pd is not None:
        df = _read_table(tmb_tsv)
        if df is not None and len(df):
            c_tmb = _pick_col(df, "tmb_estimate", "TMB", "tmb")
            c_meas = _pick_col(df, "tmb_measure", "measure")
            c_n = _pick_col(df, "n_somatic_variants", "n_variants", "n_somatic")
            row = df.iloc[0]
            tmb.update({
                "tmb_per_mb": _to_num(row.get(c_tmb)) if c_tmb else None,
                "measure": row.get(c_meas) if c_meas else None,
                "n_somatic": _to_num(row.get(c_n)) if c_n else None,
            })
    out["tmb"] = tmb

    # --- Mutational signatures ---------------------------------------------
    sig_tsv = _first(cfg.input_dir, "pcgr/*msigs.tsv.gz", "pcgr/*msigs.tsv")
    sigs: List[Dict[str, Any]] = []
    if sig_tsv and pd is not None:
        df = _read_table(sig_tsv)
        if df is not None and len(df):
            c_sig = _pick_col(df, "signature_id", "signature", "signature_aetiology")
            c_w = _pick_col(df, "prop_signature", "contribution", "weight", "exposure", "value")
            # PCGR's 'aetiology' is long prose; 'group' holds the short keyword
            c_aet = _pick_col(df, "aetiology_keyword", "group", "aetiology", "etiology",
                              "description")
            if c_sig is None or c_w is None:      # header didn't match — go by content
                v_sig, v_w = _sig_cols_by_value(df)
                c_sig = c_sig or v_sig
                c_w = c_w or v_w
                LOG.warning("PCGR msigs: header not recognised (%s); using columns "
                            "sig=%s weight=%s", list(df.columns)[:8], c_sig, c_w)

            def _w_num(v):
                """PCGR writes contribution as a percentage STRING ('15.11%')."""
                if v is None:
                    return None
                return _to_num(str(v).strip().replace("%", "").replace(",", ""))

            tmpdf = df.copy()
            if c_w:
                tmpdf["_w"] = tmpdf[c_w].map(_w_num)
                tmpdf = tmpdf.sort_values("_w", ascending=False)
            for _, r in tmpdf.iterrows():
                w = _w_num(r.get(c_w)) if c_w else None
                if w is not None and w <= 0:
                    continue
                name = str(r.get(c_sig)).strip() if c_sig else None
                aet = r.get(c_aet) if c_aet else None
                if name and " - " in name and not aet:   # 'SBS6 - MMR deficiency'
                    name, aet = [x.strip() for x in name.split(" - ", 1)]
                sigs.append({"signature": name, "weight": w, "aetiology": aet})
            # contributions may be percentages (0-100) rather than proportions
            ws = [s["weight"] for s in sigs if s["weight"] is not None]
            if ws and max(ws) > 1.5:
                for s in sigs:
                    if s["weight"] is not None:
                        s["weight"] = s["weight"] / 100.0
            # A row with no usable weight is useless for the plots — the figures
            # filter on weight>0, which previously left them silently empty.
            if not [s for s in sigs if (s["weight"] or 0) > 0]:
                LOG.warning("PCGR msigs parsed %d row(s) but no usable weights "
                            "(sig col=%s, weight col=%s, columns=%s) — trying the "
                            "PCGR HTML instead", len(sigs), c_sig, c_w,
                            list(df.columns)[:10])
                sigs = []
    if not sigs:
        # Fallback: pull them out of the PCGR HTML report's signature widget
        html = _first(cfg.input_dir, "pcgr/*.html")
        if html:
            sigs = _parse_pcgr_html_signatures(html)
            if sigs:
                LOG.info("Signatures recovered from PCGR HTML (%d)", len(sigs))
    if not sigs:
        LOG.warning("No mutational signatures parsed (tsv=%s) — signature plots "
                    "will be omitted", sig_tsv)
    out["signatures"] = {"source": os.path.basename(sig_tsv) if sig_tsv else None,
                         "items": sigs}

    # --- Gene-level CNA -----------------------------------------------------
    cna_tsv = _first(cfg.input_dir, "pcgr/*cna_gene_ann.tsv.gz", "pcgr/*cna_gene.tsv.gz")
    cnas: List[Dict[str, Any]] = []
    if cna_tsv and pd is not None:
        df = _read_table(cna_tsv)
        if df is not None and len(df):
            c_gene = _pick_col(df, "SYMBOL", "gene_symbol", "GENE")
            c_gname = _pick_col(df, "GENENAME", "gene_name", "genename")
            c_event = _pick_col(df, "VARIANT_CLASS", "variant_class", "event_type", "cna_type", "STATUS")
            c_etype = _pick_col(df, "EVENT_TYPE", "event_type", "cna_event_type")
            c_cn = _pick_col(df, "TOTAL_CN", "total_cn", "CN_TOTAL", "CN", "n_copy")
            c_log2 = _pick_col(df, "LOG2R", "log2r", "log2_ratio")
            c_band = _pick_col(df, "CYTOBAND", "cytoband", "cyto_band")
            c_tier = _pick_col(df, "TIER", "ACTIONABILITY_TIER", "AMP_LEVEL")
            c_evi = _pick_col(df, "BIOMARKER_EVIDENCE", "CLINICAL_EVIDENCE", "ACTIONABILITY")
            for _, r in df.iterrows():
                cnas.append({
                    "gene": r.get(c_gene) if c_gene else None,
                    "genename": r.get(c_gname) if c_gname else None,
                    "event": r.get(c_event) if c_event else None,
                    "event_type": r.get(c_etype) if c_etype else None,
                    "total_cn": _to_num(r.get(c_cn)) if c_cn else None,
                    "log2r": _to_num(r.get(c_log2)) if c_log2 else None,
                    "cytoband": r.get(c_band) if c_band else None,
                    "tier": r.get(c_tier) if c_tier else None,
                    "evidence": r.get(c_evi) if c_evi else None,
                })
    out["cna_genes"] = {"source": os.path.basename(cna_tsv) if cna_tsv else None,
                        "items": cnas}
    return out


def parse_cpsr(cfg: Config) -> Dict[str, Any]:
    """Germline cancer-predisposition findings from CPSR (normal sample)."""
    cdir = os.path.join(cfg.input_dir, "cpsr")
    out: Dict[str, Any] = {"available": os.path.isdir(cdir) or os.path.islink(cdir),
                           "html": None, "variants": []}
    if not out["available"]:
        return out
    out["html"] = _first(cfg.input_dir, "cpsr/*.html")
    tsv = _first(cfg.input_dir, "cpsr/*classification.tsv.gz",
                 "cpsr/*tiers.tsv.gz", "cpsr/*snvs_indels*.tsv.gz", "cpsr/*.tsv.gz")
    if tsv and pd is not None:
        df = _read_table(tsv)
        if df is not None and len(df):
            c_gene = _pick_col(df, "SYMBOL", "GENE")
            c_class = _pick_col(df, "CPSR_CLASSIFICATION", "CLINVAR_CLASSIFICATION",
                                "FINAL_CLASSIFICATION", "CLASSIFICATION")
            c_hgvsp = _pick_col(df, "HGVSp_short", "HGVSp", "PROTEIN_CHANGE")
            c_cons = _pick_col(df, "CONSEQUENCE")
            c_geno = _pick_col(df, "GENOTYPE", "GT")
            c_clinvar = _pick_col(df, "CLINVAR_TRAITS", "PHENOTYPE", "CLINVAR_PHENOTYPE")
            keep = {"PATHOGENIC", "LIKELY_PATHOGENIC", "LIKELY PATHOGENIC", "P", "LP",
                    "PATHOGENIC/LIKELY_PATHOGENIC"}
            for _, r in df.iterrows():
                cls = str(r.get(c_class, "")).upper().replace("-", "_") if c_class else ""
                is_plp = any(k in cls for k in ("PATHOGENIC",)) and "BENIGN" not in cls
                out["variants"].append({
                    "gene": r.get(c_gene) if c_gene else None,
                    "classification": r.get(c_class) if c_class else None,
                    "hgvsp": r.get(c_hgvsp) if c_hgvsp else None,
                    "consequence": r.get(c_cons) if c_cons else None,
                    "genotype": r.get(c_geno) if c_geno else None,
                    "trait": r.get(c_clinvar) if c_clinvar else None,
                    "reportable": bool(is_plp),
                })
    out["n_reportable"] = sum(1 for v in out["variants"] if v.get("reportable"))
    return out


def parse_purity_ploidy(cfg: Config) -> Dict[str, Any]:
    """Reconcile purity/ploidy across ASCAT, wakhan and SAVANA."""
    estimates: List[Dict[str, Any]] = []

    # ASCAT --------------------------------------------------------------
    a_pur = _first(cfg.input_dir, "ascat/ascatpurity.txt", "ascatpurity.txt")
    a_plo = _first(cfg.input_dir, "ascat/ascatploidy.txt", "ascatploidy.txt")
    def _scalar(path):
        try:
            with open(path) as fh:
                txt = fh.read().strip().split()
                return _to_num(txt[-1]) if txt else None
        except Exception:
            return None
    if a_pur or a_plo:
        estimates.append({"caller": "ASCAT",
                          "purity": _scalar(a_pur) if a_pur else None,
                          "ploidy": _scalar(a_plo) if a_plo else None,
                          "data_type": "allele-specific (SR-WGS)"})

    # wakhan -------------------------------------------------------------
    wak = _first(cfg.input_dir, "wakhan/solutions_ranks.tsv")
    wak_added = False
    if wak and pd is not None:
        df = _read_table(wak)
        if df is not None and len(df):
            c_p = _pick_col(df, "purity", "cellularity")
            c_pl = _pick_col(df, "ploidy")
            c_rank = _pick_col(df, "rank", "score", "confidence")
            top = df.iloc[0]
            estimates.append({"caller": "wakhan",
                              "purity": _to_num(top.get(c_p)) if c_p else None,
                              "ploidy": _to_num(top.get(c_pl)) if c_pl else None,
                              "data_type": "ONT haplotype-aware"})
            wak_added = True
    if not wak_added:
        # fall back to the solution directory name "ploidy_purity_conf"
        for d in sorted(glob.glob(os.path.join(cfg.input_dir, "wakhan", "*_*_*"))):
            base = os.path.basename(d)
            m = re.match(r"^([0-9.]+)_([0-9.]+)_([0-9.]+)$", base)
            if m:
                estimates.append({"caller": "wakhan",
                                  "purity": _to_num(m.group(2)),
                                  "ploidy": _to_num(m.group(1)),
                                  "data_type": "ONT haplotype-aware",
                                  "note": "parsed from solution dir name"})
                break

    # SAVANA -------------------------------------------------------------
    sav = _first(cfg.input_dir, "savana/*fitted_purity_ploidy.tsv")
    if sav and pd is not None:
        df = _read_table(sav)
        if df is not None and len(df):
            c_p = _pick_col(df, "purity", "cellularity")
            c_pl = _pick_col(df, "ploidy")
            top = df.iloc[0]
            estimates.append({"caller": "SAVANA",
                              "purity": _to_num(top.get(c_p)) if c_p else None,
                              "ploidy": _to_num(top.get(c_pl)) if c_pl else None,
                              "data_type": "ONT read-depth + hetSNP"})

    purs = [e["purity"] for e in estimates if e.get("purity") is not None]
    plos = [e["ploidy"] for e in estimates if e.get("ploidy") is not None]
    consensus = {
        "purity_median": (sorted(purs)[len(purs)//2] if purs else None),
        "ploidy_median": (sorted(plos)[len(plos)//2] if plos else None),
        "purity_range": [min(purs), max(purs)] if purs else None,
        "ploidy_range": [min(plos), max(plos)] if plos else None,
        "divergent": bool(purs and (max(purs) - min(purs) > 0.25)),
    }
    return {"available": bool(estimates), "estimates": estimates, "consensus": consensus}


def parse_ascat_segments(cfg: Config) -> Dict[str, Any]:
    """ASCAT allele-specific segments (Chromosome/Start/End/nMajor/nMinor) +
    the canonical ASCAT plots. The file is named inconsistently across runs
    (T9.segments.txt, DET061.segments2.txt, ascat/*.segments.txt ...), so try
    the specific names first and then any *.segments*.txt in the case dir."""
    seg = _first(cfg.input_dir,
                 "ascat/%s*.segments*.txt" % cfg.tumour_id,
                 "%s*.segments*.txt" % cfg.tumour_id,
                 "ascat/%s*.segments*.txt" % cfg.patient_id,
                 "%s*.segments*.txt" % cfg.patient_id,
                 "ascat/*.segments*.txt",
                 "*.segments*.txt")
    if seg:
        LOG.info("ASCAT segments: %s", os.path.relpath(seg, cfg.input_dir))
    else:
        LOG.warning("No ASCAT segments file found (*.segments*.txt) — the "
                    "allele-specific copy-number plot will be omitted")
    out: Dict[str, Any] = {"available": False, "n_segments": 0, "segments": [],
                           "plots": {}}
    if seg and pd is not None:
        df = _read_table(seg)
        if df is not None and len(df):
            c_chr = _pick_col(df, "chr", "chrom", "Chromosome")
            c_s = _pick_col(df, "startpos", "start", "Start")
            c_e = _pick_col(df, "endpos", "end", "End")
            c_maj = _pick_col(df, "nMajor", "nA", "major_cn")
            c_min = _pick_col(df, "nMinor", "nB", "minor_cn")
            segs = []
            for _, r in df.iterrows():
                nmaj = _to_num(r.get(c_maj)) if c_maj else None
                nmin = _to_num(r.get(c_min)) if c_min else None
                segs.append({
                    "chrom": r.get(c_chr) if c_chr else None,
                    "start": _to_num(r.get(c_s)) if c_s else None,
                    "end": _to_num(r.get(c_e)) if c_e else None,
                    "n_major": nmaj, "n_minor": nmin,
                    "total_cn": (nmaj + nmin) if (nmaj is not None and nmin is not None) else None,
                    "loh": bool(nmin == 0) if nmin is not None else None,
                })
            out.update({"available": True, "n_segments": len(segs), "segments": segs})
    # canonical ASCAT plots
    for key, pat in [("profile", "ascat/*ASCATprofile.png"),
                     ("sunrise", "ascat/*sunrise.png"),
                     ("aspcf", "ascat/*ASPCF.png"),
                     ("rawprofile", "ascat/*rawprofile.png")]:
        p = _first(cfg.input_dir, pat)
        if p:
            out["plots"][key] = _img_file_b64(p)
    return out


def _parse_pcgr_html_signatures(html_path: str) -> List[Dict[str, Any]]:
    """Fallback: read the fitted signatures out of the PCGR HTML report's
    'Key signatures identified' plotly widget. Trace names carry the signature
    AND its aetiology ('SBS6 - MMR deficiency'); y[0] is the proportion. The
    widget duplicates its traces, so de-dup by signature."""
    try:
        with open(html_path, encoding="utf-8", errors="replace") as fh:
            txt = fh.read()
    except Exception as exc:
        LOG.warning("PCGR HTML read failed (%s): %s", html_path, exc)
        return []
    i = txt.find("Key signatures identified")
    if i < 0:
        return []
    seg = txt[i:i + 500000]
    m = re.search(r'data-for="htmlwidget-[^"]+">(\{.*?\})</script>', seg, re.S)
    if not m:
        return []
    try:
        traces = json.loads(m.group(1)).get("x", {}).get("data", [])
    except Exception as exc:
        LOG.warning("PCGR HTML signature widget parse failed: %s", exc)
        return []
    seen, sigs = set(), []
    for t in traces:
        if not isinstance(t, dict):
            continue
        name = str(t.get("name") or "").strip()
        y = t.get("y")
        if not name or not isinstance(y, list) or not y:
            continue
        w = _to_num(y[0])
        if w is None or w <= 0:
            continue
        parts = name.split(" - ", 1)
        sig = parts[0].strip()
        aet = parts[1].strip() if len(parts) > 1 else None
        if sig in seen:
            continue
        seen.add(sig)
        sigs.append({"signature": sig, "weight": w, "aetiology": aet or None})
    sigs.sort(key=lambda s: -(s["weight"] or 0))
    return sigs


def _sig_cols_by_value(df) -> Tuple[Optional[str], Optional[str]]:
    """Locate signature-id and weight columns by CONTENTS, for when the header
    names don't match anything known."""
    if df is None or not len(df):
        return None, None
    sig_re = re.compile(r"^(SBS|ID|DBS|CN|SV|RefSig)\s?\d+[a-z]?$", re.I)
    c_sig = c_w = None
    for c in df.columns:
        vals = [str(v).strip() for v in df[c].head(20) if str(v).strip()]
        if vals and sum(bool(sig_re.match(v)) for v in vals) >= max(1, len(vals) // 2):
            c_sig = c
            break
    for c in df.columns:
        if c == c_sig:
            continue
        raw = [str(v) for v in df[c].head(20)]
        nums = [n for n in (_to_num(v.replace("%", "").replace(",", "").strip())
                            for v in raw) if n is not None]
        if not nums or not any(n > 0 for n in nums):
            continue
        # proportions (0-1) or percentages (0-100), but not CI ranges ('15 - 15.2%')
        if any(" - " in v for v in raw):
            continue
        if all(0 <= n <= 1.0001 for n in nums) or all(0 <= n <= 100.001 for n in nums):
            c_w = c
            break
    return c_sig, c_w


def _parse_pcgr_html_msi(html_path: str) -> Dict[str, Any]:
    """Extract the MSI classification from the PCGR HTML report. PCGR 2.x renders
    the call in a bslib value box: <... value-box-value">MSI.H (Microsatellite
    instability - high)</...> (preceded by an 'MSI status' label)."""
    try:
        with open(html_path, errors="ignore") as fh:
            txt = fh.read()
    except Exception:
        return {"available": False}
    status = None
    # bslib value boxes: pick the value whose text is an MSI/microsatellite call
    for v in re.findall(r'value-box-value[^>]*>\s*([^<]+?)\s*<', txt):
        if re.search(r"\bMSI\b|\bMSS\b|microsatellite", v, re.I):
            status = v.strip()
            break
    if status is None:
        # fallback: flat-text search with scripts stripped
        clean = re.sub(r"<script[\s\S]*?</script>", " ", txt, flags=re.I)
        flat = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", clean))
        m = re.search(r"MSI status\s*[:\-]?\s*"
                      r"(MSI[.\s\w()-]*|MSS[.\s\w()-]*|microsatellite[\s\w()-]*)", flat, re.I)
        status = m.group(1).strip() if m else None
    if not status:
        return {"available": False}
    prob = None
    pm = re.search(r"MSI[\w\s./-]*prob\w*\s*[:=]?\s*([0-9.]+)", txt, re.I)
    if pm:
        prob = _to_num(pm.group(1))
    return {"available": True, "status": status, "probability": prob,
            "source": os.path.basename(html_path)}


def parse_msi(cfg: Config) -> Dict[str, Any]:
    """MSI call taken from PCGR (its statistical MSI classifier), not MSIsensor.
    PCGR 2.x records MSI in the XLSX workbook and the HTML report; try XLSX first,
    then fall back to parsing the HTML."""
    out: Dict[str, Any] = {"available": False, "caller": "PCGR"}
    xlsx = _first(cfg.input_dir, "pcgr/*.pcgr.%s.xlsx" % cfg.genome_build,
                  "pcgr/*.pcgr.*.xlsx", "pcgr/*.xlsx")
    res = parse_pcgr_msi(xlsx) if xlsx else {"available": False}
    if not res.get("available"):
        html = _first(cfg.input_dir, "pcgr/*.pcgr.%s.html" % cfg.genome_build,
                      "pcgr/*.pcgr.*.html", "pcgr/*.html")
        if html:
            res = _parse_pcgr_html_msi(html)
    if res.get("available"):
        status = str(res.get("status", ""))
        sl = status.lower()
        norm = ("MSI" if ("unstable" in sl or "instab" in sl or "msi.h" in sl or "msi-h" in sl
                          or "msi high" in sl or sl.strip() in ("msi", "high"))
                else "MSS" if ("stable" in sl or "mss" in sl) else status)
        out.update({"available": True, "caller": "PCGR",
                    "status": norm, "status_raw": status,
                    "probability": res.get("probability"),
                    "source": res.get("source"),
                    "threshold_note": res.get("classifier",
                        "PCGR MSI classifier (TCGA exome-trained); review locally.")})
    return out


def _parse_sv_vcf(path: str, max_records: int = 5000) -> Dict[str, Any]:
    svtypes: Dict[str, int] = {}
    fusions: List[Dict[str, Any]] = []
    n = 0
    for rec in iter_vcf(path, max_records=max_records):
        chrom, pos, vid, ref, alt, qual, filt, info, fmt_keys, samples = rec
        if filt not in ("PASS", ".", ""):
            continue
        svt = info.get("SVTYPE", "BND")
        svtypes[svt] = svtypes.get(svt, 0) + 1
        # capture annotated gene fusions if present
        for key in ("FUSION", "GENE_FUSION", "GENES", "ANN", "BND_GENES"):
            if key in info and any(sep in str(info[key]) for sep in ("&", "::", "-")):
                fusions.append({"id": vid, "type": svt, "chrom": chrom,
                                "pos": pos, "annotation": info[key][:200]})
                break
        n += 1
    return {"n_pass": n, "svtype_counts": svtypes, "fusions": fusions[:50]}


def parse_structural_variants(cfg: Config) -> Dict[str, Any]:
    out: Dict[str, Any] = {"callers": {}}
    sav = _first(cfg.input_dir, "savana/*classified.somatic.vcf",
                 "savana/*classified.vcf")
    if sav:
        d = _parse_sv_vcf(sav)
        d["source"] = os.path.basename(sav)
        out["callers"]["SAVANA"] = d
    sev = _first(cfg.input_dir, "severus_out/somatic_SVs/*.vcf",
                 "severus_out/all_SVs/*.vcf", "severus_out/*.vcf")
    if sev:
        d = _parse_sv_vcf(sev)
        d["source"] = os.path.relpath(sev, cfg.input_dir)
        out["callers"]["SEVERUS"] = d
    snf = _first(cfg.input_dir, "*%s*.phased_SV.vcf" % cfg.tumour_id)
    if snf:
        d = _parse_sv_vcf(snf)
        d["source"] = os.path.basename(snf)
        out["callers"]["phased_SV"] = d
    out["available"] = bool(out["callers"])
    # union of somatic SV count (prefer SAVANA somatic)
    if "SAVANA" in out["callers"]:
        out["n_somatic_sv"] = out["callers"]["SAVANA"]["n_pass"]
    elif "SEVERUS" in out["callers"]:
        out["n_somatic_sv"] = out["callers"]["SEVERUS"]["n_pass"]

    # ---- Severus vs SAVANA breakpoint concordance ----
    sav_bedpe = _first(cfg.input_dir, "savana/*classified.somatic.bedpe")
    sev_vcf = _first(cfg.input_dir, "severus_out/somatic_SVs/*.vcf",
                     "severus_out/all_SVs/*.vcf")
    sav_bps = load_savana_bedpe_breakpoints(sav_bedpe) if sav_bedpe else []
    sev_bps = load_severus_vcf_breakpoints(sev_vcf) if sev_vcf else []
    if sav_bps or sev_bps:
        out["concordance"] = build_sv_concordance(sav_bps, sev_bps)
        out["concordance"]["sources"] = {
            "SAVANA": os.path.relpath(sav_bedpe, cfg.input_dir) if sav_bedpe else None,
            "SEVERUS": os.path.relpath(sev_vcf, cfg.input_dir) if sev_vcf else None}
    # translocations (inter-chromosomal) for the circos
    out["translocations"] = [b for b in sav_bps if b["chrom1"] != b["chrom2"]] or \
                            [b for b in sev_bps if b["chrom1"] != b["chrom2"]]

    # ---- fusions, annotated by clinical relevance ----
    raw_fusions = [f for d in out["callers"].values() for f in d.get("fusions", [])]
    out["fusions_annotated"] = annotate_fusions(raw_fusions)
    out["n_actionable_fusions"] = sum(1 for f in out["fusions_annotated"] if f.get("clinical"))
    return out


# nanoDx / methylation-classifier class-code annotations (TCGA study + paediatric)
NANODX_CLASS_MAP = {
    "ACC": "Adrenocortical carcinoma", "BLCA": "Bladder urothelial carcinoma",
    "BRCA": "Breast carcinoma", "CESC": "Cervical carcinoma", "CHOL": "Cholangiocarcinoma",
    "COAD": "Colon adenocarcinoma", "READ": "Rectal adenocarcinoma",
    "CO-READ": "Colorectal adenocarcinoma", "COADREAD": "Colorectal adenocarcinoma",
    "DLBC": "Diffuse large B-cell lymphoma", "ESCA": "Oesophageal carcinoma",
    "GBM": "Glioblastoma", "HNSC": "Head & neck squamous carcinoma",
    "KICH": "Kidney chromophobe", "KIRC": "Kidney renal clear-cell carcinoma",
    "KIRP": "Kidney renal papillary carcinoma", "LAML": "Acute myeloid leukaemia",
    "LGG": "Lower-grade glioma", "LIHC": "Liver hepatocellular carcinoma",
    "LUAD": "Lung adenocarcinoma", "LUSC": "Lung squamous carcinoma",
    "MESO": "Mesothelioma", "OV": "Ovarian serous carcinoma", "PAAD": "Pancreatic adenocarcinoma",
    "PCPG": "Phaeochromocytoma & paraganglioma", "PRAD": "Prostate adenocarcinoma",
    "SARC": "Sarcoma", "SKCM": "Cutaneous melanoma", "STAD": "Stomach adenocarcinoma",
    "TGCT": "Testicular germ-cell tumour", "THCA": "Thyroid carcinoma", "THYM": "Thymoma",
    "UCEC": "Uterine corpus endometrial carcinoma", "UCS": "Uterine carcinosarcoma",
    "UVM": "Uveal melanoma", "WT": "Wilms tumour", "NBL": "Neuroblastoma",
    "RT": "Rhabdoid tumour", "HB": "Hepatoblastoma", "OS": "Osteosarcoma",
    "EWS": "Ewing sarcoma", "RMS": "Rhabdomyosarcoma",
}


def _nanodx_desc(code: Optional[str]) -> Optional[str]:
    if not code:
        return None
    key = str(code).strip().upper()
    if key in NANODX_CLASS_MAP:
        return NANODX_CLASS_MAP[key]
    # strip common suffixes/subclass tags (e.g. "GBM_RTK2" -> "GBM")
    base = re.split(r"[_\-.\s]", key)[0]
    return NANODX_CLASS_MAP.get(base)


def parse_nanodx(cfg: Config) -> Dict[str, Any]:
    """Methylation-based classifier (nanoDx)."""
    res = _first(cfg.input_dir, "nanoDx/*nanoDx.txt", "nanoDx/*.txt")
    votes = _first(cfg.input_dir, "nanoDx/*votes.txt")
    out: Dict[str, Any] = {"available": bool(res), "top_class": None,
                           "score": None, "n_features": None, "votes": []}
    if res:
        try:
            with open(res) as fh:
                text = fh.read()
            # try table form first
            if pd is not None:
                try:
                    df = pd.read_csv(res, sep=None, engine="python")
                    c_cls = _pick_col(df, "class", "prediction", "Freq", "label")
                    c_sc = _pick_col(df, "score", "probability", "Freq", "votes")
                    if c_cls is not None:
                        df2 = df.copy()
                        if c_sc:
                            df2["_s"] = df2[c_sc].map(_to_num)
                            df2 = df2.sort_values("_s", ascending=False)
                        top = df2.iloc[0]
                        out["top_class"] = str(top.get(c_cls))
                        out["score"] = _to_num(top.get(c_sc)) if c_sc else None
                except Exception:
                    pass
            if out["top_class"] is None:
                out["raw"] = text[:2000]
        except Exception as exc:
            LOG.warning("nanoDx parse failed: %s", exc)
    if votes and pd is not None:
        try:
            dfv = pd.read_csv(votes, sep=None, engine="python")
            dfv.columns = [str(c).lstrip("#").strip() for c in dfv.columns]
            c_cls = _pick_col(dfv, "class", "label", "prediction")
            c_v = _pick_col(dfv, "votes", "Freq", "score", "probability", "cal_Freq")
            # fall back to positional (first text col, first numeric col) if unnamed
            if c_cls is None and len(dfv.columns):
                c_cls = dfv.columns[0]
            if c_v is None:
                for col in dfv.columns[1:]:
                    if dfv[col].map(_to_num).notna().any():
                        c_v = col
                        break
            if c_cls is not None and c_v is not None:
                dfv["_v"] = dfv[c_v].map(_to_num)
                dfv = dfv.sort_values("_v", ascending=False)
                out["votes"] = [{"class": r.get(c_cls), "votes": _to_num(r.get(c_v)),
                                 "class_desc": _nanodx_desc(r.get(c_cls))}
                                for _, r in dfv.head(8).iterrows()]
                # derive the top class/score from the votes if the prediction file didn't give it
                if out.get("top_class") is None and len(dfv):
                    top = dfv.iloc[0]
                    out["top_class"] = str(top.get(c_cls))
                    out["score"] = _to_num(top.get(c_v))
        except Exception as exc:
            LOG.warning("nanoDx votes parse failed: %s", exc)
    out["top_class_desc"] = _nanodx_desc(out.get("top_class"))
    return out


def parse_methylation(cfg: Config) -> Dict[str, Any]:
    """Targeted promoter methylation from the CpG bedMethyl + modkit DMR counts."""
    out: Dict[str, Any] = {"available": False, "promoters": [], "dmr": {}}
    # The CpG.bedMethyl is GRCh38 by default; only the explicitly named
    # *_hg19.cgID.* file is hg19, which we skip here.
    bm = _first(cfg.input_dir, "*%s*.CpG.bedMethyl" % cfg.tumour_id, "*.CpG.bedMethyl")
    if bm and "hg19" in os.path.basename(bm).lower():
        bm = None
    if bm and _exists(bm):
        out["bedmethyl_source"] = os.path.basename(bm)
        # Stream the (large) bedMethyl once, bucketing into panel regions.
        regions = cfg.promoter_panel
        acc = {g: {"sum": 0.0, "n": 0} for g in regions}
        try:
            with _open_text(bm) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 11:
                        continue
                    chrom = f[0]
                    try:
                        start = int(f[1])
                    except ValueError:
                        continue
                    # modkit bedMethyl: col 11 (index 10) is % modified
                    frac = _to_num(f[10])
                    if frac is None:
                        continue
                    for g, (gc, gs, ge, _v) in regions.items():
                        if chrom == gc and gs <= start <= ge:
                            acc[g]["sum"] += frac
                            acc[g]["n"] += 1
            for g, (gc, gs, ge, verified) in regions.items():
                n = acc[g]["n"]
                mean = (acc[g]["sum"] / n) if n else None
                out["promoters"].append({
                    "gene": g, "region": f"{gc}:{gs}-{ge}",
                    "mean_methylation_pct": round(mean, 1) if mean is not None else None,
                    "n_cpg": n, "coords_verified": verified,
                })
            out["available"] = any(p["n_cpg"] for p in out["promoters"])
        except Exception as exc:
            LOG.warning("bedMethyl streaming failed: %s", exc)

    # modkit differential methylation call counts
    modc = _first(cfg.input_dir, "*%s*.modcall.vcf" % cfg.tumour_id)
    if modc:
        try:
            out["dmr"] = {"modcall_source": os.path.basename(modc),
                          "n_records": count_vcf_records(modc)}
        except Exception:
            pass
    return out


_ALLELE_RE = re.compile(r"[A-Za-z][A-Za-z0-9_-]*\*[0-9A-Za-z:]+|\*[0-9]+[A-Za-z]*")


def _parse_specimmune_table(path: str) -> List[Dict[str, Any]]:
    """SpecImmune result reader (line-based, robust to ragged rows).

    Real .final.type.1.result.txt is TAB-separated LONG format:
      Locus  Chromosome  Genotype  Match_info  Reads_num  Step1_type  One_guess
    ONE ROW PER HAPLOTYPE (Chromosome 1/2). **One_guess** (the LAST field) is the
    single best call for that haplotype; Genotype is the ';'-separated candidate
    list. Some rows are ragged (a missing tab merges Match_info+Reads_num), so we
    parse by splitting on tab and take One_guess by its column index *or* the last
    field, rather than trusting pandas column alignment. Also handles a WIDE
    allele_1/allele_2 layout (IG_TR)."""
    if not path:
        return []
    try:
        raw = [l.rstrip("\n") for l in open(path, encoding="utf-8", errors="replace")]
    except Exception as exc:
        LOG.warning("SpecImmune read failed (%s): %s", path, exc)
        return []
    lines = [l for l in raw if l.strip()]
    i = 0
    while i < len(lines) and lines[i].lstrip().startswith("#") and "\t" not in lines[i]:
        i += 1
    if i >= len(lines):
        return []
    header = [h.strip().lstrip("#").lower() for h in lines[i].split("\t")]

    def idx(*names):
        for n in names:
            if n in header:
                return header.index(n)
        return None

    i_gene = idx("locus", "gene", "hla", "region")
    i_chrom = idx("chromosome", "chrom", "hap", "haplotype")
    i_geno = idx("genotype", "allele", "call")
    i_guess = idx("one_guess", "best_guess", "final_type")
    i_a1 = idx("allele_1", "allele1", "a1")
    i_a2 = idx("allele_2", "allele2", "a2")
    i_extra = idx("phenotype", "diplotype", "metabolizer", "metaboliser",
                  "function", "activity_score")
    if i_gene is None:
        i_gene = 0

    def tok(s):
        m = _ALLELE_RE.findall(str(s or ""))
        return m[0] if m else None

    body = lines[i + 1:]

    # WIDE layout (IG_TR): explicit allele columns
    if i_a1 is not None and i_a2 is not None:
        rows = []
        for l in body:
            f = l.split("\t")
            if len(f) <= max(i_a1, i_a2):
                continue
            rows.append({"gene": f[i_gene] if i_gene < len(f) else None,
                         "allele_1": f[i_a1], "allele_2": f[i_a2],
                         "extra": (f[i_extra] if i_extra is not None and i_extra < len(f) else None)})
        if any((r["allele_1"] or "").strip() or (r["allele_2"] or "").strip() for r in rows):
            return rows

    # LONG layout: One_guess (last field) is the single best allele per haplotype
    from collections import OrderedDict
    per: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
    for l in body:
        f = l.split("\t")
        if not f or not f[0].strip():
            continue
        g = f[i_gene].strip() if i_gene < len(f) else ""
        if not g or g.lower() in ("nan", "na", "."):
            continue
        best = f[i_guess] if (i_guess is not None and i_guess < len(f)) else None
        if not best or not str(best).strip() or str(best).strip().upper() == "NA":
            best = f[-1]                       # One_guess is conventionally the last field
        a = tok(best)
        if a is None and i_geno is not None and i_geno < len(f):
            a = tok(f[i_geno])                 # fall back to the Genotype candidate list
        chrom = f[i_chrom].strip() if (i_chrom is not None and i_chrom < len(f)) else ""
        d = per.setdefault(g, {"a1": None, "a2": None, "spill": [], "extra": None})
        if chrom == "1" and d["a1"] is None:
            d["a1"] = a
        elif chrom == "2" and d["a2"] is None:
            d["a2"] = a
        elif a:
            d["spill"].append(a)
        if i_extra is not None and i_extra < len(f) and not d["extra"]:
            ev = f[i_extra].strip()
            if ev and ev not in ("", ".", "NA", "nan"):
                d["extra"] = ev
    out = []
    for g, d in per.items():
        a1, a2, spill = d["a1"], d["a2"], d["spill"]
        if a1 is None and spill:
            a1 = spill.pop(0)
        if a2 is None and spill:
            a2 = spill.pop(0)
        out.append({"gene": g, "allele_1": a1, "allele_2": a2, "extra": d["extra"]})
    return out

HLA_CLASSICAL = {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRA"}

# Compact CPIC-style allele-function map for the main pharmacogenes. Values are
# activity scores (CYP2D6) or functional classes used to derive a metaboliser
# phenotype. Star-allele suffixes are stripped to the core (*4xN -> *4).
_CYP_NOFUNC = {
    "CYP2D6": {"3", "4", "5", "6", "7", "8", "11", "12", "13", "14", "15", "18", "19", "20", "36", "38", "40", "42", "44", "56", "62"},
    "CYP2C19": {"2", "3", "4", "5", "6", "7", "8"},
    "CYP2C9": {"3", "6", "8", "11", "13"},
    "DPYD": {"2A", "13"},
    "TPMT": {"2", "3A", "3B", "3C", "4"},
    "UGT1A1": {"6", "27", "28", "37"},
    "SLCO1B1": {"5", "15"},
}
_CYP_DECREASED = {
    "CYP2D6": {"9", "10", "17", "29", "41", "49", "59"},
    "CYP2C19": set(),
    "CYP2C9": {"2", "5", "9", "12"},
    "DPYD": {"9A"},
    "TPMT": set(),
    "UGT1A1": set(),
    "SLCO1B1": set(),
}
_CYP_INCREASED = {"CYP2C19": {"17"}, "CYP2D6": set()}


def _cyp_phenotype(gene: str, a1: Optional[str], a2: Optional[str]) -> Optional[str]:
    """Derive a metaboliser phenotype from a star-allele diplotype (best-effort;
    review against CPIC/PharmGKB)."""
    gene = (gene or "").upper()
    if gene not in _CYP_NOFUNC:
        return None

    def core(a):
        m = re.search(r"\*([0-9]+[A-Z]*)", str(a or ""))
        return m.group(1) if m else None

    def fclass(al):
        c = core(al)
        if c is None:
            return None
        if c in _CYP_INCREASED.get(gene, set()):
            return "inc"
        if c in _CYP_NOFUNC.get(gene, set()):
            return "none"
        if c in _CYP_DECREASED.get(gene, set()):
            return "dec"
        return "norm"          # *1 and unlisted default to normal function
    f = [fclass(a1), fclass(a2)]
    if None in f:
        return None
    none_n = f.count("none")
    dec_n = f.count("dec")
    inc_n = f.count("inc")
    if none_n == 2:
        return "Poor metaboliser"
    if none_n == 1 and dec_n == 1:
        return "Poor–intermediate metaboliser"
    if none_n == 1 or dec_n == 2:
        return "Intermediate metaboliser"
    if inc_n == 2:
        return "Ultrarapid metaboliser"
    if inc_n == 1 and none_n == 0 and dec_n == 0:
        return "Rapid metaboliser"
    if dec_n == 1:
        return "Normal–intermediate metaboliser"
    return "Normal metaboliser"


def parse_hla_neoantigen(cfg: Config) -> Dict[str, Any]:
    """SpecImmune immune-gene typing (HLA / KIR / CYP / IG-TR) + pVACseq yield."""
    out: Dict[str, Any] = {"hla": {"available": False, "alleles": []},
                           "kir": {"available": False}, "cyp": {"available": False},
                           "ig_tr": {"available": False}, "neoantigen": {"available": False}}

    # SpecImmune layout varies: HLA under SpecImmune/<sample>/ or SpecImmune/HLA/<sample>/;
    # KIR/CYP/IG_TR under SpecImmune/<FAMILY>/<sample>/.
    families = {
        "hla": ("SpecImmune/*/*.HLA.final.type.1.result.txt",
                "SpecImmune/*/*/*.HLA.final.type.1.result.txt",
                "SpecImmune/*/*.HLA.type.result.txt",
                "SpecImmune/*/*/*.HLA.type.result.txt",
                "SpecImmune/*/*.HLA.final.type.result.formatted.txt",
                "SpecImmune/*/*/*.HLA.final.type.result.formatted.txt"),
        "kir": ("SpecImmune/KIR/*/*.KIR.final.type.1.result.txt",
                "SpecImmune/KIR/*/*.KIR.type.result.txt",
                "SpecImmune/KIR/*/*.KIR.final.type.result.formatted.txt",
                "SpecImmune/*/*.KIR.type.result.txt"),
        "cyp": ("SpecImmune/CYP/*/*.CYP.final.type.1.result.txt",
                "SpecImmune/CYP/*/*.CYP.type.result.txt",
                "SpecImmune/CYP/*/*.CYP.final.type.result.formatted.txt",
                "SpecImmune/*/*.CYP.type.result.txt"),
        "ig_tr": ("SpecImmune/IG_TR/*/*.IG_TR_typing_result.txt",
                  "SpecImmune/IG_TR/*/*.IG_TR.final.type.1.result.txt",
                  "SpecImmune/*/*.IG_TR_typing_result.txt"),
    }
    for fam, pats in families.items():
        path = _first(cfg.input_dir, *pats)
        rows = _parse_specimmune_table(path) if path else []
        if fam == "hla":
            # keep only the classical loci Andrew wants (A/B/C/DRB1/DRA)
            hla_rows = [r for r in rows if str(r.get("gene") or "").upper() in HLA_CLASSICAL]
            if not hla_rows:                 # fall back to all if none matched the filter
                hla_rows = rows
            out["hla"] = {"available": bool(hla_rows),
                          "alleles": [{"locus": r["gene"], "allele_1": r["allele_1"],
                                       "allele_2": r["allele_2"]} for r in hla_rows],
                          "source": os.path.relpath(path, cfg.input_dir) if path else None}
        else:
            if fam == "cyp":
                for r in rows:               # add metaboliser phenotype where missing
                    if not r.get("extra"):
                        r["extra"] = _cyp_phenotype(r.get("gene"), r.get("allele_1"), r.get("allele_2"))
            out[fam] = {"available": bool(rows), "rows": rows, "n": len(rows),
                        "source": os.path.relpath(path, cfg.input_dir) if path else None}

    # pVACseq aggregated reports + binder summary
    # pVACseq default --binding-threshold is 500 nM (binder); <50 nM is a strong binder.
    BINDER_NM, STRONG_NM = 500.0, 50.0

    def _count_binders(mhc_dir: str) -> Dict[str, Any]:
        agg = _first(cfg.input_dir, "pvacseq/%s/*.all_epitopes.aggregated.tsv" % mhc_dir,
                     "pvacseq/%s/*aggregated.tsv" % mhc_dir,
                     "pvacseq/%s/*.filtered.tsv" % mhc_dir)
        info: Dict[str, Any] = {"available": bool(agg)}
        if agg and pd is not None:
            df = _read_table(agg)
            if df is not None and len(df):
                info["n_candidates"] = int(len(df))
                c_gene = _pick_col(df, "Gene Name", "Gene", "gene")
                c_pep = _pick_col(df, "Best Peptide", "MT Epitope Seq", "peptide")
                c_ic = _pick_col(df, "IC50 MT", "Best MT IC50 Score", "Median MT IC50 Score")
                ics = [_to_num(v) for v in df[c_ic]] if c_ic else []
                ics = [v for v in ics if v is not None]
                info["n_binders"] = sum(1 for v in ics if v <= BINDER_NM)
                info["n_strong"] = sum(1 for v in ics if v <= STRONG_NM)
                top = []
                dft = df.copy()
                if c_ic:
                    dft["_ic"] = dft[c_ic].map(_to_num)
                    dft = dft.sort_values("_ic")
                for _, r in dft.head(15).iterrows():
                    ic = _to_num(r.get(c_ic)) if c_ic else None
                    top.append({"gene": r.get(c_gene) if c_gene else None,
                                "peptide": r.get(c_pep) if c_pep else None,
                                "ic50_nM": ic,
                                "strong": bool(ic is not None and ic <= STRONG_NM)})
                info["top"] = top
                info["source"] = os.path.basename(agg)
        return info

    ci = _count_binders("MHC_Class_I")
    cii = _count_binders("MHC_Class_II")
    out["neoantigen"] = {"available": ci.get("available") or cii.get("available"),
                         "class_I": ci, "class_II": cii,
                         "binder_threshold_nM": BINDER_NM, "strong_threshold_nM": STRONG_NM,
                         "summary_html": _first(cfg.input_dir, "*_pvac.vcf_summary.html")}
    return out


def parse_wakhan_plots(cfg: Config) -> Dict[str, Any]:
    out = {"heatmap_html": _first(cfg.input_dir, "wakhan/*heatmap_ploidy_purity.html")}
    sol = _wakhan_best_solution(cfg.input_dir)
    if sol:
        out["solution"] = os.path.basename(sol)
        # interactive HTML (embedded in the HTML report)
        out["details_html"] = _first(sol, "*genome_copynumbers_details.html")
        out["subclonal_html"] = _first(sol, "*genome_copynumbers_breakpoints_subclonal.html")
        out["cn_html"] = (out["details_html"] or
                          _first(sol, "*genome_copynumbers_breakpoints.html",
                                 "*genome_copynumbers*.html"))
        # Wakhan's own PDF plots (used for the PDF report — far better than a
        # snapshot). They may sit in a DIFFERENT solution dir from the best one
        # (e.g. if PDF generation was enabled on an earlier run), so look in the
        # best solution first, then anywhere under wakhan/.
        def _pdf(pattern: str) -> Optional[str]:
            return (_first(sol, pattern)
                    or _first(cfg.input_dir, f"wakhan/*/{pattern}")
                    or _first(cfg.input_dir, f"wakhan/{pattern}"))

        out["details_pdf"] = _pdf("*genome_copynumbers_details.pdf")
        out["subclonal_pdf"] = _pdf("*genome_copynumbers_breakpoints_subclonal.pdf")
        out["breakpoints_pdf"] = _pdf("*genome_copynumbers_breakpoints.pdf")
        out["genes_pdf"] = _pdf("*genes_genome.pdf")
        for k in ("details_pdf", "subclonal_pdf"):
            if out.get(k) and os.path.dirname(out[k]) != sol:
                LOG.info("Using Wakhan PDF from %s (best solution is %s)",
                         os.path.basename(os.path.dirname(out[k])), os.path.basename(sol))
    return out


_PLOTLY_JS: Optional[str] = None
_PLOTLY_JS_TRIED = False
_PLOTLY_CDN_RE = re.compile(
    r'<script\b[^>]*\bsrc=["\']([^"\']*plotly[^"\']*)["\'][^>]*>\s*</script>', re.I)


def _plotly_bundle(cdn_url: Optional[str] = None) -> Optional[str]:
    """The plotly.js source for OFFLINE embedding. Order: (1) the installed
    plotly package (no network), (2) a local cache file, (3) a one-off CDN fetch
    that gets cached for reuse. None if all fail."""
    global _PLOTLY_JS, _PLOTLY_JS_TRIED
    if _PLOTLY_JS is not None or _PLOTLY_JS_TRIED:
        return _PLOTLY_JS
    _PLOTLY_JS_TRIED = True
    for mod_name in ("plotly.offline", "plotly.io"):
        try:
            mod = __import__(mod_name, fromlist=["get_plotlyjs"])
            if hasattr(mod, "get_plotlyjs"):
                _PLOTLY_JS = mod.get_plotlyjs()
                return _PLOTLY_JS
        except Exception:
            pass
    cache = os.path.join(os.path.expanduser("~"), ".cache", "tn_report", "plotly.min.js")
    try:
        if os.path.exists(cache):
            with open(cache, encoding="utf-8") as fh:
                _PLOTLY_JS = fh.read()
            return _PLOTLY_JS
    except Exception:
        pass
    if cdn_url:
        try:
            import urllib.request
            if cdn_url.startswith("//"):
                cdn_url = "https:" + cdn_url
            with urllib.request.urlopen(cdn_url, timeout=20) as r:
                js = r.read().decode("utf-8", "replace")
            os.makedirs(os.path.dirname(cache), exist_ok=True)
            with open(cache, "w", encoding="utf-8") as fh:
                fh.write(js)
            _PLOTLY_JS = js
            return _PLOTLY_JS
        except Exception as exc:
            LOG.warning("plotly.js offline-inline unavailable (CDN fetch failed: %s). "
                        "Install the 'plotly' package for a fully offline report.", exc)
    return None


def _inline_plotly(doc: str) -> Tuple[str, bool]:
    """Replace any CDN plotly <script src> with an inline bundle so the plot
    renders with no network. Returns (doc, offline_ok). If the doc already ships
    plotly inline, it's returned unchanged as already-offline."""
    urls = _PLOTLY_CDN_RE.findall(doc)
    if not urls:
        return doc, True
    bundle = _plotly_bundle(urls[0])
    if not bundle:
        return doc, False
    return _PLOTLY_CDN_RE.sub(lambda m: "<script>" + bundle + "</script>", doc), True


def _plotly_fig_from_html(path: str) -> Optional[Tuple[List[Dict], Dict]]:
    """Pull (data, layout) out of a plotly HTML file. Handles:
      - plotly.py write_html:  Plotly.newPlot("div", [...], {...}, {...})
      - element form:          Plotly.newPlot(document.getElementById('x'), [...], {...})
      - htmlwidget (R/plotly): <script ... data-for="htmlwidget-...">{"x":{"data":[...],"layout":{...}}}
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            doc = fh.read()
    except Exception as exc:
        LOG.warning("Plot read failed (%s): %s", path, exc)
        return None
    dec = json.JSONDecoder()

    def _try_newplot(seg: str) -> Optional[Tuple[List[Dict], Dict]]:
        """Parse the args of one Plotly.newPlot( call, given the text after '('."""
        # skip the first argument: quoted id, bare identifier, or a JS expression
        # like document.getElementById('x') — consume to the first top-level comma
        m = re.match(r"""\s*("(?:[^"\\]|\\.)*"|'(?:[^'\\]|\\.)*')\s*,\s*""", seg)
        if m:
            rest = seg[m.end():]
        else:
            depth, k = 0, 0
            while k < len(seg):
                ch = seg[k]
                if ch in "([{":
                    depth += 1
                elif ch in ")]}":
                    depth -= 1
                    if depth < 0:
                        return None
                elif ch == "," and depth == 0:
                    break
                k += 1
            rest = seg[k + 1:]
        rest = rest.lstrip()
        try:
            data, end = dec.raw_decode(rest)
        except Exception:
            return None                      # not the real call (e.g. a JS-var form)
        if not isinstance(data, list) or not data:
            return None
        tail = rest[end:].lstrip().lstrip(",").lstrip()
        try:
            layout, _ = dec.raw_decode(tail)
        except Exception:
            layout = {}
        return data, (layout if isinstance(layout, dict) else {})

    # 1) Plotly.newPlot(...) — try EVERY occurrence. When the plotly.js bundle is
    # inlined (Wakhan's default), the first hit is a doc-string example inside the
    # minified library ("Plotly.newPlot(gd, data, layout, {...})"), not the plot.
    hits = [m.end() for m in re.finditer(r"Plotly\.newPlot\(", doc)]
    for start in reversed(hits):             # the real call is normally last
        got = _try_newplot(doc[start:])
        if got:
            return got
    if hits:
        LOG.warning("Found %d Plotly.newPlot call(s) in %s but none had inline JSON args",
                    len(hits), os.path.basename(path))

    # 2) htmlwidget JSON blob
    m = re.search(r'data-for="htmlwidget-[^"]+"\s*>\s*(\{.*?\})\s*</script>', doc, re.S)
    if m:
        try:
            x = json.loads(m.group(1)).get("x", {})
            data, layout = x.get("data"), x.get("layout") or {}
            if isinstance(data, list) and data:
                return data, (layout if isinstance(layout, dict) else {})
        except Exception as exc:
            LOG.warning("htmlwidget JSON parse failed (%s): %s",
                        os.path.basename(path), exc)

    LOG.warning("No plotly figure found in %s — snapshot skipped",
                os.path.basename(path))
    return None


def _plotly_title(layout: Dict) -> str:
    t = layout.get("title")
    if isinstance(t, dict):
        return str(t.get("text") or "")
    return str(t or "")


def _trim_png_border(png_bytes: bytes) -> bytes:
    """Trim uniform white margins from a rendered page (Wakhan's PDFs sit on a
    full page, so most of it is whitespace). No-op if Pillow isn't available."""
    try:
        from PIL import Image, ImageChops
        import io as _io
        im = Image.open(_io.BytesIO(png_bytes)).convert("RGB")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bbox = ImageChops.difference(im, bg).getbbox()
        if not bbox:
            return png_bytes
        pad = 6
        box = (max(bbox[0] - pad, 0), max(bbox[1] - pad, 0),
               min(bbox[2] + pad, im.width), min(bbox[3] + pad, im.height))
        out = _io.BytesIO()
        im.crop(box).save(out, format="PNG", optimize=True)
        return out.getvalue()
    except Exception:
        return png_bytes


def render_pdf_pages(path: str, dpi: int = 170, max_pages: int = 3) -> List[str]:
    """Rasterise a PDF's pages to base64 PNG data-URIs. This is Wakhan's OWN
    rendering of its plot, so it looks exactly like the PDF it produced —
    unlike re-plotting the traces. Backends tried in order:
      1. PyMuPDF (fitz)      2. poppler's pdftoppm      3. pdf2image
    Returns [] if none is available."""
    imgs: List[str] = []

    def _add(png: bytes):
        imgs.append("data:image/png;base64," +
                    base64.b64encode(_trim_png_border(png)).decode("ascii"))

    # 1) PyMuPDF
    try:
        import fitz  # PyMuPDF
        with fitz.open(path) as doc:
            for pg in list(doc)[:max_pages]:
                _add(pg.get_pixmap(dpi=dpi).tobytes("png"))
        if imgs:
            LOG.info("Wakhan PDF rendered via PyMuPDF: %s", os.path.basename(path))
            return imgs
    except ImportError:
        pass
    except Exception as exc:
        LOG.warning("PyMuPDF render failed (%s): %s", os.path.basename(path), exc)

    # 2) poppler pdftoppm (usually present on HPC via the poppler-utils module)
    if shutil.which("pdftoppm"):
        try:
            with tempfile.TemporaryDirectory() as td:
                pre = os.path.join(td, "pg")
                subprocess.run(["pdftoppm", "-png", "-r", str(dpi),
                                "-f", "1", "-l", str(max_pages), path, pre],
                               check=True, capture_output=True, timeout=180)
                for f in sorted(glob.glob(pre + "*.png"))[:max_pages]:
                    with open(f, "rb") as fh:
                        _add(fh.read())
            if imgs:
                LOG.info("Wakhan PDF rendered via pdftoppm: %s", os.path.basename(path))
                return imgs
        except Exception as exc:
            LOG.warning("pdftoppm render failed (%s): %s", os.path.basename(path), exc)

    # 3) pdf2image
    try:
        from pdf2image import convert_from_path
        import io as _io
        for im in convert_from_path(path, dpi=dpi, first_page=1, last_page=max_pages):
            buf = _io.BytesIO()
            im.save(buf, format="PNG")
            _add(buf.getvalue())
        if imgs:
            LOG.info("Wakhan PDF rendered via pdf2image: %s", os.path.basename(path))
            return imgs
    except ImportError:
        pass
    except Exception as exc:
        LOG.warning("pdf2image render failed (%s): %s", os.path.basename(path), exc)

    if not imgs:
        LOG.warning("Could not rasterise %s — install PyMuPDF (pip install pymupdf) "
                    "or poppler-utils to embed Wakhan's PDF plots",
                    os.path.basename(path))
    return imgs


# Which Wakhan PDF plots go into the PDF report. Add the commented lines to
# include the other two Wakhan produces.
WAKHAN_PDF_PLOTS = [
    ("Haplotype-specific copy number (details)", "details_pdf", "details_html"),
    ("Subclonal copy number & breakpoints", "subclonal_pdf", "subclonal_html"),
    # ("Copy number & breakpoints", "breakpoints_pdf", None),
    # ("Gene-level copy number", "genes_pdf", None),
]


def _wakhan_pdf_figures(wk: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Images for the PDF report's Wakhan section.

    Prefers Wakhan's OWN PDF plots (rasterised at high DPI) — these are exactly
    what Wakhan drew, so they match the interactive HTML. Falls back to
    snapshotting the plotly HTML only when the PDFs aren't present (i.e. Wakhan
    was run without PDF plot generation)."""
    out: List[Tuple[str, str]] = []
    for title, pdf_key, html_key in WAKHAN_PDF_PLOTS:
        pdf = wk.get(pdf_key)
        pages = render_pdf_pages(pdf) if pdf else []
        if pages:
            if len(pages) == 1:
                out.append((title, pages[0]))
            else:
                for n, img in enumerate(pages, 1):
                    out.append((f"{title} ({n}/{len(pages)})", img))
            continue
        if pdf and not pages:
            LOG.warning("Wakhan PDF %s could not be rasterised — falling back to a "
                        "snapshot of the interactive plot", os.path.basename(pdf))
        # no usable PDF — fall back to snapshotting the interactive plot
        html = wk.get(html_key) if html_key else None
        if html:
            img = make_plot_snapshot(html, title)
            if img:
                out.append((title, img))
    return out


def make_plot_snapshot(path: str, title: str = "") -> Optional[str]:
    """Fixed snapshot PNG (base64) of a plotly HTML plot, for the PDF.
    Tries kaleido (pixel-exact) first; falls back to re-rendering the plot's own
    extracted traces with matplotlib, which needs no browser."""
    fig_args = _plotly_fig_from_html(path)
    if not fig_args:
        return None
    data, layout = fig_args
    # 1) kaleido — pixel-exact, uses the plot's own plotly styling.
    # Wakhan's HTML embeds a plotly template that can contain trace types removed
    # from newer plotly (e.g. 'heatmapgl'), which fails go.Figure() validation —
    # so go through plotly.io with validate=False, and retry without the template.
    try:
        import plotly.io as pio
        import base64 as _b64
        lay = dict(layout or {})
        lay.update(width=1150, height=520,
                   margin=dict(l=60, r=20, t=50, b=50))
        attempts = [lay]
        if "template" in lay:
            no_tpl = {k: v for k, v in lay.items() if k != "template"}
            attempts.append(no_tpl)
        last = None
        for i, L in enumerate(attempts):
            try:
                png = pio.to_image({"data": data, "layout": L}, format="png",
                                   scale=2, validate=False)
                LOG.info("Snapshot via kaleido%s: %s",
                         " (template stripped)" if i else "", os.path.basename(path))
                return "data:image/png;base64," + _b64.b64encode(png).decode()
            except Exception as exc:
                last = exc
        raise last if last else RuntimeError("kaleido produced no image")
    except ImportError as exc:
        LOG.warning("plotly not importable (%s) — falling back to matplotlib snapshot", exc)
    except Exception as exc:
        LOG.warning("kaleido snapshot failed for %s (%s) — falling back to matplotlib. "
                    "For pixel-exact plots: pip install kaleido && plotly_get_chrome",
                    os.path.basename(path), str(exc).split("\n")[0][:120])
    # 2) matplotlib fallback — replots the plot's own trace data
    if not _HAVE_MPL:
        return None
    fig, ax = plt.subplots(figsize=(11.0, 4.6))
    plotted = 0
    for tr in data:
        if not isinstance(tr, dict):
            continue
        ttype = str(tr.get("type") or "scatter").lower()
        x, y = tr.get("x"), tr.get("y")
        if y is None:
            continue
        if x is None:
            x = list(range(len(y)))
        name = tr.get("name")
        colr = None
        for k in ("marker", "line"):
            c = (tr.get(k) or {}).get("color") if isinstance(tr.get(k), dict) else None
            if isinstance(c, str):
                colr = c
                break
        try:
            if ttype.startswith("bar"):
                ax.bar(x, y, color=colr, label=name, alpha=0.9)
            elif "markers" in str(tr.get("mode") or "") and "lines" not in str(tr.get("mode") or ""):
                ax.scatter(x, y, s=6, color=colr, label=name, alpha=0.85)
            else:
                ax.plot(x, y, lw=1.0, color=colr, label=name, alpha=0.9)
            plotted += 1
        except Exception:
            continue
    if not plotted:
        plt.close(fig)
        return None
    xt = (layout.get("xaxis") or {}).get("title")
    yt = (layout.get("yaxis") or {}).get("title")
    ax.set_xlabel(xt.get("text") if isinstance(xt, dict) else (xt or ""))
    ax.set_ylabel(yt.get("text") if isinstance(yt, dict) else (yt or ""))
    ax.set_title(_plotly_title(layout) or title, fontsize=10, loc="left")
    if any(isinstance(t, dict) and t.get("name") for t in data):
        ax.legend(fontsize=7, frameon=False, ncol=4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if len(str(max((len(str(v)) for tr in data for v in (tr.get("x") or [""])), default=0))) > 0:
        plt.setp(ax.get_xticklabels(), fontsize=7, rotation=90)
    return _png_b64(fig)


# ============================================================================
# Genome-wide copy-number figure (matplotlib -> base64 PNG)
# ============================================================================
_CHROM_SIZES_GRCH38 = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555,
    "5": 181538259, "6": 170805979, "7": 159345973, "8": 145138636,
    "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345,
    "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415,
}


def make_cn_figure(ascat: Dict[str, Any]) -> Optional[str]:
    """Genome-wide allele-specific CN (ASCAT), coloured by event type:
    gain (CN>2), loss (CN<2), homozygous deletion (CN=0), copy-neutral (CN=2),
    plus an LOH track underneath."""
    if not _HAVE_MPL or not ascat.get("segments"):
        return None
    order = list(_CHROM_SIZES_GRCH38.keys())
    offsets, cum = {}, 0
    for c in order:
        offsets[c] = cum
        cum += _CHROM_SIZES_GRCH38[c]
    genome_len = cum

    C_GAIN, C_LOSS, C_HOMDEL, C_NEUT, C_LOH = (
        "#b3261e", "#1f6feb", "#0b3d91", "#7a8696", "#7a5195")
    seen = set()
    max_cn = 2
    fig, ax = plt.subplots(figsize=(11, 3.2))
    for seg in ascat["segments"]:
        chrom = str(seg.get("chrom") or "").replace("chr", "")
        if chrom not in offsets:
            continue
        s = seg.get("start"); e = seg.get("end"); tcn = seg.get("total_cn")
        if s is None or e is None or tcn is None:
            continue
        x0 = offsets[chrom] + s
        x1 = offsets[chrom] + e
        if tcn == 0:
            colour, key = C_HOMDEL, "homdel"
        elif tcn > 2:
            colour, key = C_GAIN, "gain"
        elif tcn < 2:
            colour, key = C_LOSS, "loss"
        else:
            colour, key = C_NEUT, "neutral"
        seen.add(key)
        max_cn = max(max_cn, tcn)
        ax.hlines(tcn, x0, x1, colors=colour, linewidth=3.2)
        if seg.get("loh"):
            ax.hlines(-0.25, x0, x1, colors=C_LOH, linewidth=3.0)
            seen.add("loh")
    for c in order:
        ax.axvline(offsets[c], color="#e3e8ee", linewidth=0.6, zorder=0)
    ax.set_xticks([offsets[c] + _CHROM_SIZES_GRCH38[c] / 2 for c in order])
    ax.set_xticklabels(order, fontsize=7)
    ax.set_xlim(0, genome_len)
    ax.set_ylim(-0.6, max(4, min(max_cn + 1, 12)))
    ax.set_ylabel("Total copy number")
    ax.axhline(2, color="#aab4c0", linewidth=0.8, linestyle=":")
    ax.set_title("Genome-wide allele-specific copy number (ASCAT)", fontsize=10, loc="left")
    # legend keyed by event type
    from matplotlib.lines import Line2D
    handles = []
    for key, colour, label in [
            ("gain", C_GAIN, "Gain (CN > 2)"),
            ("loss", C_LOSS, "Loss (CN < 2)"),
            ("homdel", C_HOMDEL, "Homozygous deletion (CN = 0)"),
            ("neutral", C_NEUT, "Copy-neutral (CN = 2)"),
            ("loh", C_LOH, "LOH (minor CN = 0)")]:
        if key in seen:
            handles.append(Line2D([0], [0], color=colour, lw=3.2, label=label))
    if handles:
        ax.legend(handles=handles, fontsize=7, frameon=False, ncol=len(handles),
                  loc="upper center", bbox_to_anchor=(0.5, 1.14))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return _png_b64(fig)


def _cn_track_figure(segments: List[Dict[str, Any]], title: str,
                     value: str = "total_cn") -> Optional[str]:
    """Generic genome-wide CN track (total CN or log2 ratio) from a segment list."""
    if not _HAVE_MPL or not segments:
        return None
    order = list(_CHROM_SIZES_GRCH38.keys())
    offsets, cum = {}, 0
    for c in order:
        offsets[c] = cum
        cum += _CHROM_SIZES_GRCH38[c]
    is_log2 = value == "log2r"
    fig, ax = plt.subplots(figsize=(11, 3.0))
    for seg in segments:
        chrom = str(seg.get("chrom") or "").replace("chr", "")
        if chrom not in offsets:
            continue
        s, e, val = seg.get("start"), seg.get("end"), seg.get(value)
        if s is None or e is None or val is None:
            continue
        if is_log2:
            colour = "#b3261e" if val > 0.3 else ("#1f6feb" if val < -0.3 else "#7a8696")
        else:
            colour = "#b3261e" if val > 2 else ("#1f6feb" if val < 2 else "#7a8696")
        ax.hlines(val, offsets[chrom] + s, offsets[chrom] + e, colors=colour, linewidth=3.0)
        if not is_log2 and seg.get("minor_cn") == 0:
            ax.hlines(-0.25, offsets[chrom] + s, offsets[chrom] + e, colors="#7a5195", linewidth=2.6)
    for c in order:
        ax.axvline(offsets[c], color="#e3e8ee", linewidth=0.6, zorder=0)
    ax.set_xticks([offsets[c] + _CHROM_SIZES_GRCH38[c] / 2 for c in order])
    ax.set_xticklabels(order, fontsize=7)
    ax.set_xlim(0, cum)
    ax.axhline(0 if is_log2 else 2, color="#aab4c0", linewidth=0.8, linestyle=":")
    ax.set_ylabel("log2 ratio" if is_log2 else "Total copy number")
    ax.set_title(title, fontsize=10, loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return _png_b64(fig)


def make_wakhan_cn_figure(segments: List[Dict[str, Any]]) -> Optional[str]:
    return _cn_track_figure(segments, "Genome-wide copy number (Wakhan, haplotype-summed)")


def parse_pcgr_cna_segments(cfg: Config) -> List[Dict[str, Any]]:
    """PCGR segment-level CNA (pcgr/*.cna_segment.tsv.gz) — reproduces the data
    behind PCGR's genome-wide CNA plot."""
    seg = _first(cfg.input_dir, "pcgr/*.cna_segment.tsv.gz", "pcgr/*cna_segment*.tsv*")
    if not seg or pd is None:
        return []
    df = _read_table(seg)
    if df is None or not len(df):
        return []
    c_chr = _pick_col(df, "chrom", "chromosome", "Chromosome", "chr")
    c_s = _pick_col(df, "start", "segment_start", "segstart", "loc.start", "Start")
    c_e = _pick_col(df, "end", "segment_end", "segend", "loc.end", "End")
    c_l2 = _pick_col(df, "log2r", "logr", "log_r_ratio", "LogR", "seg.mean", "logR")
    c_tot = _pick_col(df, "total_cn", "cn_total", "CN", "nCN", "copy_number")
    out = []
    for _, r in df.iterrows():
        if not c_chr or not c_s or not c_e:
            break
        out.append({"chrom": r.get(c_chr), "start": _to_num(r.get(c_s)), "end": _to_num(r.get(c_e)),
                    "log2r": _to_num(r.get(c_l2)) if c_l2 else None,
                    "total_cn": _to_num(r.get(c_tot)) if c_tot else None})
    return out


def make_pcgr_cna_figure(segments: List[Dict[str, Any]]) -> Optional[str]:
    if not segments:
        return None
    has_log2 = any(s.get("log2r") is not None for s in segments)
    return _cn_track_figure(segments, "Genome-wide copy number (PCGR segments)",
                            value="log2r" if has_log2 else "total_cn")


def make_vaf_figure(snvs: List[Dict[str, Any]]) -> Optional[str]:
    if not _HAVE_MPL:
        return None
    vafs = [v["vaf"] for v in snvs if v.get("vaf") is not None]
    if not vafs:
        return None
    fig, ax = plt.subplots(figsize=(5.2, 2.8))
    ax.hist(vafs, bins=25, color="#0e6e76", alpha=0.85, edgecolor="white")
    ax.set_xlabel("Tumour VAF")
    ax.set_ylabel("Somatic variants")
    ax.set_title("Somatic VAF distribution", fontsize=10, loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return _png_b64(fig)


def make_signature_figure(sigs: List[Dict[str, Any]]) -> Optional[str]:
    if not _HAVE_MPL or not sigs:
        return None
    items = [(s["signature"] or "?", s["weight"] or 0.0)
             for s in sigs if (s.get("weight") or 0) > 0][:10]
    if not items:
        return None
    labels = [i[0] for i in items][::-1]
    vals = [i[1] for i in items][::-1]
    fig, ax = plt.subplots(figsize=(5.2, max(2.4, 0.42 * len(items))))
    ax.barh(labels, vals, color="#243b6b", alpha=0.9)
    ax.set_xlabel("Relative contribution")
    ax.set_title("Mutational signatures (PCGR)", fontsize=10, loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return _png_b64(fig)


def make_signature_pie_figure(sigs: List[Dict[str, Any]]) -> Optional[str]:
    """Coloured pie of signature components (top 8 + 'Other')."""
    if not _HAVE_MPL or not sigs:
        return None
    items = [(s.get("signature") or "?", float(s.get("weight") or 0.0))
             for s in sigs if (s.get("weight") or 0) > 0]
    if not items:
        return None
    items.sort(key=lambda t: -t[1])
    top, rest = items[:8], items[8:]
    labels = [t[0] for t in top]
    vals = [t[1] for t in top]
    if rest:
        labels.append("Other")
        vals.append(sum(t[1] for t in rest))
    try:
        cmap = matplotlib.colormaps["tab10"]      # mpl >= 3.6
    except Exception:
        cmap = plt.get_cmap("tab10")             # older mpl

    colors = [cmap(i % 10) for i in range(len(vals))]
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    wedges, _txt, autotxt = ax.pie(
        vals, colors=colors, startangle=90, counterclock=False,
        autopct=lambda p: f"{p:.0f}%" if p >= 5 else "",
        pctdistance=0.72, wedgeprops={"edgecolor": "white", "linewidth": 1.2},
        textprops={"fontsize": 8})
    for t in autotxt:
        t.set_color("white")
        t.set_fontweight("bold")
    ax.legend(wedges, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=8, frameon=False)
    ax.set_title("Signature composition", fontsize=10, loc="left")
    ax.axis("equal")
    return _png_b64(fig)


# Approximate TCGA median coding TMB (mutations/Mb) per primary site — used to
# place the sample's TMB in context (illustrative reference, not the PCGR data).
_TCGA_TMB_MEDIAN = {
    "Skin/melanoma": 13.0, "Lung (squamous)": 9.5, "Lung (adeno)": 7.0, "Bladder": 7.0,
    "Oesophagus": 4.6, "Stomach": 4.4, "Colon/rectum": 4.0, "Endometrium": 4.0,
    "Head & neck": 3.6, "Cervix": 3.0, "Liver": 2.5, "Glioblastoma": 2.2,
    "Ovary": 1.8, "Pancreas": 1.7, "Breast": 1.5, "Kidney": 1.2, "Prostate": 1.0,
    "Thyroid": 0.4, "Paediatric (typical)": 0.2,
}


def make_tmb_context_figure(tmb_value: Optional[float]) -> Optional[str]:
    """TMB in context: a violin per tumour type showing the spread of TCGA TMB
    (log-normal approximation around the published median), with the median
    marked, and this tumour drawn as a vertical line."""
    if not _HAVE_MPL or tmb_value is None or np is None:
        return None
    items = sorted(_TCGA_TMB_MEDIAN.items(), key=lambda kv: kv[1])
    labels = [k for k, _ in items]
    meds = [v for _, v in items]
    # TCGA TMB is heavily right-skewed; approximate each cohort as log-normal
    # about its median (sigma chosen to give a realistic ~1.5-decade spread).
    rng = np.random.default_rng(0)
    dists = [np.exp(rng.normal(np.log(max(m, 0.05)), 0.85, 400)) for m in meds]
    fig, ax = plt.subplots(figsize=(7.0, max(3.2, 0.34 * len(items))))
    parts = ax.violinplot(dists, positions=range(len(items)), vert=False,
                          widths=0.85, showextrema=False, showmedians=False)
    for pc in parts["bodies"]:
        pc.set_facecolor("#c9d4e4")
        pc.set_edgecolor("#9fb0c8")
        pc.set_alpha(0.85)
    # median tick per tumour type
    ax.scatter(meds, range(len(items)), s=14, color="#3f5a80", zorder=6,
               label="TCGA median (approx.)")
    ax.axvline(tmb_value, color="#b3261e", linewidth=1.8, zorder=7)
    ax.text(tmb_value, len(items) - 0.4, f" this tumour: {tmb_value:.1f}/Mb",
            color="#b3261e", fontsize=8.5, va="top",
            ha="left" if tmb_value < max(meds) * 2 else "right")
    ax.axvline(10, color="#7a5a12", linewidth=1.0, linestyle="--", zorder=5)
    ax.text(10, -1.0, "TMB-high (10)", color="#7a5a12", fontsize=7.5, ha="center")
    ax.set_xscale("log")
    ax.set_xlim(0.05, 300)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Tumour mutational burden (mutations / Mb, log scale)")
    ax.set_title("TMB in context of tumour types", fontsize=10, loc="left")
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return _png_b64(fig)


def make_loh_figure(loh: Dict[str, Any]) -> Optional[str]:
    """Genome-wide LOH track from the Severus LOH bed."""
    regions = loh.get("regions") if isinstance(loh, dict) else None
    if not _HAVE_MPL or not regions:
        return None
    order = list(_CHROM_SIZES_GRCH38.keys())
    offsets, cum = {}, 0
    for c in order:
        offsets[c] = cum
        cum += _CHROM_SIZES_GRCH38[c]
    fig, ax = plt.subplots(figsize=(11, 1.7))
    for r in regions:
        chrom = str(r.get("chrom") or "").replace("chr", "")
        if chrom not in offsets:
            continue
        s, e = r.get("start"), r.get("end")
        if s is None or e is None:
            continue
        ax.barh(0, e - s, left=offsets[chrom] + s, height=0.6, color="#7a5195")
    for c in order:
        ax.axvline(offsets[c], color="#e3e8ee", linewidth=0.6, zorder=0)
    ax.set_xticks([offsets[c] + _CHROM_SIZES_GRCH38[c] / 2 for c in order])
    ax.set_xticklabels(order, fontsize=7)
    ax.set_xlim(0, cum)
    ax.set_yticks([])
    ax.set_title("Genome-wide loss of heterozygosity (Severus)", fontsize=10, loc="left")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    return _png_b64(fig)


# ============================================================================
# Actionability roll-up (the MTB headline)
# ============================================================================
def build_actionability(report: Dict[str, Any]) -> Dict[str, Any]:
    items: List[Dict[str, Any]] = []

    def _tier_rank(t: Optional[str]) -> int:
        if not t:
            return 9
        t = str(t).upper()
        # explicit arabic tier digit (PCGR: TIER_1 .. TIER_4)
        m = re.search(r"\b(?:TIER[_ ]?)?([1-4])\b", t)
        if m:
            return int(m.group(1))
        # roman numerals only if the token is *exactly* roman (avoids the 'I' in 'TIER')
        roman = {"IV": 4, "III": 3, "II": 2, "I": 1}
        for tok in re.split(r"[^A-Z]+", t):
            if tok in roman:
                return roman[tok]
        return 9

    for v in report.get("pcgr", {}).get("snv_indel", {}).get("variants", []):
        r = _tier_rank(v.get("tier"))
        if r <= 2:
            items.append({"category": "Somatic SNV/indel",
                          "finding": f"{v.get('gene')} {v.get('hgvsp') or v.get('consequence') or ''}".strip(),
                          "tier": v.get("tier"), "rank": r,
                          "evidence": v.get("evidence")})
    for c in report.get("pcgr", {}).get("cna_genes", {}).get("items", []):
        r = _tier_rank(c.get("tier"))
        if r <= 2 or (c.get("event") and "amp" in str(c.get("event")).lower()):
            items.append({"category": "Copy number",
                          "finding": f"{c.get('gene')} {c.get('event') or ''} (CN={c.get('total_cn')})".strip(),
                          "tier": c.get("tier"), "rank": r,
                          "evidence": c.get("evidence")})

    msi = report.get("msi", {})
    if msi.get("status") == "MSI":
        items.append({"category": "Biomarker", "finding": "MSI / microsatellite-unstable (PCGR)",
                      "tier": "Biomarker", "rank": 1,
                      "evidence": "Pembrolizumab / immune checkpoint inhibition"})
    tmb = report.get("pcgr", {}).get("tmb", {})
    if tmb.get("tmb_per_mb") is not None and tmb["tmb_per_mb"] >= 10:
        items.append({"category": "Biomarker",
                      "finding": f"TMB-High ({tmb['tmb_per_mb']:.1f}/Mb)",
                      "tier": "Biomarker", "rank": 1,
                      "evidence": "Possible ICI sensitivity (tumour-type dependent)"})
    hrd = report.get("hrd", {})
    if hrd.get("available") and hrd.get("hr_deficient"):
        _hs = hrd.get("hrd_sum")
        items.append({"category": "Biomarker",
                      "finding": "HRD / homologous-recombination deficient"
                                 + (f" (HRD-sum {_hs:.0f}, threshold \u226542)"
                                    if isinstance(_hs, (int, float)) else ""),
                      "tier": "Biomarker", "rank": 1,
                      "evidence": "Possible PARP-inhibitor / platinum sensitivity "
                                  "(tumour-type dependent; scarHRD is not a validated "
                                  "companion diagnostic — confirm with an accredited assay)"})
    for gv in report.get("cpsr", {}).get("variants", []):
        if gv.get("reportable"):
            items.append({"category": "Germline (CPSR)",
                          "finding": f"{gv.get('gene')} {gv.get('hgvsp') or ''} — {gv.get('classification')}".strip(),
                          "tier": "Germline P/LP", "rank": 1,
                          "evidence": gv.get("trait")})
    for caller, d in report.get("sv", {}).get("callers", {}).items():
        for fus in d.get("fusions", []):
            items.append({"category": "Structural / fusion",
                          "finding": f"{fus.get('annotation')}", "tier": caller,
                          "rank": 2, "evidence": "Review for actionable fusion"})

    items.sort(key=lambda x: x.get("rank", 9))
    return {"n": len(items), "items": items}


def build_qc_flags(report: Dict[str, Any], cfg: Config) -> List[Dict[str, str]]:
    flags: List[Dict[str, str]] = []
    pp = report.get("purity_ploidy", {}).get("consensus", {})
    if pp.get("divergent"):
        flags.append({"level": "warn",
                      "message": ("Purity estimates diverge >25% across callers "
                                  "(ASCAT / wakhan / SAVANA) — interpret CN-derived "
                                  "metrics with caution.")})
    cov = report.get("coverage", {})
    for role, s in cov.get("samples", {}).items():
        d = s.get("mean_depth")
        if d is not None and d < 25:
            flags.append({"level": "warn",
                          "message": f"{role.capitalize()} mean depth {d:.1f}x is below 25x; "
                                     "somatic sensitivity may be reduced."})
    if not report.get("pcgr", {}).get("available"):
        flags.append({"level": "error",
                      "message": "PCGR outputs not found — somatic annotation incomplete."})
    return flags


# ============================================================================
# Recommended additional analyses for the "ideal" T/N read-out
# ============================================================================
RECOMMENDED_ADDITIONS = [
    # Previously-recommended analyses (sample identity, coverage QC, HLA-LOH, HRD,
    # clonal architecture, orthogonal signatures) are now integrated as live
    # sections. Decoil ecDNA reconstruction intentionally omitted.
]


# ============================================================================
# Assembly
# ============================================================================

# ===========================================================================
# NEW ANALYTICS: caller concordance, genome features, PCGR-MSI, SV circos
# ===========================================================================
# --------------------------------------------------------------------------- #
GENE_COORDS: Dict[str, Tuple[str, int, int]] = {
    "WT1": ("chr11", 32_387_775, 32_435_564), "MYCN": ("chr2", 15_940_550, 15_947_007),
    "CTNNB1": ("chr3", 41_194_741, 41_260_096), "AMER1": ("chrX", 64_186_888, 64_207_180),
    "TP53": ("chr17", 7_668_421, 7_687_550), "SIX1": ("chr14", 60_180_543, 60_186_960),
    "SIX2": ("chr2", 45_010_307, 45_015_348), "DROSHA": ("chr5", 31_400_601, 31_532_181),
    "DGCR8": ("chr22", 20_080_223, 20_111_872), "DICER1": ("chr14", 95_086_228, 95_158_263),
    "MYC": ("chr8", 127_735_434, 127_742_951), "RB1": ("chr13", 48_303_751, 48_481_890),
    "PTEN": ("chr10", 87_863_113, 87_971_930), "CDKN2A": ("chr9", 21_967_752, 21_995_301),
    "EGFR": ("chr7", 54_999_240, 55_207_337), "ERBB2": ("chr17", 39_688_094, 39_728_658),
    "ALK": ("chr2", 29_192_774, 29_921_586), "MDM2": ("chr12", 68_807_330, 68_850_686),
    "CDK4": ("chr12", 57_747_727, 57_756_013), "MYCL": ("chr1", 39_894_178, 39_898_756),
    "FGFR1": ("chr8", 38_411_138, 38_468_834), "KRAS": ("chr12", 25_205_246, 25_250_936),
    "NMYC": ("chr2", 15_940_550, 15_947_007), "TERT": ("chr5", 1_253_167, 1_295_184),
}

# Curated actionable / oncogenic fusion gene-pairs (unordered), with relevance.
ACTIONABLE_FUSIONS: Dict[frozenset, Dict[str, str]] = {
    frozenset(("EWSR1", "FLI1")): {"relevance": "diagnostic", "note": "Ewing sarcoma"},
    frozenset(("EWSR1", "WT1")): {"relevance": "diagnostic", "note": "Desmoplastic small round cell tumour"},
    frozenset(("PAX3", "FOXO1")): {"relevance": "diagnostic/prognostic", "note": "Alveolar rhabdomyosarcoma"},
    frozenset(("PAX7", "FOXO1")): {"relevance": "diagnostic/prognostic", "note": "Alveolar rhabdomyosarcoma"},
    frozenset(("ETV6", "NTRK3")): {"relevance": "therapeutic", "note": "NTRK inhibitor (larotrectinib/entrectinib)"},
    frozenset(("TPM3", "NTRK1")): {"relevance": "therapeutic", "note": "NTRK inhibitor"},
    frozenset(("EML4", "ALK")): {"relevance": "therapeutic", "note": "ALK inhibitor"},
    frozenset(("KIAA1549", "BRAF")): {"relevance": "diagnostic", "note": "Pilocytic astrocytoma"},
    frozenset(("BCR", "ABL1")): {"relevance": "therapeutic", "note": "TKI (imatinib)"},
    frozenset(("CIC", "DUX4")): {"relevance": "diagnostic", "note": "CIC-rearranged sarcoma"},
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _num(x):
    try:
        return float(str(x).strip())
    except (ValueError, AttributeError, TypeError):
        return None


def _norm_chr(c: Any) -> str:
    c = str(c).strip()
    return c if c.lower().startswith("chr") else "chr" + c


def cn_at_locus(segments: List[Dict], chrom: str, start: int, end: int) -> Optional[float]:
    """Total CN of the segment with greatest overlap of [start,end] on chrom."""
    chrom = _norm_chr(chrom)
    best, best_ov = None, 0
    for s in segments:
        if _norm_chr(s["chrom"]) != chrom or s.get("total_cn") is None:
            continue
        ov = min(end, s["end"]) - max(start, s["start"])
        if ov > best_ov:
            best, best_ov = s["total_cn"], ov
    return best


# --------------------------------------------------------------------------- #
# CN segment loaders (one per caller)
# --------------------------------------------------------------------------- #
def load_ascat_segments(path: str) -> List[Dict]:
    import csv
    out = []
    with _open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        low = {c.lower(): c for c in (rd.fieldnames or [])}
        def col(*names):
            for n in names:
                if n in low: return low[n]
        c_chr, c_s, c_e = col("chr", "chromosome"), col("startpos", "start"), col("endpos", "end")
        c_maj, c_min = col("nmajor", "namajor", "major_cn"), col("nminor", "nbminor", "minor_cn")
        for r in rd:
            maj, mn = _num(r.get(c_maj)), _num(r.get(c_min))
            if None in (maj, mn): continue
            out.append({"chrom": r[c_chr], "start": int(_num(r[c_s])), "end": int(_num(r[c_e])),
                        "total_cn": maj + mn, "minor_cn": mn})
    return out


def load_savana_cn(path: str) -> List[Dict]:
    """SAVANA *_segmented_absolute_copy_number.tsv (allele-specific absolute CN)."""
    import csv
    out = []
    with _open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        low = {c.lower(): c for c in (rd.fieldnames or [])}
        def col(*names):
            for n in names:
                if n in low: return low[n]
        c_chr = col("chromosome", "chr", "chrom")
        c_s, c_e = col("start", "startpos", "start_position"), col("end", "endpos", "end_position")
        c_tot = col("copynumber", "total_cn", "absolute_copy_number", "copy_number", "cn")
        c_min = col("minor_allele_copynumber", "minorallelecopynumber", "minor_cn", "nminor")
        for r in rd:
            tot = _num(r.get(c_tot))
            if c_tot is None or tot is None:
                # some SAVANA outputs split haplotype CN: sum them
                hp = [c for c in low if "copynumber" in c]
                vals = [_num(r.get(low[c])) for c in hp]
                vals = [v for v in vals if v is not None]
                tot = sum(vals) if vals else None
            if tot is None: continue
            out.append({"chrom": r[c_chr], "start": int(_num(r[c_s])), "end": int(_num(r[c_e])),
                        "total_cn": tot, "minor_cn": _num(r.get(c_min)) if c_min else None})
    return out


def _read_wakhan_bed(path: str) -> List[Dict]:
    segs = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:                      # some Wakhan beds are space-separated
                f = line.split()
            if len(f) < 4:
                continue
            s, e = _num(f[1]), _num(f[2])        # accept int OR float coordinates
            if s is None or e is None:           # header / malformed row
                continue
            nums = [v for v in (_num(x) for x in f[3:]) if v is not None]
            if not nums:
                continue
            # CN state is an integer; prefer the first integer col (avoids coverage floats)
            ints = [v for v in nums if float(v).is_integer() and 0 <= v <= 200]
            cn = ints[0] if ints else nums[-1]
            segs.append({"chrom": f[0], "start": int(s), "end": int(e), "cn": cn})
    return segs


def _cn_point(segs, chrom, pos):
    for s in segs:
        if s["chrom"] == chrom and s["start"] <= pos < s["end"]:
            return s["cn"]
    return None


def _combine_haplotype_cn(s1: List[Dict], s2: List[Dict]) -> List[Dict]:
    """Sum haplotype-specific CN (HP_1 + HP_2) into total-CN segments."""
    out = []
    for c in sorted(set(x["chrom"] for x in s1 + s2)):
        pts = sorted(set([x["start"] for x in s1 + s2 if x["chrom"] == c] +
                         [x["end"] for x in s1 + s2 if x["chrom"] == c]))
        for a, b in zip(pts, pts[1:]):
            if b <= a:
                continue
            mid = (a + b) // 2
            cn1, cn2 = _cn_point(s1, c, mid), _cn_point(s2, c, mid)
            if cn1 is None and cn2 is None:
                continue
            out.append({"chrom": c, "start": a, "end": b,
                        "total_cn": (cn1 or 0) + (cn2 or 0),
                        "minor_cn": min(cn1 or 0, cn2 or 0)})
    return out


def _wakhan_best_solution(input_dir: str) -> Optional[str]:
    """Return the Wakhan solution directory for the best-ranked solution."""
    sols = sorted(glob.glob(os.path.join(input_dir, "wakhan", "*", "bed_output")))
    sols = [os.path.dirname(s) for s in sols]
    if not sols:
        return None
    if len(sols) == 1:
        return sols[0]
    ranks = _first(input_dir, "wakhan/solutions_ranks.tsv")
    if ranks:
        try:
            with _open(ranks) as fh:
                lines = [l for l in fh if l.strip()]
            # the top data row names the best solution (ploidy_purity_conf token)
            for l in lines[1:]:
                tok = re.search(r"[\d.]+_[\d.]+_[\d.]+", l)
                if tok:
                    for s in sols:
                        if tok.group(0) in os.path.basename(s):
                            return s
                    break
        except Exception:
            pass
    return sols[0]


def load_wakhan_cn(input_dir: str) -> List[Dict]:
    """Wakhan CN from the best solution's bed_output. Haplotype-specific CN is in
    separate *_copynumbers_segments_HP_1.bed / _HP_2.bed files, summed to a total.
    Falls back to a single combined bed if the HP split isn't present."""
    sol = _wakhan_best_solution(input_dir)
    if sol:
        hp1 = _first(sol, "bed_output/*copynumbers_segments_HP_1.bed",
                     "bed_output/*copynumbers*HP_1.bed")
        hp2 = _first(sol, "bed_output/*copynumbers_segments_HP_2.bed",
                     "bed_output/*copynumbers*HP_2.bed")
        if hp1 and hp2:
            try:
                combined = _combine_haplotype_cn(_read_wakhan_bed(hp1), _read_wakhan_bed(hp2))
                if combined:
                    return combined
                LOG.warning("Wakhan HP beds parsed to 0 segments (%s); trying fallback",
                            os.path.basename(hp1))
            except Exception as exc:
                LOG.warning("Wakhan HP combine failed: %s", exc)
    # fallback: any single CN bed, summing numeric CN columns
    cands = (glob.glob(os.path.join(input_dir, "wakhan", "**", "bed_output", "*copynumber*.bed"), recursive=True)
             or glob.glob(os.path.join(input_dir, "wakhan", "**", "*copynumber*.bed"), recursive=True))
    cands = [p for p in cands if "subclonal" not in os.path.basename(p).lower() and "loh" not in os.path.basename(p).lower()]
    out = []
    for path in sorted(cands)[:1]:
        for s in _read_wakhan_bed(path):
            out.append({"chrom": s["chrom"], "start": s["start"], "end": s["end"],
                        "total_cn": s["cn"]})
    return out


# --------------------------------------------------------------------------- #
# CN concordance across callers, per driver gene
# --------------------------------------------------------------------------- #
def build_cn_concordance(driver_genes: List[str], callers: Dict[str, List[Dict]],
                         tol: float = 0.5) -> Dict[str, Any]:
    rows = []
    n_discord = 0
    for gene in sorted(set(g for g in driver_genes if g in GENE_COORDS)):
        chrom, s, e = GENE_COORDS[gene]
        per = {}
        for name, segs in callers.items():
            if segs:
                cn = cn_at_locus(segs, chrom, s, e)
                if cn is not None:
                    per[name] = round(cn, 2)
        if len(per) < 2:
            continue
        vals = list(per.values())
        discord = (max(vals) - min(vals)) > tol
        # direction disagreement (gain vs loss vs neutral) is always flagged
        signs = {("gain" if v > 2 + tol else "loss" if v < 2 - tol else "neutral") for v in vals}
        discord = discord or len(signs) > 1
        if discord:
            n_discord += 1
        rows.append({"gene": gene, "locus": f"{chrom}:{s}-{e}", "per_caller": per,
                     "discordant": discord})
    return {"rows": rows, "n_discordant": n_discord, "callers": list(callers.keys())}


# --------------------------------------------------------------------------- #
# Genome features: WGD, 1p19q co-deletion, 1q gain
# --------------------------------------------------------------------------- #
# GRCh38 centromere midpoints (approx) for arm assignment
_CENTROMERE = {
    "chr1": 123_400_000, "chr2": 93_900_000, "chr3": 90_900_000, "chr19": 26_200_000,
}
_CHR_END = {"chr1": 248_956_422, "chr19": 58_617_616}


def compute_genome_features(ascat_segments: List[Dict], ploidy: Optional[float]) -> Dict[str, Any]:
    feats: Dict[str, Any] = {}
    # ---- WGD: fraction of autosomal genome with major allele CN >= 2 (Dentro) ----
    span_total = span_majge2 = 0.0
    for s in ascat_segments:
        chrom = _norm_chr(s["chrom"])
        if chrom in ("chrX", "chrY"):
            continue
        length = (s.get("end") or 0) - (s.get("start") or 0)
        if length <= 0:
            continue
        span_total += length
        maj = None
        if s.get("minor_cn") is not None and s.get("total_cn") is not None:
            maj = s["total_cn"] - s["minor_cn"]
        if maj is not None and maj >= 2:
            span_majge2 += length
    frac = (span_majge2 / span_total) if span_total else None
    wgd = bool(frac is not None and frac > 0.5)
    if ploidy is not None and ploidy >= 2.7:
        wgd = wgd or True
    feats["wgd"] = {"call": wgd, "frac_major_ge2": round(frac, 3) if frac is not None else None,
                    "ploidy": ploidy}

    # ---- arm-level mean CN for 1p, 1q, 19q ----
    def arm_cn(chrom, arm):
        cen = _CENTROMERE.get(chrom)
        wsum = wlen = 0.0
        for s in ascat_segments:
            if _norm_chr(s["chrom"]) != chrom or s.get("total_cn") is None:
                continue
            st, en = s.get("start") or 0, s.get("end") or 0
            if arm == "p":
                lo, hi = 0, cen
            else:
                lo, hi = cen, _CHR_END.get(chrom, en)
            ov = min(en, hi) - max(st, lo)
            if ov > 0:
                wsum += s["total_cn"] * ov
                wlen += ov
        return (wsum / wlen) if wlen else None

    base = ploidy if (ploidy and ploidy > 0) else 2.0
    cn_1p, cn_1q, cn_19q = arm_cn("chr1", "p"), arm_cn("chr1", "q"), arm_cn("chr19", "q")
    # 1q trisomy / gain: arm CN clearly above the genome baseline (Wilms prognostic)
    feats["gain_1q"] = {"call": bool(cn_1q is not None and cn_1q >= base + 0.6),
                        "arm_cn": round(cn_1q, 2) if cn_1q is not None else None}
    # 1p/19q co-deletion: both arms below baseline (oligodendroglioma marker)
    codel = bool(cn_1p is not None and cn_19q is not None
                 and cn_1p <= base - 0.6 and cn_19q <= base - 0.6)
    feats["codel_1p19q"] = {"call": codel,
                            "cn_1p": round(cn_1p, 2) if cn_1p is not None else None,
                            "cn_19q": round(cn_19q, 2) if cn_19q is not None else None}
    return feats


# --------------------------------------------------------------------------- #
# SV breakpoints + Severus/SAVANA concordance + fusion annotation
# --------------------------------------------------------------------------- #
def load_savana_bedpe_breakpoints(path: str) -> List[Dict]:
    out = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or not f[1].lstrip("-").isdigit():
                continue
            c1, p1, c2, p2 = f[0], int(_num(f[1])), f[3], int(_num(f[4]))
            name = f[6] if len(f) > 6 else ""
            out.append({"chrom1": _norm_chr(c1), "pos1": p1,
                        "chrom2": _norm_chr(c2), "pos2": p2, "name": name})
    return out


def load_severus_vcf_breakpoints(path: str) -> List[Dict]:
    """Severus somatic VCF -> breakpoint pairs. BND mates resolved via ALT/MATE."""
    out, mate_seen = [], {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, alt, info = f[0], int(_num(f[1])), f[2], f[4], f[7]
            ikv = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            svtype = ikv.get("SVTYPE", "")
            if svtype == "BND":
                m = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
                if m:
                    out.append({"chrom1": _norm_chr(chrom), "pos1": pos,
                                "chrom2": _norm_chr(m.group(1)), "pos2": int(m.group(2)),
                                "svtype": "BND"})
            else:
                end = _num(ikv.get("END"))
                out.append({"chrom1": _norm_chr(chrom), "pos1": pos,
                            "chrom2": _norm_chr(chrom), "pos2": int(end) if end else pos,
                            "svtype": svtype})
    return out


def _bp_match(a, b, window=1000):
    def same(c1, p1, c2, p2):
        return c1 == c2 and abs(p1 - p2) <= window
    direct = (same(a["chrom1"], a["pos1"], b["chrom1"], b["pos1"]) and
              same(a["chrom2"], a["pos2"], b["chrom2"], b["pos2"]))
    swap = (same(a["chrom1"], a["pos1"], b["chrom2"], b["pos2"]) and
            same(a["chrom2"], a["pos2"], b["chrom1"], b["pos1"]))
    return direct or swap


def build_sv_concordance(savana_bps, severus_bps, window=1000) -> Dict[str, Any]:
    sev_matched = [False] * len(severus_bps)
    both = sav_only = 0
    for a in savana_bps:
        hit = False
        for j, b in enumerate(severus_bps):
            if _bp_match(a, b, window):
                hit = True; sev_matched[j] = True
        both += 1 if hit else 0
        sav_only += 0 if hit else 1
    sev_only = sev_matched.count(False)
    total = len(savana_bps) + sev_only
    return {"n_savana": len(savana_bps), "n_severus": len(severus_bps),
            "n_concordant": both, "n_savana_only": sav_only, "n_severus_only": sev_only,
            "concordance_pct": round(100 * both / total, 1) if total else None, "window_bp": window}


def annotate_fusions(fusions: List[Dict]) -> List[Dict]:
    """Tag fusions by clinical relevance against the curated actionable set."""
    out = []
    for fus in fusions:
        text = str(fus.get("annotation", "") or fus.get("genes", "")).upper()
        genes = set(re.split(r"[^A-Z0-9]+", text))
        rel = None
        for pair, meta in ACTIONABLE_FUSIONS.items():
            if pair <= genes:
                rel = {"pair": "::".join(sorted(pair)), **meta}
                break
        f2 = dict(fus)
        f2["clinical"] = rel
        out.append(f2)
    return out


# --------------------------------------------------------------------------- #
# PCGR MSI from the XLSX workbook
# --------------------------------------------------------------------------- #
def parse_pcgr_msi(xlsx_path: str) -> Dict[str, Any]:
    try:
        import pandas as pd
    except Exception:
        return {"available": False}
    try:
        sheets = pd.read_excel(xlsx_path, sheet_name=None)
    except Exception:
        return {"available": False}
    out = {"available": False, "source": os.path.basename(xlsx_path)}
    for name, df in sheets.items():
        cols = {str(c).lower(): c for c in df.columns}
        msi_col = next((cols[c] for c in cols if "msi" in c and
                        ("class" in c or "status" in c or "predict" in c)), None)
        prob_col = next((cols[c] for c in cols if "msi" in c and ("prob" in c or "confid" in c)), None)
        if msi_col is not None and len(df):
            val = str(df.iloc[0][msi_col])
            out.update({"available": True, "status": val, "sheet": name,
                        "probability": (_num(df.iloc[0][prob_col]) if prob_col is not None else None),
                        "classifier": "PCGR statistical MSI classifier (TCGA-trained)"})
            break
        # some PCGR workbooks place MSI in a long key/value 'assay' sheet
        if {"setting", "value"} <= set(c.lower() for c in df.columns):
            kc = cols.get("setting"); vc = cols.get("value")
            for _, r in df.iterrows():
                if "msi" in str(r[kc]).lower():
                    out.update({"available": True, "status": str(r[vc]), "sheet": name})
                    break
            if out["available"]:
                break
    return out


def build_sv_circos(translocations: List[Dict], width: int = 360) -> Optional[str]:
    """Static SVG circos of inter-chromosomal translocations (chr1-22,X) as
    chords on a chromosome ring. Width defaults to ~50% of the page. Pure SVG so
    it renders identically in HTML and PDF (no JS)."""
    import math as _m
    CHR_LEN = {  # GRCh38
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468, "chrX": 156040895}
    order = list(CHR_LEN.keys())
    total = sum(CHR_LEN.values())
    gap = 0.012 * 2 * _m.pi          # gap between chromosomes (radians)
    ang_span = 2 * _m.pi - gap * len(order)
    # cumulative angle ranges per chromosome
    ranges, a = {}, -_m.pi / 2
    for c in order:
        sweep = ang_span * CHR_LEN[c] / total
        ranges[c] = (a, a + sweep)
        a += sweep + gap

    def ang(chrom, pos):
        a0, a1 = ranges[chrom]
        return a0 + (a1 - a0) * max(0, min(1, pos / CHR_LEN[chrom]))

    W = width; cx = cy = W / 2
    R = W * 0.42; ringw = W * 0.035
    palette = ["#3b6", "#39c", "#c63", "#96c", "#c39", "#693", "#36c", "#c93"]
    S = [f'<svg viewBox="0 0 {W} {W}" width="{W}" xmlns="http://www.w3.org/2000/svg" '
         f'font-family="-apple-system,Segoe UI,Helvetica,Arial,sans-serif">']
    # chromosome arcs + labels
    for i, c in enumerate(order):
        a0, a1 = ranges[c]
        x0, y0 = cx + R * _m.cos(a0), cy + R * _m.sin(a0)
        x1, y1 = cx + R * _m.cos(a1), cy + R * _m.sin(a1)
        large = 1 if (a1 - a0) > _m.pi else 0
        col = palette[i % len(palette)]
        S.append(f'<path d="M {x0:.1f} {y0:.1f} A {R:.1f} {R:.1f} 0 {large} 1 {x1:.1f} {y1:.1f}" '
                 f'fill="none" stroke="{col}" stroke-width="{ringw:.1f}"/>')
        am = (a0 + a1) / 2; lr = R + ringw
        lx, ly = cx + lr * _m.cos(am), cy + lr * _m.sin(am)
        S.append(f'<text x="{lx:.1f}" y="{ly:.1f}" font-size="{W*0.026:.1f}" fill="#555" '
                 f'text-anchor="middle" dominant-baseline="middle">{c[3:]}</text>')
    # chords for translocations
    rin = R - ringw / 2
    for b in translocations:
        if b["chrom1"] not in ranges or b["chrom2"] not in ranges:
            continue
        a1, a2 = ang(b["chrom1"], b["pos1"]), ang(b["chrom2"], b["pos2"])
        x1, y1 = cx + rin * _m.cos(a1), cy + rin * _m.sin(a1)
        x2, y2 = cx + rin * _m.cos(a2), cy + rin * _m.sin(a2)
        S.append(f'<path d="M {x1:.1f} {y1:.1f} Q {cx:.1f} {cy:.1f} {x2:.1f} {y2:.1f}" '
                 f'fill="none" stroke="#d2453e" stroke-width="1.1" stroke-opacity="0.55"/>')
    n = len([b for b in translocations if b["chrom1"] in ranges])
    S.append(f'<text x="{cx:.1f}" y="{cy:.1f}" font-size="{W*0.05:.1f}" fill="#333" '
             f'text-anchor="middle" dominant-baseline="middle" font-weight="700">{n}</text>')
    S.append(f'<text x="{cx:.1f}" y="{cy + W*0.045:.1f}" font-size="{W*0.026:.1f}" '
             f'fill="#777" text-anchor="middle">translocations</text>')
    S.append("</svg>")
    return "".join(S)


def parse_raita(cfg: Config) -> Dict[str, Any]:
    """Tumour-normal identity QC from RAITA (Andrew's own T/N-aware caller).
    Output is <sampleID>.RAITA.txt, a set of '#'-prefixed comment lines. We only
    surface: germline-het sites compared, whether they're concordant, and the
    PASS/FAIL verdict."""
    out: Dict[str, Any] = {"available": False}
    path = _first(cfg.input_dir,
                  "*.RAITA.txt", "Somalier/*.RAITA.txt", "somalier/*.RAITA.txt",
                  "raita/*.RAITA.txt", "RAITA/*.RAITA.txt")
    if not path:
        return out
    try:
        text = open(path, encoding="utf-8", errors="replace").read()
    except Exception as exc:
        LOG.warning("RAITA read failed (%s): %s", path, exc)
        return out
    lines = [l.lstrip("#").strip() for l in text.splitlines() if l.strip()]

    def find(pat):
        rx = re.compile(pat, re.I)
        for l in lines:
            m = rx.search(l)
            if m:
                return m
        return None

    m = find(r"germline[- ]het sites compared:\s*(\d+)")
    if m:
        out["n_sites"] = int(m.group(1))
    m = find(r"concordant[^:]*:\s*(\d+)\s*/\s*(\d+)\s*=\s*([\d.]+)\s*%")
    if m:
        out["n_concordant"] = int(m.group(1))
        out["n_compared"] = int(m.group(2))
        out["concordance_pct"] = _to_num(m.group(3))
        out.setdefault("n_sites", int(m.group(2)))
    m = find(r"discordant[^:]*:\s*(\d+)")
    if m:
        out["n_discordant"] = int(m.group(1))
    m = find(r"of which\s+het/het:\s*(\d+)\s+LoH:\s*(\d+)")
    if m:
        out["n_hethet"] = int(m.group(1))
        out["n_loh"] = int(m.group(2))
    m = find(r"VERDICT:\s*(PASS|FAIL)\b[:\s]*(.*)$")
    if m:
        out["verdict"] = m.group(1).upper()
        out["verdict_detail"] = m.group(2).strip() or None
    # concordance is "yes" when RAITA saw no foreign (discordant) alleles
    if "n_discordant" in out:
        out["concordant"] = out["n_discordant"] == 0
    elif out.get("concordance_pct") is not None:
        out["concordant"] = out["concordance_pct"] >= 99.0
    if out.get("verdict") is None and out.get("concordant") is not None:
        out["verdict"] = "PASS" if out["concordant"] else "FAIL"
    out["identity_ok"] = out.get("verdict") == "PASS"
    out["available"] = bool(out.get("verdict") or out.get("n_sites"))
    out["source"] = os.path.relpath(path, cfg.input_dir)
    return out


def parse_somalier(cfg: Config) -> Dict[str, Any]:
    """Tumour-normal sample-identity QC from somalier relate output."""
    out: Dict[str, Any] = {"available": False, "pairs": [], "samples": []}
    pairs = _first(cfg.input_dir, "*somalier.pairs.tsv", "somalier/*pairs.tsv")
    samples = _first(cfg.input_dir, "*somalier.samples.tsv", "somalier/*samples.tsv")
    if pairs and pd is not None:
        df = _read_table(pairs)
        if df is not None and len(df):
            df.columns = [str(c).lstrip("#").strip() for c in df.columns]
            c_a, c_b = _pick_col(df, "sample_a"), _pick_col(df, "sample_b")
            c_rel = _pick_col(df, "relatedness")
            c_n = _pick_col(df, "n", "shared_hets", "shared_sites")
            for _, r in df.iterrows():
                out["pairs"].append({"sample_a": r.get(c_a) if c_a else None,
                                     "sample_b": r.get(c_b) if c_b else None,
                                     "relatedness": _to_num(r.get(c_rel)) if c_rel else None,
                                     "n_sites": _to_num(r.get(c_n)) if c_n else None})
            rels = [p["relatedness"] for p in out["pairs"] if p["relatedness"] is not None]
            if rels:
                # T/N from the same individual: relatedness ~1; a swap drops it well below
                out["tn_relatedness"] = max(rels)
                out["identity_ok"] = bool(max(rels) >= 0.9)
            out["available"] = True
            out["pairs_source"] = os.path.basename(pairs)
    if samples and pd is not None:
        df = _read_table(samples)
        if df is not None and len(df):
            df.columns = [str(c).lstrip("#").strip() for c in df.columns]
            c_s = _pick_col(df, "sample_id", "sample")
            c_sex = _pick_col(df, "original_pedigree_sex", "sex", "predicted_sex")
            c_het = _pick_col(df, "n_het")
            c_dp = _pick_col(df, "depth_mean", "gt_depth_mean")
            for _, r in df.iterrows():
                out["samples"].append({"sample": r.get(c_s) if c_s else None,
                                       "sex": r.get(c_sex) if c_sex else None,
                                       "n_het": _to_num(r.get(c_het)) if c_het else None,
                                       "depth_mean": _to_num(r.get(c_dp)) if c_dp else None})
            out["available"] = True
            out["samples_source"] = os.path.basename(samples)
    return out


def parse_hla_loh(cfg: Config) -> Dict[str, Any]:
    """HLA loss-of-heterozygosity from hla_loh_phased_baf.py JSON output."""
    out: Dict[str, Any] = {"available": False, "genes": []}
    js = _first(cfg.input_dir, "hla_loh/*hla_loh*.json", "hla_loh/*.json")
    if not js:
        return out
    try:
        d = json.load(open(js))
        genes = d.get("genes", [])
        out["genes"] = [{"gene": g.get("gene"), "call": g.get("call"),
                         "cn_state": g.get("cn_state"),
                         "obs_minor_fraction": g.get("obs_minor_fraction"),
                         "n_sites_used": g.get("n_sites_used"),
                         "pvalue": g.get("pvalue")} for g in genes]
        out["n_loh"] = sum(1 for g in genes if str(g.get("call")).upper() == "LOH")
        out["n_informative"] = sum(1 for g in genes
                                   if str(g.get("call")).lower() in ("loh", "retained"))
        out["method"] = d.get("method")
        out["available"] = bool(genes)
        out["source"] = os.path.relpath(js, cfg.input_dir)
    except Exception as exc:
        LOG.warning("HLA-LOH parse failed: %s", exc)
    return out


def parse_hrd(cfg: Config) -> Dict[str, Any]:
    """HRD score (scarHRD): HRD-LOH, Telomeric AI, LST, HRD-sum."""
    out: Dict[str, Any] = {"available": False}
    f = _first(cfg.input_dir, "hrd_out/*HRDresults.txt", "*HRDresults.txt")
    if not f:
        return out
    df = _read_table(f)
    if df is not None and len(df):
        r = df.iloc[0]
        def g(*names):
            c = _pick_col(df, *names)
            return _to_num(r.get(c)) if c else None
        out.update({"available": True,
                    "hrd_loh": g("HRD", "HRD_LOH", "HRD-LOH", "HRD.LOH"),
                    "tai": g("Telomeric AI", "Telomeric.AI", "Telomeric_AI", "TAI"),
                    "lst": g("LST"),
                    "hrd_sum": g("HRD-sum", "HRD.sum", "HRD_sum", "HRDsum"),
                    "source": os.path.relpath(f, cfg.input_dir)})
        out["hr_deficient"] = bool(out["hrd_sum"] is not None and out["hrd_sum"] >= 42)
    return out


def parse_clonal(cfg: Config) -> Dict[str, Any]:
    """Clonal architecture from PyClone-VI cluster summary."""
    out: Dict[str, Any] = {"available": False, "clusters": []}
    cl = _first(cfg.input_dir, "clonal/*pyclone*clusters.tsv", "clonal/*clusters.tsv")
    if not cl:
        return out
    df = _read_table(cl)
    if df is not None and len(df):
        c_id = _pick_col(df, "cluster_id")
        c_n = _pick_col(df, "n_mutations", "size")
        c_ccf = next((c for c in df.columns if str(c).lower().endswith("mean_ccf")), None) \
            or _pick_col(df, "mean_ccf", "mean", "cellular_prevalence")
        for _, r in df.iterrows():
            out["clusters"].append({"cluster_id": r.get(c_id) if c_id else None,
                                    "n_mutations": _to_num(r.get(c_n)) if c_n else None,
                                    "mean_ccf": _to_num(r.get(c_ccf)) if c_ccf else None})
        out["clusters"].sort(key=lambda c: -(c["mean_ccf"] or 0))
        out["n_clusters"] = len(out["clusters"])
        # clonal = top CCF cluster; count subclonal clusters (mean CCF < 0.9)
        out["n_subclonal_clusters"] = sum(1 for c in out["clusters"]
                                          if (c["mean_ccf"] or 0) < 0.9)
        out["available"] = True
        out["tool"] = "PyClone-VI"
        out["source"] = os.path.relpath(cl, cfg.input_dir)
    return out


def parse_gene_fusions(cfg: Config) -> Dict[str, Any]:
    """Prioritised gene fusions from a simple_sv_annotation *.genefusions.txt.
    Robust to the file's duplicated / comma-joined records via regex + de-dup."""
    path = _first(cfg.input_dir, "*genefusions.txt", "*.gene_fusions.txt",
                  "*genefusions.tsv", "*gene_fusions.txt")
    if not path:
        return {"available": False, "rows": [], "n": 0}
    seen, rows = set(), []
    other = set()
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                for m in _FUSION_RE.finditer(line):
                    svtype, ann, genes, tx, prio, lvl = m.groups()
                    key = (svtype, ann, genes, tx, prio, lvl)
                    # Andrew: only bidirectional gene fusions are of interest
                    if "BIDIRECTIONAL_GENE_FUSION" not in ann:
                        other.add(key)
                        continue
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append({"svtype": svtype, "annotation": ann.replace("&", " & "),
                                 "genes": genes.replace("&", "-"), "genes_raw": genes,
                                 "transcript": tx, "priority": prio,
                                 "level": int(lvl), "is_fusion": "&" in genes})
    except Exception as exc:
        LOG.warning("Gene-fusion parse failed (%s): %s", path, exc)
        return {"available": False, "rows": [], "n": 0}
    # prioritised first, then by level (1 = highest), then gene name
    rows.sort(key=lambda r: (r["priority"] != "PRIORITISED", r["level"], r["genes"]))
    return {"available": bool(rows), "rows": rows, "n": len(rows),
            "n_prioritised": sum(1 for r in rows if r["priority"] == "PRIORITISED"),
            "n_pairs": sum(1 for r in rows if r["is_fusion"]),
            "n_filtered_out": len(other),
            "source": os.path.relpath(path, cfg.input_dir)}


def parse_pharmcat(cfg: Config) -> Dict[str, Any]:
    """PharmCAT calls-only TSV (<sample>.*.report.tsv). Documented columns:
    Gene / Source Diplotype / Phenotype / Activity Score / Haplotype 1 /
    Haplotype 1 Function / Haplotype 1 Activity Value / Haplotype 2 /
    Haplotype 2 Function / Haplotype 2 Activity Value / Outside Call /
    Match Score / Missing positions / Recommendation Lookup Diplotype /
    Recommendation Lookup Phenotype / Recommendation Lookup Activity Score.
    Matched by header name, so extra/reordered columns are tolerated."""
    path = _first(cfg.input_dir, "*.report.tsv", "pharmcat/*.report.tsv",
                  "*phased.report.tsv", "pharmcat/*.tsv")
    if not path:
        return {"available": False, "rows": [], "n": 0}
    try:
        lines = [l.rstrip("\n") for l in open(path, encoding="utf-8", errors="replace")]
    except Exception as exc:
        LOG.warning("PharmCAT read failed (%s): %s", path, exc)
        return {"available": False, "rows": [], "n": 0}
    # skip any preamble; the header row is the one starting with 'Gene'
    hi = next((i for i, l in enumerate(lines)
               if l.split("\t")[0].strip().lower() == "gene"), None)
    if hi is None:
        LOG.warning("PharmCAT: no 'Gene' header row in %s", path)
        return {"available": False, "rows": [], "n": 0}
    header = [h.strip().lower() for h in lines[hi].split("\t")]

    def idx(*names):
        for n in names:
            if n in header:
                return header.index(n)
        return None

    i_gene = idx("gene")
    i_dip = idx("source diplotype", "diplotype")
    i_phe = idx("phenotype")
    i_as = idx("activity score")
    i_rdip = idx("recommendation lookup diplotype")
    i_rphe = idx("recommendation lookup phenotype")
    i_ras = idx("recommendation lookup activity score")
    i_out = idx("outside call")
    i_miss = idx("missing positions")

    def get(f, i):
        if i is None or i >= len(f):
            return None
        v = f[i].strip()
        return v if v and v not in (".", "NA", "-") else None

    rows = []
    for l in lines[hi + 1:]:
        if not l.strip():
            continue
        f = l.split("\t")
        g = get(f, i_gene)
        if not g:
            continue
        # recommendation-lookup values drive prescribing; fall back to the source call
        phen = get(f, i_rphe) or get(f, i_phe)
        rows.append({
            "gene": g,
            "diplotype": get(f, i_dip),
            "phenotype": phen,
            "source_phenotype": get(f, i_phe),
            "activity_score": get(f, i_ras) or get(f, i_as),
            "lookup_diplotype": get(f, i_rdip),
            "outside_call": get(f, i_out),
            "missing_positions": get(f, i_miss),
            "uncertain": bool(get(f, i_dip) and " OR " in (get(f, i_dip) or "")),
            "actionable": bool(phen and not re.search(
                r"normal metaboliser|normal metabolizer|normal function|indeterminate",
                phen, re.I)),
        })
    rows.sort(key=lambda r: (not r["actionable"], r["gene"]))
    return {"available": bool(rows), "rows": rows, "n": len(rows),
            "n_actionable": sum(1 for r in rows if r["actionable"]),
            "source": os.path.relpath(path, cfg.input_dir)}


def parse_severus_loh(cfg: Config) -> Dict[str, Any]:
    """Severus genome-wide LOH regions (severus_out/severus_LOH.bed)."""
    out: Dict[str, Any] = {"available": False, "regions": []}
    bed = _first(cfg.input_dir, "severus_out/severus_LOH.bed", "severus_out/*LOH*.bed",
                 "*severus_LOH.bed")
    if not bed:
        return out
    try:
        total = 0
        with _open(bed) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3 or not f[1].lstrip("-").isdigit():
                    continue
                s, e = int(_num(f[1])), int(_num(f[2]))
                out["regions"].append({"chrom": f[0], "start": s, "end": e})
                total += max(0, e - s)
        out["available"] = bool(out["regions"])
        out["n"] = len(out["regions"])
        out["total_mb"] = round(total / 1e6, 1)
        out["source"] = os.path.relpath(bed, cfg.input_dir)
    except Exception as exc:
        LOG.warning("Severus LOH parse failed: %s", exc)
    return out


# simple_sv_annotation record: SVTYPE|ANNOTATION|GENE(S)|TRANSCRIPT|PRIORITY_CAT|LEVEL
# (records may be comma-joined and the whole line is often duplicated)
_FUSION_RE = re.compile(
    r"([A-Z]+)\|([A-Z_&]+)\|([^|,]+)\|([^|,]+)\|(PRIORITISED|NOT_PRIORITISED)\|([0-9])")


# --- Signature aetiologies -------------------------------------------------
# SBS: COSMIC v3.x proposed aetiologies (abbreviated).
_SBS_AETIOLOGY = {
    "SBS1": "Clock-like: spontaneous deamination of 5-methylcytosine",
    "SBS2": "APOBEC cytidine deaminase activity",
    "SBS3": "Homologous-recombination deficiency (BRCA1/2)",
    "SBS4": "Tobacco smoking",
    "SBS5": "Clock-like, unknown aetiology",
    "SBS6": "Defective DNA mismatch repair (MSI)",
    "SBS7a": "UV light exposure", "SBS7b": "UV light exposure",
    "SBS7c": "UV light exposure", "SBS7d": "UV light exposure",
    "SBS8": "Unknown; associated with HRD / NER deficiency",
    "SBS9": "Polymerase eta somatic hypermutation",
    "SBS10a": "Polymerase epsilon (POLE) proofreading deficiency",
    "SBS10b": "Polymerase epsilon (POLE) proofreading deficiency",
    "SBS10c": "Polymerase delta (POLD1) proofreading deficiency",
    "SBS10d": "Polymerase delta (POLD1) proofreading deficiency",
    "SBS11": "Temozolomide treatment",
    "SBS13": "APOBEC cytidine deaminase activity",
    "SBS14": "Defective mismatch repair + polymerase proofreading",
    "SBS15": "Defective DNA mismatch repair (MSI)",
    "SBS16": "Unknown; associated with alcohol / liver cancer",
    "SBS17a": "Unknown; associated with 5-FU chemotherapy",
    "SBS17b": "Unknown; associated with 5-FU chemotherapy",
    "SBS18": "Reactive oxygen species / oxidative damage",
    "SBS20": "Defective mismatch repair + POLD1 mutation",
    "SBS21": "Defective DNA mismatch repair (MSI)",
    "SBS22": "Aristolochic acid exposure",
    "SBS23": "Unknown", "SBS24": "Aflatoxin exposure",
    "SBS25": "Chemotherapy (Hodgkin lymphoma)",
    "SBS26": "Defective DNA mismatch repair (MSI)",
    "SBS28": "Unknown; co-occurs with POLE",
    "SBS29": "Tobacco chewing",
    "SBS30": "Base-excision repair deficiency (NTHL1)",
    "SBS31": "Platinum chemotherapy",
    "SBS32": "Azathioprine treatment",
    "SBS33": "Unknown", "SBS35": "Platinum chemotherapy",
    "SBS36": "Base-excision repair deficiency (MUTYH)",
    "SBS37": "Unknown", "SBS38": "Indirect UV damage",
    "SBS39": "Unknown", "SBS40": "Unknown; correlates with age",
    "SBS40a": "Unknown; correlates with age",
    "SBS40b": "Unknown", "SBS40c": "Unknown",
    "SBS41": "Unknown", "SBS42": "Haloalkane exposure",
    "SBS44": "Defective DNA mismatch repair (MSI)",
    "SBS84": "Activation-induced cytidine deaminase (AID)",
    "SBS85": "Indirect AID activity",
    "SBS86": "Unknown chemotherapy exposure",
    "SBS87": "Thiopurine chemotherapy",
    "SBS88": "Colibactin (pks+ E. coli) exposure",
    "SBS90": "Duocarmycin exposure",
    "SBS92": "Tobacco smoking (bladder)",
    "SBS93": "Unknown", "SBS94": "Unknown",
    "SBS27": "Possible sequencing artefact",
    "SBS43": "Possible sequencing artefact", "SBS45": "Possible sequencing artefact",
    "SBS46": "Possible sequencing artefact", "SBS47": "Possible sequencing artefact",
    "SBS48": "Possible sequencing artefact", "SBS49": "Possible sequencing artefact",
    "SBS50": "Possible sequencing artefact", "SBS51": "Possible sequencing artefact",
    "SBS52": "Possible sequencing artefact", "SBS53": "Possible sequencing artefact",
    "SBS54": "Possible sequencing artefact / germline contamination",
    "SBS55": "Possible sequencing artefact", "SBS56": "Possible sequencing artefact",
    "SBS57": "Possible sequencing artefact", "SBS58": "Possible sequencing artefact",
    "SBS59": "Possible sequencing artefact", "SBS60": "Possible sequencing artefact",
    "SBS95": "Possible sequencing artefact",
}

# Copy-number signatures — proposed aetiologies, Steele et al. Nature 2022 (COSMIC CN1-CN21)
_CN_AETIOLOGY = {
    "CN1": "Diploid, no whole-genome doubling",
    "CN2": "Tetraploid (1x whole-genome doubling)",
    "CN3": "Octoploid (2x whole-genome doubling)",
    "CN4": "Chromothripsis", "CN5": "Chromothripsis",
    "CN6": "Chromothripsis before whole-genome doubling",
    "CN7": "Chromothripsis-associated amplification",
    "CN8": "Chromothripsis-associated amplification",
    "CN9": "Chromosomal instability on a diploid background",
    "CN10": "Loss of heterozygosity with 1x whole-genome doubling",
    "CN11": "Loss of heterozygosity with 2x whole-genome doubling",
    "CN12": "Chromosomal instability with 1x whole-genome doubling",
    "CN13": "Chromosomal-scale LOH, no whole-genome doubling",
    "CN14": "Chromosomal-scale LOH with 1x whole-genome doubling",
    "CN15": "Chromosomal-scale LOH with 2x whole-genome doubling",
    "CN16": "Chromosomal LOH + heterozygous mix, 2x whole-genome doubling",
    "CN17": "Homologous-recombination deficiency (HRD)",
    "CN18": "Unknown aetiology", "CN19": "Unknown aetiology",
    "CN20": "Unknown aetiology",
    "CN21": "Unknown aetiology (associated with tandem-duplicator phenotype)",
    "CN22": "Unknown aetiology", "CN23": "Unknown aetiology",
    "CN24": "Unknown aetiology",
}

# Rearrangement / SV reference signatures (RefSig R*, Degasperi et al. 2020)
_SV_AETIOLOGY = {
    "RefSigR1": "Large (>100 kb) tandem duplications — tandem-duplicator phenotype",
    "RefSigR2": "Unclustered translocations, large deletions and inversions",
    "RefSigR3": "Short tandem duplications — associated with BRCA1 loss (HRD)",
    "RefSigR4": "Clustered interchromosomal translocations",
    "RefSigR5": "Unclustered deletions up to 100 kb — associated with BRCA1/2 loss (HRD)",
    "RefSigR6a": "Clustered rearrangements / complex focal amplification",
    "RefSigR6b": "Clustered rearrangements / complex focal amplification",
    "RefSigR7": "Unknown aetiology",
    "RefSigR8": "Short (1-10 kb) inversions",
    "RefSigR9": "Unclustered deletions up to 1-10 kb — associated with BRCA2 loss",
    "RefSigR10": "Unknown aetiology",
    "RefSigR11": "Unknown aetiology",
    "RefSigR14": "Long / dispersed tandem duplications",
    "RefSigR15": "Large deletions and inversions",
    "RefSigR16": "Unknown aetiology", "RefSigR17": "Unknown aetiology",
    "RefSigR18": "Interchromosomal translocations",
    "RefSigR19": "Unknown aetiology", "RefSigR20": "Unknown aetiology",
}


def _signature_aetiology(sig: str, kind: str) -> str:
    """What the signature is currently thought to indicate. Unknown IDs are
    labelled as such rather than guessed at."""
    s = str(sig or "").strip()
    if kind == "sbs":
        return _SBS_AETIOLOGY.get(s, "Not in the annotation table — see COSMIC")
    if kind == "cnv":
        return _CN_AETIOLOGY.get(s, "Not in the annotation table — see COSMIC")
    if kind == "sv":
        key = s.replace(".", "").replace("_", "").replace("-", "")
        if key in _SV_AETIOLOGY:
            return _SV_AETIOLOGY[key]
        m = re.match(r"(?:RefSig)?\.?R?(\d+[ab]?)$", s, re.I)
        if m:
            return _SV_AETIOLOGY.get("RefSigR" + m.group(1),
                                     "Not in the annotation table — see Signal / RefSig")
        return "Not in the annotation table — see Signal / RefSig"
    return ""


def parse_sigprofiler(cfg: Config) -> Dict[str, Any]:
    """SigProfilerAssignment refit results. Each of sbs96_signatures.tsv,
    cnv_signatures.tsv and sv_signatures.tsv is a small table:
        signature   count   proportion
    Also picks up SigProfiler's own PDF plots for the PDF report."""
    out: Dict[str, Any] = {}
    families = [
        ("sbs", "SBS96 (single base substitution)",
         ("sigprofiler/sbs96_signatures.tsv", "sigprofiler/*sbs96*signatures.tsv"),
         ("sigprofiler/sbs96/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf",
          "sigprofiler/sbs96/Assignment_Solution/Signatures/SBS_96_plots_Assignment_Solution.pdf")),
        ("cnv", "CNV48 (copy number)",
         ("sigprofiler/cnv_signatures.tsv", "sigprofiler/*cnv*signatures.tsv"),
         ("sigprofiler/cnv/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf",
          "sigprofiler/cnv/Assignment_Solution/Signatures/CNV_48_plots_Assignment_Solution.pdf")),
        ("sv", "SV32 (structural variant / rearrangement)",
         ("sigprofiler/sv_signatures.tsv", "sigprofiler/*sv*signatures.tsv"),
         ("sigprofiler/sv_refit/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf",
          "sigprofiler/sv_refit/Assignment_Solution/Signatures/SV_32_plots_Assignment_Solution.pdf")),
    ]
    for kind, label, tsv_pats, pdf_pats in families:
        path = _first(cfg.input_dir, *tsv_pats)
        rows: List[Dict[str, Any]] = []
        if path:
            try:
                lines = [l.rstrip("\n") for l in open(path, encoding="utf-8",
                                                      errors="replace") if l.strip()]
                hdr = [h.strip().lower() for h in lines[0].split("\t")]
                i_sig = hdr.index("signature") if "signature" in hdr else 0
                i_cnt = hdr.index("count") if "count" in hdr else None
                i_prop = hdr.index("proportion") if "proportion" in hdr else None
                for l in lines[1:]:
                    f = l.split("\t")
                    if len(f) <= i_sig or not f[i_sig].strip():
                        continue
                    sig = f[i_sig].strip()
                    rows.append({
                        "signature": sig,
                        "count": _to_num(f[i_cnt]) if (i_cnt is not None and i_cnt < len(f)) else None,
                        "proportion": _to_num(f[i_prop]) if (i_prop is not None and i_prop < len(f)) else None,
                        "aetiology": _signature_aetiology(sig, kind),
                    })
                rows.sort(key=lambda r: -(r["proportion"] or 0))
            except Exception as exc:
                LOG.warning("SigProfiler %s parse failed (%s): %s", kind, path, exc)
        else:
            LOG.info("SigProfiler %s: no signatures TSV found", kind)
        out[kind] = {"available": bool(rows), "label": label, "rows": rows,
                     "n": len(rows),
                     "source": os.path.relpath(path, cfg.input_dir) if path else None,
                     "plot_pdf": _first(cfg.input_dir, *pdf_pats)}
    out["available"] = any(out[k].get("available") for k in ("sbs", "cnv", "sv"))
    return out


def assemble(cfg: Config) -> Dict[str, Any]:
    LOG.info("Parsing inputs from %s", cfg.input_dir)
    report: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_by": REPORT_GENERATOR,
        "generated_at": _dt.datetime.now().isoformat(timespec="seconds"),
        "case": {
            "patient_id": cfg.patient_id,
            "display_id": cfg.display_id or cfg.patient_id,
            "tumour_sample": cfg.tumour_id,
            "normal_sample": cfg.normal_id,
            "assay": cfg.assay,
            "genome_build_somatic": cfg.genome_build,
            "input_dir": os.path.abspath(cfg.input_dir),
        },
    }
    report["coverage"] = parse_coverage(cfg)
    report["purity_ploidy"] = parse_purity_ploidy(cfg)
    report["pcgr"] = parse_pcgr(cfg)
    report["cpsr"] = parse_cpsr(cfg)
    report["ascat"] = parse_ascat_segments(cfg)
    report["msi"] = parse_msi(cfg)
    report["sv"] = parse_structural_variants(cfg)
    report["gene_fusions"] = parse_gene_fusions(cfg)
    report["pharmcat"] = parse_pharmcat(cfg)
    report["sigprofiler"] = parse_sigprofiler(cfg)
    report["nanodx"] = parse_nanodx(cfg)

    # ---- multi-caller CN concordance (ASCAT vs SAVANA vs Wakhan) ----
    ascat_segs = [{"chrom": s["chrom"], "start": s["start"], "end": s["end"],
                   "total_cn": s["total_cn"], "minor_cn": s.get("n_minor")}
                  for s in report["ascat"].get("segments", [])
                  if s.get("start") is not None and s.get("end") is not None]
    cn_callers = {"ASCAT": ascat_segs}
    sav_cn = _first(cfg.input_dir, "savana/*phased_segmented_absolute_copy_number.tsv")
    if sav_cn:
        try:
            cn_callers["SAVANA"] = load_savana_cn(sav_cn)
        except Exception as exc:
            LOG.warning("SAVANA CN load failed: %s", exc)
    try:
        wk = load_wakhan_cn(cfg.input_dir)
        if wk:
            cn_callers["Wakhan"] = wk
    except Exception as exc:
        LOG.warning("Wakhan CN load failed: %s", exc)
    driver_genes = [c.get("gene") for c in report["pcgr"].get("cna_genes", {}).get("items", [])
                    if c.get("gene")]
    driver_genes += [v.get("gene") for v in report["pcgr"].get("snv_indel", {}).get("variants", [])
                     if _tier_class(v.get("tier")) in ("t1", "t2") and v.get("gene")]
    report["cn_concordance"] = build_cn_concordance(driver_genes, cn_callers)

    # ---- genome features: WGD, 1q gain, 1p/19q co-deletion ----
    ploidy = report["purity_ploidy"].get("consensus", {}).get("ploidy_median")
    report["genome_features"] = compute_genome_features(ascat_segs, ploidy)

    report["methylation"] = parse_methylation(cfg)
    report["immune"] = parse_hla_neoantigen(cfg)
    report["hla_loh"] = parse_hla_loh(cfg)
    report["somalier"] = parse_somalier(cfg)
    report["raita"] = parse_raita(cfg)
    report["hrd"] = parse_hrd(cfg)
    report["clonal"] = parse_clonal(cfg)
    report["severus_loh"] = parse_severus_loh(cfg)
    report["wakhan"] = parse_wakhan_plots(cfg)
    report["actionability"] = build_actionability(report)
    report["qc_flags"] = build_qc_flags(report, cfg)
    report["provenance"] = parse_provenance(cfg)

    # Figures (kept out of JSON; attached to a separate object)
    figures = {
        "cn_genome": make_cn_figure(report["ascat"]),
        "vaf": make_vaf_figure(report["pcgr"].get("snv_indel", {}).get("variants", [])),
        "signatures": make_signature_figure(report["pcgr"].get("signatures", {}).get("items", [])),
        "signatures_pie": make_signature_pie_figure(report["pcgr"].get("signatures", {}).get("items", [])),
        "sigprof_figs": {
            k: (make_signature_pie_figure(
                    [{"signature": r["signature"], "weight": r.get("proportion")}
                     for r in report["sigprofiler"][k]["rows"]])
                if report.get("sigprofiler", {}).get(k, {}).get("available") else None)
            for k in ("sbs", "cnv", "sv")
        },
        "wakhan_snaps": _wakhan_pdf_figures(report.get("wakhan", {})),
        "sv_circos": build_sv_circos(report["sv"].get("translocations", [])),
        "tmb_context": make_tmb_context_figure(report["pcgr"].get("tmb", {}).get("tmb_per_mb")),
        "severus_loh": make_loh_figure(report.get("severus_loh", {})),
    }
    figures.update({f"ascat_{k}": v for k, v in report["ascat"].get("plots", {}).items()})
    report["_figures"] = figures
    # count what actually rendered — a silent zero means matplotlib is broken
    _total = _rendered = 0
    for k, v in figures.items():
        if k == "wakhan_snaps":
            for _t, img in (v or []):
                _total += 1
                _rendered += 1 if img else 0
        else:
            _total += 1
            _rendered += 1 if v else 0
    report["_figure_stats"] = {"rendered": _rendered, "total": _total}
    if _rendered == 0 and _total:
        LOG.warning("ZERO figures rendered out of %d — matplotlib is almost certainly "
                    "unusable in this environment. Run with --check-env.", _total)
    return report


# ============================================================================
# HTML / PDF rendering
# ============================================================================
CSS = r"""
:root{
  --ink:#16202c; --muted:#5a6776; --hair:#dde3ea; --canvas:#f6f8fa;
  --surface:#ffffff; --primary:#0e6e76; --primary-dark:#0a565c;
  --indigo:#243b6b; --t1:#b3261e; --t2:#c77700; --t3:#1f6feb; --t4:#6b7785;
  --gain:#b3261e; --loss:#1f6feb; --loh:#7a5195;
}
*{box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  color:var(--ink);margin:0;background:var(--canvas);font-size:16.5px;line-height:1.55}
.wrap{max-width:min(1760px,93vw);margin:0 auto;padding:0 30px 72px}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92em}
a{color:var(--primary-dark)}
.masthead{background:var(--indigo);color:#fff;padding:26px 22px;margin-bottom:0}
.masthead .wrap{padding-bottom:0}
.eyebrow{letter-spacing:.16em;text-transform:uppercase;font-size:11px;font-weight:700;
  color:#aebbe0}
.masthead h1{margin:.2em 0 .1em;font-size:30px;font-weight:700;letter-spacing:-.01em}
.masthead .sub{color:#cdd7ee;font-size:13px}
.idgrid{display:flex;flex-wrap:wrap;gap:18px 34px;margin-top:16px;font-size:13px}
.idgrid div span{display:block;color:#9fb0da;font-size:11px;text-transform:uppercase;letter-spacing:.08em}
.idgrid div b{font-weight:600;color:#fff}

/* signature strip: the molecular headline */
.strip{display:flex;flex-wrap:wrap;gap:10px;background:var(--surface);
  border:1px solid var(--hair);border-radius:10px;padding:14px;margin:-26px 0 22px;
  box-shadow:0 6px 24px rgba(20,32,44,.07);position:relative;z-index:2}
.chip{border:1px solid var(--hair);border-radius:8px;padding:9px 13px;min-width:120px;background:#fff}
.chip .k{font-size:10.5px;text-transform:uppercase;letter-spacing:.07em;color:var(--muted)}
.chip .v{font-size:18px;font-weight:700;margin-top:2px}
.chip .v small{font-size:11px;font-weight:500;color:var(--muted)}
.chip.hi{border-color:var(--t1)} .chip.hi .v{color:var(--t1)}

section{background:var(--surface);border:1px solid var(--hair);border-radius:10px;
  padding:20px 22px;margin:16px 0}
section>h2{font-size:19px;margin:0 0 4px;padding-left:12px;border-left:4px solid var(--primary);
  display:flex;align-items:baseline;gap:10px}
section>h2 .num{color:var(--muted);font-weight:600;font-size:13px}
section>.lead{color:var(--muted);margin:.1em 0 14px;padding-left:14px;font-size:15px}
h3{font-size:15.5px;margin:20px 0 7px}

table{border-collapse:collapse;width:100%;font-size:14.6px;margin:8px 0}
th,td{text-align:left;padding:7px 9px;border-bottom:1px solid var(--hair);vertical-align:top}
th{background:#f1f4f8;font-weight:600;color:#33414f;font-size:12.8px;text-transform:uppercase;
  letter-spacing:.03em}
tr:hover td{background:#fbfcfe}
.tier{display:inline-block;border-radius:4px;padding:1px 7px;font-size:11px;font-weight:700;color:#fff}
.tier.t1{background:var(--t1)} .tier.t2{background:var(--t2)} .tier.t3{background:var(--t3)}
.tier.t4{background:var(--t4)} .tier.tx{background:#9aa6b2}
.badge{display:inline-block;border-radius:4px;padding:1px 7px;font-size:11px;font-weight:600;
  background:#eef2f7;color:#33414f;border:1px solid var(--hair)}
.gain{color:var(--gain);font-weight:600}.loss{color:var(--loss);font-weight:600}
.loh{color:var(--loh);font-weight:600}
.flag{border-radius:8px;padding:10px 13px;margin:8px 0;font-size:12.8px;border:1px solid}
.flag.warn{background:#fff7e9;border-color:#f0d9a8;color:#7a5a12}
.flag.error{background:#fdecec;border-color:#f2bcbc;color:#8a1f1f}
.flag.ok{background:#eaf6ef;border-color:#bfe3cd;color:#1f6b3d}
.muted{color:var(--muted)} .small{font-size:13.8px}
.na{color:#9aa6b2;font-style:italic}
img.fig{max-width:100%;border:1px solid var(--hair);border-radius:8px;margin:8px 0;background:#fff}
iframe.plotly-embed{width:100%;height:640px;border:1px solid var(--hair);border-radius:8px;margin:8px 0 4px;background:#fff}
.figrow{display:flex;flex-wrap:wrap;gap:14px;align-items:flex-start;margin:8px 0}
img.fig.wide{width:100%;max-width:none}
.msi-high{color:#b3261e;font-weight:700}
.figrow img.fig{flex:1 1 380px;max-width:640px;margin:0}
.two{display:grid;grid-template-columns:1fr 1fr;gap:18px}
.kvlist{display:grid;grid-template-columns:auto 1fr;gap:4px 16px;font-size:13px}
.kvlist dt{color:var(--muted)} .kvlist dd{margin:0;font-weight:600}
.toc{columns:2;font-size:13px;gap:30px}
.toc a{display:block;padding:2px 0;text-decoration:none;color:var(--ink)}
.toc a:hover{color:var(--primary-dark)}
.footer{color:var(--muted);font-size:11.5px;text-align:center;margin-top:26px;line-height:1.7}
details summary{cursor:pointer;font-weight:600;color:var(--primary-dark);font-size:12.5px}
.rec{border-left:3px solid var(--primary);padding:6px 0 6px 12px;margin:8px 0}
.rec b{display:block}.rec .tool{color:var(--primary-dark);font-family:ui-monospace,monospace;font-size:12px}
.tabs{display:flex;gap:6px;margin:12px 0 0}
.tab-btn{border:1px solid var(--hair);background:#f6f8fb;border-radius:7px 7px 0 0;padding:6px 13px;cursor:pointer;font:inherit;font-size:13px;color:var(--muted)}
.tab-btn.active{background:#fff;color:var(--ink);border-bottom-color:#fff;font-weight:600}
.tab-panel{display:none;border-top:1px solid var(--hair);padding-top:10px}
.tab-panel.active{display:block}
.appendix-title{display:none}
.pager{display:flex;gap:8px;align-items:center;margin:8px 0;font-size:12.5px;color:var(--muted)}
.pager button{border:1px solid var(--hair);background:#fff;border-radius:6px;padding:3px 11px;cursor:pointer;font:inherit}
.pager button[disabled]{opacity:.4;cursor:default}
@page{size:A4;margin:11mm 8mm}
@media print{
  html,body{background:#fff;font-size:10px}
  .wrap{max-width:none;width:auto;padding:0}
  section{break-inside:avoid;box-shadow:none;border-color:#ccc;margin:0 0 6mm}
  .masthead{background:var(--indigo)!important;-webkit-print-color-adjust:exact;print-color-adjust:exact}
  .strip{box-shadow:none}
  tr:hover td{background:transparent}
  .pagebreak{break-before:page}
  .tabs{display:none}
  .tab-panel{display:block!important}
  #tab-vus{break-before:page}
  .appendix-title{display:block}
  .pager,iframe.plotly-embed{display:none}
  .paged tbody tr{display:table-row!important}
  /* keep wide tables on the page */
  img.fig.wide{width:100%;max-width:none;break-inside:avoid}
  /* keep section 6 (mutational signatures) on a single page */
  section#sig{break-inside:avoid;break-before:page;padding:8px 10px}
  section#sig img.fig{max-height:118px;width:auto;margin:1px 0}
  section#sig .figrow{gap:6px;margin:2px 0}
  section#sig table{font-size:7px;margin:1px 0}
  section#sig th{font-size:6.6px;padding:1px 3px}
  section#sig td{padding:1px 3px}
  section#sig h2{margin:0 0 4px}
  section#sig h3{margin:4px 0 1px}
  section#sig h4{margin:3px 0 0}
  section#sig p{margin:1px 0;font-size:7.2px}
  table{table-layout:fixed;width:100%;font-size:8.6px;margin:4px 0}
  th,td{overflow-wrap:anywhere;word-break:break-word;padding:2px 4px}
  th{font-size:8.2px}
  .mono,code{font-size:.9em}
  section>h2{font-size:13px}
  h3{font-size:11px}
  .chip .v{font-size:14px}
}
@media screen and (max-width:720px){.two{grid-template-columns:1fr}.toc{columns:1}}
"""

SORT_JS = r"""
document.querySelectorAll('table.sortable th').forEach(function(th,idx){
  th.style.cursor='pointer';
  th.addEventListener('click',function(){
    var tb=th.closest('table'),rows=Array.from(tb.tBodies[0].rows);
    var asc=!(th.dataset.asc==='true');th.dataset.asc=asc;
    rows.sort(function(a,b){
      var x=a.cells[idx].innerText.trim(),y=b.cells[idx].innerText.trim();
      var nx=parseFloat(x),ny=parseFloat(y);
      if(!isNaN(nx)&&!isNaN(ny))return asc?nx-ny:ny-nx;
      return asc?x.localeCompare(y):y.localeCompare(x);
    });
    rows.forEach(function(r){tb.tBodies[0].appendChild(r);});
  });
});
/* tabbed tables (Tier 1/2 vs VUS) */
document.querySelectorAll('.tab-btn').forEach(function(b){
  b.addEventListener('click',function(){
    var id=b.getAttribute('data-tab');
    b.parentNode.querySelectorAll('.tab-btn').forEach(function(x){x.classList.remove('active');});
    b.classList.add('active');
    var host=b.parentNode.parentNode;
    host.querySelectorAll('.tab-panel').forEach(function(p){p.classList.toggle('active',p.id===id);});
  });
});
/* client-side pagination for tables marked .paged (data-page-size rows/page) */
document.querySelectorAll('table.paged').forEach(function(tb){
  var size=parseInt(tb.getAttribute('data-page-size')||'10',10);
  var rows=Array.from(tb.tBodies[0].rows), pages=Math.ceil(rows.length/size), pg=0;
  if(pages<=1)return;
  var bar=document.createElement('div');bar.className='pager';
  var prev=document.createElement('button');prev.textContent='‹ Prev';
  var next=document.createElement('button');next.textContent='Next ›';
  var lab=document.createElement('span');
  bar.appendChild(prev);bar.appendChild(lab);bar.appendChild(next);
  tb.parentNode.insertBefore(bar,tb.nextSibling);
  function draw(){
    rows.forEach(function(r,i){r.style.display=(i>=pg*size&&i<(pg+1)*size)?'':'none';});
    lab.textContent='Showing '+(pg*size+1)+'–'+Math.min((pg+1)*size,rows.length)+' of '+rows.length;
    prev.disabled=pg===0;next.disabled=pg>=pages-1;
  }
  prev.onclick=function(){if(pg>0){pg--;draw();}};
  next.onclick=function(){if(pg<pages-1){pg++;draw();}};
  draw();
});
"""


def _tier_class(tier: Optional[str]) -> str:
    if not tier:
        return "tx"
    t = str(tier).upper()
    # explicit arabic tier digit (PCGR: TIER_1 .. TIER_4)
    m = re.search(r"\b(?:TIER[_ ]?)?([1-4])\b", t)
    if m:
        return "t" + m.group(1)
    # roman numerals only if a token is *exactly* roman (avoids the 'I' in 'TIER')
    roman = {"IV": "t4", "III": "t3", "II": "t2", "I": "t1"}
    for tok in re.split(r"[^A-Z]+", t):
        if tok in roman:
            return roman[tok]
    return "tx"


def _tier_rank_str(tier: Optional[str]) -> int:
    """Numeric tier rank (1-4) from a tier label; 9 for none/'.'/unranked."""
    cls = _tier_class(tier)
    return int(cls[1]) if cls[1].isdigit() else 9


def _fmt(x: Any, nd: int = 2, pct: bool = False, na: str = "—") -> str:
    if x is None or x == "":
        return f'<span class="na">{na}</span>'
    if isinstance(x, float):
        return f"{x*100:.1f}%" if pct else f"{x:.{nd}f}"
    return str(x)


def render_html(report: Dict[str, Any], cfg: Config, for_pdf: bool = False) -> str:
    fig = report.get("_figures", {})
    case = report["case"]
    pp = report["purity_ploidy"]["consensus"]
    tmb = report["pcgr"].get("tmb", {})
    msi = report.get("msi", {})
    act = report.get("actionability", {})
    sigs = report["pcgr"].get("signatures", {}).get("items", [])

    def esc(s):
        if s is None:
            return '<span class="na">—</span>'
        return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))

    H: List[str] = []
    H.append(f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>T/N report — {esc(case['patient_id'])} {esc(case['tumour_sample'])}</title>
<style>{CSS}</style></head><body>""")

    # masthead
    H.append(f"""<div class="masthead"><div class="wrap">
      <div class="eyebrow">University of <b>B</b>irmingh<b>A</b>m <b>L</b>ong read <b>T</b>umour-normal p<b>I</b>peline (<b>BALTI</b>)</div>
      <h1>Molecular tumour profile — {esc(case.get('display_id') or case['patient_id'])}</h1>
      <div class="sub">{esc(case['assay'])}</div>
      <div class="idgrid">
        <div><span>Tumour</span><b>{esc(case['tumour_sample'])}</b></div>
        <div><span>Normal</span><b>{esc(case['normal_sample'])}</b></div>
        <div><span>Somatic build</span><b>{esc(case['genome_build_somatic'])}</b></div>
        <div><span>Generated</span><b>{esc(report['generated_at'][:16].replace('T',' '))}</b></div>
        <div><span>Schema</span><b>v{esc(report['schema_version'])}</b></div>
      </div>
    </div></div>""")

    H.append('<div class="wrap">')

    # signature strip
    def chip(k, v, hi=False, sub=""):
        cls = "chip hi" if hi else "chip"
        subhtml = f" <small>{sub}</small>" if sub else ""
        return f'<div class="{cls}"><div class="k">{k}</div><div class="v">{v}{subhtml}</div></div>'

    drivers = [v for v in report["pcgr"].get("snv_indel", {}).get("variants", [])
               if _tier_class(v.get("tier")) in ("t1", "t2")]
    driver_txt = ", ".join(sorted({d.get("gene") for d in drivers if d.get("gene")})[:4]) or "—"
    H.append('<div class="strip">')
    H.append(chip("Purity (consensus)",
                  _fmt(pp.get("purity_median"), pct=True) if pp.get("purity_median") is not None else "—"))
    H.append(chip("Ploidy", _fmt(pp.get("ploidy_median")) if pp.get("ploidy_median") is not None else "—"))
    H.append(chip("TMB", _fmt(tmb.get("tmb_per_mb"), 1) if tmb.get("tmb_per_mb") is not None else "—",
                  hi=(tmb.get("tmb_per_mb") or 0) >= 10, sub="/Mb"))
    H.append(chip("MSI", msi.get("status") or "—", hi=(msi.get("status") == "MSI"),
                  sub=(f"p={msi.get('probability'):.2f}" if msi.get("probability") is not None else "")))
    H.append(chip("Somatic SNV/indel",
                  report["pcgr"].get("snv_indel", {}).get("n_total", "—")))
    H.append(chip("Somatic SV", report["sv"].get("n_somatic_sv", "—")))
    # genome-level features
    gf = report.get("genome_features", {})
    if gf:
        wgd = gf.get("wgd", {})
        H.append(chip("WGD", "Yes" if wgd.get("call") else "No", hi=bool(wgd.get("call")),
                      sub=(f"{wgd.get('frac_major_ge2'):.0%} maj≥2"
                           if wgd.get("frac_major_ge2") is not None else "")))
        if gf.get("gain_1q", {}).get("call"):
            H.append(chip("1q gain", "Yes", hi=True,
                          sub=(f"CN {gf['gain_1q'].get('arm_cn')}"
                               if gf["gain_1q"].get("arm_cn") is not None else "")))
        if gf.get("codel_1p19q", {}).get("call"):
            H.append(chip("1p/19q", "Co-deleted", hi=True))
    hrd = report.get("hrd", {})
    if hrd.get("available") and hrd.get("hrd_sum") is not None:
        hrdef = hrd.get("hr_deficient")
        H.append(chip("HRD", "Yes" if hrdef else "No", hi=bool(hrdef),
                      sub=f"sum {hrd.get('hrd_sum'):.0f}"))
    H.append(chip("Actionable", act.get("n", 0), hi=(act.get("n", 0) > 0)))
    H.append(chip("Top drivers", esc(driver_txt)))
    H.append('</div>')

    # TOC
    secs = [("1", "Clinical actionability", "act"),
            ("2", "Quality & provenance", "qc"),
            ("3", "Tumour purity & ploidy", "pp"),
            ("4", "Somatic SNV / indel", "snv"),
            ("5", "Tumour mutational burden & MSI", "tmbmsi"),
            ("6", "Mutational signatures", "sig"),
            ("7", "Copy-number alterations", "cna"),
            ("8", "Structural variants & fusions", "sv"),
            ("9", "Germline predisposition (CPSR)", "germ"),
            ("10", "Methylation & classification", "meth"),
            ("11", "HLA & neoantigens", "immune"),
            ("12", "Pharmacogenomics", "pgx")]
    H.append('<section><h2><span class="num">·</span> Contents</h2><div class="toc">')
    for n, t, a in secs:
        H.append(f'<a href="#{a}">{n}. {t}</a>')
    H.append('</div></section>')

    # ---- 1. Actionability ----
    H.append('<section id="act"><h2><span class="num">01</span> Clinical actionability</h2>')
    H.append('<p class="lead">Prioritised somatic and germline findings with therapeutic or '
             'diagnostic implications. Tiers follow PCGR/AMP-ACMG. Always confirm against the '
             'native PCGR report before sign-out.</p>')
    if act.get("items"):
        H.append('<table class="sortable"><thead><tr><th>Category</th><th>Finding</th>'
                 '<th>Tier</th><th>Evidence / implication</th></tr></thead><tbody>')
        for it in act["items"]:
            tcls = _tier_class(it.get("tier"))
            H.append(f'<tr><td>{esc(it["category"])}</td><td><b>{esc(it["finding"])}</b></td>'
                     f'<td><span class="tier {tcls}">{esc(it.get("tier"))}</span></td>'
                     f'<td>{esc(it.get("evidence"))}</td></tr>')
        H.append('</tbody></table>')
    else:
        H.append('<p class="na">No tier 1/2 actionable findings detected by automated roll-up. '
                 'Manual review of the PCGR report is still advised.</p>')
    H.append('</section>')

    # ---- 2. QC ----
    H.append('<section id="qc"><h2><span class="num">02</span> Quality &amp; provenance</h2>')
    H.append('<p class="lead">Automated QC flags. These gate clinical interpretation.</p>')
    flags = report.get("qc_flags", [])
    if not flags:
        H.append('<div class="flag ok">No QC issues detected.</div>')
    for f in flags:
        H.append(f'<div class="flag {esc(f["level"])}">{esc(f["message"])}</div>')
    cov = report.get("coverage", {})
    if cov.get("available"):
        H.append('<h3>Coverage</h3><table><thead><tr><th>Sample</th><th>Mean depth</th>'
                 '<th>Method</th></tr></thead><tbody>')
        for role, s in cov.get("samples", {}).items():
            H.append(f'<tr><td>{esc(role)}</td><td>{_fmt(s.get("mean_depth"),1)}x</td>'
                     f'<td>{esc(cov.get("method"))}</td></tr>')
        H.append('</tbody></table>')
    else:
        H.append(f'<p class="na">Coverage not computed ({esc(cov.get("note") or "no BAM/mosdepth found")}).</p>')
    # ---- tumour-normal sample identity (RAITA, else somalier) ----
    raita = report.get("raita", {})
    som = report.get("somalier", {})
    if raita.get("available"):
        H.append('<h3>Tumour-normal sample identity (RAITA)</h3>')
        ok = raita.get("verdict") == "PASS"
        cls = "flag ok" if ok else "flag error"
        conc = raita.get("concordant")
        conc_txt = ("concordant" if conc else "discordant") if conc is not None else "—"
        bits = []
        if raita.get("n_sites") is not None:
            bits.append(f'{raita["n_sites"]} germline-het sites compared')
        bits.append(f'tumour {conc_txt} with the germline')
        verdict = raita.get("verdict") or "—"
        H.append(f'<div class="{cls}"><b>{esc(verdict)}</b> — {esc("; ".join(bits))}.</div>')
        H.append(f'<p class="small muted">Source: <code>{esc(raita.get("source"))}</code>. '
                 'RAITA rules out contamination / wrong-sample but does not by itself exclude '
                 'a human-human swap between related individuals.</p>')
    elif som.get("available"):
        H.append('<h3>Tumour-normal sample identity (somalier)</h3>')
        rel = som.get("tn_relatedness")
        if rel is not None:
            # red (error) below 0.9; otherwise consistent with the same individual
            low = rel < 0.9
            cls = "flag error" if low else "flag ok"
            verdict = ("LOWER than expected for a matched pair (< 0.9) — possible sample "
                       "swap or contamination, review before trusting somatic subtraction"
                       if low else "consistent with the same individual")
            H.append(f'<div class="{cls}">Tumour/normal relatedness = {rel:.3f} — {verdict}.</div>')
        if som.get("samples"):
            H.append('<table><thead><tr><th>Sample</th><th>Sex</th><th>Het sites</th>'
                     '<th>Depth</th></tr></thead><tbody>')
            for s in som["samples"]:
                H.append(f'<tr><td>{esc(s.get("sample"))}</td><td>{esc(s.get("sex"))}</td>'
                         f'<td>{_fmt(s.get("n_het"),0)}</td><td>{_fmt(s.get("depth_mean"),1)}x</td></tr>')
            H.append('</tbody></table>')
    H.append(f'<p class="small muted">{report["provenance"]["n_files"]} input files fingerprinted '
             '(sha256) — see JSON provenance block for the full audit trail.</p>')
    H.append('</section>')

    # ---- 3. Purity/ploidy ----
    H.append('<section id="pp" class="pagebreak"><h2><span class="num">03</span> Tumour purity &amp; ploidy</h2>')
    H.append('<p class="lead">Three orthogonal callers reconciled. Divergence is itself a quality signal.</p>')
    H.append('<table><thead><tr><th>Caller</th><th>Data</th><th>Purity</th><th>Ploidy</th><th>Note</th></tr></thead><tbody>')
    for e in report["purity_ploidy"].get("estimates", []):
        H.append(f'<tr><td><b>{esc(e["caller"])}</b></td><td class="small">{esc(e.get("data_type"))}</td>'
                 f'<td>{_fmt(e.get("purity"),pct=True) if e.get("purity") is not None else "—"}</td>'
                 f'<td>{_fmt(e.get("ploidy"))}</td><td class="small muted">{esc(e.get("note"))}</td></tr>')
    H.append('</tbody></table>')
    if report["ascat"].get("plots"):
        H.append('<div class="two">')
        if fig.get("ascat_sunrise"):
            H.append(f'<img class="fig" src="{fig["ascat_sunrise"]}" alt="ASCAT sunrise">')
        if fig.get("ascat_profile"):
            H.append(f'<img class="fig" src="{fig["ascat_profile"]}" alt="ASCAT profile">')
        H.append('</div>')
    if report.get("wakhan", {}).get("heatmap_html"):
        H.append(f'<p class="small muted">Interactive wakhan purity/ploidy heatmap: '
                 f'<code>{esc(os.path.relpath(report["wakhan"]["heatmap_html"], cfg.input_dir))}</code></p>')
    H.append('</section>')

    # ---- 4. SNV/indel ----
    snv = report["pcgr"].get("snv_indel", {})
    H.append('<section id="snv" class="pagebreak"><h2><span class="num">04</span> Somatic SNV / indel</h2>')
    H.append(f'<p class="lead">{snv.get("n_total",0)} PCGR-annotated somatic variants '
             f'(source: <code>{esc(snv.get("source"))}</code>). Showing prioritised set.</p>')
    if fig.get("vaf"):
        H.append(f'<img class="fig" style="max-width:480px" src="{fig["vaf"]}" alt="VAF distribution">')
    variants = snv.get("variants", [])

    def _vtable(vs):
        rows = ['<table class="sortable"><thead><tr><th>Gene</th><th>Consequence</th>'
                '<th>Protein</th><th>Tier</th><th>Oncogenicity</th><th>VAF</th><th>DP</th>'
                '<th>Evidence</th></tr></thead><tbody>']
        for v in vs:
            rows.append(f'<tr><td><b>{esc(v.get("gene"))}</b></td><td class="small">{esc(v.get("consequence"))}</td>'
                        f'<td class="mono">{esc(v.get("hgvsp"))}</td>'
                        f'<td><span class="tier {_tier_class(v.get("tier"))}">{esc(v.get("tier"))}</span></td>'
                        f'<td class="small">{esc(v.get("oncogenicity"))}</td>'
                        f'<td>{_fmt(v.get("vaf"),pct=True) if v.get("vaf") is not None else "—"}</td>'
                        f'<td>{_fmt(v.get("dp"),0) if v.get("dp") is not None else "—"}</td>'
                        f'<td class="small">{esc(v.get("evidence"))}</td></tr>')
        rows.append('</tbody></table>')
        return "".join(rows)

    def _vsort(v):
        return -(v.get("vaf") or 0)
    actionable = sorted([v for v in variants if _tier_class(v.get("tier")) in ("t1", "t2")], key=_vsort)
    vus = sorted([v for v in variants if _tier_class(v.get("tier")) in ("t3", "t4", "tx")], key=_vsort)[:400]

    if actionable or vus:
        if for_pdf:
            # PDF: Tier 1/2 only — the Tier 3+ / VUS appendix is huge and lives in the HTML
            H.append(_vtable(actionable) if actionable
                     else '<p class="na">No Tier 1/2 (actionable) variants.</p>')
            if vus:
                H.append(f'<p class="small muted">{len(vus)} Tier 3+ / VUS variant(s) are omitted '
                         'from this PDF for length — see the HTML report or the JSON output.</p>')
        else:
            H.append('<div class="tabs" role="tablist">'
                     f'<button class="tab-btn active" data-tab="tab-tier12">Tier 1/2 ({len(actionable)})</button>'
                     f'<button class="tab-btn" data-tab="tab-vus">VUS · Tier 3+ ({len(vus)})</button></div>')
            H.append('<div id="tab-tier12" class="tab-panel active">')
            H.append(_vtable(actionable) if actionable
                     else '<p class="na">No Tier 1/2 (actionable) variants.</p>')
            H.append('</div>')
            H.append('<div id="tab-vus" class="tab-panel">')
            H.append('<h3 class="appendix-title">Appendix — variants of uncertain significance (Tier 3+)</h3>')
            H.append(_vtable(vus) if vus else '<p class="na">No Tier 3+ variants.</p>')
            H.append('</div>')
    else:
        H.append('<p class="na">No annotated somatic SNV/indel table found.</p>')
    # ---- clonal architecture (PyClone-VI) ----
    clo = report.get("clonal", {})
    if clo.get("available"):
        nsub = clo.get("n_subclonal_clusters", 0)
        H.append('<h3>Clonal architecture (PyClone-VI)</h3>')
        H.append(f'<p class="small">{clo.get("n_clusters",0)} mutation cluster(s); '
                 f'{nsub} subclonal (mean CCF &lt; 0.9). Source: '
                 f'<code>{esc(clo.get("source"))}</code>.</p>')
        H.append('<table class="sortable"><thead><tr><th>Cluster</th><th>Mutations</th>'
                 '<th>Mean CCF</th><th>Status</th></tr></thead><tbody>')
        for i, c in enumerate(clo.get("clusters", [])):
            ccf = c.get("mean_ccf")
            status = ("clonal" if i == 0 and (ccf or 0) >= 0.9
                      else "subclonal" if (ccf or 0) < 0.9 else "—")
            H.append(f'<tr><td>{esc(c.get("cluster_id"))}</td><td>{_fmt(c.get("n_mutations"),0)}</td>'
                     f'<td>{_fmt(ccf,3)}</td><td>{esc(status)}</td></tr>')
        H.append('</tbody></table>')
        H.append('<p class="small muted">Clonal vs subclonal status affects driver '
                 'interpretation — a subclonal driver is a later event.</p>')
    H.append('</section>')

    # ---- 5. TMB / MSI ----
    H.append('<section id="tmbmsi"><h2><span class="num">05</span> Tumour mutational burden &amp; MSI</h2>')
    H.append('<div class="two"><div>')
    H.append('<h3>TMB (PCGR)</h3><dl class="kvlist">')
    H.append(f'<dt>TMB</dt><dd>{_fmt(tmb.get("tmb_per_mb"),1)} /Mb</dd>')
    H.append(f'<dt>Measure</dt><dd>{esc(tmb.get("measure"))}</dd>')
    H.append(f'<dt>Somatic n</dt><dd>{_fmt(tmb.get("n_somatic"),0)}</dd></dl>')
    if fig.get("tmb_context"):
        H.append(f'<img class="fig" src="{fig["tmb_context"]}" alt="TMB in context">')
        H.append('<p class="small muted">Sample TMB (red line) against approximate TCGA median '
                 'burden per tumour type; dashed line = TMB-high (10/Mb).</p>')
    H.append('</div><div>')
    H.append('<h3>MSI</h3><dl class="kvlist">')
    _msi_high = msi.get("status") == "MSI"
    _st = esc(msi.get("status"))
    _st_html = f'<b class="msi-high">{_st}</b>' if _msi_high else _st
    H.append(f'<dt>Status</dt><dd>{_st_html}'
             f'{(" (" + esc(msi.get("status_raw")) + ")") if msi.get("status_raw") else ""}</dd>')
    H.append(f'<dt>Probability</dt><dd>{_fmt(msi.get("probability"),2)}</dd>')
    H.append(f'<dt>Source</dt><dd class="mono small">{esc(msi.get("source"))}</dd></dl>')
    if _msi_high:
        H.append('<p class="msi-high"><b>This tumour sample is predicted to be '
                 'microsatellite unstable.</b></p>')
    H.append(f'<p class="small muted">{esc(msi.get("threshold_note"))}</p>')
    H.append('</div></div></section>')

    # ---- 6. Signatures ----
    H.append('<section id="sig"><h2><span class="num">06</span> Mutational signatures</h2>')
    if sigs:
        if fig.get("signatures") or fig.get("signatures_pie"):
            H.append('<div class="figrow">')
            if fig.get("signatures"):
                H.append(f'<img class="fig" src="{fig["signatures"]}" alt="signature contributions">')
            if fig.get("signatures_pie"):
                H.append(f'<img class="fig" src="{fig["signatures_pie"]}" alt="signature composition">')
            H.append('</div>')
        # PDF: keep section 6 to one page — show the top contributors only
        sig_cap = 6 if for_pdf else 15
        H.append('<table class="sortable"><thead><tr><th>Signature</th><th>Contribution</th>'
                 '<th>Aetiology</th></tr></thead><tbody>')
        for s in sigs[:sig_cap]:
            H.append(f'<tr><td><b>{esc(s.get("signature"))}</b></td>'
                     f'<td>{_fmt(s.get("weight"),3)}</td><td class="small">{esc(s.get("aetiology"))}</td></tr>')
        H.append('</tbody></table>')
        if for_pdf and len(sigs) > sig_cap:
            H.append(f'<p class="small muted">{len(sigs) - sig_cap} further signature(s) with '
                     'smaller contributions are in the HTML report and the JSON output.</p>')
    else:
        H.append('<p class="na">No PCGR signature output found.</p>')
    # ---- SigProfilerAssignment: SBS96 / CNV48 / SV32 refits ----
    sp = report.get("sigprofiler", {})
    if sp.get("available"):
        H.append('<h3>SigProfilerAssignment refit</h3>')
        if for_pdf:
            H.append('<p class="small muted">Independent SigProfilerAssignment refits. '
                     'Aetiologies (COSMIC / Steele 2022 / RefSig) are proposed associations, '
                     'not diagnoses.</p>')
        else:
            H.append('<p class="small muted">Independent signature refits from '
                     'SigProfilerAssignment. Each signature is annotated with the process it '
                     'is currently thought to indicate — SBS per COSMIC v3, copy-number per '
                     'Steele et al. (Nature 2022), rearrangement per the RefSig / Signal '
                     'reference set. Aetiologies are proposed associations, not diagnoses.</p>')
        spf = fig.get("sigprof_figs") or {}
        for kind in ("sbs", "cnv", "sv"):
            fam = sp.get(kind, {})
            if not fam.get("available"):
                continue
            H.append(f'<h4 class="small">{esc(fam.get("label"))} — {fam["n"]} signature(s)</h4>')
            if not for_pdf:
                H.append('<div class="figrow"><div style="flex:1 1 460px">')
            H.append('<table class="sortable"><thead><tr><th>Signature</th><th>Count</th>'
                     '<th>Proportion</th><th>Indicates</th></tr></thead><tbody>')
            for r in fam["rows"]:
                H.append(f'<tr><td><b>{esc(r.get("signature"))}</b></td>'
                         f'<td>{_fmt(r.get("count"),0)}</td>'
                         f'<td>{_fmt(r.get("proportion"),pct=True) if r.get("proportion") is not None else "—"}</td>'
                         f'<td class="small">{esc(r.get("aetiology"))}</td></tr>')
            H.append('</tbody></table>')
            if not for_pdf:
                H.append('</div>')
                if spf.get(kind):
                    H.append(f'<img class="fig" src="{spf[kind]}" alt="{esc(kind)} composition">')
                H.append('</div>')
                H.append(f'<p class="small muted">Source: <code>{esc(fam.get("source"))}</code></p>')
        if sp.get("sv", {}).get("available"):
            H.append('<p class="small muted">SV reference signatures are Illumina-derived; '
                     'ONT breakpoint imprecision shifts size bins, so read as indicative.</p>')
    H.append('</section>')
    H.append('<section id="cna" class="pagebreak"><h2><span class="num">07</span> Copy-number alterations</h2>')
    if fig.get("cn_genome"):
        H.append(f'<img class="fig" src="{fig["cn_genome"]}" alt="genome-wide CN">')
    # Wakhan: inline the actual interactive plotly plots (self-contained — the
    # full plot HTML is embedded via srcdoc, so the report needs no wakhan/ dir).
    wp = report.get("wakhan", {})
    wk_plots = [("Haplotype-specific copy number (details)", wp.get("details_html")),
                ("Subclonal copy number & breakpoints", wp.get("subclonal_html"))]
    wk_plots = [(t, p) for t, p in wk_plots if p]
    wk_figs = fig.get("wakhan_snaps") or []
    if wk_plots or wk_figs:
        H.append(f'<h3>Wakhan copy number — best solution {esc(wp.get("solution"))}</h3>')
        if for_pdf:
            # PDF: Wakhan's OWN PDF plots, rasterised (falls back to a snapshot of
            # the interactive plot when Wakhan was run without PDF output)
            if wk_figs:
                from_pdf = bool(wp.get("details_pdf") or wp.get("subclonal_pdf")
                                or wp.get("breakpoints_pdf") or wp.get("genes_pdf"))
                for title, img in wk_figs:
                    H.append(f'<p class="small muted">{esc(title)}</p>')
                    H.append(f'<img class="fig wide" src="{img}" alt="{esc(title)}">')
                H.append('<p class="small muted">'
                         + ("Wakhan's own PDF plots. " if from_pdf
                            else "Static snapshots of the interactive Wakhan plots. ")
                         + 'The HTML report carries the fully interactive versions.</p>')
            else:
                for title, p in wk_plots:
                    rel = os.path.relpath(p, cfg.out_dir)
                    H.append(f'<p class="small muted">{esc(title)} — interactive plot at '
                             f'<code>{esc(rel)}</code> (open the HTML report to view).</p>')
        else:
            for title, p in wk_plots:
                doc = None
                try:
                    with open(p, encoding="utf-8", errors="replace") as fh:
                        doc = fh.read()
                except Exception as exc:
                    LOG.warning("Wakhan plot read failed (%s): %s", p, exc)
                H.append(f'<p class="small muted">{esc(title)}</p>')
                if doc:
                    doc, offline_ok = _inline_plotly(doc)   # embed plotly.js -> no network
                    if not offline_ok:
                        H.append('<p class="small muted">(interactive plot loads plotly.js from '
                                 'the CDN — install the <code>plotly</code> package to make it '
                                 'fully offline)</p>')
                    srcdoc = doc.replace("&", "&amp;").replace('"', "&quot;")
                    H.append(f'<iframe class="plotly-embed" title="{esc(title)}" '
                             f'srcdoc="{srcdoc}"></iframe>')
                else:
                    rel = os.path.relpath(p, cfg.out_dir)
                    H.append(f'<iframe class="plotly-embed" title="{esc(title)}" '
                             f'src="{esc(rel)}" loading="lazy"></iframe>')
    if fig.get("severus_loh"):
        loh = report.get("severus_loh", {})
        H.append(f'<img class="fig" src="{fig["severus_loh"]}" alt="genome-wide LOH">')
        H.append(f'<p class="small muted">{loh.get("n",0)} Severus LOH region(s), '
                 f'{loh.get("total_mb","—")} Mb total. Source: <code>{esc(loh.get("source"))}</code>.</p>')
    cnas = report["pcgr"].get("cna_genes", {}).get("items", [])
    def _is_real(v):
        return v is not None and str(v).strip() not in ("", ".", "nan", "NA", "None")
    # PCGR-style prioritisation: keep only genes with a real actionability tier
    # (TIER I/II) or a matched biomarker — not every gene overlapping a segment.
    cnas_show = [c for c in cnas if _tier_rank_str(c.get("tier")) <= 2 or _is_real(c.get("evidence"))]
    cnas_show.sort(key=lambda c: (_tier_rank_str(c.get("tier")), str(c.get("gene"))))
    cna_omitted = 0
    if for_pdf:
        # PDF: Tier 1/2 only — drop untiered biomarker-only rows to save space
        tiered = [c for c in cnas_show if _tier_rank_str(c.get("tier")) <= 2]
        cna_omitted = len(cnas_show) - len(tiered)
        cnas_show = tiered
    if cnas_show:
        paged = " paged" if len(cnas_show) > 10 else ""
        H.append(f'<h3>Actionable copy-number alterations (PCGR, AMP/ASCO/CAP) — {len(cnas_show)}</h3>')
        H.append(f'<table class="sortable{paged}" data-page-size="10"><thead><tr><th>Gene</th>'
                 '<th>Event</th><th>Type</th><th>Total CN</th><th>Cytoband</th><th>Tier</th>'
                 '<th>Biomarker evidence</th></tr></thead><tbody>')
        for c in cnas_show:
            ev = str(c.get("event") or "")
            cls = "gain" if "amp" in ev.lower() or "gain" in ev.lower() else ("loss" if "del" in ev.lower() or "loss" in ev.lower() else "")
            annot = c.get("evidence") if _is_real(c.get("evidence")) else "—"
            band = c.get("cytoband") if _is_real(c.get("cytoband")) else "—"
            etype = c.get("event_type") if _is_real(c.get("event_type")) else "—"
            H.append(f'<tr><td><b>{esc(c.get("gene"))}</b></td>'
                     f'<td class="{cls}">{esc(c.get("event"))}</td>'
                     f'<td class="small">{esc(etype)}</td>'
                     f'<td>{_fmt(c.get("total_cn"),0)}</td>'
                     f'<td class="mono small">{esc(band)}</td>'
                     f'<td><span class="tier {_tier_class(c.get("tier"))}">{esc(c.get("tier"))}</span></td>'
                     f'<td class="small">{esc(annot)}</td></tr>')
        H.append('</tbody></table>')
        if cna_omitted:
            H.append(f'<p class="small muted">{cna_omitted} further biomarker-matched but '
                     'untiered CNA(s) are omitted from this PDF for length — see the HTML '
                     'report or the JSON output.</p>')
    else:
        H.append('<p class="na">No CNAs reached AMP/ASCO/CAP Tier I/II or matched a biomarker.</p>')
    asc = report["ascat"]
    H.append(f'<p class="small muted">ASCAT: {asc.get("n_segments",0)} allele-specific segments; '
             f'{sum(1 for s in asc.get("segments",[]) if s.get("loh"))} with LOH.</p>')
    # ---- HRD (scarHRD) ----
    hrd = report.get("hrd", {})
    if hrd.get("available"):
        H.append('<h3>Homologous-recombination deficiency (scarHRD)</h3>')
        if hrd.get("hr_deficient"):
            H.append(f'<div class="flag warn">HRD-sum = {_fmt(hrd.get("hrd_sum"),0)} '
                     '(&ge;42) — HR-deficient by the common threshold; consider PARP relevance.</div>')
        H.append('<dl class="kvlist">')
        H.append(f'<dt>HRD-LOH</dt><dd>{_fmt(hrd.get("hrd_loh"),0)}</dd>')
        H.append(f'<dt>Telomeric AI</dt><dd>{_fmt(hrd.get("tai"),0)}</dd>')
        H.append(f'<dt>LST</dt><dd>{_fmt(hrd.get("lst"),0)}</dd>')
        H.append(f'<dt>HRD-sum</dt><dd><b>{_fmt(hrd.get("hrd_sum"),0)}</b> '
                 '(threshold &ge;42)</dd></dl>')
    H.append('</section>')

    # ---- 8. SV ----
    H.append('<section id="sv" class="pagebreak"><h2><span class="num">08</span> Structural variants &amp; fusions</h2>')
    svd = report.get("sv", {})
    if svd.get("callers"):
        # half-width translocation circos beside the caller table
        circ = fig.get("sv_circos")
        H.append('<div class="two"><div>')
        H.append('<table><thead><tr><th>Caller</th><th>PASS SVs</th><th>Breakdown</th>'
                 '<th>Source</th></tr></thead><tbody>')
        for caller, d in svd["callers"].items():
            bk = ", ".join(f"{k}:{v}" for k, v in sorted(d.get("svtype_counts", {}).items()))
            H.append(f'<tr><td><b>{esc(caller)}</b></td><td>{d.get("n_pass",0)}</td>'
                     f'<td class="small">{esc(bk)}</td><td class="mono small">{esc(d.get("source"))}</td></tr>')
        H.append('</tbody></table>')
        # Severus vs SAVANA concordance
        conc = svd.get("concordance")
        if conc:
            H.append('<h3>Caller concordance (SAVANA vs Severus)</h3><dl class="kvlist">')
            H.append(f'<dt>Concordant breakpoints</dt><dd>{conc.get("n_concordant",0)} '
                     f'({_fmt(conc.get("concordance_pct"),1)}%, ±{conc.get("window_bp")} bp)</dd>')
            H.append(f'<dt>SAVANA-only</dt><dd>{conc.get("n_savana_only",0)}</dd>')
            H.append(f'<dt>Severus-only</dt><dd>{conc.get("n_severus_only",0)}</dd></dl>')
            if conc.get("n_savana_only") or conc.get("n_severus_only"):
                H.append('<p class="small muted">Single-caller SVs are not corroborated by '
                         'the other caller — confirm before reporting as driver events.</p>')
        H.append('</div><div>')
        if circ:
            H.append('<h3>Translocations</h3>')
            H.append(f'<div style="width:50%;min-width:300px">{circ}</div>')
            H.append('<p class="small muted">Inter-chromosomal rearrangements drawn as chords '
                     '(SAVANA somatic BEDPE).</p>')
        H.append('</div></div>')
        # fusions annotated by clinical relevance
        # Prioritised gene fusions from simple_sv_annotation (*.genefusions.txt)
        gf = report.get("gene_fusions", {})
        if gf.get("available"):
            H.append(f'<h3>Bidirectional gene fusions (simple_sv_annotation) — {gf["n"]} calls, '
                     f'{gf.get("n_prioritised", 0)} prioritised</h3>')
            paged = " paged" if gf["n"] > 12 else ""
            H.append(f'<table class="sortable{paged}" data-page-size="12"><thead><tr>'
                     '<th>Gene(s)</th><th>SV type</th><th>Annotation</th><th>Transcript</th>'
                     '<th>Priority</th></tr></thead><tbody>')
            for r in gf["rows"]:
                if r["priority"] == "PRIORITISED":
                    badge = f'<span class="tier t{min(r["level"], 4)}">Prioritised · L{r["level"]}</span>'
                else:
                    badge = f'<span class="small muted">not prioritised · L{r["level"]}</span>'
                H.append(f'<tr><td><b>{esc(r["genes"])}</b></td>'
                         f'<td class="mono small">{esc(r["svtype"])}</td>'
                         f'<td class="small">{esc(r["annotation"])}</td>'
                         f'<td class="mono small">{esc(r["transcript"])}</td>'
                         f'<td>{badge}</td></tr>')
            H.append('</tbody></table>')
            H.append(f'<p class="small muted">Source: <code>{esc(gf.get("source"))}</code> '
                     '(AstraZeneca-NGS simple_sv_annotation); priority level 1 = highest. '
                     'Filtered to <code>BIDIRECTIONAL_GENE_FUSION</code> records only — '
                     f'{gf.get("n_filtered_out", 0)} other fusion annotations were excluded.</p>')

        fusions = svd.get("fusions_annotated", [])
        if fusions:
            H.append('<h3>Candidate gene fusions</h3><table><thead><tr><th>Type</th>'
                     '<th>Locus</th><th>Annotation</th><th>Clinical relevance</th>'
                     '</tr></thead><tbody>')
            for f in fusions[:40]:
                cl = f.get("clinical")
                clcell = (f'<span class="tier t1">{esc(cl["relevance"])}</span> '
                          f'<span class="small">{esc(cl["pair"])} — {esc(cl["note"])}</span>'
                          if cl else '<span class="small muted">not in curated set</span>')
                H.append(f'<tr><td>{esc(f.get("type"))}</td>'
                         f'<td class="mono small">{esc(f.get("chrom"))}:{esc(f.get("pos"))}</td>'
                         f'<td class="small">{esc(f.get("annotation"))}</td>'
                         f'<td>{clcell}</td></tr>')
            H.append('</tbody></table>')
            if svd.get("n_actionable_fusions"):
                H.append(f'<p class="small"><b>{svd["n_actionable_fusions"]}</b> fusion(s) match '
                         'the curated actionable/oncogenic set.</p>')
        else:
            H.append('<p class="small muted">No annotated fusions parsed from SV INFO; '
                     'run an SV gene-annotation step (e.g. AnnotSV) to populate this.</p>')
    else:
        H.append('<p class="na">No SV VCFs found.</p>')
    H.append('</section>')

    # ---- 9. Germline ----
    H.append('<section id="germ" class="pagebreak"><h2><span class="num">09</span> Germline predisposition (CPSR)</h2>')
    cpsr = report.get("cpsr", {})
    if cpsr.get("available"):
        rep = [v for v in cpsr.get("variants", []) if v.get("reportable")]
        H.append(f'<p class="lead">{cpsr.get("n_reportable",0)} pathogenic/likely-pathogenic '
                 'germline variant(s) in cancer-predisposition genes.</p>')
        show = rep if rep else cpsr.get("variants", [])[:20]
        if show:
            H.append('<table class="sortable"><thead><tr><th>Gene</th><th>Protein</th>'
                     '<th>Classification</th><th>Genotype</th><th>Trait</th></tr></thead><tbody>')
            for v in show:
                H.append(f'<tr><td><b>{esc(v.get("gene"))}</b></td>'
                         f'<td class="mono">{esc(v.get("hgvsp"))}</td>'
                         f'<td>{esc(v.get("classification"))}</td><td>{esc(v.get("genotype"))}</td>'
                         f'<td class="small">{esc(v.get("trait"))}</td></tr>')
            H.append('</tbody></table>')
        else:
            H.append('<p class="na">No germline variants parsed.</p>')
    else:
        H.append('<p class="na">CPSR output not found.</p>')
    H.append('</section>')

    # ---- 10. Methylation ----
    H.append('<section id="meth" class="pagebreak"><h2><span class="num">10</span> Methylation &amp; classification</h2>')
    nd = report.get("nanodx", {})
    H.append('<div class="two"><div><h3>nanoDx classifier</h3>')
    if nd.get("available"):
        tc = nd.get("top_class")
        desc = nd.get("top_class_desc")
        tc_disp = f'{esc(tc)}' + (f' <span class="small muted">({esc(desc)})</span>' if desc else '')
        H.append('<dl class="kvlist">')
        H.append(f'<dt>Top class</dt><dd><b>{tc_disp if tc else "—"}</b></dd>')
        H.append(f'<dt>Score</dt><dd>{_fmt(nd.get("score"),3) if nd.get("score") is not None else "—"}</dd></dl>')
        if nd.get("votes"):
            H.append('<table><thead><tr><th>Class</th><th>Meaning</th><th>Votes</th></tr></thead><tbody>')
            for v in nd["votes"]:
                H.append(f'<tr><td><b>{esc(v.get("class"))}</b></td>'
                         f'<td class="small">{esc(v.get("class_desc") or "—")}</td>'
                         f'<td>{_fmt(v.get("votes"),3)}</td></tr>')
            H.append('</tbody></table>')
    else:
        H.append('<p class="na">No nanoDx result.</p>')
    H.append('</div><div><h3>Promoter methylation</h3>')
    meth = report.get("methylation", {})
    if meth.get("promoters"):
        H.append('<table><thead><tr><th>Gene</th><th>Region</th><th>Mean %</th>'
                 '<th>CpGs</th></tr></thead><tbody>')
        for p in meth["promoters"]:
            star = "" if p.get("coords_verified") else ' <span class="badge">est.</span>'
            H.append(f'<tr><td><b>{esc(p.get("gene"))}</b>{star}</td>'
                     f'<td class="mono small">{esc(p.get("region"))}</td>'
                     f'<td>{_fmt(p.get("mean_methylation_pct"),1)}</td>'
                     f'<td>{esc(p.get("n_cpg"))}</td></tr>')
        H.append('</tbody></table>')
        if meth.get("dmr"):
            H.append(f'<p class="small muted">modkit differential calls: '
                     f'{meth["dmr"].get("n_records","—")} records.</p>')
    else:
        H.append('<p class="na">No bedMethyl promoter summary available.</p>')
    H.append('</div></div></section>')

    # ---- 11. Immune ----
    H.append('<section id="immune" class="pagebreak"><h2><span class="num">11</span> HLA &amp; neoantigens</h2>')
    imm = report.get("immune", {})
    H.append('<div class="two"><div><h3>HLA typing (SpecImmune)</h3>')
    if imm.get("hla", {}).get("available"):
        H.append('<table><thead><tr><th>Locus</th><th>Allele 1</th><th>Allele 2</th></tr></thead><tbody>')
        for a in imm["hla"]["alleles"]:
            H.append(f'<tr><td><b>{esc(a.get("locus"))}</b></td>'
                     f'<td class="mono">{esc(a.get("allele_1"))}</td>'
                     f'<td class="mono">{esc(a.get("allele_2"))}</td></tr>')
        H.append('</tbody></table>')
    else:
        H.append('<p class="na">HLA result not parsed.</p>')
    H.append('</div><div><h3>Neoantigens (pVACseq)</h3>')
    neo = imm.get("neoantigen", {})
    bt = neo.get("binder_threshold_nM", 500)
    st = neo.get("strong_threshold_nM", 50)
    # binder summary at the top of the section
    for cls_name, key in [("Class I", "class_I"), ("Class II", "class_II")]:
        d = neo.get(key, {})
        if not d.get("available"):
            continue
        nb, ns = d.get("n_binders"), d.get("n_strong")
        parts = [f'{d.get("n_candidates","—")} candidates']
        if nb is not None:
            parts.append(f'<b>{nb}</b> binders (&le;{bt:.0f} nM)')
        if ns is not None:
            parts.append(f'<b>{ns}</b> strong (&le;{st:.0f} nM)')
        H.append(f'<p><b>{cls_name}:</b> {" · ".join(parts)}</p>')
    H.append(f'<p class="small muted">Binding threshold {bt:.0f} nM (pVACseq default); '
             f'strong binders &le;{st:.0f} nM highlighted.</p>')
    ci_top = neo.get("class_I", {}).get("top", [])
    if ci_top:
        H.append('<table><thead><tr><th>Gene</th><th>Peptide</th><th>IC50 (nM)</th></tr></thead><tbody>')
        for t in ci_top:
            rowcls = ' class="gain"' if t.get("strong") else ""
            H.append(f'<tr{rowcls}><td>{esc(t.get("gene"))}</td>'
                     f'<td class="mono small">{esc(t.get("peptide"))}</td>'
                     f'<td>{_fmt(t.get("ic50_nM"),1)}</td></tr>')
        H.append('</tbody></table>')
    H.append('<p class="small muted">Note: presentation depends on HLA allele '
             'retention — see the HLA-LOH table below.</p>')
    H.append('</div></div>')
    # ---- KIR / IG-TR typing (CYP now lives in section 12: Pharmacogenomics) ----
    for fam, title, note in [("kir", "KIR typing (SpecImmune)", ""),
                             ("ig_tr", "IG / TCR typing (SpecImmune)",
                              "Germline immunoglobulin / T-cell receptor gene alleles.")]:
        fd = imm.get(fam, {})
        if not fd.get("available"):
            continue
        rows = fd.get("rows", [])
        # PDF: IG/TCR is a long, low-yield table — summarise the loci in prose instead
        if for_pdf and fam == "ig_tr":
            genes = [str(r.get("gene")) for r in rows if r.get("gene")]
            fams = sorted({re.match(r"(IGH|IGK|IGL|TRA|TRB|TRG|TRD)", g).group(1)
                           for g in genes if re.match(r"(IGH|IGK|IGL|TRA|TRB|TRG|TRD)", g)})
            H.append(f'<h3>{title}</h3>')
            H.append(f'<p class="small">{len(genes)} immunoglobulin / T-cell receptor genes were '
                     f'sequenced and typed{(" across the " + ", ".join(fams) + " loci") if fams else ""}: '
                     f'{esc(", ".join(genes[:25]))}'
                     f'{", …" if len(genes) > 25 else ""}. '
                     'Per-gene allele calls are available in the HTML report and the JSON output.</p>')
            continue
        H.append(f'<h3>{title}</h3>')
        has_extra = any(r.get("extra") for r in rows)
        paged = " paged" if fam == "ig_tr" and len(rows) > 10 else ""
        H.append(f'<table class="sortable{paged}" data-page-size="10"><thead><tr><th>Gene</th>'
                 '<th>Allele 1</th><th>Allele 2</th>' +
                 ('<th>Annotation</th>' if has_extra else '') + '</tr></thead><tbody>')
        for r in rows:
            extra = f'<td class="small">{esc(r.get("extra"))}</td>' if has_extra else ''
            H.append(f'<tr><td><b>{esc(r.get("gene"))}</b></td>'
                     f'<td class="mono">{esc(r.get("allele_1"))}</td>'
                     f'<td class="mono">{esc(r.get("allele_2"))}</td>{extra}</tr>')
        H.append('</tbody></table>')
        if note:
            H.append(f'<p class="small muted">{note}</p>')
    # ---- HLA loss of heterozygosity ----
    hl = report.get("hla_loh", {})
    if hl.get("available"):
        H.append('<h3>HLA loss of heterozygosity (phased BAF)</h3>')
        if hl.get("n_loh"):
            H.append(f'<div class="flag warn">{hl["n_loh"]} HLA gene(s) with LOH — predicted '
                     'neoantigens at those loci may not be presentable.</div>')
        H.append('<table class="sortable"><thead><tr><th>Gene</th><th>Call</th>'
                 '<th>CN state</th><th>Minor fraction</th><th>Sites used</th>'
                 '<th>p-value</th></tr></thead><tbody>')
        for g in hl.get("genes", []):
            call = str(g.get("call") or "")
            cls = "loss" if call.upper() == "LOH" else ""
            H.append(f'<tr><td><b>{esc(g.get("gene"))}</b></td>'
                     f'<td class="{cls}">{esc(call)}</td>'
                     f'<td class="mono small">{esc(g.get("cn_state"))}</td>'
                     f'<td>{_fmt(g.get("obs_minor_fraction"),3)}</td>'
                     f'<td>{_fmt(g.get("n_sites_used"),0)}</td>'
                     f'<td class="small">{esc(g.get("pvalue"))}</td></tr>')
        H.append('</tbody></table>')
        H.append(f'<p class="small muted">Method: {esc(hl.get("method"))}. '
                 f'Source: <code>{esc(hl.get("source"))}</code>.</p>')
    H.append('</section>')

    # ---- 12. Pharmacogenomics (PharmCAT + SpecImmune CYP) ----
    pc = report.get("pharmcat", {})
    cyp = report.get("immune", {}).get("cyp", {})
    if pc.get("available") or cyp.get("available"):
        H.append('<section id="pgx" class="pagebreak"><h2><span class="num">12</span> '
                 'Pharmacogenomics</h2>')
        H.append('<p class="lead">Germline pharmacogenomic calls. PharmCAT provides '
                 'CPIC/DPWG-aligned diplotypes and prescribing phenotypes; SpecImmune provides '
                 'long-read CYP star-allele calls as an orthogonal check.</p>')
        if pc.get("available"):
            H.append(f'<h3>PharmCAT — {pc["n"]} genes, '
                     f'{pc.get("n_actionable", 0)} with a non-normal phenotype</h3>')
            paged = " paged" if pc["n"] > 12 else ""
            H.append(f'<table class="sortable{paged}" data-page-size="12"><thead><tr>'
                     '<th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Activity score</th>'
                     '<th>Flags</th></tr></thead><tbody>')
            for r in pc["rows"]:
                phen = r.get("phenotype")
                pcell = (f'<b class="loss">{esc(phen)}</b>' if r.get("actionable")
                         else f'{esc(phen) if phen else "—"}')
                flags = []
                if r.get("uncertain"):
                    flags.append('<span class="tier tx">ambiguous call</span>')
                if r.get("outside_call") and str(r["outside_call"]).lower() == "yes":
                    flags.append('<span class="small muted">outside call</span>')
                if r.get("missing_positions"):
                    flags.append('<span class="small muted">missing positions</span>')
                H.append(f'<tr><td><b>{esc(r.get("gene"))}</b></td>'
                         f'<td class="mono small">{esc(r.get("diplotype"))}</td>'
                         f'<td>{pcell}</td>'
                         f'<td>{esc(r.get("activity_score")) if r.get("activity_score") else "—"}</td>'
                         f'<td class="small">{" ".join(flags) if flags else "—"}</td></tr>')
            H.append('</tbody></table>')
            H.append(f'<p class="small muted">Source: <code>{esc(pc.get("source"))}</code>. '
                     'Phenotype/activity score are the recommendation-lookup values where '
                     'PharmCAT provides them. Ambiguous (" OR ") diplotypes arise from unphased '
                     'or incomplete data and should be treated as uncalled.</p>')
        if cyp.get("available"):
            H.append('<h3>CYP star alleles (SpecImmune, long-read)</h3>')
            H.append('<table class="sortable"><thead><tr><th>Gene</th><th>Allele 1</th>'
                     '<th>Allele 2</th><th>Metaboliser (derived)</th><th>PharmCAT phenotype</th>'
                     '</tr></thead><tbody>')
            pc_by_gene = {str(r.get("gene", "")).upper(): r for r in pc.get("rows", [])}
            for r in cyp.get("rows", []):
                pcr = pc_by_gene.get(str(r.get("gene") or "").upper())
                pcp = pcr.get("phenotype") if pcr else None
                agree = ""
                if pcp and r.get("extra"):
                    same = pcp.strip().lower()[:6] == str(r["extra"]).strip().lower()[:6]
                    agree = ('<span class="small muted"> · concordant</span>' if same else
                             '<span class="small loss"> · discordant</span>')
                H.append(f'<tr><td><b>{esc(r.get("gene"))}</b></td>'
                         f'<td class="mono">{esc(r.get("allele_1"))}</td>'
                         f'<td class="mono">{esc(r.get("allele_2"))}</td>'
                         f'<td class="small">{esc(r.get("extra")) if r.get("extra") else "—"}</td>'
                         f'<td class="small">{(esc(pcp) + agree) if pcp else "—"}</td></tr>')
            H.append('</tbody></table>')
            H.append('<p class="small muted">SpecImmune metaboliser status is derived from the '
                     'star-allele diplotype; PharmCAT is the authority for prescribing. '
                     'Review any discordance and confirm against CPIC/PharmGKB before clinical '
                     'use. Neither call is a validated pharmacogenomic diagnostic.</p>')
        H.append('</section>')

    # footer
    H.append(f'<div class="footer">Generated by {esc(REPORT_GENERATOR)} '
             f'(schema v{esc(report["schema_version"])}) · {esc(report["generated_at"])}<br>'
             'Research/clinical-decision-support output — not a substitute for an accredited '
             'diagnostic report. Confirm all findings against native PCGR/CPSR outputs and an '
             'orthogonal assay where required.<br>'
             'The BALTI pipeline is &copy; Andrew Beggs 2026.</div>')

    H.append('</div>')  # wrap
    if not for_pdf:
        H.append(f'<script>{SORT_JS}</script>')
    H.append('</body></html>')
    return "".join(H)


def write_pdf(html: str, out_pdf: str) -> bool:
    """Try WeasyPrint, then xhtml2pdf, then give up gracefully."""
    # WeasyPrint logs benign 'unknown property' warnings (box-shadow,
    # print-color-adjust) via its own logger; pin it to ERROR and stop it
    # propagating to the root handler so the console stays clean.
    for _ln in ("weasyprint", "weasyprint.css", "fontTools", "fontTools.subset"):
        _lg = logging.getLogger(_ln)
        _lg.setLevel(logging.ERROR)
        _lg.propagate = False
    try:
        from weasyprint import HTML as _WHTML
        _WHTML(string=html).write_pdf(out_pdf)
        LOG.info("PDF written via WeasyPrint: %s", out_pdf)
        return True
    except Exception as exc:
        LOG.warning("WeasyPrint unavailable/failed (%s); trying xhtml2pdf.", exc)
    try:
        from xhtml2pdf import pisa
        with open(out_pdf, "wb") as fh:
            pisa.CreatePDF(src=html, dest=fh)
        LOG.info("PDF written via xhtml2pdf: %s", out_pdf)
        return True
    except Exception as exc:
        LOG.warning("xhtml2pdf also failed (%s). Skipping PDF.", exc)
        return False


# ============================================================================
# JSON (figures stripped out)
# ============================================================================
def write_json(report: Dict[str, Any], out_json: str) -> None:
    clean = {k: v for k, v in report.items() if k != "_figures"}
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump(clean, fh, indent=2, default=str)
    LOG.info("JSON written: %s", out_json)


# ============================================================================
# CLI
# ============================================================================
def check_environment(verbose: bool = True) -> Dict[str, Any]:
    """Report which optional dependencies are usable. Missing matplotlib is the
    single most common cause of 'the plots aren't rendering' — every figure in
    this report (signatures, purity/ploidy, TMB, CN tracks, VAF, and the Wakhan
    snapshot fallback) is a matplotlib PNG."""
    status: Dict[str, Any] = {}
    status["pandas"] = "ok" if pd is not None else f"MISSING ({_PD_ERR})"
    status["matplotlib"] = "ok" if _HAVE_MPL else f"MISSING ({_MPL_ERR})"
    try:
        import plotly  # noqa: F401
        status["plotly"] = "ok"
    except Exception as exc:
        status["plotly"] = f"missing ({exc})"
    try:
        import kaleido  # noqa: F401
        status["kaleido"] = "ok"
    except Exception as exc:
        status["kaleido"] = f"missing ({exc})"
    try:
        import weasyprint  # noqa: F401
        status["weasyprint"] = "ok"
    except Exception as exc:
        status["weasyprint"] = f"missing ({exc})"
    try:
        import fitz  # noqa: F401  (PyMuPDF — renders Wakhan's PDF plots)
        status["pymupdf"] = "ok"
    except Exception as exc:
        status["pymupdf"] = (f"missing ({exc}) — needed to embed Wakhan's PDF plots "
                             "unless pdftoppm is available")
    status["pdftoppm"] = "ok" if shutil.which("pdftoppm") else "missing (poppler-utils)"
    if verbose:
        for k, v in status.items():
            print(f"  {k:<12} {v}")
    return status


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Generate a tumour-normal molecular report from the DETERMINE "
                    "ONT pipeline working directory.")
    ap.add_argument("-i", "--input-dir", required=False, help="Per-case Tissue/ directory")
    ap.add_argument("-o", "--out-dir", default=None, help="Output dir (default: input-dir)")
    ap.add_argument("--patient-id", default="PATIENT05")
    ap.add_argument("--sample-id", default=None,
                    help="ID shown in the report title (defaults to --patient-id)")
    ap.add_argument("--tumour-id", default="T9")
    ap.add_argument("--normal-id", default="K6")
    ap.add_argument("--genome-build", default="grch38")
    ap.add_argument("--skip-coverage", action="store_true",
                    help="Do not compute BAM coverage (faster)")
    ap.add_argument("--no-pdf", action="store_true", help="Skip PDF generation")
    ap.add_argument("--check-env", action="store_true",
                    help="Print dependency status and exit (diagnoses missing plots)")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)s %(message)s")
    # third-party PDF/font libraries are extremely chatty at INFO; pin them up
    for _noisy in ("fontTools", "fontTools.subset", "fontTools.ttLib", "PIL"):
        logging.getLogger(_noisy).setLevel(logging.WARNING)
    # WeasyPrint emits benign CSS-support warnings (box-shadow, screen @media,
    # print-color-adjust) that don't affect the rendered PDF — keep them quiet
    logging.getLogger("weasyprint").setLevel(logging.ERROR)
    for _noisy in ("kaleido", "choreographer", "logistro"):
        logging.getLogger(_noisy).setLevel(logging.ERROR)

    if args.check_env:
        print("Dependency status:")
        check_environment()
        return 0
    if not args.input_dir:
        ap.error("-i/--input-dir is required")

    # LOUD failure: a silently-missing matplotlib means every plot disappears
    if not _HAVE_MPL:
        LOG.warning("matplotlib is NOT usable (%s) — ALL figures (signatures, "
                    "purity/ploidy, TMB, copy-number, VAF, Wakhan snapshots) will be "
                    "omitted. Fix with: pip install matplotlib   (run "
                    "`%s --check-env` for full dependency status)",
                    _MPL_ERR, os.path.basename(sys.argv[0]))
    if pd is None:
        LOG.warning("pandas is NOT usable (%s) — most tables will be empty", _PD_ERR)

    cfg = Config(
        input_dir=args.input_dir,
        out_dir=args.out_dir or args.input_dir,
        patient_id=args.patient_id,
        display_id=args.sample_id or args.patient_id,
        tumour_id=args.tumour_id,
        normal_id=args.normal_id,
        genome_build=args.genome_build,
        skip_coverage=args.skip_coverage,
        no_pdf=args.no_pdf,
        samtools=args.samtools,
    )
    if not os.path.isdir(cfg.input_dir):
        ap.error(f"input dir not found: {cfg.input_dir}")
    os.makedirs(cfg.out_dir, exist_ok=True)

    report = assemble(cfg)

    stem = f"{cfg.patient_id}_{cfg.tumour_id}_{cfg.normal_id}.tn_report"
    out_html = os.path.join(cfg.out_dir, stem + ".html")
    out_json = os.path.join(cfg.out_dir, stem + ".json")
    out_pdf = os.path.join(cfg.out_dir, stem + ".pdf")

    html = render_html(report, cfg, for_pdf=False)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(html)
    LOG.info("HTML written: %s", out_html)

    write_json(report, out_json)

    pdf_ok = False
    if not cfg.no_pdf:
        pdf_html = render_html(report, cfg, for_pdf=True)
        pdf_ok = write_pdf(pdf_html, out_pdf)

    print("Report outputs:")
    print(f"  HTML : {out_html}")
    print(f"  JSON : {out_json}")
    print(f"  PDF  : {out_pdf if pdf_ok else '(not generated — install weasyprint)'}")
    print(f"  Actionable findings: {report['actionability']['n']}")
    _figs = report.get("_figure_stats") or {}
    if _figs:
        print(f"  Figures rendered: {_figs.get('rendered', 0)}/{_figs.get('total', 0)}")
        if not _figs.get("rendered"):
            print("  WARNING: no figures rendered — run with --check-env "
                  "(matplotlib is almost certainly missing)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
