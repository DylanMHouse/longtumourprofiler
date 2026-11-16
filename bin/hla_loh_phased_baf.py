#!/usr/bin/env python3
"""
hla_loh_phased_baf.py
=====================

Detect loss of heterozygosity across the HLA / MHC region from Oxford Nanopore
tumour-normal data using a *phased B-allele-frequency* approach, rather than the
per-base mismatch model used by LOHHLA.

Why this design for ONT
-----------------------
LOHHLA assigns each read to one of the two homologous HLA alleles by counting
mismatches at allele-distinguishing bases. That depends on Illumina-grade
per-base accuracy and breaks on native ONT reads, where the per-read error rate
swamps the handful of distinguishing positions and pulls allele-specific
coverage toward 0.5 (biasing *away* from real LOH).

This tool instead relies on *phasing*, which is where ONT is strong:

  1. The germline (normal) phased VCF defines which heterozygous SNVs are real
     and, for each, which allele lies on haplotype 1 vs haplotype 2 (the GT
     "0|1" / "1|0" and the PS phase-set tag from WhatsHap / longphase).
  2. In the tumour BAM we simply *count* reads supporting each allele at those
     same sites — no realignment, no mismatch model.
  3. Because phase gives a consistent orientation within a phase set, we
     aggregate haplotype-1 vs haplotype-2 support across many SNVs *coherently*.
     Under no LOH the aggregate sits at 0.5; under LOH one haplotype is
     systematically depleted. Aggregating dozens of sites averages out the
     per-site ONT substitution noise (central-limit), so the call rests on
     phasing accuracy and site count, not per-base accuracy.

Caveats (honest ones)
---------------------
* Het SNVs *inside* the most polymorphic HLA exons are unreliably mapped to the
  linear reference. We therefore (a) use every reliably-callable het SNV the
  germline caller phased across each gene window, including the better-behaved
  intronic/flanking sites that long reads phase into the same block, and
  (b) optionally use the *normal* BAM to drop sites whose germline allele
  balance is itself skewed (a mapping/CN artefact rather than a clean het).
* The H1/H2 labels are arbitrary *per phase set*; we never assume orientation is
  shared across phase sets. Multiple phase sets in one gene are tested
  separately and combined by Fisher's method, which avoids the folding bias you
  get from pooling "minor" counts across blocks.
* This measures allelic imbalance at the locus. With purity+ploidy we translate
  that into the most likely allele-specific copy-number state and an explicit
  LOH call, but a copy-neutral LOH and a single-copy loss both present as
  imbalance — the fitted CN state tells you which.

Outputs
-------
  <prefix>.hla_loh.tsv         per gene (and MHC-wide) summary + call
  <prefix>.hla_loh.sites.tsv   per-site detail (counts, BAFs, phase set)
  <prefix>.hla_loh.json        machine-readable, ready for the T/N report
  <prefix>.hla_loh.png         per-gene tumour vs normal BAF (if matplotlib)

Author: built for the DETERMINE ONT tumour-normal pipeline.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import sys
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

LOG = logging.getLogger("hla_loh")

# --------------------------------------------------------------------------- #
# HLA / MHC gene windows, GRCh38 (1-based, inclusive). Padded modestly so that
# long ONT reads pull flanking/intronic phased het SNVs into the same block.
# Override the whole set with --regions-bed if your build differs.
# --------------------------------------------------------------------------- #
HLA_GENES_GRCH38: Dict[str, Tuple[int, int]] = {
    "HLA-A":    (29_941_260, 29_945_884),
    "HLA-C":    (31_268_749, 31_272_130),
    "HLA-B":    (31_353_872, 31_357_188),
    "HLA-DRB1": (32_577_902, 32_589_848),
    "HLA-DQA1": (32_628_179, 32_647_062),
    "HLA-DQB1": (32_659_467, 32_668_383),
    "HLA-DPB1": (33_075_990, 33_089_696),
}
# MHC-wide window (classical + extended core) for a single high-powered call.
MHC_WINDOW_GRCH38: Tuple[int, int] = (29_900_000, 33_200_000)


@dataclass
class Config:
    normal_vcf: str
    tumor_bam: str
    out_prefix: str
    normal_bam: Optional[str] = None
    sample: Optional[str] = None          # normal sample name in the VCF
    purity: Optional[float] = None
    ploidy: Optional[float] = None
    flank: int = 1000
    min_base_qual: int = 7                 # ONT-appropriate (Q7 ~ default pass)
    min_mapq: int = 20
    min_site_depth: int = 10               # tumour depth to use a site
    min_normal_depth: int = 10
    min_sites: int = 8                     # min informative sites to call a gene
    max_normal_skew: float = 0.15          # drop het sites with |normalBAF-0.5|>this
    chrom: Optional[str] = None            # auto-detected if None
    regions_bed: Optional[str] = None
    mhc_only: bool = False


# --------------------------------------------------------------------------- #
# Small stats helpers (scipy if present, else pure-python fallbacks)
# --------------------------------------------------------------------------- #
def _binom_two_sided(k: int, n: int, p: float = 0.5) -> float:
    if n == 0:
        return 1.0
    try:
        from scipy.stats import binomtest
        return float(binomtest(k, n, p, alternative="two-sided").pvalue)
    except Exception:
        # normal approximation fallback
        mu = n * p
        sd = math.sqrt(n * p * (1 - p)) or 1e-9
        z = abs(k - mu) / sd
        # two-sided tail via erfc
        return math.erfc(z / math.sqrt(2))


def _fisher_combine(pvals: List[float]) -> float:
    pvals = [min(max(p, 1e-300), 1.0) for p in pvals if p is not None]
    if not pvals:
        return 1.0
    if len(pvals) == 1:
        return pvals[0]
    stat = -2.0 * sum(math.log(p) for p in pvals)
    df = 2 * len(pvals)
    try:
        from scipy.stats import chi2
        return float(chi2.sf(stat, df))
    except Exception:
        # Wilson-Hilferty approximation to the chi-square survival function
        x = stat / df
        z = (x ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
        return 0.5 * math.erfc(z / math.sqrt(2))


def expected_minor_fraction(purity: float, n_major: int, n_minor: int) -> float:
    """Expected observed minor-allele fraction for an allele-specific CN state
    (n_major, n_minor) in a tumour of given purity against a diploid 1:1 normal."""
    tum = purity * (n_major + n_minor)
    norm = (1 - purity) * 2
    minor = purity * n_minor + (1 - purity) * 1
    denom = tum + norm
    return minor / denom if denom > 0 else 0.5


def fit_cn_state(obs_minor: float, purity: float,
                 max_cn: int = 4) -> Tuple[Tuple[int, int], float]:
    """Pick the integer allele-specific CN state whose expected minor fraction is
    closest to the observed one. Returns ((n_major, n_minor), expected_frac)."""
    best, best_d, best_e = (1, 1), 1e9, 0.5
    for tot in range(1, max_cn + 1):
        for n_minor in range(0, tot // 2 + 1):
            n_major = tot - n_minor
            e = expected_minor_fraction(purity, n_major, n_minor)
            d = abs(e - obs_minor)
            if d < best_d:
                best, best_d, best_e = (n_major, n_minor), d, e
    return best, best_e


# --------------------------------------------------------------------------- #
# VCF parsing: phased heterozygous SNVs from the normal
# --------------------------------------------------------------------------- #
@dataclass
class HetSite:
    chrom: str
    pos: int            # 1-based
    ref: str
    alt: str
    h1: str             # base on haplotype 1
    h2: str             # base on haplotype 2
    ps: str             # phase set id (string; "." if none)


def _open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load_phased_hets(cfg: Config, chrom: str,
                     start: int, end: int) -> List[HetSite]:
    """Return phased het biallelic SNVs from the normal VCF within [start,end]."""
    sites: List[HetSite] = []
    sample_idx: Optional[int] = None
    with _open_text(cfg.normal_vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                samples = cols[9:]
                if cfg.sample:
                    if cfg.sample not in samples:
                        raise SystemExit(
                            f"Sample {cfg.sample!r} not in VCF; have {samples}")
                    sample_idx = samples.index(cfg.sample)
                else:
                    sample_idx = 0
                    if len(samples) > 1:
                        LOG.warning("VCF has %d samples; using first (%s). "
                                    "Set --sample to be explicit.",
                                    len(samples), samples[0])
                continue
            if sample_idx is None:
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            pos = int(f[1])
            if pos < start or pos > end:
                continue
            ref, alt = f[3], f[4]
            # biallelic SNV only
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            fmt = f[8].split(":")
            smp = f[9 + sample_idx].split(":")
            ff = dict(zip(fmt, smp))
            gt = ff.get("GT", "")
            if "|" not in gt:          # must be phased
                continue
            a, b = gt.split("|")[:2]
            if {a, b} != {"0", "1"}:   # must be het 0/1
                continue
            h1 = ref if a == "0" else alt
            h2 = ref if b == "0" else alt
            ps = ff.get("PS", ".") or "."
            sites.append(HetSite(chrom, pos, ref, alt, h1, h2, ps))
    return sites


# --------------------------------------------------------------------------- #
# Allele counting from a BAM at a single position
# --------------------------------------------------------------------------- #
def count_alleles(bam, chrom: str, pos: int, ref: str, alt: str,
                  min_bq: int, min_mapq: int) -> Tuple[int, int, int]:
    """Return (ref_count, alt_count, other_count) at 1-based pos."""
    p0 = pos - 1
    rc = ac = oc = 0
    for col in bam.pileup(chrom, p0, p0 + 1, truncate=True,
                          stepper="samtools", min_base_quality=min_bq,
                          max_depth=1_000_000):
        if col.reference_pos != p0:
            continue
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            if pr.alignment.mapping_quality < min_mapq:
                continue
            base = pr.alignment.query_sequence[pr.query_position].upper()
            if base == ref:
                rc += 1
            elif base == alt:
                ac += 1
            else:
                oc += 1
    return rc, ac, oc


# --------------------------------------------------------------------------- #
# Core per-region analysis
# --------------------------------------------------------------------------- #
@dataclass
class SiteResult:
    gene: str
    chrom: str
    pos: int
    ps: str
    h1: str
    h2: str
    t_h1: int
    t_h2: int
    t_other: int
    t_baf_h1: Optional[float]
    n_h1: Optional[int] = None
    n_h2: Optional[int] = None
    n_baf_h1: Optional[float] = None
    used: bool = True
    drop_reason: str = ""


@dataclass
class GeneResult:
    gene: str
    chrom: str
    start: int
    end: int
    n_sites_total: int = 0
    n_sites_used: int = 0
    tumor_h1: int = 0
    tumor_h2: int = 0
    informative_reads: int = 0
    obs_minor_fraction: Optional[float] = None
    median_site_minor: Optional[float] = None
    pvalue: Optional[float] = None
    n_phase_sets: int = 0
    expected_minor_if_loss: Optional[float] = None
    cn_state: Optional[str] = None
    cn_expected_minor: Optional[float] = None
    call: str = "uninformative"
    notes: str = ""


def analyse_region(cfg: Config, tbam, nbam, gene: str, chrom: str,
                   start: int, end: int) -> Tuple[GeneResult, List[SiteResult]]:
    hets = load_phased_hets(cfg, chrom, start, end)
    gr = GeneResult(gene=gene, chrom=chrom, start=start, end=end,
                    n_sites_total=len(hets))
    site_results: List[SiteResult] = []
    # per phase set aggregation: ps -> [h1_total, h2_total]
    ps_tot: Dict[str, List[int]] = {}
    per_site_minor: List[float] = []

    for h in hets:
        t_rc, t_ac, t_oc = count_alleles(tbam, h.chrom, h.pos, h.ref, h.alt,
                                         cfg.min_base_qual, cfg.min_mapq)
        # map ref/alt counts onto haplotype counts
        t_h1 = t_rc if h.h1 == h.ref else t_ac
        t_h2 = t_rc if h.h2 == h.ref else t_ac
        depth = t_h1 + t_h2
        sr = SiteResult(gene=gene, chrom=h.chrom, pos=h.pos, ps=h.ps,
                        h1=h.h1, h2=h.h2, t_h1=t_h1, t_h2=t_h2, t_other=t_oc,
                        t_baf_h1=(t_h1 / depth) if depth else None)

        # optional normal-BAM sanity filter
        if nbam is not None:
            n_rc, n_ac, n_oc = count_alleles(nbam, h.chrom, h.pos, h.ref, h.alt,
                                             cfg.min_base_qual, cfg.min_mapq)
            n_h1 = n_rc if h.h1 == h.ref else n_ac
            n_h2 = n_rc if h.h2 == h.ref else n_ac
            n_depth = n_h1 + n_h2
            sr.n_h1, sr.n_h2 = n_h1, n_h2
            sr.n_baf_h1 = (n_h1 / n_depth) if n_depth else None
            if n_depth < cfg.min_normal_depth:
                sr.used = False
                sr.drop_reason = "low_normal_depth"
            elif sr.n_baf_h1 is not None and abs(sr.n_baf_h1 - 0.5) > cfg.max_normal_skew:
                sr.used = False
                sr.drop_reason = "normal_allele_skew"

        if depth < cfg.min_site_depth:
            sr.used = False
            sr.drop_reason = sr.drop_reason or "low_tumor_depth"

        if sr.used:
            ps_tot.setdefault(h.ps, [0, 0])
            ps_tot[h.ps][0] += t_h1
            ps_tot[h.ps][1] += t_h2
            gr.tumor_h1 += t_h1
            gr.tumor_h2 += t_h2
            per_site_minor.append(min(t_h1, t_h2) / depth)
        site_results.append(sr)

    gr.n_sites_used = sum(1 for s in site_results if s.used)
    gr.informative_reads = gr.tumor_h1 + gr.tumor_h2
    gr.n_phase_sets = len(ps_tot)

    if gr.n_sites_used < cfg.min_sites or gr.informative_reads == 0:
        gr.call = "uninformative"
        gr.notes = (f"only {gr.n_sites_used} usable sites / "
                    f"{gr.informative_reads} reads")
        return gr, site_results

    # aggregate minor fraction: fold WITHIN each phase set (H1/H2 labels are only
    # consistent inside a phase set) then pool across phase sets — label-safe even
    # when a region spans several blocks (e.g. the MHC-wide window).
    minor_sum = sum(min(a, b) for a, b in ps_tot.values())
    total = sum(a + b for a, b in ps_tot.values())
    gr.obs_minor_fraction = minor_sum / total if total else None
    gr.median_site_minor = sorted(per_site_minor)[len(per_site_minor) // 2]

    # per-phase-set two-sided binomial test, combined by Fisher
    pvals = []
    n_ps_sig = 0
    for ps, (a, b) in ps_tot.items():
        n = a + b
        if n == 0:
            continue
        p_ps = _binom_two_sided(min(a, b), n, 0.5)
        pvals.append(p_ps)
        if p_ps < 0.01:
            n_ps_sig += 1
    gr.pvalue = _fisher_combine(pvals)

    # interpret with purity if available
    if cfg.purity:
        gr.expected_minor_if_loss = expected_minor_fraction(cfg.purity, 1, 0)
        (nM, nm), exp = fit_cn_state(gr.obs_minor_fraction, cfg.purity)
        gr.cn_state = f"{nM}:{nm}"
        gr.cn_expected_minor = exp
        loss_like = (nm == 0)
        # midpoint threshold (on the minor-fraction scale) between balanced and loss
        thr = (0.5 + gr.expected_minor_if_loss) / 2.0
        sig = gr.pvalue is not None and gr.pvalue < 0.01
        if sig and gr.obs_minor_fraction < thr and loss_like:
            gr.call = "LOH"
        elif sig and gr.obs_minor_fraction < thr:
            gr.call = "allelic_imbalance"   # imbalance but fitted state retains both
        else:
            gr.call = "retained"
        # multi-block region whose pooled fraction looks balanced but which
        # contains an individually-significant phase set: focal LOH is being
        # diluted — point the reader at the per-gene rows.
        if gr.n_phase_sets > 1 and n_ps_sig and gr.call in ("retained", "allelic_imbalance"):
            gr.notes = (f"{n_ps_sig}/{gr.n_phase_sets} phase sets individually "
                        f"imbalanced (p<0.01) — focal event diluted by pooling; "
                        f"see per-gene calls")
            if gr.call == "retained":
                gr.call = "heterogeneous"
    else:
        sig = gr.pvalue is not None and gr.pvalue < 0.01
        gr.call = "imbalance_no_purity" if sig else "retained"
        gr.notes = "supply --purity for an explicit LOH / CN call"

    return gr, site_results


# --------------------------------------------------------------------------- #
# Plot
# --------------------------------------------------------------------------- #
def make_plot(prefix: str, genes: List[GeneResult],
              sites: List[SiteResult]) -> Optional[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        LOG.warning("matplotlib unavailable; skipping plot")
        return None
    genes = [g for g in genes if g.gene != "MHC-wide" and g.n_sites_used]
    if not genes:
        return None
    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(genes)), 4.2))
    order = {g.gene: i for i, g in enumerate(genes)}
    for s in sites:
        if s.gene not in order or not s.used or s.t_baf_h1 is None:
            continue
        x = order[s.gene] + (hash((s.pos,)) % 100 - 50) / 320.0
        ax.scatter(x, s.t_baf_h1, s=10, alpha=0.45, color="#3a6ea5",
                   edgecolors="none", zorder=2)
    for g in genes:
        i = order[g.gene]
        if g.obs_minor_fraction is not None:
            ax.scatter(i, g.obs_minor_fraction, marker="_", s=900,
                       color="#c0392b", zorder=3, linewidths=2.5)
        if g.expected_minor_if_loss is not None:
            ax.scatter(i, g.expected_minor_if_loss, marker="x", s=40,
                       color="#7f8c8d", zorder=3)
    ax.axhline(0.5, color="#999", lw=1, ls="--", zorder=1)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([g.gene for g in genes], rotation=30, ha="right")
    ax.set_ylabel("Tumour haplotype-1 BAF (per het site)")
    ax.set_ylim(0, 1)
    ax.set_title("Phased allelic balance across HLA genes", fontsize=11, pad=22)
    ax.text(0.5, 1.04,
            "red — aggregate minor fraction   |   grey × — expected if single-copy loss",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
            color="#555")
    fig.tight_layout()
    out = f"{prefix}.hla_loh.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def detect_chrom(tbam, want: str = "6") -> str:
    refs = set(tbam.references)
    for cand in (f"chr{want}", want):
        if cand in refs:
            return cand
    raise SystemExit("Could not find chr6/6 in BAM header; use --chrom")


def load_regions(cfg: Config, chrom: str) -> List[Tuple[str, str, int, int]]:
    out: List[Tuple[str, str, int, int]] = []
    if cfg.regions_bed:
        with open(cfg.regions_bed) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                c, s, e, *rest = line.split()
                name = rest[0] if rest else f"{c}:{s}-{e}"
                out.append((name, c, int(s) + 1, int(e)))
        return out
    if not cfg.mhc_only:
        for g, (s, e) in HLA_GENES_GRCH38.items():
            out.append((g, chrom, s - cfg.flank, e + cfg.flank))
    out.append(("MHC-wide", chrom, MHC_WINDOW_GRCH38[0], MHC_WINDOW_GRCH38[1]))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Phased-BAF HLA/MHC LOH detection for ONT tumour-normal data")
    ap.add_argument("-N", "--normal-vcf", required=True,
                    help="phased germline VCF (WhatsHap/longphase) for the normal")
    ap.add_argument("-T", "--tumor-bam", required=True, help="tumour BAM (indexed)")
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--normal-bam", help="normal BAM (indexed) for het-site QC filter")
    ap.add_argument("--sample", help="normal sample name in the VCF (default: first)")
    ap.add_argument("--purity", type=float, help="tumour purity (e.g. ASCAT 0.91)")
    ap.add_argument("--ploidy", type=float, help="tumour ploidy (e.g. ASCAT 2.64)")
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--min-base-qual", type=int, default=7)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-site-depth", type=int, default=10)
    ap.add_argument("--min-normal-depth", type=int, default=10)
    ap.add_argument("--min-sites", type=int, default=8)
    ap.add_argument("--max-normal-skew", type=float, default=0.15)
    ap.add_argument("--chrom", help="override chr6 contig name")
    ap.add_argument("--regions-bed", help="custom regions (BED) instead of built-in HLA genes")
    ap.add_argument("--mhc-only", action="store_true",
                    help="skip per-gene windows, do the single MHC-wide call only")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)s %(message)s")

    cfg = Config(
        normal_vcf=args.normal_vcf, tumor_bam=args.tumor_bam,
        out_prefix=args.out_prefix, normal_bam=args.normal_bam,
        sample=args.sample, purity=args.purity, ploidy=args.ploidy,
        flank=args.flank, min_base_qual=args.min_base_qual,
        min_mapq=args.min_mapq, min_site_depth=args.min_site_depth,
        min_normal_depth=args.min_normal_depth, min_sites=args.min_sites,
        max_normal_skew=args.max_normal_skew, chrom=args.chrom,
        regions_bed=args.regions_bed, mhc_only=args.mhc_only)

    import pysam
    tbam = pysam.AlignmentFile(cfg.tumor_bam, "rb")
    nbam = pysam.AlignmentFile(cfg.normal_bam, "rb") if cfg.normal_bam else None

    chrom = cfg.chrom or detect_chrom(tbam)
    LOG.info("Using contig %s", chrom)
    regions = load_regions(cfg, chrom)

    genes: List[GeneResult] = []
    all_sites: List[SiteResult] = []
    for name, c, s, e in regions:
        LOG.info("Analysing %s %s:%d-%d", name, c, s, e)
        gr, srs = analyse_region(cfg, tbam, nbam, name, c, s, e)
        genes.append(gr)
        all_sites.extend(srs)
        LOG.info("  %s: %d/%d sites used, call=%s minorFrac=%s p=%s",
                 name, gr.n_sites_used, gr.n_sites_total, gr.call,
                 f"{gr.obs_minor_fraction:.3f}" if gr.obs_minor_fraction else "NA",
                 f"{gr.pvalue:.2e}" if gr.pvalue else "NA")

    # ---- write outputs ----
    tsv = f"{cfg.out_prefix}.hla_loh.tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join([
            "gene", "chrom", "start", "end", "n_sites_total", "n_sites_used",
            "n_phase_sets", "tumor_h1", "tumor_h2", "obs_minor_fraction",
            "median_site_minor", "expected_minor_if_loss", "cn_state",
            "cn_expected_minor", "pvalue", "call", "notes"]) + "\n")
        for g in genes:
            fh.write("\t".join(str(x) if x is not None else "NA" for x in [
                g.gene, g.chrom, g.start, g.end, g.n_sites_total, g.n_sites_used,
                g.n_phase_sets, g.tumor_h1, g.tumor_h2,
                _r(g.obs_minor_fraction), _r(g.median_site_minor),
                _r(g.expected_minor_if_loss), g.cn_state,
                _r(g.cn_expected_minor),
                f"{g.pvalue:.3e}" if g.pvalue is not None else "NA",
                g.call, g.notes]) + "\n")

    sites_tsv = f"{cfg.out_prefix}.hla_loh.sites.tsv"
    with open(sites_tsv, "w") as fh:
        fh.write("\t".join([
            "gene", "chrom", "pos", "phase_set", "hap1_base", "hap2_base",
            "tumor_hap1", "tumor_hap2", "tumor_other", "tumor_baf_hap1",
            "normal_hap1", "normal_hap2", "normal_baf_hap1",
            "used", "drop_reason"]) + "\n")
        for s in all_sites:
            fh.write("\t".join(str(x) if x is not None else "NA" for x in [
                s.gene, s.chrom, s.pos, s.ps, s.h1, s.h2,
                s.t_h1, s.t_h2, s.t_other, _r(s.t_baf_h1),
                s.n_h1, s.n_h2, _r(s.n_baf_h1),
                int(s.used), s.drop_reason]) + "\n")

    png = make_plot(cfg.out_prefix, genes, all_sites)

    js = {
        "tool": "hla_loh_phased_baf",
        "method": "phased B-allele-frequency (normal-phased het SNVs, tumour allele counts)",
        "purity": cfg.purity, "ploidy": cfg.ploidy, "chrom": chrom,
        "genes": [asdict(g) for g in genes],
        "params": {k: getattr(cfg, k) for k in (
            "min_base_qual", "min_mapq", "min_site_depth", "min_sites",
            "max_normal_skew", "flank")},
    }
    json_path = f"{cfg.out_prefix}.hla_loh.json"
    with open(json_path, "w") as fh:
        json.dump(js, fh, indent=2)

    # ---- console summary ----
    print("\nHLA / MHC phased-BAF LOH summary")
    print(f"  purity={cfg.purity}  ploidy={cfg.ploidy}  contig={chrom}")
    if cfg.purity:
        print(f"  expected minor fraction if single-copy loss: "
              f"{expected_minor_fraction(cfg.purity,1,0):.3f}")
    print(f"  {'gene':10s} {'sites':>7s} {'minorFrac':>9s} {'CN':>5s} "
          f"{'p':>10s}  call")
    for g in genes:
        print(f"  {g.gene:10s} {g.n_sites_used:>3d}/{g.n_sites_total:<3d} "
              f"{(f'{g.obs_minor_fraction:.3f}' if g.obs_minor_fraction is not None else 'NA'):>9s} "
              f"{(g.cn_state or 'NA'):>5s} "
              f"{(f'{g.pvalue:.1e}' if g.pvalue is not None else 'NA'):>10s}  {g.call}")
    print(f"\n  TSV : {tsv}\n  sites: {sites_tsv}\n  JSON: {json_path}"
          + (f"\n  plot: {png}" if png else ""))
    return 0


def _r(x, nd=4):
    return round(x, nd) if isinstance(x, float) else x


if __name__ == "__main__":
    sys.exit(main())
