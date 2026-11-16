#!/usr/bin/env python3
"""
run_sigprofiler_assignment.py
=============================

Refit COSMIC mutational signatures to a DETERMINE somatic call set with
SigProfilerAssignment. Two complementary read-outs in one run:

  * SBS / ID / DBS signatures  from the ClairS VCF        (--vcf)
  * Copy-number (CNV48) signatures from ASCAT segments    (--ascat-segments)

SigProfilerAssignment is a *Python* tool (AlexandrovLab). For SNV input it uses
SigProfilerMatrixGenerator to build the 96-context matrix (needs the reference
genome installed, ~3 GB one-time). For copy-number input it reads the ASCAT
allele-specific segments directly (no genome needed) and refits the COSMIC CN
signature catalogue (Steele et al. 2022).

At least one of --vcf / --ascat-segments must be given; supply both to get the
SBS and CN read-outs side by side.

Install (Python, into the env you run this from):
    pip install SigProfilerAssignment    # pulls SigProfilerMatrixGenerator + SigProfilerPlotting
    # one-time reference genome, only needed for --vcf:
    python -c "from SigProfilerMatrixGenerator import install as g; g.install('GRCh38')"

Outputs (under --output):
    sbs/   full SigProfilerAssignment output for the SBS refit (if --vcf)
    cnv/   full output for the CN refit (if --ascat-segments)
    signature_summary.json   combined SBS + CN summary
    sbs_signatures.tsv / cnv_signatures.tsv   per-signature counts + proportions

Usage:
    python run_sigprofiler_assignment.py \
        --vcf ClairS/snv.vcf.gz --ascat-segments T9.segments.txt \
        --sample-id PATIENT05_T9 --genome GRCh38 --output sigprofiler_out -v
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import sys
from typing import Dict, List, Optional, Tuple

LOG = logging.getLogger("sigprofiler_assignment")

_VALID_CONTEXTS = ("96", "288", "1536", "DINUC", "ID")
_CONTEXT_LABEL = {"96": "sbs96", "288": "sbs288", "1536": "sbs1536",
                  "DINUC": "dbs78", "ID": "id83"}


# --------------------------------------------------------------------------- #
# Input preparation
# --------------------------------------------------------------------------- #
def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _classify(ref: str, alt: str) -> str:
    if "," in alt:
        return "multiallelic"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    return "indel"


def prepare_vcf(in_vcfs: List[str], out_vcf: str, sample_id: str,
                include_indels: bool, strip_chr: bool) -> Tuple[int, int, int]:
    """Write a clean single-sample VCF (PASS SNVs, optionally indels) that
    SigProfilerMatrixGenerator parses reliably. Returns (n_snv, n_indel, n_skip)."""
    n_snv = n_indel = n_skip = 0
    with open(out_vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + sample_id + "\n")
        for vcf in in_vcfs:
            with _open(vcf) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 8:
                        continue
                    chrom, pos, ref, alt, filt = f[0], f[1], f[3], f[4], f[6]
                    if filt not in ("PASS", "."):
                        n_skip += 1
                        continue
                    kind = _classify(ref, alt)
                    if kind == "multiallelic":
                        n_skip += 1
                        continue
                    if kind == "indel" and not include_indels:
                        n_skip += 1
                        continue
                    if strip_chr and chrom.lower().startswith("chr"):
                        chrom = chrom[3:]
                    out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t"
                              f"GT\t0/1\n")
                    if kind == "snv":
                        n_snv += 1
                    else:
                        n_indel += 1
    return n_snv, n_indel, n_skip


def prepare_ascat_seg(in_seg: str, out_seg: str, sample_id: str) -> int:
    """Normalise an ASCAT *.segments.txt to the exact columns the
    SigProfilerMatrixGenerator 'seg:ASCAT' parser reads (sample, chr, startpos,
    endpos, nMajor, nMinor). Returns the number of segments written.

    The CNV48 classifier uses only sample / nMajor / nMinor / startpos / endpos,
    so chromosome naming (chr-prefixed or not) is irrelevant here."""
    with _open(in_seg) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        low = [h.lower() for h in header]

        def col(*names, required=True):
            for n in names:
                if n in low:
                    return low.index(n)
            if required:
                raise SystemExit(f"ASCAT segments missing column {names}; have {header}")
            return None

        i_samp = col("sample", "sampleid", "sample_id", required=False)
        i_chr = col("chr", "chromosome", required=False)
        i_s = col("startpos", "start", "start_position")
        i_e = col("endpos", "end", "end_position")
        i_maj = col("nmajor", "nmaj", "major_cn", "a_cn")
        i_min = col("nminor", "nmin", "minor_cn", "b_cn")

        n = 0
        with open(out_seg, "w") as out:
            out.write("sample\tchr\tstartpos\tendpos\tnMajor\tnMinor\n")
            for line in fh:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                try:
                    maj, mn = int(float(f[i_maj])), int(float(f[i_min]))
                    s, e = int(float(f[i_s])), int(float(f[i_e]))
                except (ValueError, IndexError):
                    continue
                samp = f[i_samp] if i_samp is not None else sample_id
                chrom = f[i_chr] if i_chr is not None else "NA"
                out.write(f"{samp}\t{chrom}\t{s}\t{e}\t{maj}\t{mn}\n")
                n += 1
    return n


# --------------------------------------------------------------------------- #
# Result parsing (identical layout for SBS and CN runs)
# --------------------------------------------------------------------------- #
def parse_activities(output_dir: str, sample_id: str) -> Optional[Dict]:
    act = os.path.join(output_dir, "Assignment_Solution", "Activities",
                       "Assignment_Solution_Activities.txt")
    if not os.path.exists(act):
        LOG.warning("Activities file not found at %s", act)
        return None
    rows = list(csv.DictReader(open(act), delimiter="\t"))
    row = next((r for r in rows if r.get("Samples") == sample_id), None)
    if row is None and rows:
        row = rows[0]
    if row is None:
        return None
    sigs = {k: int(float(v)) for k, v in row.items()
            if k != "Samples" and float(v) > 0}
    total = sum(sigs.values()) or 1
    ranked = sorted(sigs.items(), key=lambda kv: -kv[1])
    return {"sample": row.get("Samples"), "total_assigned": total,
            "signatures": [{"signature": s, "count": c,
                            "proportion": round(c / total, 4)} for s, c in ranked]}


def parse_recon_stats(output_dir: str, sample_id: str) -> Optional[float]:
    stats = os.path.join(output_dir, "Assignment_Solution", "Solution_Stats",
                         "Assignment_Solution_Samples_Stats.txt")
    if not os.path.exists(stats):
        return None
    try:
        rows = list(csv.DictReader(open(stats), delimiter="\t"))
        key = next((k for k in rows[0] if "cosine" in k.lower()), None)
        if not key:
            return None
        r = next((r for r in rows if r.get(list(r)[0]) == sample_id), rows[0])
        return float(r[key])
    except Exception:
        return None


def _write_sig_tsv(path: str, summary: Dict):
    with open(path, "w") as fh:
        fh.write("signature\tcount\tproportion\n")
        for s in summary["signatures"]:
            fh.write(f"{s['signature']}\t{s['count']}\t{s['proportion']}\n")


def _print_summary(label: str, summary: Dict, cosine: Optional[float]):
    cos = f"{cosine:.4f}" if cosine is not None else "n/a"
    print(f"\n  [{label}]  {summary['total_assigned']} events; reconstruction cosine = {cos}")
    print(f"  {'signature':>10s} {'count':>7s} {'prop':>7s}")
    for s in summary["signatures"]:
        print(f"  {s['signature']:>10s} {s['count']:>7d} {s['proportion']:>7.1%}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Refit COSMIC SBS (VCF) and CN (ASCAT) signatures with SigProfilerAssignment")
    ap.add_argument("-V", "--vcf", help="ClairS somatic SNV VCF (SBS/ID/DBS signatures)")
    ap.add_argument("-A", "--ascat-segments",
                    help="ASCAT allele-specific segments (CNV48 signatures)")
    ap.add_argument("-o", "--output", required=True, help="output directory")
    ap.add_argument("--sample-id", default="PATIENT05_T9")
    ap.add_argument("--genome", default="GRCh38",
                    choices=["GRCh37", "GRCh38", "mm9", "mm10", "rn6"])
    ap.add_argument("--cosmic-version", type=float, default=3.5)
    ap.add_argument("--sbs-contexts", default="96",
                    help="comma-separated mutation contexts from the VCF: any of "
                         "96,288,1536 (SBS), DINUC (DBS78), ID (indels/ID83). "
                         "e.g. --sbs-contexts 96,ID,DINUC")
    ap.add_argument("--indel-vcf", help="optional ClairS indel VCF to include")
    ap.add_argument("--include-indels", action="store_true")
    ap.add_argument("--per-mutation", action="store_true",
                    help="export per-mutation SBS signature probabilities")
    ap.add_argument("--exclude-artifacts", action="store_true",
                    help="exclude possible-sequencing-artifact COSMIC subgroups (SBS)")
    ap.add_argument("--strip-chr", action="store_true",
                    help="strip 'chr' when writing the prepared VCF")
    ap.add_argument("--install-genome", action="store_true",
                    help="install the reference genome first (one-time ~3GB; SBS only)")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)s %(message)s")

    if not args.vcf and not args.ascat_segments:
        raise SystemExit("Provide at least one of --vcf / --ascat-segments.")

    os.makedirs(args.output, exist_ok=True)
    from SigProfilerAssignment import Analyzer as Analyze

    combined: Dict = {"sample": args.sample_id, "genome": args.genome,
                      "cosmic_version": args.cosmic_version}

    # ================== SBS / ID / DBS from the VCF (per context) ========== #
    if args.vcf:
        contexts = [c.strip() for c in args.sbs_contexts.split(",") if c.strip()]
        bad = [c for c in contexts if c not in _VALID_CONTEXTS]
        if bad:
            raise SystemExit(f"invalid context(s) {bad}; choose from {_VALID_CONTEXTS}")
        want_indels = args.include_indels or ("ID" in contexts)

        # one prepared VCF (SNVs always; indels too if any ID context requested)
        # is shared by every context run.
        shared_in = os.path.join(args.output, "spa_vcf_input")
        os.makedirs(shared_in, exist_ok=True)
        in_vcfs = [args.vcf] + ([args.indel_vcf] if args.indel_vcf else [])
        prepared = os.path.join(shared_in, f"{args.sample_id}.vcf")
        n_snv, n_indel, n_skip = prepare_vcf(in_vcfs, prepared, args.sample_id,
                                             want_indels, args.strip_chr)
        LOG.info("VCF prepared: %d SNV, %d indel, %d skipped", n_snv, n_indel, n_skip)
        if n_snv + n_indel == 0:
            raise SystemExit("No PASS variants written — check FILTER values.")
        if "ID" in contexts and n_indel == 0:
            LOG.warning("ID context requested but no indels were written — "
                        "pass --indel-vcf or check the VCF.")
        if args.install_genome:
            LOG.info("Installing reference genome %s (one-time)...", args.genome)
            from SigProfilerMatrixGenerator import install as genInstall
            genInstall.install(args.genome)

        for ctx in contexts:
            label = _CONTEXT_LABEL[ctx]
            ctx_out = os.path.join(args.output, label)
            os.makedirs(ctx_out, exist_ok=True)
            LOG.info("cosmic_fit %s (context=%s, genome=%s, COSMIC v%s)",
                     label, ctx, args.genome, args.cosmic_version)
            Analyze.cosmic_fit(
                samples=shared_in, output=ctx_out, input_type="vcf",
                context_type=ctx, collapse_to_SBS96=(ctx == "96"),
                cosmic_version=args.cosmic_version, genome_build=args.genome,
                make_plots=not args.no_plots,
                sample_reconstruction_plots=("pdf" if not args.no_plots else False),
                export_probabilities=True,
                export_probabilities_per_mutation=args.per_mutation,
                exclude_signature_subgroups=(["Sequencing_Artifact_signatures"]
                                             if args.exclude_artifacts else None),
                verbose=args.verbose)
            res = parse_activities(ctx_out, args.sample_id)
            if res:
                res["reconstruction_cosine"] = parse_recon_stats(ctx_out, args.sample_id)
                res["context"] = ctx
                combined[label] = res
                _write_sig_tsv(os.path.join(args.output, f"{label}_signatures.tsv"), res)

    # ===================== Copy-number (CNV48) from ASCAT ================== #
    if args.ascat_segments:
        cnv_out = os.path.join(args.output, "cnv")
        os.makedirs(cnv_out, exist_ok=True)
        seg_clean = os.path.join(cnv_out, f"{args.sample_id}.ascat.seg.txt")
        n_seg = prepare_ascat_seg(args.ascat_segments, seg_clean, args.sample_id)
        LOG.info("ASCAT segments prepared: %d segments", n_seg)
        if n_seg < 20:
            LOG.warning("Only %d segments — CN signature fit may be unstable.", n_seg)
        LOG.info("cosmic_fit CN (CNV48, seg:ASCAT, COSMIC v%s)", args.cosmic_version)
        Analyze.cosmic_fit(
            samples=seg_clean, output=cnv_out, input_type="seg:ASCAT",
            context_type="CNV48", collapse_to_SBS96=False,
            cosmic_version=args.cosmic_version, genome_build=args.genome,
            make_plots=not args.no_plots, verbose=args.verbose)
        cnv = parse_activities(cnv_out, args.sample_id)
        if cnv:
            cnv["reconstruction_cosine"] = parse_recon_stats(cnv_out, args.sample_id)
            cnv["n_segments"] = n_seg
            combined["cnv"] = cnv
            _write_sig_tsv(os.path.join(args.output, "cnv_signatures.tsv"), cnv)

    # =============================== summary ============================== #
    with open(os.path.join(args.output, "signature_summary.json"), "w") as fh:
        json.dump(combined, fh, indent=2)

    print(f"\nSigProfilerAssignment — {args.sample_id} "
          f"(COSMIC v{args.cosmic_version}, {args.genome})")
    _LABEL_TITLE = {"sbs96": "SBS96", "sbs288": "SBS288", "sbs1536": "SBS1536",
                    "dbs78": "DBS78 (doublets)", "id83": "ID83 (indels)",
                    "cnv": "CNV48 (copy number)"}
    for label in ("sbs96", "sbs288", "sbs1536", "dbs78", "id83", "cnv"):
        if label in combined:
            _print_summary(_LABEL_TITLE[label], combined[label],
                           combined[label].get("reconstruction_cosine"))
    print(f"\n  combined summary : {os.path.join(args.output, 'signature_summary.json')}")
    for label in ("sbs96", "sbs288", "sbs1536", "dbs78", "id83", "cnv"):
        if label in combined:
            print(f"  {label:>8s} output : {os.path.join(args.output, label, 'Assignment_Solution')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
