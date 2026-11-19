#!/usr/bin/env python3
"""
run_sv_signatures.py
====================

Structural-variant (rearrangement) signature analysis for the DETERMINE ONT
pipeline, from SAVANA / Severus SV calls.

WHY THIS IS SEPARATE FROM run_sigprofiler_assignment.py
-------------------------------------------------------
SigProfilerAssignment's cosmic_fit refits the COSMIC reference catalogues, which
cover SBS, DBS, ID and CN signatures only — there are NO COSMIC reference SV
signatures. SV signatures therefore take a different route:

  1. SigProfilerMatrixGenerator.generateSVMatrix : BEDPE -> SV32 matrix
     (clustered/non-clustered x {del, tds, inv, trans} x size bins). No genome
     download needed; classification is coordinate/strand based.
  2. (optional) refit that SV32 matrix against a REFERENCE rearrangement
     signature set you supply via --sv-signatures, using cosmic_fit with a
     custom signature_database. The usual reference is the RefSig / Signal
     rearrangement signatures (Degasperi et al.; signal.mutographs.org) — these
     are not bundled with SigProfiler and are not COSMIC, so you provide the tsv.

Without --sv-signatures the script stops at the SV32 matrix, which is itself the
input for de novo extraction (SigProfilerExtractor) or a later refit.

ONT CAVEAT
----------
Reference SV signatures were derived from Illumina short-read cohorts. ONT SV
calls (SAVANA/Severus) differ in breakpoint precision and in sensitivity within
repeats/segmental duplications, so the SV32 spectrum — and any refit against
short-read reference signatures — should be read as indicative, not definitive.
Breakpoint-coordinate accuracy mostly affects the size-bin assignment, not the
SV-type channel.

Input BEDPE columns required by generateSVMatrix:
    sample, chrom1, start1, end1, chrom2, start2, end2,
    and EITHER svclass OR (strand1 AND strand2)   [svclass computed from strands]

Usage:
    python run_sv_signatures.py --bedpe savana.classified.somatic.bedpe \
        --sample-id PATIENT05_T9 --output sv_out
    # with a refit:
    python run_sv_signatures.py --bedpe savana....bedpe --sample-id PATIENT05_T9 \
        --output sv_out --sv-signatures RefSigv2_Rearr.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import shutil
import sys
from typing import Dict, List, Optional

LOG = logging.getLogger("sv_signatures")

# map common caller SV-type vocabularies -> SigProfiler full names
_SVCLASS_NORM = {
    "del": "deletion", "deletion": "deletion", "dele": "deletion",
    "dup": "tandem-duplication", "tandem-duplication": "tandem-duplication",
    "tds": "tandem-duplication", "duplication": "tandem-duplication",
    "inv": "inversion", "inversion": "inversion",
    "h2hinv": "inversion", "t2tinv": "inversion",
    "tra": "translocation", "trans": "translocation",
    "translocation": "translocation", "bnd": "translocation", "ctx": "translocation",
}

# SAVANA / standard-BEDPE breakpoint orientation -> SV class. NB this is the
# bedpe convention (+- del, -+ dup, ++/-- inv), which is NOT SigProfiler's
# Sanger strand convention — so we translate to svclass and never pass these
# through as strand1/strand2.
_ORIENT_TO_SVCLASS = {
    "+-": "deletion", "-+": "tandem-duplication",
    "++": "inversion", "--": "inversion",
}


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _find(cols_low: List[str], *names) -> Optional[int]:
    for n in names:
        if n in cols_low:
            return cols_low.index(n)
    return None


def _orient_from_name(name: str) -> Optional[str]:
    """Pull the trailing breakpoint-orientation token (++ +- -+ --) out of a
    SAVANA name field like 'ID_42030|0bp|TUMOUR_6|--'."""
    for tok in reversed(name.split("|")):
        tok = tok.strip()
        if tok in _ORIENT_TO_SVCLASS:
            return tok
    return None


def convert_bedpe(in_bedpe: str, out_bedpe: str, sample_id: str) -> int:
    """Read a SAVANA/Severus/generic BEDPE and write the column set
    generateSVMatrix needs (svclass-only; see note below). Returns count written.

    SAVANA's classified.somatic.bedpe is HEADERLESS: columns are
    chrom1 start1 end1 chrom2 start2 end2 name, where `name` ends in the
    breakpoint orientation (++ +- -+ --). SAVANA/standard-BEDPE orientation maps
    +- -> deletion, -+ -> tandem-duplication, ++/-- -> inversion, and any
    inter-chromosomal pair -> translocation.

    IMPORTANT: SigProfiler's generateSVMatrix derives the SV class from
    strand1/strand2 *whenever present*, using a DIFFERENT (Sanger) convention,
    and ignores an svclass column. So we always emit svclass ONLY (no strands),
    translating each caller's orientation/type into svclass ourselves — this
    keeps the classification correct regardless of caller convention."""
    rows_out = []
    with _open(in_bedpe) as fh:
        first = None
        for line in fh:
            if line.startswith("##") or not line.strip():
                continue
            first = line.rstrip("\n").split("\t")
            break
        if first is None:
            raise SystemExit("empty BEDPE")

        # headerless? -> SAVANA-style: cols 1 and 4 are integer start positions
        headerless = len(first) >= 6 and first[1].lstrip("-").isdigit() \
            and first[4].lstrip("-").isdigit()

        def emit(c1, s1, e1, c2, s2, e2, svclass):
            rows_out.append((sample_id if True else sample_id,
                             c1, s1, e1, c2, s2, e2, svclass))

        def handle_positional(f):
            try:
                c1, s1, e1 = f[0], int(float(f[1])), int(float(f[2]))
                c2, s2, e2 = f[3], int(float(f[4])), int(float(f[5]))
            except (ValueError, IndexError):
                return
            if c1 != c2:
                svclass = "translocation"
            else:
                name = f[6] if len(f) > 6 else ""
                orient = _orient_from_name(name)
                svclass = _ORIENT_TO_SVCLASS.get(orient, "") if orient else ""
                if not svclass:
                    return                      # insertion / unparseable -> skip
            emit(c1, s1, e1, c2, s2, e2, svclass)

        if headerless:
            handle_positional(first)
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                handle_positional(line.rstrip("\n").split("\t"))
        else:
            # header present: fuzzy column lookup (generic / Severus-style BEDPE)
            low = [h.strip().lower() for h in [c.lstrip("#") for c in first]]
            i_c1 = _find(low, "chrom1", "chr1", "chr_1")
            i_s1 = _find(low, "start1", "pos1", "start_1")
            i_e1 = _find(low, "end1", "end_1")
            i_c2 = _find(low, "chrom2", "chr2", "chr_2")
            i_s2 = _find(low, "start2", "pos2", "start_2")
            i_e2 = _find(low, "end2", "end_2")
            i_cls = _find(low, "svclass", "svtype", "type", "sv_type", "class")
            i_orient = _find(low, "orientation", "bp_notation", "strands", "notation")
            i_samp = _find(low, "sample", "sample_id", "sampleid")
            if None in (i_c1, i_s1, i_c2, i_s2):
                raise SystemExit(f"BEDPE missing breakpoint columns; header={first}")
            if i_cls is None and i_orient is None:
                raise SystemExit("BEDPE has no svclass/svtype or orientation/strands "
                                 "column to classify rearrangements")
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                try:
                    c1, s1 = f[i_c1], int(float(f[i_s1]))
                    c2, s2 = f[i_c2], int(float(f[i_s2]))
                except (ValueError, IndexError):
                    continue
                e1 = int(float(f[i_e1])) if i_e1 is not None and f[i_e1] not in ("", ".") else s1 + 1
                e2 = int(float(f[i_e2])) if i_e2 is not None and f[i_e2] not in ("", ".") else s2 + 1
                if c1 != c2:
                    svclass = "translocation"
                elif i_cls is not None and i_cls < len(f):
                    svclass = _SVCLASS_NORM.get(f[i_cls].strip().lower(), "")
                elif i_orient is not None and i_orient < len(f):
                    svclass = _ORIENT_TO_SVCLASS.get(f[i_orient].strip(), "")
                else:
                    svclass = ""
                if not svclass:
                    continue
                samp = f[i_samp] if i_samp is not None and i_samp < len(f) else sample_id
                rows_out.append((samp, c1, s1, e1, c2, s2, e2, svclass))

    with open(out_bedpe, "w") as out:
        out.write("sample\tchrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsvclass\n")
        for r in rows_out:
            out.write("\t".join(map(str, r)) + "\n")
    return len(rows_out)


def parse_activities(output_dir: str, sample_id: str) -> Optional[Dict]:
    act = os.path.join(output_dir, "Assignment_Solution", "Activities",
                       "Assignment_Solution_Activities.txt")
    if not os.path.exists(act):
        return None
    rows = list(csv.DictReader(open(act), delimiter="\t"))
    row = next((r for r in rows if r.get("Samples") == sample_id), rows[0] if rows else None)
    if row is None:
        return None
    sigs = {k: int(float(v)) for k, v in row.items()
            if k != "Samples" and float(v) > 0}
    total = sum(sigs.values()) or 1
    ranked = sorted(sigs.items(), key=lambda kv: -kv[1])
    return {"sample": row.get("Samples"), "total_assigned": total,
            "signatures": [{"signature": s, "count": c,
                            "proportion": round(c / total, 4)} for s, c in ranked]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Structural-variant (SV32) signatures from SAVANA/Severus calls")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--bedpe", help="SAVANA/Severus/generic BEDPE of somatic SVs")
    src.add_argument("--sv-vcf", help="SV VCF (converted to BEDPE by SigProfiler)")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--sample-id", default="PATIENT05_T9")
    ap.add_argument("--project", default="DETERMINE_SV")
    ap.add_argument("--sv-signatures",
                    help="reference rearrangement signature tsv (e.g. RefSig/Signal) "
                         "to refit the SV32 matrix against; omit to stop at the matrix")
    ap.add_argument("--cosmic-version", type=float, default=3.5)
    ap.add_argument("--genome", default="GRCh38")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)s %(message)s")
    os.makedirs(args.output, exist_ok=True)

    # ---- 1. assemble an input dir with one BEDPE (or VCF) for the sample ----
    sv_in = os.path.join(args.output, "sv_input")
    os.makedirs(sv_in, exist_ok=True)
    if args.bedpe:
        out_bedpe = os.path.join(sv_in, f"{args.sample_id}.bedpe")
        n = convert_bedpe(args.bedpe, out_bedpe, args.sample_id)
        LOG.info("Converted %d rearrangements -> %s", n, out_bedpe)
        if n < 10:
            LOG.warning("Only %d SVs — the SV32 spectrum will be sparse.", n)
    else:
        # let SigProfiler convert the VCF; it keys the sample off the filename
        dst = os.path.join(sv_in, f"{args.sample_id}.vcf")
        shutil.copy(args.sv_vcf, dst)
        LOG.info("Copied SV VCF -> %s (SigProfiler will convert to BEDPE)", dst)

    # ---- 2. generate the SV32 matrix ----
    # SigProfilerMatrixGenerator's SV generator still calls np.in1d, which was
    # removed in numpy 2.0; np.isin is the drop-in replacement. Shim it so the
    # matrix step works regardless of the installed numpy version.
    import numpy as _np
    if not hasattr(_np, "in1d"):
        _np.in1d = _np.isin
    from SigProfilerMatrixGenerator.scripts import SVMatrixGenerator as svmg
    mat_dir = os.path.join(args.output, "sv_matrix")
    LOG.info("Generating SV32 matrix (project=%s)", args.project)
    svmg.generateSVMatrix(sv_in, args.project, mat_dir)
    sv32 = os.path.join(mat_dir, f"{args.project}.SV32.matrix.tsv")
    if not os.path.exists(sv32):
        raise SystemExit(f"SV32 matrix not produced at {sv32}")
    LOG.info("SV32 matrix: %s", sv32)

    result: Dict = {"tool": "SV32 / SigProfilerMatrixGenerator",
                    "sample": args.sample_id, "sv32_matrix": sv32}

    # ---- 3. optional refit against a reference SV signature set ----
    if args.sv_signatures:
        from SigProfilerAssignment import Analyzer as Analyze
        refit_out = os.path.join(args.output, "sv_refit")
        LOG.info("Refitting SV32 against %s", args.sv_signatures)
        Analyze.cosmic_fit(
            samples=sv32, output=refit_out, input_type="matrix",
            signatures=None, signature_database=args.sv_signatures,
            genome_build=args.genome, cosmic_version=args.cosmic_version,
            collapse_to_SBS96=False, make_plots=not args.no_plots,
            verbose=args.verbose)
        summ = parse_activities(refit_out, args.sample_id)
        if summ:
            result["refit"] = summ
            with open(os.path.join(args.output, "sv_signatures.tsv"), "w") as fh:
                fh.write("signature\tcount\tproportion\n")
                for s in summ["signatures"]:
                    fh.write(f"{s['signature']}\t{s['count']}\t{s['proportion']}\n")
    else:
        LOG.info("No --sv-signatures given; stopping at the SV32 matrix.")

    with open(os.path.join(args.output, "sv_summary.json"), "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"\nSV signature analysis — {args.sample_id}")
    print(f"  SV32 matrix : {sv32}")
    if "refit" in result:
        r = result["refit"]
        print(f"  refit against {os.path.basename(args.sv_signatures)} "
              f"({r['total_assigned']} SVs):")
        for s in r["signatures"]:
            print(f"    {s['signature']:>10s} {s['count']:>5d} {s['proportion']:>7.1%}")
    else:
        print("  (no reference SV signatures supplied — matrix only; feed it to "
              "SigProfilerExtractor for de novo, or re-run with --sv-signatures)")
    print(f"  summary JSON: {os.path.join(args.output, 'sv_summary.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
