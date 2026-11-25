#!/usr/bin/env python3
"""
Build a PyClone-VI input table from the ClairS/DeepSomatic consensus VCF and ASCAT output.

Every somatic SNV is matched to the ASCAT segment it sits on and written as one row:

    mutation_id  sample_id  ref_counts  alt_counts  normal_cn  major_cn  minor_cn  tumour_content  chrom

PyClone-VI ignores chrom. PhyClone reads the same table and uses chrom to find clusters that
may have been lost.

Inputs
    --vcf           Consensus VCF (*.consensus.final.vcf.gz). It has one sample, the tumour,
                    whose FORMAT/AD holds the ref,alt read counts.
    --segments      ASCAT segments.txt with the columns: sample chr startpos endpos nMajor nMinor
    --purityploidy  ASCAT purityploidy.txt with the columns: AberrantCellFraction Ploidy

Only autosomes are kept, because the normal copy number is written as 2.
"""

import argparse
import sys

import numpy as np
import pandas as pd
from cyvcf2 import VCF

# One row per VCF record, as read by read_consensus_vcf()
VARIANT_COLUMNS = {
    "CHROM": str,
    "POS": np.int64,
    "REF": str,
    "ALT": str,
    "passed": bool,
    "n_alt_alleles": np.int64,
    "is_snv": bool,
    "ref_counts": np.int64,
    "alt_counts": np.int64,
}
SEGMENT_COLUMNS = ["chr", "startpos", "endpos", "nMajor", "nMinor"]
PYCLONE_COLUMNS = [
    "mutation_id",
    "sample_id",
    "ref_counts",
    "alt_counts",
    "normal_cn",
    "major_cn",
    "minor_cn",
    "tumour_content",
    "chrom",
]
AUTOSOMES = [str(number) for number in range(1, 23)]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, help="Consensus somatic VCF")
    parser.add_argument("--segments", required=True, help="ASCAT segments.txt")
    parser.add_argument("--purityploidy", required=True, help="ASCAT purityploidy.txt")
    parser.add_argument("--sample-id", required=True, help="Value for the sample_id column")
    parser.add_argument("--output", required=True, help="PyClone-VI input TSV to write")
    parser.add_argument("--min-depth", type=int, default=1, help="Minimum ref + alt reads (default: 1)")
    parser.add_argument("--min-alt", type=int, default=1, help="Minimum alt reads (default: 1)")
    parser.add_argument(
        "--include-indels",
        action="store_true",
        help="Keep indels too. Off by default because ONT indel allele fractions are noisy",
    )
    return parser.parse_args()


# --------------------------------------------------------------------------------------
# Reading the inputs
# --------------------------------------------------------------------------------------


def read_consensus_vcf(path):
    """Read every record of the consensus VCF into a table with one row per variant."""
    records = []
    for variant in VCF(path):
        ref_counts, alt_counts = tumour_read_counts(variant)
        records.append(
            {
                "CHROM": variant.CHROM,
                "POS": variant.POS,
                "REF": variant.REF,
                "ALT": ",".join(variant.ALT),
                "passed": variant.FILTER is None,
                "n_alt_alleles": len(variant.ALT),
                "is_snv": variant.is_snp,
                "ref_counts": ref_counts,
                "alt_counts": alt_counts,
            }
        )

    variants = pd.DataFrame(records, columns=list(VARIANT_COLUMNS)).astype(VARIANT_COLUMNS)
    # The VCF writes 'chr1' while ASCAT writes '1', so both are compared without the prefix
    variants["chrom"] = variants["CHROM"].str.removeprefix("chr")
    return variants


def tumour_read_counts(variant):
    """Return the tumour's (ref, alt) reads from FORMAT/AD, or (-1, -1) when AD is not usable.

    AD is usable when it holds exactly two counts, ref then alt. The consensus splits a
    multi-allelic site into one record per ALT but keeps the full AD (ref,alt1,alt2) on each,
    so it is not possible to tell which count belongs to the record's ALT.
    """
    allelic_depths = variant.format("AD")
    if allelic_depths is None or allelic_depths.shape[1] != 2:
        return -1, -1

    ref_counts, alt_counts = allelic_depths[0]
    if ref_counts < 0 or alt_counts < 0:
        return -1, -1
    return int(ref_counts), int(alt_counts)


def read_ascat_segments(path):
    """Read ASCAT segments, with major_cn always the larger of the two allele copy numbers."""
    segments = pd.read_csv(path, sep="\t", dtype={"chr": str})

    missing = [column for column in SEGMENT_COLUMNS if column not in segments.columns]
    if missing:
        sys.exit(f"ERROR: {path} is missing the columns {', '.join(missing)}")

    segments = segments.dropna(subset=SEGMENT_COLUMNS)
    if segments.empty:
        sys.exit(f"ERROR: {path} has no copy number segments")

    return pd.DataFrame(
        {
            "chrom": segments["chr"].str.removeprefix("chr"),
            "startpos": segments["startpos"].astype(np.int64),
            "endpos": segments["endpos"].astype(np.int64),
            "major_cn": np.maximum(segments["nMajor"], segments["nMinor"]).astype(np.int64),
            "minor_cn": np.minimum(segments["nMajor"], segments["nMinor"]).astype(np.int64),
        }
    )


def read_ascat_purity(path):
    """Read the tumour purity (AberrantCellFraction) from ASCAT purityploidy.txt."""
    purity_ploidy = pd.read_csv(path, sep="\t")
    if "AberrantCellFraction" not in purity_ploidy.columns or purity_ploidy.empty:
        sys.exit(f"ERROR: {path} has no AberrantCellFraction value")

    purity = float(purity_ploidy["AberrantCellFraction"].iloc[0])
    # ASCAT writes a purity of 0 when it finds no solution
    if not 0 < purity <= 1:
        sys.exit(f"ERROR: the ASCAT purity in {path} is {purity}, so there is no usable solution")
    return purity


# --------------------------------------------------------------------------------------
# Adding copy number to each variant
# --------------------------------------------------------------------------------------


def add_copy_number(variants, segments):
    """Add major_cn and minor_cn from the ASCAT segment each variant falls in.

    merge_asof pairs every variant with the last segment that starts at or before it on the
    same chromosome. A variant past that segment's end is in a gap and gets NaN.
    """
    merged = pd.merge_asof(
        variants.rename_axis("vcf_row").reset_index().sort_values("POS"),
        segments.sort_values("startpos"),
        left_on="POS",
        right_on="startpos",
        by="chrom",
        direction="backward",
    )
    in_gap = merged["POS"] > merged["endpos"]
    merged.loc[in_gap, ["major_cn", "minor_cn"]] = np.nan

    # merge_asof needs the rows sorted by position, so put them back in VCF order
    return merged.sort_values("vcf_row")


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------


def keep(table, wanted, reason, dropped):
    """Keep the rows where `wanted` is True and record how many were dropped and why."""
    dropped[reason] = int((~wanted).sum())
    return table[wanted]


def main():
    args = parse_args()

    variants = read_consensus_vcf(args.vcf)
    segments = read_ascat_segments(args.segments)
    purity = read_ascat_purity(args.purityploidy)
    total_records = len(variants)
    dropped = {}

    variants = keep(variants, variants["passed"], "not PASS", dropped)
    variants = keep(variants, variants["n_alt_alleles"] == 1, "more than one ALT allele", dropped)
    if not args.include_indels:
        variants = keep(variants, variants["is_snv"], "not an SNV", dropped)
    variants = keep(variants, variants["chrom"].isin(AUTOSOMES), "not on an autosome", dropped)
    variants = keep(variants, variants["alt_counts"] >= 0, "no usable ref,alt FORMAT/AD", dropped)

    total_reads = variants["ref_counts"] + variants["alt_counts"]
    enough_reads = (variants["alt_counts"] >= args.min_alt) & (total_reads >= args.min_depth)
    variants = keep(variants, enough_reads, "too few reads", dropped)

    variants = add_copy_number(variants, segments)
    variants = keep(variants, variants["major_cn"].notna(), "not on an ASCAT segment", dropped)
    # A major copy number of 0 means both alleles are lost, so the mutation cannot be present
    variants = keep(variants, variants["major_cn"] > 0, "copy number 0", dropped)

    mutation_id = variants["CHROM"] + ":" + variants["POS"].astype(str) + ":" + variants["REF"] + ":" + variants["ALT"]
    pyclone_input = pd.DataFrame(
        {
            "mutation_id": mutation_id,
            "sample_id": args.sample_id,
            "ref_counts": variants["ref_counts"],
            "alt_counts": variants["alt_counts"],
            "normal_cn": 2,
            "major_cn": variants["major_cn"].astype(np.int64),
            "minor_cn": variants["minor_cn"].astype(np.int64),
            "tumour_content": purity,
            "chrom": variants["CHROM"],
        },
        columns=PYCLONE_COLUMNS,
    )
    pyclone_input = keep(pyclone_input, ~pyclone_input["mutation_id"].duplicated(), "duplicate record", dropped)
    pyclone_input.to_csv(args.output, sep="\t", index=False)

    print(f"{args.sample_id}: kept {len(pyclone_input)} of {total_records} VCF records", file=sys.stderr)
    for reason, count in dropped.items():
        if count:
            print(f"    dropped {count:>7}  {reason}", file=sys.stderr)
    if pyclone_input.empty:
        print(f"WARNING: no mutations left for {args.sample_id}, so there is nothing to cluster", file=sys.stderr)


if __name__ == "__main__":
    main()
