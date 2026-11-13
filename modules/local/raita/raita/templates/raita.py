#!/usr/bin/env python3
"""RAITA - Relatedness Assessment via IBS and Tumour Alleles.

Germline/tumour SNP identity comparison, anchored on germline-het sites.
Accompanies the BALTI nanopore tumour/normal pipeline.

PRIMARY MODE (default): compares the two BAMs directly at sites where the
GERMLINE is heterozygous. These sites are robust to somatic point mutations
(a heterozygote already carries both common alleles, so a somatic SNV cannot
manufacture a spurious "mismatch" the way it can at a homozygous site). For
each germline-het site the germline's own two alleles are identified, then
the tumour is described relative to them:

    het            tumour keeps both germline alleles
    LoH>ALLELE     tumour lost one germline allele (copy-number / cn-LoH)
    Foreign        tumour is dominated by a NON-germline allele (red flag)

Concordant = the tumour still consists of germline alleles (het or LoH).
Only 'Foreign' is evidence against identity. Depths are shown so individual
calls can be vetted.

A SECONDARY block reports germline-homozygous sites and the opposite-
homozygote (IBS0) count. Treat these with caution: in a mutated/aneuploid
tumour, somatic SNVs and germline allelic dropout corrupt homozygous calls
and can inflate IBS0 above the unrelated expectation.

--all-sites restores the older pooled-allele concordance table + parametric
identity likelihood over every callable site.

This file is a Nextflow template. Nextflow interpolates dollar-prefixed
variables and consumes one level of backslash escaping before Python runs it.
Options are read from task.ext.args.

BED columns: chrom  start  end  label  [maf]
"""

import argparse
import contextlib
import math
import platform
import shlex
import sys
from collections import Counter

import pysam

ACGT = ("A", "C", "G", "T")
NT_ORDER = ("A", "C", "T", "G")
HET_LO, HET_HI = 0.25, 0.75


BANNER = r"""
============================================================
 ____      _    ___ _____  _
|  _ \\    / \\  |_ _|_   _|/ \\
| |_) |  / _ \\  | |  | | / _ \\
|  _ <  / ___ \\ | |  | |/ ___ \\
|_| \\_\\/_/   \\_\\___| |_/_/   \\_\\

  Relatedness Assessment via IBS and Tumour Alleles
  germline/tumour SNP identity QC   (accompanies BALTI)
============================================================
"""


def make_emit(report_path):
    """Return an emit(line) that writes console output to stderr and,
    if report_path is set, tees it to that text file. Returns (emit, fh)."""
    fh = open(report_path, "w") if report_path else None

    def emit(line="", to_file=False):
        print(line, file=sys.stderr)
        if fh is not None and to_file:
            print(line, file=fh)
    return emit, fh


def base_counts(bam, chrom, pos1, min_bq, min_mq, max_depth):
    pos0 = pos1 - 1
    c = Counter()
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True,
                          min_base_quality=min_bq, min_mapping_quality=min_mq,
                          max_depth=max_depth, stepper="samtools"):
        for read in col.pileups:
            if read.is_del or read.is_refskip:
                continue
            b = read.alignment.query_sequence[read.query_position].upper()
            if b in ACGT:
                c[b] += 1
    return c


def iter_bed(path):
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            label = f[3] if len(f) > 3 else ""
            maf = float(f[4]) if len(f) > 4 and f[4] not in ("", ".") else None
            for pos0 in range(start, end):
                pos1 = pos0 + 1
                yield chrom, pos1, (label or f"{chrom}:{pos1}"), maf


def hwe_freqs(maf):
    p = maf if maf <= 0.5 else 1.0 - maf
    q = 1.0 - p
    return [q * q, 2 * p * q, p * p]


# --------------------------------------------------------------------------
# PRIMARY MODE: germline-het direct comparison
# --------------------------------------------------------------------------
def run_het_mode(args, bam1, bam2, emit):
    print("\\t".join(["label", "chrom", "pos", "g_dp", "t_dp",
                     "germline", "g_vaf", "tum_state", "t_vaf",
                     "foreign_pct", "concordant"]))
    n_het = concordant = hethet = loh = foreign = 0
    n_germ_nocall = n_tum_nocall = 0
    vaf_absdiff = []
    ibs0_obs = 0
    e_ibs0 = 0.0
    n_hom = 0

    for chrom, pos1, label, maf in iter_bed(args.bed):
        gc = base_counts(bam1, chrom, pos1, args.min_bq, args.min_mq, args.max_depth)
        tc = base_counts(bam2, chrom, pos1, args.min_bq, args.min_mq, args.max_depth)
        g_dp, t_dp = sum(gc.values()), sum(tc.values())

        g_ranked = sorted(ACGT, key=lambda b: (-gc[b], b))
        a1, a2 = g_ranked[0], g_ranked[1]         # germline's two commonest alleles
        g_inform = gc[a1] + gc[a2]
        if g_inform < args.min_depth:
            n_germ_nocall += 1
            continue
        g_minor_frac = gc[a2] / g_inform
        germ_het = HET_LO <= g_minor_frac <= HET_HI

        # ---- secondary IBS0 bookkeeping over homozygous germline sites ----
        if not germ_het:
            n_hom += 1
            if maf is not None:
                f = hwe_freqs(maf)
                e_ibs0 += 2 * f[0] * f[2]
            t_ranked = sorted(ACGT, key=lambda b: (-tc[b], b))
            t_inform = tc[t_ranked[0]] + tc[t_ranked[1]]
            if t_inform >= args.min_depth:
                t_minor = tc[t_ranked[1]] / t_inform
                tum_hom = t_minor < HET_LO
                if tum_hom and t_ranked[0] != a1:
                    ibs0_obs += 1
            continue

        # ---- germline-het site: the primary comparison ----
        n_het += 1
        t_germ = tc[a1] + tc[a2]
        t_foreign = t_dp - t_germ
        if t_dp < args.min_depth:
            n_tum_nocall += 1
            print("\\t".join([label, chrom, str(pos1), str(g_dp), str(t_dp),
                             f"{a1}/{a2}", f"{g_minor_frac*100:.0f}", "lowdepth",
                             "NA", "NA", "NA"]))
            continue
        foreign_frac = t_foreign / t_dp
        germ_frac = t_germ / t_dp

        if germ_frac < 0.5:
            state, is_conc = "Foreign", False
            foreign += 1
            t_vaf_str = "NA"
        else:
            t_minor_frac = tc[a2] / t_germ if t_germ else 0.0
            if t_minor_frac < HET_LO:
                state = f"LoH>{a1}"
                loh += 1
            elif t_minor_frac > HET_HI:
                state = f"LoH>{a2}"
                loh += 1
            else:
                state = "het"
                hethet += 1
                vaf_absdiff.append(abs(g_minor_frac - t_minor_frac) * 100)
            is_conc = True
            t_vaf_str = f"{t_minor_frac*100:.0f}"
            concordant += 1

        print("\\t".join([label, chrom, str(pos1), str(g_dp), str(t_dp),
                         f"{a1}/{a2}", f"{g_minor_frac*100:.0f}", state,
                         t_vaf_str, f"{foreign_frac*100:.0f}",
                         "Y" if is_conc else "N"]))

    # ---------------- summary ----------------
    emit(f"\\n# ===== PRIMARY: germline-heterozygous sites =====", to_file=True)
    emit(f"# germline-het sites compared: {n_het}"
        + (f"   (germline NoCall: {n_germ_nocall}, tumour NoCall: {n_tum_nocall})"
           if (n_germ_nocall or n_tum_nocall) else ""), to_file=True)
    if n_het:
        pct = 100 * concordant / n_het
        emit(f"# concordant (retain a germline allele): {concordant}/{n_het} = {pct:.1f}%", to_file=True)
        emit(f"#   of which het/het: {hethet}   LoH: {loh}", to_file=True)
        emit(f"# discordant (foreign non-germline allele): {foreign}  <- the only swap signal here", to_file=True)
        if vaf_absdiff:
            md = sum(vaf_absdiff) / len(vaf_absdiff)
            emit(f"# het/het VAF agreement: mean |Δ| = {md:.1f}% over {len(vaf_absdiff)} sites", to_file=True)
        if foreign == 0:
            verdict = ("PASS: no foreign alleles; tumour consistent with this germline + LoH. "
                       "Rules out contamination/wrong-sample, but does NOT by itself exclude a "
                       "human-human swap (unrelated people share these common alleles) -> see IBS0.")
        else:
            verdict = (f"{foreign} foreign-allele site(s): investigate contamination / wrong sample")
        emit(f"# VERDICT: {verdict}", to_file=True)
    emit(f"\\n# ----- SECONDARY: germline-homozygous sites (artifact-prone) -----")
    emit(f"# germline-hom sites: {n_hom}")
    emit(f"# opposite-homozygote (IBS0): observed {ibs0_obs}"
        + (f"  vs {e_ibs0:.2f} expected if unrelated" if e_ibs0 else ""))
    emit(f"# NB: somatic SNVs and germline allelic dropout corrupt these calls in a")
    emit(f"#     mutated tumour and can push IBS0 above the unrelated line. Vet with")
    emit(f"#     depth (columns above), strand balance and mappability before using.")


# --------------------------------------------------------------------------
# --all-sites MODE: pooled-allele concordance + parametric LR (legacy)
# --------------------------------------------------------------------------
def snap(frac):
    return 0 if frac < HET_LO else (50 if frac <= HET_HI else 100)


def geno_index(pct):
    return {0: 0, 50: 1, 100: 2}[pct]


def transition_matrix(eps, loh):
    e2, e1 = 2 * eps / 3, eps / 3
    t = [[0.0] * 3 for _ in range(3)]
    t[0] = [1 - eps, e2, e1]
    t[2] = [e1, e2, 1 - eps]
    stay = (1 - loh) * (1 - eps)
    side = loh / 2 + (1 - loh) * eps / 2
    t[1] = [side, stay, side]
    return t


def nt_percentages(c):
    tot = sum(c[b] for b in ACGT)
    return {b: (100.0 * c[b] / tot if tot else None) for b in NT_ORDER}


def fmt(x):
    return f"{x:.1f}" if x is not None else "NA"


def run_all_sites(args, bam1, bam2, emit):
    T = transition_matrix(args.error, args.loh)
    cols = ["label", "chrom", "pos", "bam1_pct", "bam2_pct"]
    cols += [f"bam1_{b}" for b in NT_ORDER] + [f"bam2_{b}" for b in NT_ORDER]
    cols += ["result", "LoH", "log10_LR"]
    print("\\t".join(cols))
    sum_log10lr = 0.0
    n_used = n_nocall = n_nomaf = 0
    for chrom, pos1, label, maf in iter_bed(args.bed):
        c1 = base_counts(bam1, chrom, pos1, args.min_bq, args.min_mq, args.max_depth)
        c2 = base_counts(bam2, chrom, pos1, args.min_bq, args.min_mq, args.max_depth)
        pooled = c1 + c2
        ranked = sorted(ACGT, key=lambda b: (-pooled[b], b))
        a_allele, b_allele = ranked[0], ranked[1]

        def call(c):
            a, b = c[a_allele], c[b_allele]
            n = a + b
            return (snap(b / n), None) if n >= args.min_depth else (None, None)

        pct1, _ = call(c1)
        pct2, _ = call(c2)
        nt1, nt2 = nt_percentages(c1), nt_percentages(c2)
        loh = ""
        lr_str = "NA"
        if pct1 is None or pct2 is None:
            result, p1s, p2s = "NoCall", "NA", "NA"
            n_nocall += 1
        else:
            p1s, p2s = str(pct1), str(pct2)
            if pct1 == pct2:
                result = "Match"
            elif pct1 == 50 and pct2 in (0, 100):
                result, loh = "Match", "LoH"
            else:
                result = "No match"
            if maf is None:
                n_nomaf += 1
            else:
                f = hwe_freqs(maf)
                g1, g2 = geno_index(pct1), geno_index(pct2)
                if f[g2] > 0:
                    lr = (1 - args.gamma) * T[g1][g2] / f[g2] + args.gamma
                    sum_log10lr += math.log10(lr)
                    n_used += 1
                    lr_str = f"{math.log10(lr):.3f}"
        out = [label, chrom, str(pos1), p1s, p2s]
        out += [fmt(nt1[b]) for b in NT_ORDER] + [fmt(nt2[b]) for b in NT_ORDER]
        out += [result, loh, lr_str]
        print("\\t".join(out))

    emit(f"\\n# all-sites identity likelihood  (loh={args.loh}, somatic_noise={args.gamma})", to_file=True)
    emit(f"# sites used for LR: {n_used}   NoCall: {n_nocall}   no-MAF: {n_nomaf}", to_file=True)
    if n_used:
        total = 10 ** sum_log10lr
        emit(f"# total log10 LR: {sum_log10lr:.3f}   LR: {total:.3g}", to_file=True)
        po = (args.prior / (1 - args.prior)) * total
        emit(f"# posterior P(same individual): {po/(1+po):.10f}", to_file=True)


def main():
    ap = argparse.ArgumentParser(description="Germline/tumour SNP identity comparison.")
    ap.add_argument("--min-bq", type=int, default=13)
    ap.add_argument("--min-mq", type=int, default=20)
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--max-depth", type=int, default=1_000_000)
    ap.add_argument("--all-sites", action="store_true",
                    help="Legacy: pooled-allele concordance + parametric LR over all sites")
    ap.add_argument("--error", type=float, default=0.02, help="(all-sites) genotyping error")
    ap.add_argument("--loh", type=float, default=0.3, help="(all-sites) P(LoH at het)")
    ap.add_argument("--somatic-noise", type=float, default=0.1, dest="gamma",
                    help="(all-sites) per-site LR floor")
    ap.add_argument("--prior", type=float, default=0.5, help="(all-sites) prior P(same)")
    args = ap.parse_args(shlex.split('''${task.ext.args ?: ''}'''))
    args.bed = "${bed}"
    prefix = "${task.ext.prefix ?: meta.id}"

    emit, report_fh = make_emit(f"{prefix}.raita_report.txt")
    emit(BANNER)

    bam1 = pysam.AlignmentFile("${germline_bam}")
    bam2 = pysam.AlignmentFile("${tumour_bam}")
    with open(f"{prefix}.raita.tsv", "w") as tsv, contextlib.redirect_stdout(tsv):
        if args.all_sites:
            run_all_sites(args, bam1, bam2, emit)
        else:
            run_het_mode(args, bam1, bam2, emit)
    bam1.close()
    bam2.close()
    report_fh.close()

    with open("versions.yml", "w") as fh:
        fh.write('"${task.process}":\\n')
        fh.write(f"    python: {platform.python_version()}\\n")
        fh.write(f"    pysam: {pysam.__version__}\\n")


if __name__ == "__main__":
    main()
