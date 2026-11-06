#!/usr/bin/env bash
#
# ont_somatic_consensus.sh
# -----------------------------------------------------------------------------
# Stratified consensus of ClairS + DeepSomatic somatic small-variant calls
# (matched tumour/normal, Oxford Nanopore) into a single high-accuracy VCF.
#
# Implements stages 1-5:
#   1. Normalise both VCFs (left-align, split multiallelics, dedup) vs same REF
#   2. Split SNV / indel streams (PASS only) and intersect each separately
#   3. Stratified consensus rule:
#        - shared calls           -> accept
#        - SNV private  (high VAF) -> rescue either caller at base QUAL floor
#        - SNV private  (low  VAF) -> DeepSomatic moderate; ClairS strict+min-alt
#        - indel private          -> DeepSomatic moderate; ClairS strict+AF floor
#                                     (both pass homopolymer/STR mask)
#   4. ONT-specific precision filters: homopolymer/STR mask (private indels),
#        optional panel-of-normals, tumour depth/alt floors, user passthrough
#   5. Recombine, sort, emit consensus; optional stratified som.py eval vs truth
#
# IMPORTANT: every QUAL/VAF/depth threshold below is a PLACEHOLDER default.
#   The accuracy of the output depends entirely on tuning these against a truth
#   set (CASTLE / COLO829 / HCC1395) for YOUR coverage & purity regime. Run with
#   --truth-vcf/--highconf-bed to score, sweep the thresholds, then lock them in.
#
# Dependencies: bcftools (>=1.18 recommended), bgzip, tabix.
#               Optional: som.py (hap.py suite) for the eval stage.
# -----------------------------------------------------------------------------

set -euo pipefail

# ----------------------------- defaults --------------------------------------
REF=""
CLAIRS_VCF=""
DS_VCF=""
HP_STR_BED=""                # homopolymer + tandem-repeat BED to mask (required)
PON_VCF=""                   # optional panel-of-normals VCF (sites excluded)
OUTDIR="consensus_out"
PREFIX="consensus"
THREADS=4
KEEP_TMP=0

# Sample handling: name of the TUMOUR sample column in each caller's VCF.
# "auto" -> use the first sample (index 0), which is the usual case for these
# somatic output VCFs. Override if your VCFs carry both T and N columns.
CLAIRS_TUMOR="auto"
DS_TUMOR="auto"

# FORMAT tag names carrying tumour allele fraction / depth in each caller.
CLAIRS_AF="AF"               # ClairS reports AF in FORMAT
CLAIRS_DP="DP"
DS_AF="VAF"                  # DeepSomatic (DeepVariant lineage) reports VAF
DS_DP="DP"

# Stratification knob.
VAF_SPLIT="0.2"              # boundary between "high" and "low" VAF SNV regimes

# Generic floors (applied to every rescued, non-shared call).
MIN_DP_T=12                  # min tumour depth
MIN_ALT=4                    # min tumour ALT-supporting reads (computed AF*DP)

# Tier-specific QUAL / AF floors  (PLACEHOLDERS - TUNE AGAINST TRUTH).
DS_SNV_QUAL=2                # DeepSomatic-private SNV  (reliable at low VAF)
CLAIRS_SNV_QUAL_HI=2         # ClairS-private SNV, VAF >= VAF_SPLIT
CLAIRS_SNV_QUAL_LO=12        # ClairS-private SNV, VAF <  VAF_SPLIT (strict)
DS_INDEL_QUAL=4              # DeepSomatic-private indel (high precision arm)
CLAIRS_INDEL_QUAL=12         # ClairS-private indel (recall arm, strict)
CLAIRS_INDEL_AF=0.10         # ClairS-private indel AF floor

# Optional caller-agnostic extra filter expressions (e.g. strand bias /
# read-position tags that vary by build). Applied with `bcftools view -i`.
EXTRA_SNV_FILTER=""
EXTRA_INDEL_FILTER=""

# Optional evaluation.
TRUTH_VCF=""
HIGHCONF_BED=""

# ----------------------------- usage -----------------------------------------
usage() {
  cat <<EOF
Usage: $(basename "$0") --ref REF.fa --clairs CLAIRS.vcf.gz --deepsomatic DS.vcf.gz \\
                        --mask homopolymer_str.bed [options]

Required:
  --ref FILE              Reference FASTA (must match both VCFs, .fai present)
  --clairs FILE           ClairS VCF (bgzipped, indexed)
  --deepsomatic FILE      DeepSomatic VCF (bgzipped, indexed)
  --mask FILE             Homopolymer + STR BED for indel masking

Output:
  --outdir DIR            Output directory                 (default: $OUTDIR)
  --prefix STR            Output file prefix               (default: $PREFIX)
  --threads N             Threads for bcftools             (default: $THREADS)
  --keep-tmp              Keep intermediate files

Samples / tags:
  --clairs-tumor NAME     Tumour sample in ClairS VCF      (default: $CLAIRS_TUMOR)
  --ds-tumor NAME         Tumour sample in DeepSomatic VCF (default: $DS_TUMOR)
  --clairs-af TAG         ClairS FORMAT AF tag             (default: $CLAIRS_AF)
  --ds-af TAG             DeepSomatic FORMAT AF tag        (default: $DS_AF)
  --clairs-dp TAG / --ds-dp TAG   FORMAT depth tags        (default: DP / DP)

Thresholds (PLACEHOLDERS - tune against truth):
  --vaf-split F           High/low VAF boundary            (default: $VAF_SPLIT)
  --min-dp N              Min tumour depth                 (default: $MIN_DP_T)
  --min-alt N             Min tumour ALT reads (AF*DP)     (default: $MIN_ALT)
  --ds-snv-qual F         DeepSomatic-private SNV QUAL     (default: $DS_SNV_QUAL)
  --clairs-snv-qual-hi F  ClairS-private SNV QUAL, hiVAF   (default: $CLAIRS_SNV_QUAL_HI)
  --clairs-snv-qual-lo F  ClairS-private SNV QUAL, loVAF   (default: $CLAIRS_SNV_QUAL_LO)
  --ds-indel-qual F       DeepSomatic-private indel QUAL   (default: $DS_INDEL_QUAL)
  --clairs-indel-qual F   ClairS-private indel QUAL        (default: $CLAIRS_INDEL_QUAL)
  --clairs-indel-af F     ClairS-private indel AF floor    (default: $CLAIRS_INDEL_AF)
  --extra-snv-filter EXPR    bcftools -i expr added to all SNV rescues
  --extra-indel-filter EXPR  bcftools -i expr added to all indel rescues

Optional precision / eval:
  --pon FILE              Panel-of-normals VCF (sites excluded)
  --truth-vcf FILE        Truth VCF for som.py eval
  --highconf-bed FILE     High-confidence BED for som.py eval

  -h, --help              Show this help
EOF
}

# ----------------------------- arg parsing -----------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2;;
    --clairs) CLAIRS_VCF="$2"; shift 2;;
    --deepsomatic) DS_VCF="$2"; shift 2;;
    --mask) HP_STR_BED="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --prefix) PREFIX="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --keep-tmp) KEEP_TMP=1; shift;;
    --clairs-tumor) CLAIRS_TUMOR="$2"; shift 2;;
    --ds-tumor) DS_TUMOR="$2"; shift 2;;
    --clairs-af) CLAIRS_AF="$2"; shift 2;;
    --ds-af) DS_AF="$2"; shift 2;;
    --clairs-dp) CLAIRS_DP="$2"; shift 2;;
    --ds-dp) DS_DP="$2"; shift 2;;
    --vaf-split) VAF_SPLIT="$2"; shift 2;;
    --min-dp) MIN_DP_T="$2"; shift 2;;
    --min-alt) MIN_ALT="$2"; shift 2;;
    --ds-snv-qual) DS_SNV_QUAL="$2"; shift 2;;
    --clairs-snv-qual-hi) CLAIRS_SNV_QUAL_HI="$2"; shift 2;;
    --clairs-snv-qual-lo) CLAIRS_SNV_QUAL_LO="$2"; shift 2;;
    --ds-indel-qual) DS_INDEL_QUAL="$2"; shift 2;;
    --clairs-indel-qual) CLAIRS_INDEL_QUAL="$2"; shift 2;;
    --clairs-indel-af) CLAIRS_INDEL_AF="$2"; shift 2;;
    --extra-snv-filter) EXTRA_SNV_FILTER="$2"; shift 2;;
    --extra-indel-filter) EXTRA_INDEL_FILTER="$2"; shift 2;;
    --pon) PON_VCF="$2"; shift 2;;
    --truth-vcf) TRUTH_VCF="$2"; shift 2;;
    --highconf-bed) HIGHCONF_BED="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1;;
  esac
done

# ----------------------------- checks ----------------------------------------
log()  { echo "[$(date +%H:%M:%S)] $*" >&2; }
die()  { echo "ERROR: $*" >&2; exit 1; }

command -v bcftools >/dev/null || die "bcftools not found on PATH"
command -v bgzip   >/dev/null || die "bgzip not found on PATH"
command -v tabix   >/dev/null || die "tabix not found on PATH"

[[ -n "$REF"         ]] || { usage; die "--ref is required"; }
[[ -n "$CLAIRS_VCF"  ]] || { usage; die "--clairs is required"; }
[[ -n "$DS_VCF"      ]] || { usage; die "--deepsomatic is required"; }
[[ -n "$HP_STR_BED"  ]] || { usage; die "--mask is required"; }
[[ -f "$REF"         ]] || die "REF not found: $REF"
[[ -f "${REF}.fai"   ]] || die "REF index ${REF}.fai not found (run: samtools faidx $REF)"
[[ -f "$CLAIRS_VCF"  ]] || die "ClairS VCF not found: $CLAIRS_VCF"
[[ -f "$DS_VCF"      ]] || die "DeepSomatic VCF not found: $DS_VCF"
[[ -f "$HP_STR_BED"  ]] || die "Mask BED not found: $HP_STR_BED"
[[ -z "$PON_VCF"     || -f "$PON_VCF" ]] || die "PoN VCF not found: $PON_VCF"

mkdir -p "$OUTDIR"
TMP="$(mktemp -d "${OUTDIR}/tmp.XXXXXX")"
cleanup() { [[ "$KEEP_TMP" -eq 0 ]] && rm -rf "$TMP" || log "Kept tmp: $TMP"; }
trap cleanup EXIT

BT="bcftools"
TH="--threads $THREADS"

# Resolve a tumour sample NAME to a 0-based column index within a VCF.
# "auto" or not-found -> 0.
resolve_idx() {
  local vcf="$1" name="$2" line
  if [[ "$name" == "auto" ]]; then echo 0; return; fi
  line=$($BT query -l "$vcf" | grep -nxF "$name" | head -n1 | cut -d: -f1 || true)
  if [[ -z "$line" ]]; then
    log "WARN: tumour sample '$name' not in $(basename "$vcf"); using index 0"
    echo 0
  else
    echo $(( line - 1 ))
  fi
}

# ============================ STAGE 1: normalise =============================
# Left-align + split multiallelics + collapse exact duplicates against REF.
# Without this, indel concordance is artificially low and consensus is junk.
log "STAGE 1  sanitising headers + normalising VCFs"

# Header self-heal: some caller outputs (esp. merged/concatenated ClairS, or
# DeepSomatic's normal-sample NAF) declare a FORMAT field as Number=A or R but
# write a mismatching value count, which makes `bcftools norm` abort on read
# (e.g. "wrong number of fields in FMT/NAF ... expected 3, found 2"). We relax
# every Number=A/R FORMAT field to Number=. so norm stops enforcing arity and
# carries it through unsplit -- EXCEPT the AF tag we actually filter on, which
# is left strict so -m-any subsets it correctly per allele and any genuine
# arity problem there still surfaces rather than being silently masked.
sanitise_header() {
  local in="$1" keep_af="$2" out="$3"
  $BT view -h "$in" | awk -v keep="$keep_af" '
    /^##FORMAT=<ID=/ {
      id=$0; sub(/.*##FORMAT=<ID=/,"",id); sub(/,.*/,"",id)
      if (id != keep && ($0 ~ /Number=A/ || $0 ~ /Number=R/))
        sub(/Number=[AR]/, "Number=.")
    }
    { print }
  ' > "$out.hdr"
  $BT reheader -h "$out.hdr" -o "$out" "$in"
  $BT index $TH -t "$out"
}

#Old norm 
#
#norm() {
#  local in="$1" out="$2"
#  $BT norm $TH -m-any -f "$REF" -Ou "$in" \
#    | $BT norm $TH -d exact -Oz -o "$out"
#  $BT index $TH -t "$out"
#}

norm() {
  local in="$1" out="$2" nbad
  nbad=$($BT view -H -i 'REF~"[^ACGTNacgtn]"' "$in" | wc -l)
  if (( nbad )); then
    echo "[$(date +%T)]   $(basename "$in"): dropping $nbad record(s) with non-ACGTN REF" >&2
  fi
  $BT view $TH -e 'REF~"[^ACGTNacgtn]"' -Ou "$in" \
    | $BT norm $TH -m-any -f "$REF" -Ou \
    | $BT norm $TH -d exact -Oz -o "$out"
  $BT index $TH -t "$out"
}

sanitise_header "$CLAIRS_VCF" "$CLAIRS_AF" "$TMP/clairs.sanitised.vcf.gz"
sanitise_header "$DS_VCF"     "$DS_AF"     "$TMP/ds.sanitised.vcf.gz"
norm "$TMP/clairs.sanitised.vcf.gz" "$TMP/clairs.norm.vcf.gz"
norm "$TMP/ds.sanitised.vcf.gz"     "$TMP/ds.norm.vcf.gz"

CI=$(resolve_idx "$TMP/clairs.norm.vcf.gz" "$CLAIRS_TUMOR")
DI=$(resolve_idx "$TMP/ds.norm.vcf.gz"     "$DS_TUMOR")
log "tumour indices -> ClairS:$CI  DeepSomatic:$DI"

# ====================== STAGE 2: split + intersect ==========================
# PASS-only, SNV and indel streams handled separately so we can apply
# class-specific consensus rules in stage 3.
log "STAGE 2  splitting SNV/indel streams and intersecting"
split_pass() {
  local in="$1" pre="$2" type="$3"
  $BT view $TH -f PASS -v "$type" "$in" -Oz -o "${pre}.vcf.gz"
  $BT index $TH -t "${pre}.vcf.gz"
}
split_pass "$TMP/clairs.norm.vcf.gz" "$TMP/clairs.snv"   snps
split_pass "$TMP/clairs.norm.vcf.gz" "$TMP/clairs.indel" indels
split_pass "$TMP/ds.norm.vcf.gz"     "$TMP/ds.snv"       snps
split_pass "$TMP/ds.norm.vcf.gz"     "$TMP/ds.indel"     indels

# isec layout per class:  0000 = ClairS-private, 0001 = DS-private,
#                         0002 = shared (ClairS repr), 0003 = shared (DS repr)
$BT isec $TH -p "$TMP/isec_snv"   "$TMP/clairs.snv.vcf.gz"   "$TMP/ds.snv.vcf.gz"
$BT isec $TH -p "$TMP/isec_indel" "$TMP/clairs.indel.vcf.gz" "$TMP/ds.indel.vcf.gz"

# ====================== STAGE 3: stratified consensus =======================
# Helper: filter a partition, keep tumour sample only, reheader to TUMOR,
# rename caller AF tag -> CAF, stamp provenance INFO, emit bgzipped+indexed.
# Args: in_vcf out_vcf filter_expr src_label class_label tumour_idx af_tag
harmonise() {
  local in="$1" out="$2" expr="$3" src="$4" cls="$5" idx="$6" aftag="$7"
  local sample rmap hdr tmpv rh tab
  sample=$($BT query -l "$in" | sed -n "$((idx+1))p")

  # rename caller AF tag (ClairS AF / DeepSomatic VAF) -> CAF
  rmap="$TMP/rename.${src}.${cls}.txt"
  echo "FORMAT/${aftag} CAF" > "$rmap"

  # provenance INFO header definitions (written once, reused)
  hdr="$TMP/consinfo.hdr"
  if [[ ! -f "$hdr" ]]; then
    {
      echo '##INFO=<ID=CONS_SRC,Number=1,Type=String,Description="Consensus provenance: SHARED|CLAIRS|DEEPSOMATIC">'
      echo '##INFO=<ID=CONS_CLASS,Number=1,Type=String,Description="Variant class: SNV|INDEL">'
    } > "$hdr"
  fi

  # filter -> keep tumour sample only -> drop unused normal-sample NAF
  # -> reheader sample to TUMOR -> rename AF/VAF to CAF.
  tmpv="$TMP/harm.${src}.${cls}.vcf.gz"
  $BT view $TH ${expr:+-i "$expr"} "$in" -Ou \
    | $BT view $TH -s "$sample" -Ou \
    | $BT annotate $TH -x FORMAT/NAF -Ou 2>/dev/null \
    | $BT reheader -s <(echo "TUMOR") \
    | $BT annotate $TH --rename-annots "$rmap" -Oz -o "$tmpv"
  $BT index $TH -t "$tmpv"

  # canonicalise CAF header (post -m-any every record is biallelic -> Number=1)
  # so concat does not warn about combining CAF definitions of different lengths
  rh="${tmpv%.vcf.gz}.rh.vcf.gz"
  $BT view -h "$tmpv" \
    | sed -E 's|##FORMAT=<ID=CAF,[^>]*>|##FORMAT=<ID=CAF,Number=1,Type=Float,Description="Tumour allele fraction (consensus)">|' \
    > "$TMP/harm.${src}.${cls}.canon.hdr"
  $BT reheader -h "$TMP/harm.${src}.${cls}.canon.hdr" -o "$rh" "$tmpv"
  $BT index $TH -t "$rh"

  # empty partition: just add the INFO defs so headers stay concat-compatible
  if [[ "$($BT view -H "$rh" | head -1 | wc -l)" -eq 0 ]]; then
    $BT annotate $TH -h "$hdr" "$rh" -Oz -o "$out"
    $BT index $TH -t "$out"
    return
  fi

  # stamp CONS_SRC / CONS_CLASS as constant String INFO via a tab annotation
  # (--mark-sites only sets flag tags; a tab file is the reliable way to set a
  #  String value on every record)
  tab="$TMP/ann.${src}.${cls}.tab"
  $BT query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$rh" \
    | awk -v s="$src" -v c="$cls" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,s,c}' > "$tab"
  bgzip -f "$tab"; tabix -s1 -b2 -e2 "$tab.gz"
  $BT annotate $TH -a "$tab.gz" -c CHROM,POS,REF,ALT,CONS_SRC,CONS_CLASS -h "$hdr" \
    "$rh" -Oz -o "$out"
  $BT index $TH -t "$out"
}

log "STAGE 3  applying stratified consensus rules"

# --- SNV tiers ---------------------------------------------------------------
GEN_SNV="FORMAT/${CLAIRS_DP}[${CI}:0] >= ${MIN_DP_T}"        # ClairS generic
GEN_SNV_DS="FORMAT/${DS_DP}[${DI}:0] >= ${MIN_DP_T}"         # DS generic
ALT_C="(FORMAT/${CLAIRS_AF}[${CI}:0]*FORMAT/${CLAIRS_DP}[${CI}:0]) >= ${MIN_ALT}"
EXSNV="${EXTRA_SNV_FILTER:+ && ($EXTRA_SNV_FILTER)}"

# ClairS-private SNV: high-VAF lenient, low-VAF strict + min-alt
C_SNV_EXPR="${GEN_SNV} && ${ALT_C} && ( \
  (FORMAT/${CLAIRS_AF}[${CI}:0] >= ${VAF_SPLIT} && QUAL >= ${CLAIRS_SNV_QUAL_HI}) || \
  (FORMAT/${CLAIRS_AF}[${CI}:0] <  ${VAF_SPLIT} && QUAL >= ${CLAIRS_SNV_QUAL_LO}) )${EXSNV}"

# DeepSomatic-private SNV: single moderate floor (reliable at low VAF)
D_SNV_EXPR="${GEN_SNV_DS} && QUAL >= ${DS_SNV_QUAL}${EXSNV}"

harmonise "$TMP/isec_snv/0000.vcf" "$TMP/snv.c_priv.vcf.gz" "$C_SNV_EXPR" CLAIRS      SNV "$CI" "$CLAIRS_AF"
harmonise "$TMP/isec_snv/0001.vcf" "$TMP/snv.d_priv.vcf.gz" "$D_SNV_EXPR" DEEPSOMATIC SNV "$DI" "$DS_AF"
# Shared SNV: accept as-is (take ClairS representation, 0002)
harmonise "$TMP/isec_snv/0002.vcf" "$TMP/snv.shared.vcf.gz" ""           SHARED      SNV "$CI" "$CLAIRS_AF"

# --- indel tiers -------------------------------------------------------------
ALT_CI="(FORMAT/${CLAIRS_AF}[${CI}:0]*FORMAT/${CLAIRS_DP}[${CI}:0]) >= ${MIN_ALT}"
EXIND="${EXTRA_INDEL_FILTER:+ && ($EXTRA_INDEL_FILTER)}"

# ClairS-private indel: recall arm, strict QUAL + AF floor + min-alt
C_IND_EXPR="FORMAT/${CLAIRS_DP}[${CI}:0] >= ${MIN_DP_T} && ${ALT_CI} && \
  QUAL >= ${CLAIRS_INDEL_QUAL} && FORMAT/${CLAIRS_AF}[${CI}:0] >= ${CLAIRS_INDEL_AF}${EXIND}"

# DeepSomatic-private indel: precision arm, moderate floor
D_IND_EXPR="FORMAT/${DS_DP}[${DI}:0] >= ${MIN_DP_T} && QUAL >= ${DS_INDEL_QUAL}${EXIND}"

harmonise "$TMP/isec_indel/0000.vcf" "$TMP/indel.c_priv.vcf.gz" "$C_IND_EXPR" CLAIRS      INDEL "$CI" "$CLAIRS_AF"
harmonise "$TMP/isec_indel/0001.vcf" "$TMP/indel.d_priv.vcf.gz" "$D_IND_EXPR" DEEPSOMATIC INDEL "$DI" "$DS_AF"
# Shared indel: accept (take DeepSomatic representation, 0003 - better indel anno)
harmonise "$TMP/isec_indel/0003.vcf" "$TMP/indel.shared.vcf.gz" ""           SHARED      INDEL "$DI" "$DS_AF"

# ================= STAGE 4: ONT-specific precision filters ===================
# Homopolymer/STR mask applies to PRIVATE indels only; shared indels are exempt
# (both callers agreeing in a repeat is far more trustworthy than one).
log "STAGE 4  masking private indels in homopolymer/STR regions + optional PoN"
mask_private_indel() {
  local in="$1" out="$2"
  $BT view $TH -T "^${HP_STR_BED}" "$in" -Oz -o "$out"
  $BT index $TH -t "$out"
}
mask_private_indel "$TMP/indel.c_priv.vcf.gz" "$TMP/indel.c_priv.masked.vcf.gz"
mask_private_indel "$TMP/indel.d_priv.vcf.gz" "$TMP/indel.d_priv.masked.vcf.gz"

# ===================== STAGE 5: recombine + (eval) ===========================
log "STAGE 5  concatenating tiers"
RAW="$TMP/${PREFIX}.raw.vcf.gz"
$BT concat $TH -a \
  "$TMP/snv.shared.vcf.gz" "$TMP/snv.c_priv.vcf.gz" "$TMP/snv.d_priv.vcf.gz" \
  "$TMP/indel.shared.vcf.gz" "$TMP/indel.c_priv.masked.vcf.gz" "$TMP/indel.d_priv.masked.vcf.gz" \
  | $BT sort -Oz -o "$RAW"
$BT index $TH -t "$RAW"

FINAL="$OUTDIR/${PREFIX}.final.vcf.gz"
if [[ -n "$PON_VCF" ]]; then
  log "       excluding panel-of-normals sites"
  $BT index -t "$PON_VCF" 2>/dev/null || true
  # drop any consensus site present in the PoN
  $BT isec $TH -p "$TMP/pon_isec" -C "$RAW" "$PON_VCF"
  $BT view $TH "$TMP/pon_isec/0000.vcf" -Oz -o "$FINAL"
else
  cp "$RAW" "$FINAL"
fi
$BT index $TH -t "$FINAL"

# Summary counts by tier/class.
log "consensus written: $FINAL"
{
  echo "=== consensus summary ($PREFIX) ==="
  printf "%-22s %s\n" "total records:" "$($BT view -H "$FINAL" | wc -l)"
  for s in SHARED CLAIRS DEEPSOMATIC; do
    for c in SNV INDEL; do
      n=$($BT view -H -i "INFO/CONS_SRC=\"$s\" && INFO/CONS_CLASS=\"$c\"" "$FINAL" 2>/dev/null | wc -l)
      printf "%-22s %s\n" "${s} ${c}:" "$n"
    done
  done
} | tee "$OUTDIR/${PREFIX}.summary.txt" >&2

# ---- optional stratified evaluation vs truth --------------------------------
if [[ -n "$TRUTH_VCF" && -n "$HIGHCONF_BED" ]]; then
  if command -v som.py >/dev/null; then
    log "EVAL    running som.py (overall + SNV-only + indel-only)"
    EVAL="$OUTDIR/eval"; mkdir -p "$EVAL"
    som.py "$TRUTH_VCF" "$FINAL" -f "$HIGHCONF_BED" -r "$REF" \
      -o "$EVAL/${PREFIX}.all" >/dev/null 2>&1 || log "WARN: som.py (all) failed"
    $BT view $TH -v snps   "$FINAL" -Oz -o "$TMP/final.snv.vcf.gz";   $BT index -t "$TMP/final.snv.vcf.gz"
    $BT view $TH -v indels "$FINAL" -Oz -o "$TMP/final.indel.vcf.gz"; $BT index -t "$TMP/final.indel.vcf.gz"
    som.py "$TRUTH_VCF" "$TMP/final.snv.vcf.gz"   -f "$HIGHCONF_BED" -r "$REF" \
      -o "$EVAL/${PREFIX}.snv"   >/dev/null 2>&1 || log "WARN: som.py (snv) failed"
    som.py "$TRUTH_VCF" "$TMP/final.indel.vcf.gz" -f "$HIGHCONF_BED" -r "$REF" \
      -o "$EVAL/${PREFIX}.indel" >/dev/null 2>&1 || log "WARN: som.py (indel) failed"
    log "eval stats in: $EVAL/  (see *.stats.csv: precision/recall/F1)"
    for f in "$EVAL/${PREFIX}".{all,snv,indel}.stats.csv; do
      [[ -f "$f" ]] && { echo "--- $(basename "$f") ---" >&2; column -t -s, "$f" >&2 || cat "$f" >&2; }
    done
  else
    log "WARN: --truth-vcf given but som.py not on PATH; skipping eval"
  fi
fi

log "done."
