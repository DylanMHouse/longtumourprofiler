#!/usr/bin/env Rscript
# ===========================================================================
# run_scarhrd.R
# ---------------------------------------------------------------------------
# Compute HRD scores (HRD-LOH, telomeric AI, LST, and their sum) with scarHRD,
# starting from ASCAT allele-specific segments + ploidy.
#
# scarHRD wants a tab-delimited allele-specific segmentation file with columns:
#   SampleID  Chromosome  Start_position  End_position  total_cn  A_cn  B_cn  ploidy
# where A_cn = major allele CN, B_cn = minor allele CN (the package README gives
# exactly this ASCAT mapping: A_cn=nMajor, B_cn=nMinor, total_cn=nMajor+nMinor).
#
# This handler reads the standard ASCAT *.segments.txt (sample, chr, startpos,
# endpos, nMajor, nMinor), pulls ploidy from ascatploidy.txt (or --ploidy),
# builds the scarHRD input, and calls scar_score().
#
# Usage:
#   Rscript run_scarhrd.R --segments T9.segments.txt --ploidy_file ascatploidy.txt \
#       --sample PATIENT05_T9 --reference grch38 --outdir hrd_out
#   # or pass the number directly: --ploidy 2.64
# ===========================================================================

.libPaths("/rds/projects/b/beggsa-clinicalnanopore/adb/Rpackages/")
options(bitmapType='cairo')
options(stringsAsFactors = FALSE)
suppressMessages(library(optparse))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--segments",    type = "character", default = NULL,
              help = "ASCAT segments file (sample, chr, startpos, endpos, nMajor, nMinor)"),
  make_option("--ploidy",      type = "double",    default = NA,
              help = "tumour ploidy (e.g. 2.64); overrides --ploidy_file"),
  make_option("--ploidy_file", type = "character", default = NULL,
              help = "ASCAT ploidy file (ascatploidy.txt) if --ploidy not given"),
  make_option("--sample",      type = "character", default = NULL,
              help = "sample id (default: taken from the segments file)"),
  make_option("--reference",   type = "character", default = "grch38",
              help = "grch38 or grch37 [grch38]"),
  make_option("--outdir",      type = "character", default = "hrd_out"),
  make_option("--keep_sex",    action = "store_true", default = FALSE,
              help = "keep sex chromosomes (default: autosomes only)")
)))

if (is.null(opt$segments)) stop("--segments is required")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# --- read ASCAT segments, match columns by name (robust to ordering) -------
seg <- read.table(opt$segments, header = TRUE, sep = "\t", check.names = FALSE)
low <- tolower(colnames(seg))
pick <- function(...) {
  for (n in c(...)) if (n %in% low) return(which(low == n)[1])
  stop(paste("segments file missing one of:", paste(c(...), collapse = "/")))
}
i_samp <- tryCatch(pick("sample", "sampleid", "sample_id"), error = function(e) NA)
i_chr  <- pick("chr", "chromosome")
i_s    <- pick("startpos", "start", "start_position")
i_e    <- pick("endpos", "end", "end_position")
i_maj  <- pick("nmajor", "nmaj", "a_cn", "major_cn")
i_min  <- pick("nminor", "nmin", "b_cn", "minor_cn")

# --- ploidy ----------------------------------------------------------------
if (!is.na(opt$ploidy)) {
  ploidy <- opt$ploidy
} else if (!is.null(opt$ploidy_file)) {
  raw <- readLines(opt$ploidy_file, warn = FALSE)
  nums <- suppressWarnings(as.numeric(gsub('[^0-9.eE+-]', '', unlist(strsplit(raw, "[ \t,]+")))))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0) stop("could not parse a ploidy value from --ploidy_file")
  ploidy <- nums[length(nums)]   # last numeric token (skips any header label)
} else {
  stop("provide --ploidy or --ploidy_file")
}

sampleid <- if (!is.null(opt$sample)) opt$sample else
            if (!is.na(i_samp)) as.character(seg[1, i_samp]) else "TUMOUR"

# --- build scarHRD input ---------------------------------------------------
chrom <- as.character(seg[[i_chr]])
chrom <- ifelse(grepl("^chr", chrom, ignore.case = TRUE), chrom, paste0("chr", chrom))

scar <- data.frame(
  SampleID       = sampleid,
  Chromosome     = chrom,
  Start_position = as.integer(seg[[i_s]]),
  End_position   = as.integer(seg[[i_e]]),
  total_cn       = as.integer(seg[[i_maj]]) + as.integer(seg[[i_min]]),
  A_cn           = as.integer(seg[[i_maj]]),
  B_cn           = as.integer(seg[[i_min]]),
  ploidy         = ploidy
)

# drop sex chromosomes by default (HRD scores are defined on autosomes; sex
# chromosomes confound LOH especially in a male sample)
if (!opt$keep_sex) {
  scar <- scar[!scar$Chromosome %in% c("chrX", "chrY", "chr23", "chr24"), ]
}
# drop any segment with missing CN
scar <- scar[!is.na(scar$A_cn) & !is.na(scar$B_cn) & !is.na(scar$total_cn), ]

cat(sprintf("scarHRD input: sample=%s | %d segments | ploidy=%.3f | ref=%s\n",
            sampleid, nrow(scar), ploidy, opt$reference))
if (nrow(scar) < 20)
  cat("WARNING: few segments — HRD scores may be unreliable on a quiet genome.\n")

infile <- file.path(opt$outdir, paste0(sampleid, ".scarHRD_input.txt"))
write.table(scar, infile, row.names = FALSE, sep = "\t", quote = FALSE)

# --- run scarHRD -----------------------------------------------------------
suppressMessages(library(scarHRD))
res <- scar_score(infile, reference = opt$reference, seqz = FALSE,
                  ploidy = ploidy, chr.in.names = TRUE)

# scar_score returns a 1-row matrix: HRD (LOH), Telomeric AI, LST, HRD-sum
res_df <- as.data.frame(res)
res_df$SampleID <- sampleid
outfile <- file.path(opt$outdir, paste0(sampleid, ".HRDresults.txt"))
write.table(res_df, outfile, row.names = FALSE, sep = "\t", quote = FALSE)

cat("\nHRD scores:\n")
print(res)
cat(sprintf("\n  HRD-LOH      : %s\n", res[1, 1]))
cat(sprintf("  Telomeric AI : %s\n", res[1, 2]))
cat(sprintf("  LST          : %s\n", res[1, 3]))
cat(sprintf("  HRD-sum      : %s   (Myriad-style HRD; >=42 is the common HR-deficient threshold)\n", res[1, 4]))
cat(sprintf("\n  input  : %s\n  result : %s\n", infile, outfile))
