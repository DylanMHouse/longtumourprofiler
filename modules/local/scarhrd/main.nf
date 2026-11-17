process SCARHRD {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/r-base_r-sequenza_r-data.table_r-remotes:211c1ac371d8aaf5'
        : 'community.wave.seqera.io/library/r-base_r-sequenza_r-data.table_r-remotes:004b4e24b3e5b8b0'}"

    input:
    tuple val(meta), path(segments), path(purityploidy)

    output:
    tuple val(meta), path("${prefix}.scarhrd.tsv"),       emit: scores
    tuple val(meta), path("${prefix}.scarhrd_input.tsv"), emit: segments
    // scarHRD has no CLI version flag; 0.1.1 is the version at the commit pinned below
    tuple val("${task.process}"), val('scarHRD'), val('0.1.1'), topic: versions, emit: versions_scarhrd

    when:
    task.ext.when == null || task.ext.when

    script:
    // scarHRD is not on bioconda, so it is installed from this pinned GitHub commit (package version 0.1.1)
    def scarhrd_commit = 'c98f8bc42ed0810393a98677d415114360616725'
    // Settings arrive as a map in ext.args, the same way the ASCAT module takes them
    def args = task.ext.args ?: [:]
    prefix = task.ext.prefix ?: "${meta.id}"

    def reference      = args.reference      ?: 'grch38'
    def size_limit_loh = args.sizelimitLOH   ?: '15e6'
    def ploidy_manual  = args.ploidy         ?: 'NULL'
    // NULL means detect a 'chr' prefix from the segments themselves
    def chr_in_names   = args.chr_in_names != null ? (args.chr_in_names ? 'TRUE' : 'FALSE') : 'NULL'
    """
    Rscript - <<'EOF'
    options(bitmapType = 'cairo')
    options(stringsAsFactors = FALSE)
    # The container has no tar binary, so unpack the GitHub tarball with R's own
    Sys.setenv(TAR = "internal")

    lib <- file.path(getwd(), "R_libs")
    dir.create(lib, showWarnings = FALSE)
    .libPaths(c(lib, .libPaths()))
    if (!requireNamespace("scarHRD", quietly = TRUE)) {
        remotes::install_github("sztup/scarHRD@${scarhrd_commit}", lib = lib, upgrade = "never", dependencies = FALSE)
    }
    suppressMessages(library(scarHRD))

    # ASCAT segments: sample, chr, startpos, endpos, nMajor, nMinor
    seg <- read.table("${segments}", header = TRUE, sep = "\\t", check.names = FALSE)
    required <- c("chr", "startpos", "endpos", "nMajor", "nMinor")
    missing_cols <- setdiff(required, colnames(seg))
    if (length(missing_cols) > 0) {
        stop(paste("Segments file is missing columns:", paste(missing_cols, collapse = ", "),
                   "- scarHRD expects the ASCAT segments.txt output."))
    }
    if (nrow(seg) == 0) stop("Segments file has no rows.")

    # Ploidy comes from ASCAT purityploidy.txt unless set in ext.args
    ploidy <- ${ploidy_manual}
    if (is.null(ploidy)) {
        pp <- read.table("${purityploidy}", header = TRUE, sep = "\\t")
        ploidy <- as.numeric(pp[["Ploidy"]][1])
    }
    # ASCAT writes 0 when it finds no solution
    if (is.na(ploidy) || ploidy <= 0) {
        stop("No valid ploidy: ASCAT found no solution for this sample. Set ploidy in ext.args to override.")
    }

    chr_in_names <- ${chr_in_names}
    if (is.null(chr_in_names)) chr_in_names <- any(grepl("^chr", seg[["chr"]]))

    # scarHRD reads columns by position: sample, chromosome, start, end,
    # total_cn, A_cn, B_cn, ploidy
    scar_input <- data.frame(
        SampleID       = "${prefix}",
        Chromosome     = seg[["chr"]],
        Start_position = seg[["startpos"]],
        End_position   = seg[["endpos"]],
        total_cn       = seg[["nMajor"]] + seg[["nMinor"]],
        A_cn           = seg[["nMajor"]],
        B_cn           = seg[["nMinor"]],
        ploidy         = ploidy,
        check.names    = FALSE
    )
    write.table(scar_input, "${prefix}.scarhrd_input.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    # scar_score also writes its own results file; keep it out of the task directory
    res <- scar_score(
        "${prefix}.scarhrd_input.tsv",
        reference    = "${reference}",
        seqz         = FALSE,
        chr.in.names = chr_in_names,
        ploidy       = ploidy,
        sizelimitLOH = ${size_limit_loh},
        outputdir    = tempdir()
    )

    out <- data.frame(sample = "${prefix}", res, ploidy = ploidy, check.names = FALSE)
    write.table(out, "${prefix}.scarhrd.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)
    EOF
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\tHRD\\tTelomeric AI\\tLST\\tHRD-sum\\tploidy\\n' > ${prefix}.scarhrd.tsv
    touch ${prefix}.scarhrd_input.tsv
    """
}
