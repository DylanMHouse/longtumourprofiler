process WAKHAN_ALL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/21/21575f17544aed81bcea4f3da16ffbbe58d962a50430b3142a73daf8fbab4941/data'
        : 'community.wave.seqera.io/library/wakhan:0.5.0--6fc755c715dbafb7'}"

    input:
    tuple val(meta), path(tumour_bam), path(tumour_bai), path(normal_vcf), path(tumour_vcf), path(breakpoints)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(cancer_genes)
    tuple val(meta4), path(centromere_bed)

    output:
    tuple val(meta), path("${prefix}"), emit: outdir
    tuple val(meta), path("${prefix}/solutions_ranks.tsv"), emit: ranks
    tuple val(meta), path("${prefix}/solution_rank_1/integer_profile.bed"), emit: bed
    tuple val(meta), path("${prefix}/solution_rank_1/integer_profile.vcf"), emit: vcf
    tuple val(meta), path("${prefix}/solution_rank_1/genes/genes_copynumber_states.bed"), emit: genes
    tuple val(meta), path("${prefix}/solution_rank_1/HiScanner_plots_data.zip"), emit: hiscanner
    tuple val(meta), path("${prefix}/phasing_output/rephased.vcf.gz"), path("${prefix}/phasing_output/rephased.vcf.gz.csi"), emit: rephased_vcf, optional: true
    tuple val("${task.process}"), val('wakhan'), eval("wakhan --version | sed 's/^wakhan //'"), topic: versions, emit: versions_wakhan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cancer_genes_arg = cancer_genes ? "--cancer-genes \"${cancer_genes}\"" : ''
    def centromere_arg = centromere_bed ? "--centromere-bed \"${centromere_bed}\"" : ''
    if (normal_vcf && tumour_vcf) {
        error("WAKHAN_ALL: provide either normal_vcf (tumour-normal) or tumour_vcf (tumour-only), not both, for ${meta.id}")
    }
    if (!normal_vcf && !tumour_vcf) {
        error("WAKHAN_ALL: one of normal_vcf or tumour_vcf is required for ${meta.id}")
    }
    if (!breakpoints) {
        error("WAKHAN_ALL: breakpoints VCF is required for ${meta.id}")
    }
    def phased_vcf_arg = tumour_vcf ? "--tumor-phased-vcf \"${tumour_vcf}\"" : "--normal-phased-vcf \"${normal_vcf}\""

    """
    wakhan \\
        all \\
        --threads ${task.cpus} \\
        --reference "${fasta}" \\
        --target-bam "${tumour_bam}" \\
        ${phased_vcf_arg} \\
        --breakpoints "${breakpoints}" \\
        --genome-name "${prefix}" \\
        --out-dir-plots "${prefix}" \\
        ${cancer_genes_arg} \\
        ${centromere_arg} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def solution = "${prefix}/solution_2.0_1.0_1.0"
    def mode_files = tumour_vcf
        ? "mkdir -p \"${prefix}/snps_loh_plots\"\n    touch \"${prefix}/snps_loh_plots/SNPs_LOH_INDEX.html\" \"${prefix}/coverage_data/${prefix}_loh_segments.csv\""
        : "touch \"${prefix}/coverage_data/pileup_SNPs.csv\""
    """
    echo ${args}

    mkdir -p "${solution}/genes" "${prefix}/coverage_data" "${prefix}/coverage_plots" "${prefix}/phasing_output"
    touch "${prefix}/solutions_ranks.tsv"
    touch "${prefix}/COPY_NUMBERS.html" "${prefix}/${prefix}_genome_coverage.html"
    touch "${solution}/integer_profile.bed" "${solution}/integer_profile.vcf" "${solution}/integer_profile.html"
    touch "${solution}/genes/genes_copynumber_states.bed"
    python -c "import zipfile; zipfile.ZipFile('${solution}/HiScanner_plots_data.zip', 'w').close()"
    ln -s "\$PWD/${solution}" "${prefix}/solution_rank_1"
    touch "${prefix}/coverage_data/coverage.csv" "${prefix}/coverage_data/coverage_ps.csv" "${prefix}/coverage_data/phase_corrected_coverage.csv" "${prefix}/coverage_data/baf.csv" "${prefix}/coverage_data/cancer_genes_coverage.csv"
    touch "${prefix}/coverage_plots/COVERAGE_INDEX.html" "${prefix}/phasing_output/PHASE_CORRECTION_INDEX.html"
    echo "" | bgzip > "${prefix}/phasing_output/rephased.vcf.gz"
    touch "${prefix}/phasing_output/rephased.vcf.gz.csi"
    ${mode_files}
    """
}
