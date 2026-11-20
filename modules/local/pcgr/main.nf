process PCGR {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/sigven/pcgr:2.2.5.singularity'
        : 'nf-core/pcgr:2.2.5'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna), path(cpsr_classification), path(cpsr_yaml), val(purity), val(ploidy)
    path pcgr_dir
    path vep_cache

    output:
    tuple val(meta), path("${prefix}"), emit: pcgr_reports
    tuple val(meta), path("${prefix}/*.pcgr.*.html"), emit: html, optional: true
    tuple val(meta), path("${prefix}/*.pcgr.*.xlsx"), emit: xlsx, optional: true
    tuple val(meta), path("${prefix}/*.pcgr.*.snv_indel_ann.tsv.gz"), emit: snv_indel_tsv, optional: true
    tuple val("${task.process}"), val('pcgr'), eval("pcgr --version | sed 's/pcgr //g'"), topic: versions, emit: versions_pcgr

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome = task.ext.genome ?: ''
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def cna_arg = cna ? "--input_cna ${cna}" : ''
    // ASCAT writes 0 when it finds no solution, which also counts as no estimate
    def purity_arg = purity ? "--tumor_purity ${purity}" : ''
    def ploidy_arg = ploidy ? "--tumor_ploidy ${ploidy}" : ''
    // PCGR only accepts the CPSR classification together with the CPSR config
    def cpsr_arg = cpsr_classification && cpsr_yaml ? "--input_cpsr ${cpsr_classification} --input_cpsr_yaml ${cpsr_yaml}" : ''
    """
    export XDG_CACHE_HOME=/tmp
    export XDG_DATA_HOME=/tmp
    export QUARTO_PRINT_STACK=true

    mkdir -p ${prefix}

    pcgr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache} \\
        --refdata_dir ${pcgr_dir} \\
        --output_dir ${prefix} \\
        --genome_assembly ${genome} \\
        --sample_id ${prefix} \\
        ${cna_arg} \\
        ${purity_arg} \\
        ${ploidy_arg} \\
        ${cpsr_arg} \\
        ${args}
    """

    stub:
    def genome = task.ext.genome ?: 'grch38'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.pcgr.${genome}.html
    touch ${prefix}/${prefix}.pcgr.${genome}.xlsx
    echo "" | gzip > ${prefix}/${prefix}.pcgr.${genome}.snv_indel_ann.tsv.gz
    """
}
