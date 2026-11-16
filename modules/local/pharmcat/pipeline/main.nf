process PHARMCAT_PIPELINE {
    tag "${meta.id}"
    label 'process_low'

    // The official image bundles the GRCh38 reference that pharmcat_pipeline needs.
    // The bioconda package downloads it into the environment on the first run.
    conda "${moduleDir}/environment.yml"
    container 'docker.io/pgkb/pharmcat:3.4.0'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.report.tsv")          , emit: report_tsv      , optional: true
    tuple val(meta), path("*.report.html")         , emit: report_html     , optional: true
    tuple val(meta), path("*.report.json")         , emit: report_json     , optional: true
    tuple val(meta), path("*.match.json")          , emit: match_json
    tuple val(meta), path("*.phenotype.json")      , emit: phenotype_json
    tuple val(meta), path("*.preprocessed.vcf.bgz"), emit: preprocessed_vcf
    tuple val(meta), path("*.missing_pgx_var.vcf") , emit: missing_pgx_var , optional: true
    tuple val("${task.process}"), val('pharmcat'), eval("pharmcat_pipeline --version | sed 's/PharmCAT v//'"), topic: versions, emit: versions_pharmcat

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.pharmcat"
    def max_memory = task.memory ? "--max-memory ${task.memory.toGiga()}G" : ''
    """
    pharmcat_pipeline \\
        ${args} \\
        --output-dir . \\
        --base-filename ${prefix} \\
        --max-concurrent-processes ${task.cpus} \\
        ${max_memory} \\
        ${vcf}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.pharmcat"
    """
    echo ${args}

    touch ${prefix}.report.tsv
    touch ${prefix}.match.json
    touch ${prefix}.phenotype.json
    echo "" | gzip > ${prefix}.preprocessed.vcf.bgz
    """
}
