process RAITA_RAITA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0':
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(germline_bam), path(germline_bai), path(tumour_bam), path(tumour_bai)
    path bed

    output:
    tuple val(meta), path("*.raita.tsv")       , emit: tsv
    tuple val(meta), path("*.raita_report.txt"), emit: report
    path "versions.yml"                        , emit: versions_raita, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'raita.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raita.tsv
    touch ${prefix}.raita_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pysam: \$(python3 -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
}
