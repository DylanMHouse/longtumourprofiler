process SIMPLESVANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/05/05da49ea791f47ba15d5e0e498760cefbef0ce77fe387523e3c293db34c07bbf/data'
        : 'community.wave.seqera.io/library/python_pyvcf:71c2dc57bdf865ed'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.simpleann.vcf"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    simple_sv_annotation.py \\
        ${args} \\
        --output ${prefix}.simpleann.vcf \\
        ${vcf}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.simpleann.vcf
    """
}
