process SIMPLESVANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/51/516ef42cf72f6bbc3393ad30c3a100581789a9d434f5c46cdd3c7f6e2317f806/data'
        : 'community.wave.seqera.io/library/python_pyvcf:a216efba3d3414de'}"

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
