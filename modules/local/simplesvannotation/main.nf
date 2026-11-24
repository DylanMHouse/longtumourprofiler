process SIMPLESVANNOTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0d/0de903c2323d2b699be590e018b1294f65b18385ee3f8a2f46ba6fd544e133e2/data'
        : 'community.wave.seqera.io/library/python_pyvcf:36f0a163c22db62a'}"

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
