process LIFTOVER_LIFTOVER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/15/15744931def04387cf095376c9d8c322bf21ca542d64946a7258c15b40b12b82/data'
        : 'community.wave.seqera.io/library/ucsc-liftover:482--f8639164c0f3b5bf'}"

    input:
    tuple val(meta), path(methyl)
    val bedPlus
    tuple val(meta2), path(chain)

    output:
    tuple val(meta), path("${prefix}.bedMethyl"), emit: mapped
    // liftOver has no version flag, so the version is pinned to the container build
    tuple val("${task.process}"), val('ucsc'), val('482'), topic: versions, emit: versions_ucsc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    liftOver \\
        ${args} \\
        -bedPlus=${bedPlus} \\
        ${methyl} \\
        ${chain} \\
        ${prefix}.bedMethyl \\
        unmapped.bed

    rm unmapped.bed
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedMethyl
    """
}
