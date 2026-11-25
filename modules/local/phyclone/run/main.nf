process PHYCLONE_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5c798b6bd1461dcca776693c5b41bc5f26ce2d1d0fbd5a2094d5878ccb49fd1/data'
        : 'community.wave.seqera.io/library/phyclone:0.8.0--4b0e1bf122b4b674'}"

    input:
    tuple val(meta), path(tsv), path(clusters)

    output:
    tuple val(meta), path("${prefix}.phyclone_trace.h5"), emit: trace
    tuple val("${task.process}"), val('phyclone'), eval("phyclone --version | sed 's/.*version //'"), topic: versions, emit: versions_phyclone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    phyclone run \\
        --in-file ${tsv} \\
        --cluster-file ${clusters} \\
        --out-file ${prefix}.phyclone_trace.h5 \\
        --num-chains ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.phyclone_trace.h5
    """
}
