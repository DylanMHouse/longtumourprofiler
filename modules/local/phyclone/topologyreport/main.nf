process PHYCLONE_TOPOLOGYREPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5c798b6bd1461dcca776693c5b41bc5f26ce2d1d0fbd5a2094d5878ccb49fd1/data'
        : 'community.wave.seqera.io/library/phyclone:0.8.0--4b0e1bf122b4b674'}"

    input:
    tuple val(meta), path(trace)

    output:
    tuple val(meta), path("${prefix}.phyclone_topology_report.tsv"), emit: report
    tuple val(meta), path("${prefix}.phyclone_topologies.tar.gz"), emit: archive
    tuple val("${task.process}"), val('phyclone'), eval("phyclone --version | sed 's/.*version //'"), topic: versions, emit: versions_phyclone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    phyclone topology-report \\
        --in-file ${trace} \\
        --out-file ${prefix}.phyclone_topology_report.tsv \\
        --topologies-archive ${prefix}.phyclone_topologies.tar.gz \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'topology_id\\talpha\\titer\\tlog_p\\tnum_nodes\\tnum_outliers\\tnum_roots\\tchain_num\\tcount\\n' > ${prefix}.phyclone_topology_report.tsv
    echo "" | gzip > ${prefix}.phyclone_topologies.tar.gz
    """
}
