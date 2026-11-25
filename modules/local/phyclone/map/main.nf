process PHYCLONE_MAP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d5/d5c798b6bd1461dcca776693c5b41bc5f26ce2d1d0fbd5a2094d5878ccb49fd1/data'
        : 'community.wave.seqera.io/library/phyclone:0.8.0--4b0e1bf122b4b674'}"

    input:
    tuple val(meta), path(trace)

    output:
    tuple val(meta), path("${prefix}.phyclone_map.tsv"), emit: table
    tuple val(meta), path("${prefix}.phyclone_map.nwk"), emit: tree
    tuple val(meta), path("${prefix}.phyclone_map_prevalence.tsv"), emit: prevalence
    tuple val("${task.process}"), val('phyclone'), eval("phyclone --version | sed 's/.*version //'"), topic: versions, emit: versions_phyclone

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    phyclone map \\
        --in-file ${trace} \\
        --out-table-file ${prefix}.phyclone_map.tsv \\
        --out-tree-file ${prefix}.phyclone_map.nwk \\
        --out-sample-prev-table ${prefix}.phyclone_map_prevalence.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'mutation_id\\tclone_id\\tcluster_id\\tsample_id\\tccf\\tclonal_prev\\n' > ${prefix}.phyclone_map.tsv
    echo '(0)root;' > ${prefix}.phyclone_map.nwk
    printf 'sample_id\\tclone_id\\tccf\\tclonal_prev\\n' > ${prefix}.phyclone_map_prevalence.tsv
    """
}
