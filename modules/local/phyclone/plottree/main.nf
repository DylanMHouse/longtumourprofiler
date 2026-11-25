process PHYCLONE_PLOTTREE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/88837692f8fc42ed5ad1acfde3554d1c0abca8599a20a594e58a57df4595de7f/data'
        : 'community.wave.seqera.io/library/colorcet_graphviz_networkx_pandas_pruned:3417a3b632a5429c'}"

    input:
    tuple val(meta), path(tree), path(table)

    output:
    tuple val(meta), path("${prefix}.phyclone_map_tree.png"), emit: png
    tuple val("${task.process}"), val('graphviz'), eval("dot -V 2>&1 | sed 's/.*version \\([^ ]*\\).*/\\1/'"), topic: versions, emit: versions_graphviz
    tuple val("${task.process}"), val('scikit-bio'), eval("python -c 'import skbio; print(skbio.__version__)'"), topic: versions, emit: versions_scikit_bio

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_phyclone_tree.py \\
        --newick ${tree} \\
        --results-table ${table} \\
        --output ${prefix}.phyclone_map_tree.png
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.phyclone_map_tree.png
    """
}
