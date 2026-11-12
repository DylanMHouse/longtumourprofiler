process NANODX_CLASSIFY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/919c1ab167d2f8f4786a6f1a0191b404468baa02a0ae9c9dac2aace80a8cf625/data'
        : 'community.wave.seqera.io/library/python_numpy_pandas_scikit-learn_pytorch:86a697850ca72dd1'}"

    input:
    tuple val(meta), path(bedmethyl)
    tuple val(meta2), path(model)

    output:
    tuple val(meta), path("${prefix}.votes.tsv"), emit: votes
    tuple val(meta), path("${prefix}.summary.txt"), emit: summary

    // The nanoDx classifier script (bin/classify_NN_bedMethyl_commandline.py) has no version flag,
    // so the version of the underlying pytorch runtime is reported instead
    tuple val("${task.process}"), val('pytorch'), eval("python -c 'import torch; print(torch.__version__)'"), topic: versions, emit: versions_pytorch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    classify_NN_bedMethyl_commandline.py \\
        ${args} \\
        --model ${model} \\
        --bed ${bedmethyl} \\
        --output_votes ${prefix}.votes.tsv \\
        --output_txt ${prefix}.summary.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.votes.tsv
    touch ${prefix}.summary.txt
    """
}
