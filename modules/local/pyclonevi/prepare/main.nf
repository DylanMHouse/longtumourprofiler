process PYCLONEVI_PREPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/84/84073edfa808b06614ca5278a21d4e4b20992245934144df530e41a61b488488/data'
        : 'community.wave.seqera.io/library/cyvcf2_numpy_pandas:699767bb316a199c'}"

    input:
    tuple val(meta), path(vcf), path(segments), path(purityploidy)

    output:
    tuple val(meta), path("${prefix}.pyclonevi_input.tsv"), emit: tsv
    tuple val("${task.process}"), val('cyvcf2'), eval("python -c 'import cyvcf2; print(cyvcf2.__version__)'"), topic: versions, emit: versions_cyvcf2
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pyclonevi_prepare.py \\
        --vcf ${vcf} \\
        --segments ${segments} \\
        --purityploidy ${purityploidy} \\
        --sample-id ${meta.id} \\
        --output ${prefix}.pyclonevi_input.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'mutation_id\\tsample_id\\tref_counts\\talt_counts\\tnormal_cn\\tmajor_cn\\tminor_cn\\ttumour_content\\n' > ${prefix}.pyclonevi_input.tsv
    """
}
