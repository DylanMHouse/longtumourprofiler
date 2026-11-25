process PYCLONEVI_FIT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pyclone-vi:0.2.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/pyclone-vi:0.2.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${prefix}.pyclonevi.h5"), emit: h5
    tuple val(meta), path("${prefix}.pyclonevi.tsv"), emit: results
    tuple val("${task.process}"), val('pyclone-vi'), eval("pyclone-vi --version | sed 's/.*version //'"), topic: versions, emit: versions_pyclonevi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pyclone-vi fit \\
        --in-file ${tsv} \\
        --out-file ${prefix}.pyclonevi.h5 \\
        --num-threads ${task.cpus} \\
        ${args}

    pyclone-vi write-results-file \\
        --in-file ${prefix}.pyclonevi.h5 \\
        --out-file ${prefix}.pyclonevi.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pyclonevi.h5
    printf 'mutation_id\\tsample_id\\tcluster_id\\tcellular_prevalence\\tcellular_prevalence_std\\tcluster_assignment_prob\\n' > ${prefix}.pyclonevi.tsv
    """
}
