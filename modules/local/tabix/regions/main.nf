process TABIX_REGIONS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/914b9b38bbfa22c458ceee6c578fd3d70ee1bffc33d3a064b0663ed4944e7872/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.25--130aa7603bf3bf5c'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path(bed)

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf
    tuple val("${task.process}"), val('tabix'), eval("tabix --version | sed '1!d; s/^.*(htslib) //'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-fh'
    prefix = task.ext.prefix ?: "${meta.id}.regions"
    """
    tabix \\
        ${args} \\
        -R ${bed} \\
        ${vcf} \\
        > ${prefix}.vcf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.regions"
    """
    touch ${prefix}.vcf
    """
}
