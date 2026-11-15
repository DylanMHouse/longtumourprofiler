process ONTSOMATICCONSENSUS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/914b9b38bbfa22c458ceee6c578fd3d70ee1bffc33d3a064b0663ed4944e7872/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.25--130aa7603bf3bf5c'}"

    input:
    tuple val(meta), path(deepsomatic_vcf), path(deepsomatic_tbi), path(clairs_vcf), path(clairs_tbi)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(mask)

    output:
    tuple val(meta), path("${prefix}.final.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.final.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${prefix}.summary.txt"), emit: summary
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), emit: versions_bcftools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ont_somatic_consensus.sh \\
        --ref ${fasta} \\
        --deepsomatic ${deepsomatic_vcf} \\
        --clairs ${clairs_vcf} \\
        --mask ${mask} \\
        --prefix ${prefix} \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.final.vcf.gz
    touch ${prefix}.final.vcf.gz.tbi
    touch ${prefix}.summary.txt
    """
}
