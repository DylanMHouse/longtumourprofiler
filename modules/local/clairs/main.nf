process CLAIRS {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/hkubal/clairs:v0.5.1"

    input:
    tuple val(meta), path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai), val(platform)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("${prefix}.indel.vcf.gz"), emit: indel_vcf, optional: true
    tuple val(meta), path("${prefix}.indel.vcf.gz.tbi"), emit: indel_tbi, optional: true
    tuple val("${task.process}"), val('clairs'), eval("run_clairs --version | sed 's/.* //'"), emit: versions_clairs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CLAIRS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clairs \\
        --tumor_bam_fn ${tumour_bam} \\
        --normal_bam_fn ${normal_bam} \\
        --ref_fn ${fasta} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir . \\
        --output_prefix ${prefix} \\
        --indel_output_prefix ${prefix}.indel \\
        --sample_name ${meta.id} \\
        ${args}
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CLAIRS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def indel_stub = args.contains('--enable_indel_calling')
        ? """
        echo "" | gzip > ${prefix}.indel.vcf.gz
        touch ${prefix}.indel.vcf.gz.tbi
        """
        : ''
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${indel_stub}
    """
}
