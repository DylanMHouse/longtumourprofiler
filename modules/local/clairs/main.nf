process CLAIRS {
    tag "${meta.id}"
    label 'process_high'

    container "docker.io/hkubal/clairs:v0.5.1"

    input:
    tuple val(meta), path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai), val(platform)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path bed
    path pileup_model
    path full_alignment_model

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('clairs'), eval("run_clairs --version | sed 's/.* //'"), emit: versions_clairs, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CLAIRS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_arg = bed ? "--bed_fn ${bed}" : ''
    def pileup_model_arg = pileup_model ? "--pileup_model_path ${pileup_model}" : ''
    def full_alignment_model_arg = full_alignment_model ? "--full_alignment_model_path ${full_alignment_model}" : ''
    """
    run_clairs \\
        --tumor_bam_fn ${tumour_bam} \\
        --normal_bam_fn ${normal_bam} \\
        --ref_fn ${fasta} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir . \\
        --output_prefix ${prefix} \\
        --sample_name ${meta.id} \\
        ${bed_arg} \\
        ${pileup_model_arg} \\
        ${full_alignment_model_arg} \\
        ${args}
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CLAIRS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
