process PCGR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/sigven/pcgr:2.2.5.singularity'
        : 'nf-core/pcgr:2.2.5'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('pcgr'), eval("pcgr --version"), topic: versions, emit: versions_pcgr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pcgr \\
        $args \\
        --input_vcf <placeholder>\\
        --output_dir <placeholder>\\
        --sample_id <placeholder>\\
        --genome_assembly <placeholder>\\
        --input_cna <placeholder>\\
        --tumor_ploidy <placeholder>\\
        --tumor_purity <placeholder>\\
        --refdata_dir <placeholder>\\
        --input_cpsr <placeholder>\\
        --input_cpsr_yaml <placeholder>\\
        --estimate_signatures \\
        --all_reference_signatures \\
        --include_artefact_signatures \\
        --estimate_tmb \\
        --estimate_msi \\
        --force_overwrite

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    touch ${prefix}.bam
    """
}
