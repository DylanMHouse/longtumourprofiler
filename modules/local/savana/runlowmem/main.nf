process SAVANA_RUN_LOWMEM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/632dd754c4bdc8b682b7594323684ddc435a51775c420fd11986dd83686d51c8/data'
        : 'community.wave.seqera.io/library/savana:1.3.8--9ae333d36c906f2e'}"

    input:
    tuple val(meta), path(tumour), path(tumour_index), path(normal), path(normal_index)
    tuple val(meta2), path(ref), path(ref_index)
    tuple val(meta3), path(contigs)

    output:
    tuple val(meta), path("${prefix}.sv_breakpoints.vcf"), emit: sv_breakpoints_vcf
    tuple val(meta), path("${prefix}.sv_breakpoints.bedpe"), emit: sv_breakpoints_bedpe
    tuple val(meta), path("${prefix}.sv_breakpoints_read_support.tsv"), emit: sv_breakpoints_read_support
    tuple val(meta), path("${prefix}.inserted_sequences.fa"), emit: inserted_sequences
    tuple val("${task.process}"), val("savana"), eval("python -c \"import importlib.metadata as m; print(m.version('savana'))\""), emit: versions_savana, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def contigs_arg = contigs ? "--contigs ${contigs}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    savana_lowmem.py \\
        --tumour ${tumour} \\
        --normal ${normal} \\
        --ref ${ref} \\
        --ref_index ${ref_index} \\
        --outdir ./outdir \\
        --sample ${prefix} \\
        --threads ${task.cpus} \\
        --spill_dir "\${TMPDIR:-/tmp}/savana_spill_${prefix}_${task.attempt}" \\
        ${contigs_arg} \\
        ${args}

    mv ./outdir/* .
    rmdir ./outdir
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.sv_breakpoints.vcf"
    touch "${prefix}.sv_breakpoints.bedpe"
    touch "${prefix}.sv_breakpoints_read_support.tsv"
    touch "${prefix}.inserted_sequences.fa"
    """
}
