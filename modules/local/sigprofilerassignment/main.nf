process SIGPROFILERASSIGNMENT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sigprofilerassignment:1.1.5--pyh106432d_0'
        : 'quay.io/biocontainers/sigprofilerassignment:1.1.5--pyh106432d_0'}"

    input:
    tuple val(meta), path(vcf), path(segments)
    path(genome_dir)

    output:
    tuple val(meta), path("${prefix}/signature_summary.json"), emit: summary
    tuple val(meta), path("${prefix}/*_signatures.tsv"), emit: signatures, optional: true
    tuple val(meta), path("${prefix}"), emit: outdir
    tuple val("${task.process}"), val('sigprofilerassignment'), eval("python -c 'import SigProfilerAssignment; print(SigProfilerAssignment.__version__)'"), topic: versions, emit: versions_sigprofilerassignment

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if (!vcf && !segments) {
        error("SIGPROFILERASSIGNMENT needs a somatic VCF, ASCAT segments, or both")
    }
    def vcf_arg = vcf ? "--vcf ${vcf}" : ''
    def segments_arg = segments ? "--ascat-segments ${segments}" : ''
    // The container cannot hold the ~3 GB SBS reference, so point SigProfilerMatrixGenerator at a staged
    // install made with install(genome, volume=dir). Only the VCF read-out needs it.
    def genome_volume = genome_dir ? "export SIGPROFILERMATRIXGENERATOR_VOLUME=\$PWD/${genome_dir}" : ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    ${genome_volume}

    run_sigprofiler_assignment.py \\
        ${vcf_arg} \\
        ${segments_arg} \\
        --sample-id ${meta.id} \\
        --output ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_sbs = vcf ? "printf 'signature\\tcount\\tproportion\\n' > ${prefix}/sbs96_signatures.tsv" : ''
    def touch_cnv = segments ? "printf 'signature\\tcount\\tproportion\\n' > ${prefix}/cnv_signatures.tsv" : ''
    """
    mkdir -p ${prefix}
    echo '{}' > ${prefix}/signature_summary.json
    ${touch_sbs}
    ${touch_cnv}
    """
}
