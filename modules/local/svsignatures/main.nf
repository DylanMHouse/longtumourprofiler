process SVSIGNATURES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sigprofilerassignment:1.1.5--pyh106432d_0'
        : 'quay.io/biocontainers/sigprofilerassignment:1.1.5--pyh106432d_0'}"

    input:
    tuple val(meta), path(bedpe)
    path(sv_signatures)

    output:
    tuple val(meta), path("${prefix}/sv_summary.json"), emit: summary
    tuple val(meta), path("${prefix}/sv_matrix/*.SV32.matrix.tsv"), emit: matrix
    tuple val(meta), path("${prefix}/sv_signatures.tsv"), emit: signatures, optional: true
    tuple val(meta), path("${prefix}"), emit: outdir
    tuple val("${task.process}"), val('sigprofilermatrixgenerator'), eval("python -c 'import SigProfilerMatrixGenerator; print(SigProfilerMatrixGenerator.__version__)'"), topic: versions, emit: versions_sigprofilermatrixgenerator
    tuple val("${task.process}"), val('sigprofilerassignment'), eval("python -c 'import SigProfilerAssignment; print(SigProfilerAssignment.__version__)'"), topic: versions, emit: versions_sigprofilerassignment

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def refit = sv_signatures ? "--sv-signatures ${sv_signatures}" : ''
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    run_sv_signatures.py \\
        --bedpe ${bedpe} \\
        --sample-id ${meta.id} \\
        --output ${prefix} \\
        ${refit} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def touch_refit = sv_signatures ? "printf 'signature\\tcount\\tproportion\\n' > ${prefix}/sv_signatures.tsv" : ''
    """
    mkdir -p ${prefix}/sv_matrix
    echo '{}' > ${prefix}/sv_summary.json
    touch ${prefix}/sv_matrix/DETERMINE_SV.SV32.matrix.tsv
    ${touch_refit}
    """
}
