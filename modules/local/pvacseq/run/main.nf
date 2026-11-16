process PVACSEQ_RUN {
    tag "${meta.id}"
    label 'process_high'

    // The official image bundles the IEDB tools (/opt/iedb) and the MHCflurry models.
    // The bioconda package has neither, so the module does not support Conda.
    container 'docker.io/griffithlab/pvactools:7.1.5'

    input:
    tuple val(meta), path(vcf), val(hla_alleles)
    val algorithms

    output:
    tuple val(meta), path("MHC_Class_*"), emit: outdir
    tuple val(meta), path("MHC_Class_*/*.all_epitopes.tsv"), emit: all_epitopes
    tuple val(meta), path("MHC_Class_*/*.all_epitopes.aggregated.tsv"), emit: aggregated, optional: true
    tuple val(meta), path("MHC_Class_*/*.filtered.tsv"), emit: filtered, optional: true
    tuple val("${task.process}"), val('pvactools'), eval("python -c 'from importlib.metadata import version; print(version(\"pvactools\"))'"), emit: versions_pvactools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("PVACSEQ_RUN module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    // The sample name must match the tumour sample column of the VCF
    def sample_name = task.ext.prefix ?: "${meta.id}"
    def alleles = hla_alleles instanceof List ? hla_alleles.join(',') : hla_alleles
    def algorithm_list = algorithms instanceof List ? algorithms.join(' ') : algorithms
    """
    pvacseq run \\
        ${vcf} \\
        ${sample_name} \\
        '${alleles}' \\
        ${algorithm_list} \\
        . \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def sample_name = task.ext.prefix ?: "${meta.id}"
    """
    mkdir MHC_Class_I
    touch MHC_Class_I/${sample_name}.all_epitopes.tsv
    touch MHC_Class_I/${sample_name}.all_epitopes.aggregated.tsv
    touch MHC_Class_I/${sample_name}.filtered.tsv
    """
}
