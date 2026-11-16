process HLALOH {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f4/f49ac2eb585fa3baef887e0a440a590957a25b9a6a2698cf7d33def50440a087/data'
        : 'community.wave.seqera.io/library/pysam_scipy_matplotlib-base:73605de772cf1390'}"

    input:
    tuple val(meta), path(tumour_bam), path(tumour_bai), path(normal_bam), path(normal_bai), path(normal_vcf), val(purity), val(ploidy)
    tuple val(meta2), path(regions_bed)

    output:
    tuple val(meta), path("${prefix}.hla_loh.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.hla_loh.sites.tsv"), emit: sites
    tuple val(meta), path("${prefix}.hla_loh.json"), emit: json
    tuple val(meta), path("${prefix}.hla_loh.png"), emit: png, optional: true
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), emit: versions_pysam, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // The normal BAM only filters skewed het sites, so the script runs without it
    def normal_bam_arg = normal_bam ? "--normal-bam ${normal_bam}" : ''
    // Without a purity the script reports allelic imbalance instead of an LOH call.
    // ASCAT writes 0 when it finds no solution, which also counts as no purity.
    def purity_arg = purity ? "--purity ${purity}" : ''
    def ploidy_arg = ploidy ? "--ploidy ${ploidy}" : ''
    def regions_arg = regions_bed ? "--regions-bed ${regions_bed}" : ''
    """
    hla_loh_phased_baf.py \\
        --normal-vcf ${normal_vcf} \\
        --tumor-bam ${tumour_bam} \\
        ${normal_bam_arg} \\
        ${purity_arg} \\
        ${ploidy_arg} \\
        ${regions_arg} \\
        --out-prefix ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hla_loh.tsv
    touch ${prefix}.hla_loh.sites.tsv
    touch ${prefix}.hla_loh.json
    touch ${prefix}.hla_loh.png
    """
}
