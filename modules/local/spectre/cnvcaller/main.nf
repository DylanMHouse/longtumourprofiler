process SPECTRE_CNVCALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker.io/ontresearch/spectre:sha42472d37a5a992c3ee27894a23dce5e2fff66d27':
        'docker.io/ontresearch/spectre:sha42472d37a5a992c3ee27894a23dce5e2fff66d27' }"

    input:
    tuple val(meta), path(summary_txt), path(regions_bed), path(regions_csi), path(snv)
    tuple val(meta2), path(fasta)
    path metadata  // optional
    path blacklist // optional

    output:
    tuple val(meta), path("${prefix}_cnv.vcf")                , emit: vcf
    tuple val(meta), path("${prefix}_cnv.bed")                , emit: bed
    tuple val(meta), path("${prefix}.spc.gz")                 , emit: spc
    tuple val(meta), path("${prefix}_predicted_karyotype.txt"), emit: karyotype
    tuple val("${task.process}"), val("spectre"), eval("spectre version 2>&1 | sed 's/.*version: //'"), emit: versions_spectre, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_spectre"
    def metadata_arg = metadata ? "--metadata ${metadata}" : ''
    def blacklist_arg = blacklist ? "--blacklist ${blacklist}" : ''
    def bin_size = params.spectre_bin_size
    """
    mkdir -p coverage
    ln -s ../${summary_txt} coverage/${prefix}.mosdepth.summary.txt
    ln -s ../${regions_bed} coverage/${prefix}.regions.bed.gz
    ln -s ../${regions_csi} coverage/${prefix}.regions.bed.gz.csi

    spectre \\
        CNVCaller \\
        --bin-size $bin_size \\
        --coverage coverage/ \\
        --sample-id ${prefix} \\
        --output-dir spectre/ \\
        --reference ${fasta} \\
        --snv ${snv} \\
        --threads $task.cpus \\
        $metadata_arg \\
        $blacklist_arg \\
        $args

    mv spectre/${prefix}.vcf ${prefix}_cnv.vcf
    mv spectre/${prefix}_cnv.bed ${prefix}_cnv.bed
    mv spectre/${prefix}.spc.gz ${prefix}.spc.gz
    mv spectre/predicted_karyotype.txt ${prefix}_predicted_karyotype.txt

    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_spectre"
    """
    echo $args
    touch ${prefix}_cnv.vcf
    touch ${prefix}_cnv.bed
    echo "" | gzip > ${prefix}.spc.gz
    touch ${prefix}_predicted_karyotype.txt
    """
}
