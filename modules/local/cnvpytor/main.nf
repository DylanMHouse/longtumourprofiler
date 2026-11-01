process CNVPYTOR {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/dylanmhouse/cnvpytor:1.3.2'

    input:
    tuple val(meta), path(bam), path(bai), path(snv)
    tuple val(meta2), path(fasta)
    val(bin_size)

    output:
    tuple val(meta), path("${prefix}.pytor")        , emit: pytor
    tuple val(meta), path("${prefix}.cnvpytor.vcf") , emit: vcf
    tuple val(meta), path("${prefix}.cnvpytor.tsv") , emit: tsv
    tuple val("${task.process}"), val('cnvpytor'), eval("cnvpytor --version 2>&1 | sed -n 's/.*CNVpytor //p'"), emit: versions_cnvpytor, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$(pwd)/.mplconfig

    cnvpytor -root ${prefix}.pytor -j ${task.cpus} -snp ${snv}
    cnvpytor -root ${prefix}.pytor -j ${task.cpus} -rd ${bam} -T ${fasta} ${args}
    cnvpytor -root ${prefix}.pytor -j ${task.cpus} -mask_snps
    cnvpytor -root ${prefix}.pytor -j ${task.cpus} -baf ${bin_size}
    cnvpytor -root ${prefix}.pytor -j ${task.cpus} -his ${bin_size}
    cnvpytor -root ${prefix}.pytor -j ${task.cpus} -call combined ${bin_size} > ${prefix}.cnvpytor.tsv
    cnvpytor -root ${prefix}.pytor -view ${bin_size}

    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    touch ${prefix}.pytor
    touch ${prefix}.cnvpytor.vcf
    touch ${prefix}.cnvpytor.tsv
    """
}
