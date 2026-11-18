process MAFSMITH {
    tag "$meta.id"
    label 'process_medium'

    container 'ghcr.io/dylanmhouse/mafsmith:0.1.0-e009e98'

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta), path(fai)
    path gff3

    output:
    tuple val(meta), path("${prefix}.maf"), emit: maf
    tuple val("${task.process}"), val('mafsmith'), eval("mafsmith --version | sed 's/^mafsmith //'"), topic: versions, emit: versions_mafsmith

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mafsmith vcf2maf \\
        --input-vcf ${vcf} \\
        --output-maf ${prefix}.maf \\
        --ref-fasta ${fasta} \\
        --gff3 ${gff3} \\
        --tumor-id ${meta.id} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf
    """
}
