process MAFSMITH_FETCH {
    tag "${genome}"
    label 'process_single'

    container 'ghcr.io/dylanmhouse/mafsmith:0.1.0-e009e98'

    input:
    val genome

    output:
    path "mafsmith_data/*/genes.gff3.gz", emit: gff3
    tuple val("${task.process}"), val('mafsmith'), eval("mafsmith --version | sed 's/^mafsmith //'"), topic: versions, emit: versions_mafsmith

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def build = [grch38: 'GRCh38', grch37: 'GRCh37', grcm39: 'GRCm39'][genome]
    if (!build) {
        error("MAFSMITH_FETCH: unsupported genome '${genome}', expected grch38, grch37 or grcm39")
    }

    """
    mkdir -p mafsmith_data/${build}
    touch mafsmith_data/${build}/reference.fa.gz mafsmith_data/${build}/reference.fa

    mafsmith fetch \\
        --genome ${genome} \\
        --data-dir mafsmith_data \\
        --skip-fastvep \\
        ${args}

    rm mafsmith_data/${build}/reference.fa.gz mafsmith_data/${build}/reference.fa
    """

    stub:
    def build = [grch38: 'GRCh38', grch37: 'GRCh37', grcm39: 'GRCm39'][genome]
    """
    mkdir -p mafsmith_data/${build}
    echo "" | gzip > mafsmith_data/${build}/genes.gff3.gz
    """
}
