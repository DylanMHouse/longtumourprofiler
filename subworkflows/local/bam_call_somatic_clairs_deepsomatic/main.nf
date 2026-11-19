include { CLAIRS } from '../../../modules/local/clairs/main'
include { DEEPSOMATIC } from '../../../modules/nf-core/deepsomatic/main'
include { PARABRICKS_DEEPVARIANT } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { BCFTOOLS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as DEEPSOMATIC_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'

// Split the reference into n_chunks BED chunks of near equal length, cutting contigs only at chunk boundaries.
// Returns [ [chunk_name, bed_text], ... ].
def splitFaiIntoChunks(fai, n_chunks) {
    def contigs = fai
        .readLines()
        .collect { line -> line.tokenize('\t') }
        .collect { fields -> [fields[0], fields[1].toLong()] }
        .findAll { _name, length -> length > 0 }
    def total = contigs.sum { _name, length -> length }
    def target = Math.max(1L, (total + n_chunks - 1).intdiv(n_chunks))
    def offset = 0L
    def pieces = contigs.collectMany { name, length ->
        def contig_start = offset
        def contig_end = offset + length
        offset = contig_end
        (contig_start.intdiv(target)..(contig_end - 1).intdiv(target)).collect { k ->
            def start = Math.max(contig_start, k * target) - contig_start
            def end = Math.min(contig_end, (k + 1) * target) - contig_start
            [k, "${name}\t${start}\t${end}"]
        }
    }
    return pieces
        .groupBy { k, _line -> k }
        .collect { k, chunk_pieces -> [String.format('chunk%03d', k as Integer), chunk_pieces.collect { _k, line -> line }.join('\n') + '\n'] }
}

workflow BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC {
    take:
    ch_tumour_normal // channel: [ val(tumour_meta), path(tumour_bam), path(tumour_bai), val(normal_meta), path(normal_bam), path(normal_bai) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    deepsomatic_chunks // integer: number of regions DeepSomatic is scattered over

    main:

    ch_clairs_input = ch_tumour_normal.map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
        [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, 'ont_r10_dorado_sup_5khz_ssrs']
    }

    CLAIRS(
        ch_clairs_input,
        ch_fasta_fai,
    )

    // Concatenate the ClairS SNVs and indels into one VCF with a single tumour sample column
    ch_clairs_concat_input = CLAIRS.out.vcf
        .join(CLAIRS.out.tbi, by: 0)
        .join(CLAIRS.out.indel_vcf, by: 0)
        .join(CLAIRS.out.indel_tbi, by: 0)
        .map { meta, snv_vcf, snv_tbi, indel_vcf, indel_tbi ->
            [meta, [snv_vcf, indel_vcf], [snv_tbi, indel_tbi]]
        }

    BCFTOOLS_CONCAT(ch_clairs_concat_input)

    // DeepSomatic is slow on a whole genome, so scatter it over reference chunks and merge the calls
    def ch_deepsomatic_beds = ch_fasta_fai
        .flatMap { _meta, _fasta, fai -> splitFaiIntoChunks(fai, deepsomatic_chunks) }
        .collectFile { chunk_name, bed_text -> ["${chunk_name}.bed", bed_text] }
    def ch_deepsomatic_n_chunks = ch_deepsomatic_beds.count()

    ch_deepsomatic_input = ch_tumour_normal
        .combine(ch_deepsomatic_beds)
        .combine(ch_deepsomatic_n_chunks)
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai, bed, n_chunks ->
            [tumour_meta + [chunk: bed.baseName, n_chunks: n_chunks], normal_bam, normal_bai, tumour_bam, tumour_bai, bed]
        }

    DEEPSOMATIC(
        ch_deepsomatic_input,
        ch_fasta_fai,
        channel.value([[], []]),
    )

    // groupKey lets each sample merge as soon as its own chunks finish
    ch_deepsomatic_concat_input = DEEPSOMATIC.out.vcf
        .join(DEEPSOMATIC.out.vcf_tbi, by: 0)
        .map { meta, vcf, tbi ->
            [groupKey(meta - meta.subMap(['chunk', 'n_chunks']), meta.n_chunks), vcf, tbi]
        }
        .groupTuple()
        .map { key, vcfs, tbis -> [key.getGroupTarget(), vcfs, tbis] }

    DEEPSOMATIC_CONCAT(ch_deepsomatic_concat_input)

    def ch_deepvariant_vcf = channel.empty()
    if ('gpu' in workflow.profile.tokenize(',')) {
        ch_deepvariant_input = ch_tumour_normal.map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, _normal_bam, _normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, []]
        }

        PARABRICKS_DEEPVARIANT(
            ch_deepvariant_input,
            ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        )
        ch_deepvariant_vcf = PARABRICKS_DEEPVARIANT.out.vcf
    }

    emit:
    clairs_vcf = BCFTOOLS_CONCAT.out.vcf // channel: [ val(meta), path(vcf.gz) ]
    clairs_tbi = BCFTOOLS_CONCAT.out.index // channel: [ val(meta), path(vcf.gz.tbi) ]
    deepsomatic_vcf = DEEPSOMATIC_CONCAT.out.vcf // channel: [ val(meta), path(vcf.gz) ]
    deepsomatic_tbi = DEEPSOMATIC_CONCAT.out.index // channel: [ val(meta), path(vcf.gz.tbi) ]
    deepvariant_vcf = ch_deepvariant_vcf
}
