include { CLAIRS                 } from '../../../modules/local/clairs/main'
include { DEEPSOMATIC            } from '../../../modules/nf-core/deepsomatic/main'
include { PARABRICKS_DEEPVARIANT } from '../../../modules/nf-core/parabricks/deepvariant/main'
include { BCFTOOLS_MERGE         } from '../../../modules/nf-core/bcftools/merge/main'

workflow BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC {

    take:
    ch_tumour_normal // channel: [ val(tumour_meta), path(tumour_bam), path(tumour_bai), val(normal_meta), path(normal_bam), path(normal_bai) ]
    ch_fasta_fai        // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_clairs_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, 'ont_r10_dorado_sup_5khz_ssrs']
        }

    CLAIRS(
        ch_clairs_input,
        ch_fasta_fai
    )

    // Merge SNV & indels from Clairs
    ch_clairs_merge_input = CLAIRS.out.vcf
        .join(CLAIRS.out.tbi, by: 0)
        .join(CLAIRS.out.indel_vcf, by: 0)
        .join(CLAIRS.out.indel_tbi, by: 0)
        .map { meta, snv_vcf, snv_tbi, indel_vcf, indel_tbi ->
            [meta, [snv_vcf, indel_vcf], [snv_tbi, indel_tbi], []]
        }

    BCFTOOLS_MERGE(
        ch_clairs_merge_input,
        [[], [], []],
    )

    ch_deepsomatic_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, normal_bam, normal_bai, tumour_bam, tumour_bai, []]
        }

    DEEPSOMATIC(
        ch_deepsomatic_input,
        ch_fasta_fai,
        channel.value([[], []]),
    )
    def ch_deepvariant_vcf = channel.empty()
    if ('gpu' in workflow.profile.tokenize(',')) {
    ch_deepvariant_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, _normal_bam, _normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, []]
        }

    PARABRICKS_DEEPVARIANT(
        ch_deepvariant_input,
        ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] }
    )
    ch_deepvariant_vcf = PARABRICKS_DEEPVARIANT.out.vcf
    }

    emit:
    clairs_vcf        = BCFTOOLS_MERGE.out.vcf             // channel: [ val(meta), path(vcf.gz) ]
    clairs_tbi        = BCFTOOLS_MERGE.out.index           // channel: [ val(meta), path(vcf.gz.tbi) ]
    deepsomatic_vcf   = DEEPSOMATIC.out.vcf                   // channel: [ val(meta), path(vcf.gz) ]
    deepsomatic_tbi   = DEEPSOMATIC.out.vcf_tbi               // channel: [ val(meta), path(vcf.gz.tbi) ]
    deepvariant_vcf   = ch_deepvariant_vcf
}
