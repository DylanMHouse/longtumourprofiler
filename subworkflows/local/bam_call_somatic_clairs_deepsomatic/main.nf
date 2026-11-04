include { CLAIRS                } from '../../../modules/local/clairs/main'
include { DEEPSOMATIC           } from '../../../modules/nf-core/deepsomatic/main'
include { PARABRICKS_DEEPVARIANT } from '../../../modules/nf-core/parabricks/deepvariant/main'

workflow BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC {

    take:
    ch_tumour_normal // channel: [ val(tumour_meta), path(tumour_bam), path(tumour_bai), val(normal_meta), path(normal_bam), path(normal_bai) ]
    ch_fasta_fai        // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_clairs_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, 'ont']
        }

    CLAIRS(
        ch_clairs_input,
        ch_fasta_fai
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

    // TODO add the DeepVariant model through ext.args (--model) in conf/modules.config
    ch_deepvariant_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, _normal_bam, _normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, []]
        }

    PARABRICKS_DEEPVARIANT(
        ch_deepvariant_input,
        ch_fasta_fai,
    )

    emit:
    clairs_vcf        = CLAIRS.out.vcf                        // channel: [ val(meta), path(vcf.gz) ]
    clairs_tbi        = CLAIRS.out.tbi                        // channel: [ val(meta), path(vcf.gz.tbi) ]
    deepsomatic_vcf   = DEEPSOMATIC.out.vcf                   // channel: [ val(meta), path(vcf.gz) ]
    deepsomatic_tbi   = DEEPSOMATIC.out.vcf_tbi               // channel: [ val(meta), path(vcf.gz.tbi) ]
    deepvariant_vcf   = PARABRICKS_DEEPVARIANT.out.vcf        // channel: [ val(meta), path(vcf.gz) ]
}
