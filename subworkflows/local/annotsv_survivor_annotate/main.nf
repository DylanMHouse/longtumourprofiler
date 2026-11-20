// Merge the somatic SVs from Severus and SAVANA with SURVIVOR, then annotate the merged calls with AnnotSV

include { BCFTOOLS_VIEW as SAVANA_PASS_BCFTOOLS } from '../../../modules/nf-core/bcftools/view/main'
include { SURVIVOR_MERGE                        } from '../../../modules/nf-core/survivor/merge/main'
include { BCFTOOLS_SORT as SAVANA_SEVERUS_SORT  } from '../../../modules/nf-core/bcftools/sort/main'
include { ANNOTSV_ANNOTSV                       } from '../../../modules/nf-core/annotsv/annotsv/main'

workflow ANNOTSV_SURVIVOR_ANNOTATE {

    take:
    ch_savana_severus // channel: [ val(meta), path(savana_vcf), path(severus_vcf) ]
    ch_annotations    // channel: [ val(meta2), path(annotsv_annotations_dir) ]

    main:

    SAVANA_PASS_BCFTOOLS(
        ch_savana_severus.map { meta, savana_vcf, _severus_vcf -> [meta, savana_vcf, []] },
        [],
        [],
        [],
    )

    ch_survivor_input = ch_savana_severus
        .map { meta, _savana_vcf, severus_vcf -> [meta, severus_vcf] }
        .join(SAVANA_PASS_BCFTOOLS.out.vcf, by: 0)
        .map { meta, severus_vcf, savana_pass_vcf -> [meta, [severus_vcf, savana_pass_vcf]] }

    SURVIVOR_MERGE(
        ch_survivor_input,
        1000,  // max distance between breakpoints (bp)
        2,     // min number of callers supporting an SV
        1,     // SV type must agree
        0,     // ignore strand
        0,     // do not estimate distance from SV size
        30,    // min SV size (bp)
    )

    SAVANA_SEVERUS_SORT(SURVIVOR_MERGE.out.vcf)

    ANNOTSV_ANNOTSV(
        SAVANA_SEVERUS_SORT.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        ch_annotations,
        [[], []],
        [[], []],
        [[], []],
    )

    emit:
    savana_pass_vcf = SAVANA_PASS_BCFTOOLS.out.vcf // channel: [ val(meta), path(vcf) ]
    merged_vcf      = SAVANA_SEVERUS_SORT.out.vcf  // channel: [ val(meta), path(vcf) ]
    annotsv_tsv     = ANNOTSV_ANNOTSV.out.tsv      // channel: [ val(meta), path(tsv) ]
    annotsv_vcf     = ANNOTSV_ANNOTSV.out.vcf      // channel: [ val(meta), path(vcf) ]
}
