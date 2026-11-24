// Large CNV calling (>100kb) with Spectre from mosdepth binned coverage and SNV calls

include { MOSDEPTH as MOSDEPTH_SPECTRE } from '../../../modules/nf-core/mosdepth/main'
include { SPECTRE_CNVCALLER } from '../../../modules/local/spectre/cnvcaller/main'

workflow BAM_CALL_CNV_SPECTRE {
    take:
    ch_input // channel: [ val(meta), path(bam), path(bai) ]
    ch_snv // channel: [ val(meta), path(vcf) ] small variant calls for the same samples
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_spectre_metadata // channel: path(mdr) or [] (built from the fasta when missing)
    ch_spectre_blacklist // channel: path(bed) or []

    main:

    // The window size is set with ext.args (--by) in conf/modules.config
    MOSDEPTH_SPECTRE(
        ch_input.map { meta, bam, bai -> [meta, bam, bai, []] },
        ch_fasta,
        [],
    )

    ch_spectre_input = MOSDEPTH_SPECTRE.out.summary_txt
        .join(MOSDEPTH_SPECTRE.out.regions_bed, by: 0)
        .join(MOSDEPTH_SPECTRE.out.regions_csi, by: 0)
        .join(ch_snv, by: 0)

    SPECTRE_CNVCALLER(
        ch_spectre_input,
        ch_fasta,
        ch_spectre_metadata,
        ch_spectre_blacklist,
    )

    emit:
    cnv_vcf = SPECTRE_CNVCALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    cnv_bed = SPECTRE_CNVCALLER.out.bed // channel: [ val(meta), path(bed) ]
    karyotype = SPECTRE_CNVCALLER.out.karyotype // channel: [ val(meta), path(txt) ]
}
