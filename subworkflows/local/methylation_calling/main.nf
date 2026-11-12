include { LIFTOVER } from '../../../modules/local/liftover/liftover/main'
include { BEDTOOLS_INTERSECT } from '../../../modules/nf-core/bedtools/intersect/main'
include { NANODX_CLASSIFY } from '../../../modules/local/nanodx/classify/main'

workflow METHYLATION_CALLING {

    take:
    ch_methylation
    bedMethyl
    ch_bedchain

    main:

    LIFTOVER(
        ch_methylation,
        bedMethyl,
        ch_bedchain
    )
    ch_intersect = LIFTOVER.out.mapped

    emit:

}
