include { LIFTOVER_LIFTOVER } from '../../../modules/local/liftover/liftover/main'
include { BEDTOOLS_INTERSECT } from '../../../modules/nf-core/bedtools/intersect/main'
include { NANODX_CLASSIFY } from '../../../modules/local/nanodx/classify/main'

workflow METHYLATION_CALLING {
    take:
    ch_bedmethyl // channel: [ val(meta), path(bedmethyl) ] modkit pileup on the assembly the reads were aligned to
    ch_chain     // channel: [ val(meta), path(chain) ] aligned assembly to probe assembly
    ch_probes    // channel: [ val(meta), path(bed) ] BED6 of array probes, probe ID in the name column
    ch_model     // channel: [ val(meta), path(pkl) ] nanoDx model

    main:

    LIFTOVER_LIFTOVER(ch_bedmethyl, 9, ch_chain)

    BEDTOOLS_INTERSECT(
        LIFTOVER_LIFTOVER.out.mapped.combine(ch_probes.map { _meta, bed -> bed }),
        [[:], []],
        'bed',
    )

    NANODX_CLASSIFY(BEDTOOLS_INTERSECT.out.intersect, ch_model)

    emit:
    bedmethyl = BEDTOOLS_INTERSECT.out.intersect // channel: [ val(meta), path(bed) ]
    votes     = NANODX_CLASSIFY.out.votes        // channel: [ val(meta), path(tsv) ]
    summary   = NANODX_CLASSIFY.out.summary      // channel: [ val(meta), path(txt) ]
}
