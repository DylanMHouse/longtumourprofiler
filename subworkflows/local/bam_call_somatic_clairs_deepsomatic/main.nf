include { CLAIRS             } from '../../../modules/local/clairs/main'
include { DEEPSOMATIC        } from '../../../modules/nf-core/deepsomatic/main'

workflow BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC {

    take:
    ch_normal_bam,
    ch_somatic_bam,


    main:
    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
}
