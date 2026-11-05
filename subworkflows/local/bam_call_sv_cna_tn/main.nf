include { SAVANA_RUN_LOWMEM } from '../../../modules/local/savana/runlowmem/main'
include { SAVANA_CLASSIFY } from '../../../modules/nf-core/savana/classify/main'
include { SAVANA_CNA } from '../../../modules/nf-core/savana/cna/main'
include { SEVERUS_LOWMEM } from '../../../modules/local/severuslowmem/main'
include { WAKHAN_ALL } from '../../../modules/nf-core/wakhan/main'
include { ASCAT } from '../../../modules/nf-core/ascat/main'

workflow BAM_CALL_SV_CNA_TN {

    take:
    main:
    emit:

    }
