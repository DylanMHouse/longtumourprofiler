include { SAVANA_RUN_LOWMEM } from '../../../modules/local/savana/runlowmem/main'
include { SAVANA_CLASSIFY } from '../../../modules/nf-core/savana/classify/main'
include { SAVANA_CNA } from '../../../modules/nf-core/savana/cna/main'
include { SEVERUS_LOWMEM } from '../../../modules/local/severuslowmem/main'

workflow BAM_CALL_SV_CNA_TN {

    take:
    ch_tumour_normal // channel: [ val(tumour_meta), path(tumour_bam), path(tumour_bai), val(normal_meta), path(normal_bam), path(normal_bai) ]
    ch_fasta_fai     // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_savana_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai]
        }

    SAVANA_RUN_LOWMEM(
        ch_savana_input,
        ch_fasta_fai,
        [[], []],
    )

    SAVANA_CLASSIFY(
        SAVANA_RUN_LOWMEM.out.sv_breakpoints_vcf,
    )

    ch_savana_cna_input = ch_tumour_normal
        .join(SAVANA_RUN_LOWMEM.out.sv_breakpoints_vcf, by: 0)
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai, breakpoints ->
            [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, [], [], breakpoints]
        }

    SAVANA_CNA(
        ch_savana_cna_input,
        ch_fasta_fai,
        [],
        [],
        '',
    )

    ch_severus_input = ch_tumour_normal
        .map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, []]
        }

    SEVERUS_LOWMEM(
        ch_severus_input,
        [[], []],
    )

    emit:
    savana_vcf                    = SAVANA_RUN_LOWMEM.out.sv_breakpoints_vcf
    savana_bedpe                  = SAVANA_RUN_LOWMEM.out.sv_breakpoints_bedpe
    savana_read_support           = SAVANA_RUN_LOWMEM.out.sv_breakpoints_read_support
    savana_inserted_sequences     = SAVANA_RUN_LOWMEM.out.inserted_sequences
    savana_classified_vcf         = SAVANA_CLASSIFY.out.classified_vcf
    savana_somatic_vcf            = SAVANA_CLASSIFY.out.somatic_vcf
    savana_germline_vcf           = SAVANA_CLASSIFY.out.germline_vcf
    savana_somatic_bedpe          = SAVANA_CLASSIFY.out.somatic_bedpe
    savana_cna                    = SAVANA_CNA.out.cna
    savana_cna_segments           = SAVANA_CNA.out.segmented_log2r
    savana_cna_raw_counts         = SAVANA_CNA.out.raw_read_counts
    savana_cna_binned_ref         = SAVANA_CNA.out.binned_ref
    savana_cna_ranked_solutions   = SAVANA_CNA.out.ranked_solutions
    savana_cna_purity_ploidy      = SAVANA_CNA.out.fitted_purity_ploidy
    savana_cna_no_fit             = SAVANA_CNA.out.no_fit
    savana_cna_allele_counts      = SAVANA_CNA.out.allele_counts
    severus_all_vcf               = SEVERUS_LOWMEM.out.all_vcf
    severus_somatic_vcf           = SEVERUS_LOWMEM.out.somatic_vcf
    severus_all_clusters          = SEVERUS_LOWMEM.out.all_breakpoints_clusters
    severus_all_clusters_list     = SEVERUS_LOWMEM.out.all_breakpoints_clusters_list
    severus_somatic_clusters      = SEVERUS_LOWMEM.out.somatic_breakpoints_clusters
    severus_somatic_clusters_list = SEVERUS_LOWMEM.out.somatic_breakpoints_clusters_list
    severus_log                   = SEVERUS_LOWMEM.out.log
    severus_breakpoints_double    = SEVERUS_LOWMEM.out.breakpoints_double
    severus_read_ids              = SEVERUS_LOWMEM.out.read_ids
    severus_collapsed_dup         = SEVERUS_LOWMEM.out.collapsed_dup
    severus_loh                   = SEVERUS_LOWMEM.out.loh
    severus_all_plots             = SEVERUS_LOWMEM.out.all_plots
    severus_somatic_plots         = SEVERUS_LOWMEM.out.somatic_plots
}
