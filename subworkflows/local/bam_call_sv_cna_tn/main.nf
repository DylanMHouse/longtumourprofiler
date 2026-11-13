include { SAVANA_RUN_LOWMEM } from '../../../modules/local/savana/runlowmem/main'
include { SAVANA_CLASSIFY } from '../../../modules/nf-core/savana/classify/main'
include { SAVANA_CNA } from '../../../modules/nf-core/savana/cna/main'
include { SEVERUS_LOWMEM } from '../../../modules/local/severuslowmem/main'

// Copy Number Alterations modules
include { WAKHAN_ALL } from '../../../modules/local/wakhan/all/main'
include { ASCAT } from '../../../modules/nf-core/ascat/main'

workflow BAM_CALL_SV_CNA_TN {
    take:
    ch_tumour_normal // channel: [ val(tumour_meta), path(tumour_bam), path(tumour_bai), val(normal_meta), path(normal_bam), path(normal_bai) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    ch_normal_phased_vcf // channel: [ val(normal_meta), path(vcf) ] LongPhase phased germline SNVs
    ch_wakhan_cancer_genes // channel: [ val(meta), path(tsv) ] or [ [:], [] ]
    ch_wakhan_centromere_bed // channel: [ val(meta), path(bed) ] or [ [:], [] ]
    ch_ascat_alleles // channel: path(allele_files)
    ch_ascat_loci // channel: path(loci_files)
    ch_ascat_gc // channel: path(gc_file) or []
    ch_ascat_rt // channel: path(rt_file) or []

    main:

    ch_savana_input = ch_tumour_normal.map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
        [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai]
    }

    def ch_savana_contigs = channel.fromPath(params.savana_contigs, checkIfExists: true)
        .map { savana_contigs -> tuple([id: "contigs"], savana_contigs) }


    SAVANA_RUN_LOWMEM(
        ch_savana_input,
        ch_fasta_fai,
        ch_savana_contigs,
    )

    SAVANA_CLASSIFY(
        SAVANA_RUN_LOWMEM.out.sv_breakpoints_vcf
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

    ch_severus_input = ch_tumour_normal.map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
        [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, []]
    }

    SEVERUS_LOWMEM(
        ch_severus_input,
        [[], []],
    )

    def ch_wakhan_outdir = channel.empty()
    def ch_wakhan_ranks = channel.empty()
    def ch_wakhan_bed = channel.empty()
    def ch_wakhan_vcf = channel.empty()
    def ch_wakhan_genes = channel.empty()
    def ch_wakhan_rephased_vcf = channel.empty()
    if (params.wakhan) {
        def ch_wakhan_input = ch_tumour_normal
            .map { tumour_meta, tumour_bam, tumour_bai, normal_meta, _normal_bam, _normal_bai ->
                [normal_meta, tumour_meta, tumour_bam, tumour_bai]
            }
            .combine(ch_normal_phased_vcf, by: 0)
            .map { _normal_meta, tumour_meta, tumour_bam, tumour_bai, normal_vcf ->
                [tumour_meta, tumour_bam, tumour_bai, normal_vcf]
            }
            .join(SEVERUS_LOWMEM.out.somatic_vcf, by: 0)
            .map { tumour_meta, tumour_bam, tumour_bai, normal_vcf, breakpoints ->
                [tumour_meta, tumour_bam, tumour_bai, normal_vcf, [], breakpoints]
            }

        WAKHAN_ALL(
            ch_wakhan_input,
            ch_fasta_fai,
            ch_wakhan_cancer_genes,
            ch_wakhan_centromere_bed,
        )
        ch_wakhan_outdir = WAKHAN_ALL.out.outdir
        ch_wakhan_ranks = WAKHAN_ALL.out.ranks
        ch_wakhan_bed = WAKHAN_ALL.out.bed
        ch_wakhan_vcf = WAKHAN_ALL.out.vcf
        ch_wakhan_genes = WAKHAN_ALL.out.genes
        ch_wakhan_rephased_vcf = WAKHAN_ALL.out.rephased_vcf
    }

    def ch_ascat_segments = channel.empty()
    def ch_ascat_cnvs = channel.empty()
    def ch_ascat_purityploidy = channel.empty()
    def ch_ascat_purity = channel.empty()
    def ch_ascat_ploidy = channel.empty()
    def ch_ascat_metrics = channel.empty()
    def ch_ascat_bafs = channel.empty()
    def ch_ascat_logrs = channel.empty()
    def ch_ascat_png = channel.empty()
    if (params.ascat) {
        def ch_ascat_input = ch_tumour_normal.map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, normal_bam, normal_bai, tumour_bam, tumour_bai]
        }

        ASCAT(
            ch_ascat_input,
            ch_ascat_alleles,
            ch_ascat_loci,
            [],
            [],
            ch_ascat_gc,
            ch_ascat_rt,
        )
        ch_ascat_segments = ASCAT.out.segments
        ch_ascat_cnvs = ASCAT.out.cnvs
        ch_ascat_purityploidy = ASCAT.out.purityploidy

        // ASCAT writes 0 for both values when it finds no solution
        def ch_ascat_purityploidy_values = ASCAT.out.purityploidy
            .splitCsv(header: true, sep: '\t')
            .multiMap { meta, row ->
                purity: [meta, row.AberrantCellFraction as Double]
                ploidy: [meta, row.Ploidy as Double]
            }
        ch_ascat_purity = ch_ascat_purityploidy_values.purity
        ch_ascat_ploidy = ch_ascat_purityploidy_values.ploidy
        ch_ascat_metrics = ASCAT.out.metrics
        ch_ascat_bafs = ASCAT.out.bafs
        ch_ascat_logrs = ASCAT.out.logrs
        ch_ascat_png = ASCAT.out.png
    }

    emit:
    savana_vcf = SAVANA_RUN_LOWMEM.out.sv_breakpoints_vcf
    savana_bedpe = SAVANA_RUN_LOWMEM.out.sv_breakpoints_bedpe
    savana_read_support = SAVANA_RUN_LOWMEM.out.sv_breakpoints_read_support
    savana_inserted_sequences = SAVANA_RUN_LOWMEM.out.inserted_sequences
    savana_classified_vcf = SAVANA_CLASSIFY.out.classified_vcf
    savana_somatic_vcf = SAVANA_CLASSIFY.out.somatic_vcf
    savana_germline_vcf = SAVANA_CLASSIFY.out.germline_vcf
    savana_somatic_bedpe = SAVANA_CLASSIFY.out.somatic_bedpe
    savana_cna = SAVANA_CNA.out.cna
    savana_cna_segments = SAVANA_CNA.out.segmented_log2r
    savana_cna_raw_counts = SAVANA_CNA.out.raw_read_counts
    savana_cna_binned_ref = SAVANA_CNA.out.binned_ref
    savana_cna_ranked_solutions = SAVANA_CNA.out.ranked_solutions
    savana_cna_purity_ploidy = SAVANA_CNA.out.fitted_purity_ploidy
    savana_cna_no_fit = SAVANA_CNA.out.no_fit
    savana_cna_allele_counts = SAVANA_CNA.out.allele_counts
    severus_all_vcf = SEVERUS_LOWMEM.out.all_vcf
    severus_somatic_vcf = SEVERUS_LOWMEM.out.somatic_vcf
    severus_all_clusters = SEVERUS_LOWMEM.out.all_breakpoints_clusters
    severus_all_clusters_list = SEVERUS_LOWMEM.out.all_breakpoints_clusters_list
    severus_somatic_clusters = SEVERUS_LOWMEM.out.somatic_breakpoints_clusters
    severus_somatic_clusters_list = SEVERUS_LOWMEM.out.somatic_breakpoints_clusters_list
    severus_log = SEVERUS_LOWMEM.out.log
    severus_breakpoints_double = SEVERUS_LOWMEM.out.breakpoints_double
    severus_read_ids = SEVERUS_LOWMEM.out.read_ids
    severus_collapsed_dup = SEVERUS_LOWMEM.out.collapsed_dup
    severus_loh = SEVERUS_LOWMEM.out.loh
    severus_all_plots = SEVERUS_LOWMEM.out.all_plots
    severus_somatic_plots = SEVERUS_LOWMEM.out.somatic_plots
    wakhan_outdir = ch_wakhan_outdir
    wakhan_ranks = ch_wakhan_ranks
    wakhan_bed = ch_wakhan_bed
    wakhan_vcf = ch_wakhan_vcf
    wakhan_genes = ch_wakhan_genes
    wakhan_rephased_vcf = ch_wakhan_rephased_vcf
    ascat_segments = ch_ascat_segments
    ascat_cnvs = ch_ascat_cnvs
    ascat_purityploidy = ch_ascat_purityploidy
    ascat_purity = ch_ascat_purity // channel: [ val(tumour_meta), val(purity) ]
    ascat_ploidy = ch_ascat_ploidy // channel: [ val(tumour_meta), val(ploidy) ]
    ascat_metrics = ch_ascat_metrics
    ascat_bafs = ch_ascat_bafs
    ascat_logrs = ch_ascat_logrs
    ascat_png = ch_ascat_png
}
