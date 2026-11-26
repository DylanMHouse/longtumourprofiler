/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
include { NANOPLOT } from '../modules/nf-core/nanoplot/main'

include { MOSDEPTH } from '../modules/nf-core/mosdepth/main'
include { BAM_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_stats_samtools'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PREPARE_SV } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_INPUT } from '../modules/nf-core/samtools/index/main'
include { FASTQ_ALIGN_MINIMAP2 } from '../subworkflows/local/fastq_align_minimap2'
include { BAM_CONVERT_FASTQ_QC } from '../subworkflows/local/bam_convert_fastq_qc'
include { BAM_CALL_SNV_CLAIR3 } from '../subworkflows/local/bam_call_snv_clair3'
include { SNIFFLES } from '../modules/nf-core/sniffles/main'
include { BAM_CALL_CNV_SPECTRE } from '../subworkflows/local/bam_call_cnv_spectre'
include { BAM_PHASE_HAPLOTAG_LONGPHASE } from '../subworkflows/local/bam_phase_haplotag_longphase'
include { BAM_TYPE_IMMUNE_SPECIMMUNE } from '../subworkflows/local/bam_type_immune_specimmune'
include { CNVPYTOR } from '../modules/local/cnvpytor/main'
include { STRDUST } from '../modules/nf-core/strdust/main'
include { MODKIT_PILEUP } from '../modules/nf-core/modkit/pileup/main'
include { METHYLATION_CALLING } from '../subworkflows/local/methylation_calling'
include { VCF_REPORT_CPSR } from '../subworkflows/local/vcf_report_cpsr'
include { VCF_REPORT_PCGR } from '../subworkflows/local/vcf_report_pcgr'
include { PHARMCAT_PIPELINE } from '../modules/local/pharmcat/pipeline/main'
include { HLALOH } from '../modules/local/hlaloh/main'
include { SCARHRD } from '../modules/local/scarhrd/main'
include { PYCLONEVI_PREPARE } from '../modules/local/pyclonevi/prepare/main'
include { PYCLONEVI_FIT } from '../modules/local/pyclonevi/fit/main'
include { PHYCLONE_RUN } from '../modules/local/phyclone/run/main'
include { PHYCLONE_MAP } from '../modules/local/phyclone/map/main'
include { PHYCLONE_TOPOLOGYREPORT } from '../modules/local/phyclone/topologyreport/main'
include { PHYCLONE_PLOTTREE } from '../modules/local/phyclone/plottree/main'
include { VCF_NEOANTIGEN_VEP_PVACSEQ } from '../subworkflows/local/vcf_neoantigen_vep_pvacseq/main'

include { BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC } from '../subworkflows/local/bam_call_somatic_clairs_deepsomatic/main'
include { BAM_CALL_SV_CNA_TN } from '../subworkflows/local/bam_call_sv_cna_tn/main'
include { SOMALIER_RAITA_IDENTITY } from '../subworkflows/local/somalier_raita_identity/main'
include { ONTCONSENSUS } from '../subworkflows/local/ontconsensus/main'
include { MUTATIONAL_SIGNATURES_SIGPROFILER } from '../subworkflows/local/mutational_signatures_sigprofiler/main'
include { ANNOTSV_SURVIVOR_ANNOTATE } from '../subworkflows/local/annotsv_survivor_annotate/main'
include { VCF_ANNOTATE_FUSIONS_SIMPLESV } from '../subworkflows/local/vcf_annotate_fusions_simplesv/main'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow LONGTUMOURPROFILER {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    outdir

    main:



    def ch_versions = channel.empty()
    def ch_input_bam = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_' + 'longtumourprofiler_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCES & OPTIONAL INPUT CONFIGURATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // --fasta takes priority over --genome if provided of course.
    def fasta_path = params.fasta ?: params.genome
        ? params.genomes.containsKey(params.genome)
            ? params.genomes[params.genome].fasta
            : error("Genome '${params.genome}' not found in iGenomes config. Check --genome or set --igenomes_ignore false")
        : error("Provide a reference with --fasta or --genome")

    def ch_fasta = channel.fromPath(fasta_path, checkIfExists: true)
        .map { fasta -> tuple([id: "ref"], fasta) }
        .first()

    // Indexing the reference to .fai is necessary here as it is mandatory for some steps.
    SAMTOOLS_FAIDX(ch_fasta.map { meta, fasta -> [meta, fasta, []] }, false)
    def ch_fasta_fai = ch_fasta
        .join(SAMTOOLS_FAIDX.out.fai, by: 0)
        .first()

    // Optional inputs: an empty list stands in for a missing file.
    def ch_intervals = params.intervals
        ? channel.fromPath(params.intervals, checkIfExists: true).first()
        : channel.value([])

    // minimap2 accepts a prebuilt .mmi in place of the fasta; parabricks always needs the fasta
    def ch_minimap2_ref = params.minimap2_index && params.alignment_software == 'minimap_nogpu'
        ? channel.fromPath(params.minimap2_index, checkIfExists: true).map { mmi -> tuple([id: "ref"], mmi) }.first()
        : ch_fasta


    // Clair3 SNVs feed CNVpytor, longphase and CPSR; Spectre is skipped without them
    if (!params.clair3) {
        def needs_snvs = ['cnvpytor', 'phasing', 'cpsr'].findAll { step -> params[step] }
        if (needs_snvs) {
            error("--${needs_snvs.join(', --')} need the Clair3 SNV calls, so they need --clair3")
        }
    }
    if (!params.minimap2 && params.input_type != 'aligned_bam') {
        error("--minimap2 false needs --input_type aligned_bam")
    }

    // Clair3 model: a packaged model name (in the container) OR a user model directory, not both
    if (params.clair3) {
        if (params.clair3_model && params.clair3_home_model) {
            error("Use either --clair3_model (packaged model name) or --clair3_home_model (model directory), not both")
        }
        if (!params.clair3_model && !params.clair3_home_model) {
            error("Clair3 needs a model: set --clair3_model or --clair3_home_model")
        }
    }
    // channel: [ packaged_model, user_model, platform ]
    def ch_clair3_model = params.clair3_home_model
        ? channel.fromPath(params.clair3_home_model, type: 'dir', checkIfExists: true).map { model_dir -> [[], model_dir, 'ont'] }.first()
        : channel.value([params.clair3_model, [], 'ont'])

    def specimmune_genes = params.specimmune_genes.tokenize(',')*.trim()
    if (params.specimmune) {
        if (!params.specimmune_db) {
            error("--specimmune needs --specimmune_db, the folder built by SpecImmune's make_db.py")
        }
        def unknown_genes = specimmune_genes - ['HLA', 'KIR', 'CYP', 'IG_TR']
        if (unknown_genes) {
            error("--specimmune_genes takes HLA, KIR, CYP and IG_TR, not: ${unknown_genes.join(', ')}")
        }
        if ('IG_TR' in specimmune_genes && !params.specimmune_hg38) {
            error("IG_TR typing needs --specimmune_hg38, a no-alt hg38 fasta (bwa-indexed, with .fai)")
        }
    }

    def ch_specimmune_db = params.specimmune_db
        ? channel.fromPath(params.specimmune_db, type: 'dir', checkIfExists: true).first()
        : channel.value([])

    // Stage the fasta with its .fai and bwa index files, otherwise SpecImmune rebuilds the index in every task
    def ch_specimmune_hg38 = params.specimmune_hg38
        ? channel.value(['', '.fai', '.amb', '.ann', '.bwt', '.pac', '.sa'].collect { suffix -> file("${params.specimmune_hg38}${suffix}", checkIfExists: true) })
        : channel.value([])

    // STRdust optional BED file of STR regions; the module takes [ meta, bed ]
    def ch_str_bed = params.strbed
        ? channel.fromPath(params.strbed, checkIfExists: true).map { bed -> tuple([id: "strbed"], bed) }.first()
        : channel.value([[:], []])

    def ch_methyl_bed = params.methylbed
        ? channel.fromPath(params.methylbed, checkIfExists: true).map { bed -> tuple([id: "methylbed"], bed) }.first()
        : channel.value([[:], []])

    if (params.nanodx) {
        if (!params.methyl) {
            error("--nanodx classifies the modkit pileup, so it needs --methyl")
        }
        if (!(params.nanodx_model && params.nanodx_probes && params.nanodx_chain)) {
            error("--nanodx needs --nanodx_model, --nanodx_probes and --nanodx_chain")
        }
    }

    def ch_nanodx_model = params.nanodx_model
        ? channel.fromPath(params.nanodx_model, checkIfExists: true).map { model -> tuple([id: "nanodx_model"], model) }.first()
        : channel.value([[:], []])

    def ch_nanodx_probes = params.nanodx_probes
        ? channel.fromPath(params.nanodx_probes, checkIfExists: true).map { bed -> tuple([id: "nanodx_probes"], bed) }.first()
        : channel.value([[:], []])

    def ch_nanodx_chain = params.nanodx_chain
        ? channel.fromPath(params.nanodx_chain, checkIfExists: true).map { chain -> tuple([id: "nanodx_chain"], chain) }.first()
        : channel.value([[:], []])

    def ch_spectre_metadata = params.spectre_metadata
        ? channel.fromPath(params.spectre_metadata, checkIfExists: true).first()
        : channel.value([])

    def ch_spectre_blacklist = params.spectre_blacklist
        ? channel.fromPath(params.spectre_blacklist, checkIfExists: true).first()
        : channel.value([])

    def ch_tandem_repeats = params.sniffles_tandem_repeats
        ? channel.fromPath(params.sniffles_tandem_repeats, checkIfExists: true).map { bed -> tuple([id: "tandem_repeats"], bed) }.first()
        : channel.value([[:], []])

    def ch_cpsr_regions = params.cpsr_regions
        ? channel.fromPath(params.cpsr_regions, checkIfExists: true).first()
        : channel.value([])

    // CPSR needs the PCGR reference bundle and a VEP cache
    if (params.cpsr && !(params.cpsr_refdata && params.cpsr_vep_cache)) {
        error("--cpsr needs --cpsr_refdata (PCGR reference data) and --cpsr_vep_cache (VEP cache directory)")
    }
    if (params.pcgr && !(params.pcgr_refdata && params.pcgr_vep_cache)) {
        error("--pcgr needs --pcgr_refdata (PCGR reference data) and --pcgr_vep_cache (VEP cache directory)")
    }
    // PCGR is run on the exome-cut consensus so the TMB and MSI estimates match --assay WES
    if (params.pcgr && !params.exome_bed) {
        error("--pcgr reports the exome-cut consensus VCF, so it needs --exome_bed")
    }

    if (params.annotsv && !params.annotsv_annotations) {
        error("--annotsv needs --annotsv_annotations (AnnotSV annotations directory)")
    }

    if (params.wakhan && !params.phasing) {
        error("--wakhan needs the LongPhase phased germline SNVs, so it needs --phasing")
    }
    if (params.pharmcat && !params.phasing) {
        error("--pharmcat needs the LongPhase phased germline SNVs, so it needs --phasing")
    }
    if (params.pharmcat && params.genome == 'GRCh37') {
        error("--pharmcat only supports GRCh38")
    }
    if (params.pvacseq) {
        if (!(params.specimmune && 'HLA' in specimmune_genes)) {
            error("--pvacseq takes its HLA alleles from SpecImmune, so it needs --specimmune with HLA in --specimmune_genes")
        }
        if (!(params.pvacseq_vep_cache && params.pvacseq_vep_cache_version && params.pvacseq_vep_plugins)) {
            error("--pvacseq needs --pvacseq_vep_cache, --pvacseq_vep_cache_version and --pvacseq_vep_plugins")
        }
    }
    if (params.hla_loh && !params.phasing) {
        error("--hla_loh needs the LongPhase phased germline SNVs, so it needs --phasing")
    }
    // The built-in HLA gene windows are GRCh38 coordinates
    if (params.hla_loh && params.genome == 'GRCh37' && !params.hla_loh_regions) {
        error("--hla_loh has GRCh38 HLA coordinates built in, so GRCh37 needs --hla_loh_regions")
    }
    if (params.ascat && !(params.ascat_alleles && params.ascat_loci)) {
        error("--ascat needs --ascat_alleles and --ascat_loci")
    }
    if (params.scarhrd && !params.ascat) {
        error("--scarhrd scores the ASCAT allele-specific copy number segments, so it needs --ascat")
    }
    if (params.pyclonevi && !params.ascat) {
        error("--pyclonevi takes its copy number and purity from ASCAT, so it needs --ascat")
    }
    if (params.phyclone && !params.pyclonevi) {
        error("--phyclone builds trees from the PyClone-VI clusters, so it needs --pyclonevi")
    }

    if ((params.somalier_sites as Boolean) != (params.raita_bed as Boolean)) {
        error("The tumour-normal identity check needs both --somalier_sites and --raita_bed")
    }

    def ch_somalier_sites = params.somalier_sites
        ? channel.fromPath(params.somalier_sites, checkIfExists: true).map { sites -> tuple([id: "somalier_sites"], sites) }.first()
        : channel.value([[:], []])

    def ch_raita_bed = params.raita_bed
        ? channel.fromPath(params.raita_bed, checkIfExists: true).first()
        : channel.value([])

    def ch_wakhan_cancer_genes = params.wakhan_cancer_genes
        ? channel.fromPath(params.wakhan_cancer_genes, checkIfExists: true).map { genes -> tuple([id: "cancer_genes"], genes) }.first()
        : channel.value([[:], []])

    def ch_wakhan_centromere_bed = params.wakhan_centromere_bed
        ? channel.fromPath(params.wakhan_centromere_bed, checkIfExists: true).map { bed -> tuple([id: "centromeres"], bed) }.first()
        : channel.value([[:], []])

    // ASCAT takes each reference as a directory or a single file
    def ch_ascat_alleles = params.ascat_alleles
        ? channel.fromPath(params.ascat_alleles, type: 'any', checkIfExists: true).first()
        : channel.value([])

    def ch_ascat_loci = params.ascat_loci
        ? channel.fromPath(params.ascat_loci, type: 'any', checkIfExists: true).first()
        : channel.value([])

    def ch_ascat_gc = params.ascat_gc
        ? channel.fromPath(params.ascat_gc, type: 'any', checkIfExists: true).first()
        : channel.value([])

    def ch_pvacseq_vep_cache = params.pvacseq_vep_cache
        ? channel.fromPath(params.pvacseq_vep_cache, type: 'dir', checkIfExists: true).map { cache -> tuple([id: "vep_cache"], cache) }.first()
        : channel.value([[:], []])

    // Folder holding Frameshift.pm and Wildtype.pm
    def ch_pvacseq_vep_plugins = params.pvacseq_vep_plugins
        ? channel.fromPath(params.pvacseq_vep_plugins, type: 'dir', checkIfExists: true).first()
        : channel.value([])

    def ch_hla_loh_regions = params.hla_loh_regions
        ? channel.fromPath(params.hla_loh_regions, checkIfExists: true).map { bed -> tuple([id: "hla_loh_regions"], bed) }.first()
        : channel.value([[:], []])

    def ch_ascat_rt = params.ascat_rt
        ? channel.fromPath(params.ascat_rt, type: 'any', checkIfExists: true).first()
        : channel.value([])




    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    READ QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.qc && params.input_type == 'fastq') {
        NANOPLOT(ch_samplesheet)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDEXING,  ALIGNMENT & QC CHECKS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.input_type == 'aligned_bam') {
        // Reuse an existing .bai next to the BAM, otherwise index it
        def ch_bam_index = ch_samplesheet
            .map { meta, bams -> [meta, bams[0]] }
            .branch { _meta, bam ->
                indexed: file("${bam}.bai").exists()
                unindexed: true
            }
        SAMTOOLS_INDEX_INPUT(ch_bam_index.unindexed)
        ch_input_bam = ch_bam_index.indexed
            .map { meta, bam -> [meta, bam, file("${bam}.bai")] }
            .mix(ch_bam_index.unindexed.join(SAMTOOLS_INDEX_INPUT.out.index, by: 0))
    }
    else if (params.minimap2) {

        if (params.input_type == 'fastq') {

            FASTQ_ALIGN_MINIMAP2(ch_samplesheet, ch_minimap2_ref)
        }
        else if (params.input_type == 'bam' || params.input_type == 'single_ubam') {
            BAM_CONVERT_FASTQ_QC(ch_samplesheet)

            FASTQ_ALIGN_MINIMAP2(BAM_CONVERT_FASTQ_QC.out.fastq, ch_minimap2_ref)
        }

        // Keep the full meta (id, patient, type) so tumour and germline can be told apart downstream
        ch_input_bam = FASTQ_ALIGN_MINIMAP2.out.bam.join(FASTQ_ALIGN_MINIMAP2.out.bai, by: 0)
    }

    // Generate some statistics of aligned reads
    BAM_STATS_SAMTOOLS(
        ch_input_bam,
        ch_fasta_fai,
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COVERAGE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.coverage) {
        MOSDEPTH(
            // wrap the bed so an empty [] still fills the bed slot: [ meta, bam, bai, bed ]
            ch_input_bam.combine(ch_intervals.map { bed -> [bed] }),
            ch_fasta,
            [],
        )
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SNV calls
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    def ch_snv_vcf = channel.empty()
    if (params.clair3) {
        BAM_CALL_SNV_CLAIR3(
            ch_input_bam.combine(ch_clair3_model),
            ch_fasta,
            SAMTOOLS_FAIDX.out.fai.first(),
        )
        ch_snv_vcf = BAM_CALL_SNV_CLAIR3.out.vcf
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GERMLINE VARIANT REPORTING WITH CPSR
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // CPSR reports cancer predisposition variants, so it only takes the germline calls
    def ch_cpsr_classification = channel.empty()
    def ch_cpsr_conf_yaml = channel.empty()
    if (params.cpsr) {
        VCF_REPORT_CPSR(
            ch_snv_vcf.filter { meta, _vcf -> meta.type == 'germline' },
            channel.fromPath(params.cpsr_refdata, type: 'dir', checkIfExists: true).first(),
            channel.fromPath(params.cpsr_vep_cache, type: 'dir', checkIfExists: true).first(),
            ch_cpsr_regions,
        )
        ch_cpsr_classification = VCF_REPORT_CPSR.out.classification
        ch_cpsr_conf_yaml = VCF_REPORT_CPSR.out.conf_yaml
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SV CALLS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Sniffles runs on the unphased BAM because LongPhase co-phases its SVs
    SNIFFLES(
        ch_input_bam.map { meta, bam, bai -> [meta, bam, bai, [], []] },
        ch_fasta,
        ch_tandem_repeats,
        true,
        false,
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHASING & HAPLOTAGGING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Every step from here on reads the haplotagged BAM, or the alignment BAM with --phasing false
    def ch_phased_bam = ch_input_bam
    def ch_phased_snv_vcf = channel.empty()
    if (params.phasing) {
        // Co-phase each sample's SNVs and Sniffles SVs (and 5mC calls with --longphase_modcall)
        def ch_phase_input = ch_input_bam
            .join(ch_snv_vcf, by: 0)
            .join(SNIFFLES.out.vcf, by: 0)

        BAM_PHASE_HAPLOTAG_LONGPHASE(
            ch_phase_input,
            ch_fasta,
            SAMTOOLS_FAIDX.out.fai.first(),
            params.longphase_modcall,
        )

        ch_phased_bam = BAM_PHASE_HAPLOTAG_LONGPHASE.out.bam.join(BAM_PHASE_HAPLOTAG_LONGPHASE.out.bai, by: 0)
        ch_phased_snv_vcf = BAM_PHASE_HAPLOTAG_LONGPHASE.out.snv_vcf
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNV CALLS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Spectre uses the Clair3 SNVs of the same sample to estimate diploid coverage
    BAM_CALL_CNV_SPECTRE(
        ch_phased_bam,
        ch_snv_vcf,
        ch_fasta,
        ch_spectre_metadata,
        ch_spectre_blacklist,
    )

    // CNVpytor combines read depth with B-allele frequencies from the Clair3 SNVs
    if (params.cnvpytor) {
        CNVPYTOR(
            ch_phased_bam.join(ch_snv_vcf, by: 0),
            ch_fasta,
            params.cnvpytor_bin_size,
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHARMACOGENOMICS WITH PHARMCAT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // PharmCAT calls star alleles from the phased germline SNVs
    if (params.pharmcat) {
        PHARMCAT_PIPELINE(
            ch_phased_snv_vcf.filter { meta, _vcf -> meta.type == 'germline' }
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMMUNE GENE TYPING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    def ch_specimmune_outdir = channel.empty()
    if (params.specimmune) {
        BAM_TYPE_IMMUNE_SPECIMMUNE(
            ch_phased_bam,
            ch_fasta,
            channel.fromList(specimmune_genes),
            ch_specimmune_db,
            ch_specimmune_hg38,
        )
        ch_specimmune_outdir = BAM_TYPE_IMMUNE_SPECIMMUNE.out.outdir
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    METHYLATION CALLS WITH MODKIT_PILEUP, LIFTOVER, NANODX, AND BEDTOOLS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Reads keep their MM/ML tags through SAMTOOLS_FASTQ (-T MM,ML,mv) and minimap2 (-y)
    if (params.methyl) {
        MODKIT_PILEUP(
            ch_phased_bam,
            ch_fasta_fai,
            ch_methyl_bed,
        )

        // nanoDx models classify tumours, so the germline pileup is left out
        if (params.nanodx) {
            METHYLATION_CALLING(
                MODKIT_PILEUP.out.bedgz.filter { meta, _bed -> meta.type == 'tumour' },
                ch_nanodx_chain,
                ch_nanodx_probes,
                ch_nanodx_model,
            )
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STR ANALYSIS with STRDust
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Runs on the haplotagged BAM when --phasing, as in the shell workflow
    if (params.strdust) {
        STRDUST(
            ch_phased_bam,
            ch_fasta,
            SAMTOOLS_FAIDX.out.fai.first(),
            ch_str_bed,
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Tumour-Normal profiler Workflow.
    Note: If the user wishes to run only profiling, please ensure you have phased BAM, coverage


    SOMATIC CALLS WITH CLAIRS, DEEPSOMATIC AND DEEPVARIANT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    def ch_samples_by_type = ch_phased_bam.branch { meta, _bam, _bai ->
        tumour: meta.type == 'tumour'
        normal: meta.type == 'germline'
    }

    // Create two independent channels to separate normal vs tumour.
    def ch_tumour_keyed = ch_samples_by_type.tumour.map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
    def ch_normal_keyed = ch_samples_by_type.normal.map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

    // Pair each tumour with its germline sample: [ tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai ]
    def ch_tumour_normal_input = ch_tumour_keyed
        .join(ch_normal_keyed)
        .map { _patient, tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai]
        }

    BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC(
        ch_tumour_normal_input,
        ch_fasta_fai,
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NEOANTIGEN PREDICTION WITH PVACSEQ
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // The germline HLA type is used because the tumour can lose an HLA haplotype
    if (params.pvacseq) {
        VCF_NEOANTIGEN_VEP_PVACSEQ(
            BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC.out.clairs_vcf,
            ch_specimmune_outdir.filter { meta, gene_class, _hla_dir -> gene_class == 'HLA' && meta.type == 'germline' }.map { meta, _gene_class, hla_dir -> [meta, hla_dir] },
            ch_fasta,
            params.genome == 'GRCh37' ? 'GRCh37' : 'GRCh38',
            params.pvacseq_vep_cache_version,
            ch_pvacseq_vep_cache,
            ch_pvacseq_vep_plugins,
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TUMOUR-NORMAL IDENTITY CHECK WITH SOMALIER AND RAITA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.somalier_sites && params.raita_bed) {
        SOMALIER_RAITA_IDENTITY(
            ch_tumour_normal_input,
            ch_fasta_fai,
            ch_somalier_sites,
            ch_raita_bed,
        )
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    STRUCTURAL VARIANTS AND COPY NUMBER ALTERATIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Preparing BAM input for SV calls. Removes tags
    SAMTOOLS_VIEW_PREPARE_SV(
        ch_phased_bam,
        ch_fasta_fai,
        [[], []],
        [[], []],
        'bai',
    )

    def ch_prepare_bamSV = SAMTOOLS_VIEW_PREPARE_SV.out.bam.join(SAMTOOLS_VIEW_PREPARE_SV.out.bai, by: 0)

    def ch_prepare_bamSV_by_type = ch_prepare_bamSV.branch { meta, _bam, _bai ->
        tumour: meta.type == 'tumour'
        normal: meta.type == 'germline'
    }

    def ch_prepare_bamSV_tumour = ch_prepare_bamSV_by_type.tumour.map { meta, bam, bai -> [meta.patient, meta, bam, bai] }
    def ch_prepare_bamSV_normal = ch_prepare_bamSV_by_type.normal.map { meta, bam, bai -> [meta.patient, meta, bam, bai] }

    def ch_tumour_normal_forSV = ch_prepare_bamSV_tumour
        .join(ch_prepare_bamSV_normal)
        .map { _patient, tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai ->
            [tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai]
        }

    BAM_CALL_SV_CNA_TN(
        ch_tumour_normal_forSV,
        ch_fasta_fai,
        ch_phased_snv_vcf.filter { meta, _vcf -> meta.type == 'germline' },
        ch_wakhan_cancer_genes,
        ch_wakhan_centromere_bed,
        ch_ascat_alleles,
        ch_ascat_loci,
        ch_ascat_gc,
        ch_ascat_rt,
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC SV MERGE (SURVIVOR) AND ANNOTATION (ANNOTSV)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.annotsv) {
        ANNOTSV_SURVIVOR_ANNOTATE(
            BAM_CALL_SV_CNA_TN.out.savana_somatic_vcf.join(BAM_CALL_SV_CNA_TN.out.severus_somatic_vcf, by: 0),
            channel.fromPath(params.annotsv_annotations, type: 'dir', checkIfExists: true).map { annotations -> [[id: 'annotsv_annotations'], annotations] }.first(),
        )

        // Gene fusion calls from the merged SVs: snpEff + simple_sv_annotation
        VCF_ANNOTATE_FUSIONS_SIMPLESV(ANNOTSV_SURVIVOR_ANNOTATE.out.merged_vcf, params.snpeff_db)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HLA LOSS OF HETEROZYGOSITY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.hla_loh) {
        def ch_hla_loh_pairs = ch_tumour_normal_input
            .map { tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai ->
                [normal_meta, tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai]
            }
            .combine(ch_phased_snv_vcf.filter { meta, _vcf -> meta.type == 'germline' }, by: 0)
            .map { _normal_meta, tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, normal_vcf ->
                [tumour_meta, tumour_bam, tumour_bai, normal_bam, normal_bai, normal_vcf]
            }


        def ch_hla_loh_input = params.ascat
            ? ch_hla_loh_pairs.join(BAM_CALL_SV_CNA_TN.out.ascat_purity, by: 0).join(BAM_CALL_SV_CNA_TN.out.ascat_ploidy, by: 0)
            : ch_hla_loh_pairs.map { pair -> pair + [[], []] }

        HLALOH(
            ch_hla_loh_input,
            ch_hla_loh_regions,
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCARHRD: HOMOLOGOUS RECOMBINATION DEFICIENCY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.scarhrd) {
        // scarHRD scores the ASCAT allele-specific segments using the ASCAT ploidy
        SCARHRD(
            BAM_CALL_SV_CNA_TN.out.ascat_segments.join(BAM_CALL_SV_CNA_TN.out.ascat_purityploidy, by: 0)
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ONTCONSENSUS: GATHER CONSENSUS VARIANTS FROM CLAIRS AND DEEPSOMATIC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    // Both callers emit with the tumour meta: [ meta, deepsomatic_vcf, deepsomatic_tbi, clairs_vcf, clairs_tbi ]
    def ch_ontconsenus = BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC.out.deepsomatic_vcf
        .join(BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC.out.deepsomatic_tbi, by: 0)
        .join(BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC.out.clairs_vcf, by: 0)
        .join(BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC.out.clairs_tbi, by: 0)

    // ont_somatic_consensus.sh exits without a mask, so fail before any calling starts
    if (!params.consensus_mask) {
        error("The ClairS/DeepSomatic consensus needs --consensus_mask, a homopolymer and STR BED")
    }
    def ch_mask = channel.fromPath(params.consensus_mask, checkIfExists: true)
        .map { mask -> tuple([id: "mask"], mask) }
        .first()

    def ch_exome_bed = params.exome_bed
        ? channel.fromPath(params.exome_bed, checkIfExists: true).first()
        : channel.empty()

    ONTCONSENSUS(
        ch_ontconsenus,
        ch_fasta_fai,
        ch_mask,
        ch_exome_bed,
        params.genome.toLowerCase(),
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLONAL DECONVOLUTION WITH PYCLONE-VI
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.pyclonevi) {
        def ch_ascat_solved = BAM_CALL_SV_CNA_TN.out.ascat_purity.filter { meta, purity ->
            if (purity <= 0) {
                log.warn("Skipping PyClone-VI for ${meta.id}: ASCAT found no purity solution")
            }
            purity > 0
        }

        PYCLONEVI_PREPARE(
            ONTCONSENSUS.out.consensus_vcf.join(BAM_CALL_SV_CNA_TN.out.ascat_segments, by: 0).join(BAM_CALL_SV_CNA_TN.out.ascat_purityploidy, by: 0).join(ch_ascat_solved, by: 0).map { meta, vcf, segments, purityploidy, _purity -> [meta, vcf, segments, purityploidy] }
        )

        PYCLONEVI_FIT(
            PYCLONEVI_PREPARE.out.tsv.filter { meta, tsv ->
                def has_mutations = tsv.countLines() > 1
                if (!has_mutations) {
                    log.warn("Skipping PyClone-VI for ${meta.id}: no SNVs left on ASCAT copy number segments")
                }
                has_mutations
            }
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLONE TREES WITH PHYCLONE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.phyclone) {
        PHYCLONE_RUN(
            PYCLONEVI_PREPARE.out.tsv.join(PYCLONEVI_FIT.out.results, by: 0)
        )

        PHYCLONE_MAP(PHYCLONE_RUN.out.trace)
        PHYCLONE_TOPOLOGYREPORT(PHYCLONE_RUN.out.trace)

        PHYCLONE_PLOTTREE(
            PHYCLONE_MAP.out.tree.join(PHYCLONE_MAP.out.table, by: 0)
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC VARIANT REPORTING WITH PCGR
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.pcgr) {
        // CPSR runs on the germline sample, so its results are rekeyed to the paired tumour
        def ch_cpsr_by_tumour = ch_tumour_normal_input
            .map { tumour_meta, _tumour_bam, _tumour_bai, normal_meta, _normal_bam, _normal_bai -> [normal_meta, tumour_meta] }
            .combine(ch_cpsr_classification.join(ch_cpsr_conf_yaml, by: 0), by: 0)
            .map { _normal_meta, tumour_meta, classification, conf_yaml -> [tumour_meta, classification, conf_yaml] }

        VCF_REPORT_PCGR(
            ONTCONSENSUS.out.exome_vcf.join(ONTCONSENSUS.out.exome_tbi, by: 0),
            params.ascat ? BAM_CALL_SV_CNA_TN.out.ascat_segments : channel.empty(),
            params.ascat ? BAM_CALL_SV_CNA_TN.out.ascat_purity : channel.empty(),
            params.ascat ? BAM_CALL_SV_CNA_TN.out.ascat_ploidy : channel.empty(),
            ch_cpsr_by_tumour,
            channel.fromPath(params.pcgr_refdata, type: 'dir', checkIfExists: true).first(),
            channel.fromPath(params.pcgr_vep_cache, type: 'dir', checkIfExists: true).first(),
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MUTATIONAL, COPY NUMBER AND SV SIGNATURES WITH SIGPROFILER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (params.sigprofiler) {
        MUTATIONAL_SIGNATURES_SIGPROFILER(
            params.sigprofiler_genome ? BAM_CALL_SOMATIC_CLAIRS_DEEPSOMATIC.out.clairs_vcf : channel.empty(),
            params.ascat ? BAM_CALL_SV_CNA_TN.out.ascat_segments : channel.empty(),
            BAM_CALL_SV_CNA_TN.out.savana_somatic_bedpe,
            params.sigprofiler_genome ? channel.fromPath(params.sigprofiler_genome, type: 'dir', checkIfExists: true).first() : channel.value([]),
            params.sv_signatures ? channel.fromPath(params.sv_signatures, checkIfExists: true).first() : channel.value([]),
        )
    }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
