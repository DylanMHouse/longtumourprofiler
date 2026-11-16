// Predict neoantigens with pVACseq from somatic calls annotated by VEP and the HLA alleles typed by SpecImmune

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'
include { PVACSEQ_RUN as PVACSEQ_RUN_CLASSI } from '../../../modules/local/pvacseq/run/main'
include { PVACSEQ_RUN as PVACSEQ_RUN_CLASSII } from '../../../modules/local/pvacseq/run/main'

workflow VCF_NEOANTIGEN_VEP_PVACSEQ {
    take:
    ch_vcf // channel: [ val(tumour_meta), path(vcf) ] somatic calls, the tumour sample column is named meta.id
    ch_hla_dir // channel: [ val(meta), path(dir) ] SpecImmune HLA output folder of the sample to take the alleles from
    ch_fasta // channel: [ val(meta), path(fasta) ]
    vep_genome // val: GRCh38 or GRCh37
    vep_cache_version // val: version of the VEP cache
    ch_vep_cache // channel: [ val(meta), path(cache) ]
    ch_vep_plugins // channel: path(dir) folder holding Frameshift.pm and Wildtype.pm

    main:

    // pVACseq reads the wildtype and frameshift protein sequences that the VEP plugins write
    ENSEMBLVEP_VEP(
        ch_vcf.map { meta, vcf -> [meta, vcf, []] },
        vep_genome,
        'homo_sapiens',
        vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_plugins,
        [[], []],
    )

    def ch_hla_alleles = ch_hla_dir.map { meta, hla_dir ->
        def alleles = specimmuneHlaAlleles(hla_dir.resolve("${meta.id}/${meta.id}.HLA.final.type.1.result.txt"))
        [meta.patient, alleles.class_i, alleles.class_ii]
    }

    // Pair each tumour VCF with the alleles of the same patient
    def ch_pvacseq_input = ENSEMBLVEP_VEP.out.vcf
        .map { meta, vcf -> [meta.patient, meta, vcf] }
        .combine(ch_hla_alleles, by: 0)
        .multiMap { _patient, meta, vcf, class_i, class_ii ->
            class_i: [meta, vcf, class_i]
            class_ii: [meta, vcf, class_ii]
        }

    // A locus SpecImmune could not type leaves that class without alleles
    PVACSEQ_RUN_CLASSI(
        ch_pvacseq_input.class_i.filter { _meta, _vcf, alleles -> alleles },
        'NetMHCpanEL NetMHCpan MHCflurry',
    )
    PVACSEQ_RUN_CLASSII(
        ch_pvacseq_input.class_ii.filter { _meta, _vcf, alleles -> alleles },
        'NetMHCIIpanEL NetMHCIIpan',
    )

    emit:
    vep_vcf = ENSEMBLVEP_VEP.out.vcf // channel: [ val(meta), path(vcf.gz) ]
    vep_tbi = ENSEMBLVEP_VEP.out.tbi // channel: [ val(meta), path(vcf.gz.tbi) ]
    class_i_outdir = PVACSEQ_RUN_CLASSI.out.outdir // channel: [ val(meta), path(MHC_Class_I) ]
    class_i_all_epitopes = PVACSEQ_RUN_CLASSI.out.all_epitopes // channel: [ val(meta), path(tsv) ]
    class_i_aggregated = PVACSEQ_RUN_CLASSI.out.aggregated // channel: [ val(meta), path(tsv) ]
    class_i_filtered = PVACSEQ_RUN_CLASSI.out.filtered // channel: [ val(meta), path(tsv) ]
    class_ii_outdir = PVACSEQ_RUN_CLASSII.out.outdir // channel: [ val(meta), path(MHC_Class_II) ]
    class_ii_all_epitopes = PVACSEQ_RUN_CLASSII.out.all_epitopes // channel: [ val(meta), path(tsv) ]
    class_ii_aggregated = PVACSEQ_RUN_CLASSII.out.aggregated // channel: [ val(meta), path(tsv) ]
    class_ii_filtered = PVACSEQ_RUN_CLASSII.out.filtered // channel: [ val(meta), path(tsv) ]
}

// Two-field HLA alleles from SpecImmune's HLA result table (locus in column 1, allele in column 7).
// pVACseq writes class I alleles as HLA-A*02:01 and class II DRB1 alleles as DRB1*15:01.
def specimmuneHlaAlleles(result) {
    def typed = result.readLines()
        .collect { line -> line.split('\t') }
        .findAll { fields -> fields.size() > 6 && fields[6] != 'NA' }
    def class_i = typed
        .findAll { fields -> fields[0] ==~ /HLA-[ABC]/ }
        .collect { fields -> hlaTwoField(fields[6]) }
        .unique()
    def class_ii = typed
        .findAll { fields -> fields[0] == 'HLA-DRB1' }
        .collect { fields -> hlaTwoField(fields[6]).replaceFirst(/^HLA-/, '') }
        .unique()
    return [class_i: class_i.join(','), class_ii: class_ii.join(',')]
}

def hlaTwoField(allele) {
    return allele.split(':').take(2).join(':')
}
