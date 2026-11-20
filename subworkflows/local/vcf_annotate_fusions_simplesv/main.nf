// Annotate the merged Severus/SAVANA SVs with snpEff, simplify the annotation with simple_sv_annotation and list the gene fusions

include { SNPEFF_SNPEFF                          } from '../../../modules/nf-core/snpeff/snpeff/main'
include { SIMPLESVANNOTATION                     } from '../../../modules/local/simplesvannotation/main'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_FUSIONS } from '../../../modules/nf-core/bcftools/query/main'

workflow VCF_ANNOTATE_FUSIONS_SIMPLESV {

    take:
    ch_sv_vcf

    main:

    SNPEFF_SNPEFF(
        ch_sv_vcf,
        'hg38',
        [[], []],
    )

    SIMPLESVANNOTATION(SNPEFF_SNPEFF.out.vcf)

    BCFTOOLS_QUERY_FUSIONS(
        SIMPLESVANNOTATION.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        [],
        [],
        [],
    )

    emit:
    snpeff_vcf   = SNPEFF_SNPEFF.out.vcf          // channel: [ val(meta), path(vcf) ]
    simpleann_vcf = SIMPLESVANNOTATION.out.vcf    // channel: [ val(meta), path(vcf) ]
    gene_fusions = BCFTOOLS_QUERY_FUSIONS.out.output // channel: [ val(meta), path(txt) ]
}
