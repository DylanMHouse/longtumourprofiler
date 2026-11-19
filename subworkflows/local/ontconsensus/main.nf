include { ONTSOMATICCONSENSUS } from '../../../modules/local/ontsomaticconsensus/main'
include { MAFSMITH } from '../../../modules/local/mafsmith/main'
include { MAFSMITH as MAFSMITH_EXOME } from '../../../modules/local/mafsmith/main'
include { TABIX_REGIONS as TABIX_REGIONS_EXOME } from '../../../modules/local/tabix/regions/main'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_EXOME } from '../../../modules/nf-core/bcftools/sort/main'

workflow ONTCONSENSUS {

    take:
    deepsomatic_clairs_vcfs // channel: [ val(meta), path(deepsomatic_vcf), path(deepsomatic_tbi), path(clairs_vcf), path(clairs_tbi) ]
    fasta_fai               // channel: [ val(meta), path(fasta), path(fai) ]
    mask                    // channel: [ val(meta), path(bed) ]
    exome_bed               // channel: path(bed), or empty to skip the exome VCF and MAF

    main:

    ONTSOMATICCONSENSUS(
        deepsomatic_clairs_vcfs,
        fasta_fai,
        mask
    )

    consensus_vcf = ONTSOMATICCONSENSUS.out.vcf

    MAFSMITH(
        consensus_vcf,
        fasta_fai
    )

    TABIX_REGIONS_EXOME(
        consensus_vcf.join(ONTSOMATICCONSENSUS.out.tbi, by: 0),
        exome_bed
    )

    BCFTOOLS_SORT_EXOME(TABIX_REGIONS_EXOME.out.vcf)

    MAFSMITH_EXOME(
        BCFTOOLS_SORT_EXOME.out.vcf,
        fasta_fai
    )

    emit:
    consensus_vcf = consensus_vcf
    consensus_maf = MAFSMITH.out.maf
    exome_vcf     = BCFTOOLS_SORT_EXOME.out.vcf   // channel: [ val(meta), path(vcf.gz) ]
    exome_tbi     = BCFTOOLS_SORT_EXOME.out.index // channel: [ val(meta), path(vcf.gz.tbi) ]
    exome_maf     = MAFSMITH_EXOME.out.maf
}
