include { ONTSOMATICCONSENSUS } from '../../../modules/local/ontsomaticconsensus/main'
include { MAFSMITH } from '../../../modules/local/mafsmith/main'

workflow ONTCONSENSUS {

    take:
    deepsomatic_clairs_vcfs
    fasta_fai
    mask

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


    emit:
    consensus_vcf = consensus_vcf
    consensus_maf = MAFSMITH.out.maf
}
