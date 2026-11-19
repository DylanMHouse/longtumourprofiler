include { SIGPROFILERASSIGNMENT } from '../../../modules/local/sigprofilerassignment/main'
include { SVSIGNATURES } from '../../../modules/local/svsignatures/main'

workflow MUTATIONAL_SIGNATURES_SIGPROFILER {
    take:
    ch_vcf // channel: [ val(meta), path(vcf) ] somatic SNVs and indels, or empty to skip the SBS refit
    ch_segments // channel: [ val(meta), path(segments) ] ASCAT segments, or empty to skip the CN refit
    ch_bedpe // channel: [ val(meta), path(bedpe) ] SAVANA classified somatic BEDPE
    ch_genome_dir // channel: path(genome_dir) SigProfilerMatrixGenerator reference install, or []
    ch_sv_signatures // channel: path(tsv) reference SV32 signatures, or [] for the matrix only

    main:

    def ch_sigprofiler_input = ch_vcf
        .join(ch_segments, by: 0, remainder: true)
        .map { meta, vcf, segments -> [meta, vcf ?: [], segments ?: []] }

    SIGPROFILERASSIGNMENT(
        ch_sigprofiler_input,
        ch_genome_dir,
    )

    SVSIGNATURES(
        ch_bedpe,
        ch_sv_signatures,
    )

    emit:
    summary = SIGPROFILERASSIGNMENT.out.summary // channel: [ val(meta), path(json) ]
    signatures = SIGPROFILERASSIGNMENT.out.signatures // channel: [ val(meta), path(tsv) ]
    outdir = SIGPROFILERASSIGNMENT.out.outdir // channel: [ val(meta), path(dir) ]
    sv_summary = SVSIGNATURES.out.summary // channel: [ val(meta), path(json) ]
    sv_matrix = SVSIGNATURES.out.matrix // channel: [ val(meta), path(tsv) ]
    sv_signatures = SVSIGNATURES.out.signatures // channel: [ val(meta), path(tsv) ]
    sv_outdir = SVSIGNATURES.out.outdir // channel: [ val(meta), path(dir) ]
}
