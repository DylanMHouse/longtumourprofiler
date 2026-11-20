include { GAWK as GAWK_PCGR_CNA } from '../../../modules/nf-core/gawk/main'
include { PCGR } from '../../../modules/local/pcgr/main'

workflow VCF_REPORT_PCGR {
    take:
    ch_vcf // channel: [ val(meta), path(vcf), path(tbi) ] somatic SNV/indel calls
    ch_segments // channel: [ val(meta), path(segments) ] ASCAT segments, or empty
    ch_purity // channel: [ val(meta), val(purity) ], or empty
    ch_ploidy // channel: [ val(meta), val(ploidy) ], or empty
    ch_cpsr // channel: [ val(meta), path(classification), path(conf_yaml) ] keyed by tumour meta, or empty
    ch_pcgr_dir // channel: path(pcgr_refdata)
    ch_vep_cache // channel: path(vep_cache)

    main:

    // ASCAT segments (sample, chr, startpos, endpos, nMajor, nMinor) to the PCGR CNA layout
    GAWK_PCGR_CNA(
        ch_segments.map { meta, segments -> [meta, segments, 'tsv'] },
        [],
        false,
    )

    // Every optional input is packed into one value so a missing one pads with a single null
    def ch_pcgr_input = ch_vcf
        .map { meta, vcf, tbi -> [meta, [vcf, tbi]] }
        .join(GAWK_PCGR_CNA.out.output, by: 0, remainder: true)
        .join(ch_cpsr.map { meta, classification, conf_yaml -> [meta, [classification, conf_yaml]] }, by: 0, remainder: true)
        .join(ch_purity, by: 0, remainder: true)
        .join(ch_ploidy, by: 0, remainder: true)
        .filter { _meta, vcf_tbi, _cna, _cpsr, _purity, _ploidy -> vcf_tbi != null }
        .map { meta, vcf_tbi, cna, cpsr, purity, ploidy ->
            [
                meta,
                vcf_tbi[0],
                vcf_tbi[1],
                cna ?: [],
                cpsr ? cpsr[0] : [],
                cpsr ? cpsr[1] : [],
                purity ?: [],
                ploidy ?: [],
            ]
        }

    PCGR(
        ch_pcgr_input,
        ch_pcgr_dir,
        ch_vep_cache,
    )

    emit:
    cna = GAWK_PCGR_CNA.out.output // channel: [ val(meta), path(tsv) ]
    reports = PCGR.out.pcgr_reports // channel: [ val(meta), path(dir) ]
    html = PCGR.out.html // channel: [ val(meta), path(html) ]
}
