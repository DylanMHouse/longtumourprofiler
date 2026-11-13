include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'
include { RAITA_RAITA      } from '../../../modules/local/raita/raita/main'

workflow SOMALIER_RAITA_IDENTITY {

    take:
    ch_tumour_normal  // channel: [ val(tumour_meta), path(tumour_bam), path(tumour_bai), val(normal_meta), path(normal_bam), path(normal_bai) ]
    ch_fasta_fai      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_somalier_sites // channel: [ val(meta), path(sites_vcf) ]
    ch_raita_bed      // channel: path(bed)

    main:

    // A germline sample shared by several tumours is extracted once
    def ch_extract_input = ch_tumour_normal
        .flatMap { tumour_meta, tumour_bam, tumour_bai, normal_meta, normal_bam, normal_bai ->
            [[tumour_meta, tumour_bam, tumour_bai], [normal_meta, normal_bam, normal_bai]]
        }
        .unique()

    SOMALIER_EXTRACT(
        ch_extract_input,
        ch_fasta_fai.map { meta, fasta, _fai -> [meta, fasta] },
        ch_fasta_fai.map { meta, _fasta, fai -> [meta, fai] },
        ch_somalier_sites,
    )

    // channel: [ val(tumour_meta), path(tumour_somalier), path(normal_somalier) ]
    def ch_extract_pairs = ch_tumour_normal
        .map { tumour_meta, _tumour_bam, _tumour_bai, normal_meta, _normal_bam, _normal_bai -> [tumour_meta, normal_meta] }
        .join(SOMALIER_EXTRACT.out.extract, by: 0)
        .map { tumour_meta, normal_meta, tumour_somalier -> [normal_meta, tumour_meta, tumour_somalier] }
        .combine(SOMALIER_EXTRACT.out.extract, by: 0)
        .map { _normal_meta, tumour_meta, tumour_somalier, normal_somalier -> [tumour_meta, tumour_somalier, normal_somalier] }

    // somalier names each digest after the sample name stored in it, which comes from
    // the BAM read group and can differ from meta.id
    def ch_groups = ch_extract_pairs
        .collectFile { tumour_meta, tumour_somalier, normal_somalier ->
            ["${tumour_meta.id}.groups.csv", "${tumour_somalier.baseName},${normal_somalier.baseName}\n"]
        }
        .map { groups -> [groups.name - '.groups.csv', groups] }

    def ch_relate_input = ch_extract_pairs
        .map { tumour_meta, tumour_somalier, normal_somalier -> [tumour_meta.id, tumour_meta, [tumour_somalier, normal_somalier]] }
        .join(ch_groups, by: 0)
        .multiMap { _id, tumour_meta, extracts, groups ->
            extract: [tumour_meta, extracts, []]
            groups: groups
        }

    SOMALIER_RELATE(
        ch_relate_input.extract,
        ch_relate_input.groups,
    )

    RAITA_RAITA(
        ch_tumour_normal.map { tumour_meta, tumour_bam, tumour_bai, _normal_meta, normal_bam, normal_bai ->
            [tumour_meta, normal_bam, normal_bai, tumour_bam, tumour_bai]
        },
        ch_raita_bed,
    )

    emit:
    somalier_extract = SOMALIER_EXTRACT.out.extract    // channel: [ val(meta), path(somalier) ]
    somalier_html    = SOMALIER_RELATE.out.html        // channel: [ val(tumour_meta), path(html) ]
    somalier_pairs   = SOMALIER_RELATE.out.pairs_tsv   // channel: [ val(tumour_meta), path(pairs.tsv) ]
    somalier_samples = SOMALIER_RELATE.out.samples_tsv // channel: [ val(tumour_meta), path(samples.tsv) ]
    raita_tsv        = RAITA_RAITA.out.tsv             // channel: [ val(tumour_meta), path(raita.tsv) ]
    raita_report     = RAITA_RAITA.out.report          // channel: [ val(tumour_meta), path(raita_report.txt) ]
}
