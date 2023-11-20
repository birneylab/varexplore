include { SAMTOOLS_VIEW as FILTER_CRAM                } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE                              } from '../../modules/nf-core/samtools/merge'
include { BCFTOOLS_INDEX                              } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_QUERY as EXTRACT_SAMPLES_BY_GT     } from '../../modules/nf-core/bcftools/query'
include { BCFTOOLS_QUERY as EXTRACT_POSSIBLE_GT       } from '../../modules/nf-core/bcftools/query'

workflow PREPROCESSING {
    take:
    reads // channel: [mandatory] [meta, cram, crai]
    vars  // channel: [mandatory] [meta, chr, start, end]
    vcf   // value  : [mandatory] vcf_file
    fasta // value  : [mandatory] reference fasta file

    main:
    versions = Channel.empty()

    reads
    .combine ( vars )
    .map {
        meta1, cram, crai, meta2, variant_id, chr, start, end ->
        def meta = [:]
        meta.sample  = meta1.id
        meta.variant = meta2.id
        meta.group   = meta1.group
        meta.id      = "variant_${meta.variant}_sample_${meta.sample}"
        meta.region  = "${chr}:${start}-${end}"
        [ meta, cram, crai ]
    }
    .set { filter_cram_in_ch }
    FILTER_CRAM ( filter_cram_in_ch, [ [ id: null ], fasta ], [] )
    FILTER_CRAM.out.cram
    .map {
        meta, cram ->
        def new_meta     = [:]
        new_meta.id      = "variant_${meta.variant}_group_${meta.group}"
        new_meta.variant = meta.variant
        new_meta.group   = meta.group
        new_meta.region  = meta.region
        [ new_meta, meta.sample, cram ]
    }
    .groupTuple ()
    .map {
        meta, samples, crams ->
        def new_meta     = meta.clone()
        new_meta.samples = samples
        [ new_meta, crams ]
    }
    .set { grouped_crams }
    SAMTOOLS_MERGE ( grouped_crams, [ [ id: null ], fasta], [ [ id: null ], [] ] )

    Channel.fromPath( vcf )
    .map { [ [ id: null ], it ] }
    .set { vcf_ch }
    BCFTOOLS_INDEX ( vcf_ch )
    vcf_ch
    .join ( BCFTOOLS_INDEX.out.csi, by: 0, failOnDuplicate: true, failOnMismatch: true )
    .combine ( grouped_crams )
    .first()
    .map {
        meta1, vcf, csi, meta2, crams ->
        [ meta2, vcf, csi ]
    }
    .set { vcf_groups_ch }
    EXTRACT_POSSIBLE_GT ( vcf_groups_ch, [], [], [] )
    EXTRACT_POSSIBLE_GT.out.output
    .splitText ( elem: 1 )
    .unique ()
    .map { meta, gt -> [ meta, gt.trim() ]}
    .set { possible_gt }
    vcf_groups_ch
    .join ( possible_gt, by: 0, failOnMismatch: true )
    .view()
    //EXTRACT_SAMPLES_BY_GT ( vcf_groups_ch, [], [], [] )

    //versions = versions.mix( FILTER_CRAM              .out.versions )
    //versions = versions.mix( SAMTOOLS_MERGE           .out.versions )

    emit:
    versions              // channel: [ versions.yml ]
}
