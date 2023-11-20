include { SAMTOOLS_VIEW as FILTER_CRAM                } from '../../modules/nf-core/samtools/view'                                                                                                                                                                                              
include { SAMTOOLS_INDEX as INDEX_SINGLE_SAMPLE_CRAMS } from '../../modules/nf-core/samtools/index'                                                                                                                                                                                              

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
        meta1, cram, crai, meta2, chr, start, end ->
        meta = [:]
        meta.sample  = meta1.id
        meta.variant = meta2.id
        meta.group   = meta1.group
        meta.id      = "variant_${meta.variant}_sample_${meta.sample}"
        meta.region  = "${chr}:${start}-${end}"
        [ meta, cram, crai ]
    }
    .set { filter_cram_in_ch }
    FILTER_CRAM ( filter_cram_in_ch, [ [ id:null ], fasta], [] )
    INDEX_SINGLE_SAMPLE_CRAMS ( FILTER_CRAM.out.cram )
    FILTER_CRAM.out.cram
    .join ( INDEX_SINGLE_SAMPLE_CRAMS.out.crai , by: 0, failOnDuplicate: true, failOnMismatch: true)
    .map { meta, cram, crai -> [ [ group : meta.group, variant: meta.variant, region: meta.region ], cram, crai ] }
    .groupTuple ()
    .view()
    .set { grouped_crams }

    versions = versions.mix( FILTER_CRAM              .out.versions )
    versions = versions.mix( INDEX_SINGLE_SAMPLE_CRAMS.out.versions )

    emit:
    versions              // channel: [ versions.yml ]
}
