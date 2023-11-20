include { SAMTOOLS_VIEW as FILTER_CRAM } from '../../modules/nf-core/samtools/view'                                                                                                                                                                                              

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
    .filter { meta.sample == "s1" }
    .set { filter_cram_in_ch }
    FILTER_CRAM ( filter_cram_in_ch, [ [ id:null ], fasta], [] )
    FILTER_CRAM.out.cram
    .view()

    reads
    .map { meta, cram, crai -> [ [ group : meta.group ], cram, crai ] }
    .groupTuple ()
    .set { grouped_crams }

    //versions = versions.mix( INPUT_CHECK.out.versions )

    emit:
    versions              // channel: [ versions.yml ]
}
