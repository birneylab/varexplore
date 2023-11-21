include { SAMTOOLS_VIEW as FILTER_CRAM                } from '../../modules/nf-core/samtools/view'

workflow VARIANT_CALLING {
    take:
    merged_crams // channel: [mandatory] [meta, cram, crai]
    fasta        // value  : [mandatory] reference fasta file

    main:
    versions = Channel.empty()

    merged_crams.view()

    //versions = versions.mix( BCFTOOLS_INDEX       .out.versions )

    emit:

    versions     // channel: [ versions.yml ]
}
