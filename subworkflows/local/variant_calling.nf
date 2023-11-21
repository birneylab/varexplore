include { GATK4_HAPLOTYPECALLER                                  } from '../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS                                        } from '../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_GENOMICSDBIMPORT                                 } from '../../modules/nf-core/gatk4/genomicsdbimport'

workflow VARIANT_CALLING {
    take:
    merged_crams // channel: [mandatory] [meta, cram, crai]
    fasta        // value  : [mandatory] reference fasta file
    fa_idx       // value  : [mandatory] reference fasta file fai index (and if fasta is compressed also gzi, as a single emission)
    fa_dict      // value  : [mandatory] reference fasta file dictionary

    main:
    versions = Channel.empty()

    merged_crams
    .map { meta, cram, crai -> [ meta, cram, crai, meta.region, [] ] }
    .set { haplotypecaller_in }
    GATK4_HAPLOTYPECALLER ( haplotypecaller_in, fasta, fa_idx, fa_dict, [], [] )
    GATK4_HAPLOTYPECALLER.out.vcf
    .join ( GATK4_HAPLOTYPECALLER.out.tbi )
    .map {
        meta, vcf, tbi ->
        def new_meta = [
            id: "variant_${meta.variant}",
            variant: meta.variant,
            chr    : meta.chr,
            start  : meta.start,
            end    : meta.end,
            region : meta.region
        ]
        [ new_meta, vcf, tbi ]
    }
    .groupTuple ()
    .map { meta, vcfs, tbis -> [ meta, vcfs, tbis, [], meta.region, [] ] }
    .set { genomicsdbimport_in }
    GATK4_GENOMICSDBIMPORT ( genomicsdbimport_in, false, false, false )
    GATK4_GENOMICSDBIMPORT.out.genomicsdb
    .view()

    //versions = versions.mix( BCFTOOLS_INDEX       .out.versions )

    emit:

    versions     // channel: [ versions.yml ]
}