include { GATK4_HAPLOTYPECALLER  } from '../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_GENOMICSDBIMPORT } from '../../modules/nf-core/gatk4/genomicsdbimport'
include { GATK4_GENOTYPEGVCFS    } from '../../modules/nf-core/gatk4/genotypegvcfs'

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
    .map { meta, genomicsdb -> [ meta, genomicsdb, [], meta.region, [] ] }
    .set { genotype_input }
    GATK4_GENOTYPEGVCFS(genotype_input, fasta, fa_idx, fa_dict, [], [] )
    GATK4_GENOTYPEGVCFS.out.vcf
    .join ( GATK4_GENOTYPEGVCFS.out.tbi, by: 0, failOnDuplicate: true, failOnMismatch: true )
    .set { vcf }

    versions = versions.mix( GATK4_HAPLOTYPECALLER .out.versions )
    versions = versions.mix( GATK4_GENOMICSDBIMPORT.out.versions )
    versions = versions.mix( GATK4_GENOTYPEGVCFS   .out.versions )

    emit:
    vcf       // channel: [ meta, vcf, tbi ]

    versions  // channel: [ versions.yml ]
}