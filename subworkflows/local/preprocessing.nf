include { SAMTOOLS_VIEW as FILTER_CRAM                } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE                              } from '../../modules/nf-core/samtools/merge'
include { BCFTOOLS_INDEX                              } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_QUERY as EXTRACT_SAMPLES_BY_GT     } from '../../modules/nf-core/bcftools/query'
include { BCFTOOLS_QUERY as EXTRACT_POSSIBLE_GT       } from '../../modules/nf-core/bcftools/query'
include { SAMTOOLS_REHEADER as RENAME_SM_TAG          } from '../../modules/nf-core/samtools/reheader'  

workflow PREPROCESSING {
    take:
    reads // channel: [mandatory] [meta, cram, crai]
    vars  // channel: [mandatory] [meta, chr, start, end]
    vcf   // value  : [mandatory] vcf_file
    fasta // value  : [mandatory] reference fasta file

    main:
    versions = Channel.empty()

    reads
    .map { meta, cram, crai -> [ meta.group, meta.id ] }
    .groupTuple ()
    .map { group, samples -> [ group, samples.join(",") ] }
    .set { samples_per_group }

    Channel.fromPath( vcf )
    .map { [ [ id: null ], it ] }
    .set { vcf_ch }
    BCFTOOLS_INDEX ( vcf_ch )
    vcf_ch
    .join ( BCFTOOLS_INDEX.out.csi, by: 0, failOnDuplicate: true, failOnMismatch: true )
    .combine ( samples_per_group )
    .combine ( vars )
    .map {
        meta1, vcf, csi, group, samples, meta2, variant_id, chr, start, end ->
        def meta = [
            id     : "group_${group}_variant_${variant_id}",
            group  : group,
            samples: samples,
            variant: variant_id,
            chr    : chr,
            start  : start,
            end    : end,
            region : "${chr}:${start}-${end}",
        ]
        [ meta, vcf, csi ]
    }
    .set { vcf_groups }
    EXTRACT_POSSIBLE_GT ( vcf_groups, [], [], [] )
    EXTRACT_POSSIBLE_GT.out.output
    .splitText ( elem: 1 )
    .unique ()
    .map {
        meta, gt ->
        [ meta, gt.trim() ]
    }
    .set { possible_gt }
    vcf_groups
    .combine ( possible_gt, by: 0 )
    .map {
        meta, vcf, csi, gt ->
        def new_meta    = meta.clone()
        new_meta.gt     = gt
        new_meta.id     = "group_${meta.group}_variant_${meta.variant}_gt_${gt.replaceAll('/', '-').replaceAll('\\.', 'miss')}"
        [ new_meta, vcf, csi ]
    }
    .set { vcf_groups_gts }
    EXTRACT_SAMPLES_BY_GT ( vcf_groups_gts, [], [], [] )

    reads
    .combine ( vars )
    .map {
        meta1, cram, crai, meta2, variant_id, chr, start, end ->
        def meta = [
            id     : "sample_${meta1.id}_variant_${variant_id}",
            group  : meta1.group,
            sample : meta1.id,
            variant: variant_id,
            chr    : chr,
            start  : start,
            end    : end,
            region: "${chr}:${start}-${end}",
        ]
        [ meta, cram, crai ]
    }
    .set { filter_cram_in_ch }
    FILTER_CRAM ( filter_cram_in_ch, [ [ id: null ], fasta ], [] )
    FILTER_CRAM.out.cram
    .map { 
        meta, cram ->
        def merge_col = [ meta.group, meta.variant ]
        [ merge_col, meta, cram ] 
    }
    .set { filtered_crams }

    EXTRACT_SAMPLES_BY_GT.out.output
    .splitText ( elem: 1 )
    .map { meta, sample -> [ meta, sample.trim() ] }
    .groupTuple ()
    .map {
        meta, samples ->
        def new_meta = meta.clone()
        // before this contained all the samples in group, now overwrite with only the samples
        // that have the right GT
        new_meta.samples = samples
        def merge_col = [ meta.group, meta.variant ]
        [ merge_col, new_meta ]
    }
    .combine ( filtered_crams, by: 0 )
    .map {
        merge_col, meta1, meta2, cram ->
        assert meta1.group   == meta2.group
        assert meta1.variant == meta2.variant
        assert meta1.chr     == meta2.chr
        assert meta1.start   == meta2.start
        assert meta1.end     == meta2.end
        assert meta1.region  == meta2.region
        [ meta1, cram ]
    }
    .groupTuple ()
    .set { crams_to_merge }
    SAMTOOLS_MERGE ( crams_to_merge, [ [ id: null ], fasta ], [ [ id: null ], [] ] )
    RENAME_SM_TAG ( SAMTOOLS_MERGE.out.cram )
    RENAME_SM_TAG.out.cram.view()


    //versions = versions.mix( FILTER_CRAM              .out.versions )
    //versions = versions.mix( SAMTOOLS_MERGE           .out.versions )

    emit:
    versions              // channel: [ versions.yml ]
}
