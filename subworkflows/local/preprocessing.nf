include { SAMTOOLS_VIEW as FILTER_CRAM            } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE                          } from '../../modules/nf-core/samtools/merge'
include { BCFTOOLS_INDEX                          } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_QUERY as EXTRACT_SAMPLES_BY_GT } from '../../modules/nf-core/bcftools/query'
include { BCFTOOLS_QUERY as EXTRACT_POSSIBLE_GT   } from '../../modules/nf-core/bcftools/query'
include { SAMTOOLS_REHEADER as RENAME_SM_TAG      } from '../../modules/nf-core/samtools/reheader'  
include { SAMTOOLS_INDEX                          } from '../../modules/nf-core/samtools/index'
include { SAMTOOLS_FAIDX                          } from '../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_DICT                           } from '../../modules/nf-core/samtools/dict'

workflow PREPROCESSING {
    take:
    reads // channel: [mandatory] [meta, cram, crai]
    vars  // channel: [mandatory] [meta, chr, start, end]
    vcf   // value  : [mandatory] [meta, vcf_file]
    fasta // value  : [mandatory] [meta, fasta]

    main:
    versions = Channel.empty()

    reads
    .map { meta, cram, crai -> [ meta.group, meta.id ] }
    .groupTuple ()
    .map { group, samples -> [ group, samples.join(",") ] }
    .set { samples_per_group }

    BCFTOOLS_INDEX ( vcf )
    Channel.value ( vcf )
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
    FILTER_CRAM ( filter_cram_in_ch, fasta, [] )
    FILTER_CRAM.out.cram
    .map { 
        meta, cram ->
        def merge_col = [ meta.group, meta.sample, meta.variant ]
        [ merge_col, meta, cram ] 
    }
    .set { filtered_crams }

    EXTRACT_SAMPLES_BY_GT.out.output
    .splitText ( elem: 1 )
    .map {
        meta, sample ->
        def new_meta = meta.clone()
        new_meta.remove ( "samples" ) // this was referring to the whole group, not GT specific
        def sample_trimmed = sample.trim()
        def merge_col = [ meta.group, sample_trimmed, meta.variant ]
        [ merge_col, sample_trimmed, new_meta ]
    }
    .combine ( filtered_crams, by: 0 )
    .map { merge_col, sample, meta1, meta2, cram -> [ meta1, sample, cram ] }
    .groupTuple ()
    .map {
        meta, samples, crams ->
        def new_meta = meta.clone()
        new_meta.samples = samples.sort()
        [ new_meta, crams ]
    }
    .set { crams_to_merge }
    SAMTOOLS_MERGE ( crams_to_merge, fasta, [ [ id: null ], [] ] )
    RENAME_SM_TAG ( SAMTOOLS_MERGE.out.cram )
    RENAME_SM_TAG.out.cram
    .map {
        meta, cram ->
        def new_meta = meta.clone()
        new_meta.remove ( "samples" )
        [ new_meta, cram ]
    }
    .set { renamed_crams }
    SAMTOOLS_INDEX ( renamed_crams )
    renamed_crams
    .join( SAMTOOLS_INDEX.out.crai, by: 0, failOnDuplicate: true, failOnMismatch: true )
    .set { merged_crams }

    SAMTOOLS_FAIDX ( fasta, [[ id: null ], [] ] )
    SAMTOOLS_DICT ( fasta )
    SAMTOOLS_FAIDX.out.fai .map { meta, fai  -> fai  }.set { fa_fai  }
    SAMTOOLS_FAIDX.out.gzi .map { meta, gzi  -> gzi  }.set { fa_gzi  }
    fa_fai.mix ( fa_gzi ).collect().set { fa_idx }
    SAMTOOLS_DICT .out.dict.map { meta, dict -> dict }.set { fa_dict }

    versions = versions.mix( BCFTOOLS_INDEX       .out.versions )
    versions = versions.mix( EXTRACT_POSSIBLE_GT  .out.versions )
    versions = versions.mix( EXTRACT_SAMPLES_BY_GT.out.versions )
    versions = versions.mix( FILTER_CRAM          .out.versions )
    versions = versions.mix( SAMTOOLS_MERGE       .out.versions )
    versions = versions.mix( RENAME_SM_TAG        .out.versions )
    versions = versions.mix( SAMTOOLS_INDEX       .out.versions )

    emit:
    merged_crams // channel: [ meta, cram, crai ]
    fa_idx       // value  : fasta index (fai and optionally gzi as a single emission)
    fa_dict      // value  : fasta dict 

    versions     // channel: [ versions.yml ]
}