include { ENSEMBLVEP_DOWNLOAD } from '../../modules/nf-core/ensemblvep/download'
include { ENSEMBLVEP_VEP      } from '../../modules/nf-core/ensemblvep/vep'
include { FORMAT_VEP_OUTPUT   } from '../../modules/local/format_vep_output'

workflow PREDICT_EFFECTS {
    take:
    vcf             // channel: [mandatory] [meta, vcf, tbi]
    ensemblvep_info // value  : [mandatory] [meta, assembly, species, chache_version]
    fasta           // value  : [mandatory] [meta, fasta]

    main:
    versions = Channel.empty()

    vcf
    .map{ meta, vcf, tbi -> [ meta, vcf, [] ] }
    .set { vcf_for_vep }
    ENSEMBLVEP_DOWNLOAD ( ensemblvep_info )
    ENSEMBLVEP_DOWNLOAD.out.cache
    .map{ meta, cache -> [ cache ] }
    .set { vep_cache }
    ENSEMBLVEP_VEP (
        vcf_for_vep,
        ensemblvep_info[1],
        ensemblvep_info[2],
        ensemblvep_info[3],
        vep_cache,
        fasta,
        [],
    )
    FORMAT_VEP_OUTPUT ( ENSEMBLVEP_VEP.out.tab )

    versions = versions.mix( ENSEMBLVEP_DOWNLOAD.out.versions )
    versions = versions.mix( ENSEMBLVEP_VEP     .out.versions )
    versions = versions.mix( FORMAT_VEP_OUTPUT  .out.versions )

    emit:

    versions  // channel: [ versions.yml ]
}