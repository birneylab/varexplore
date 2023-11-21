include { ENSEMBLVEP_DOWNLOAD } from '../../modules/nf-core/ensemblvep/download'

workflow PREDICT_EFFECTS {
    take:
    vcf             // channel: [mandatory] [meta, vcf, tbi]
    ensemblvep_info // value  : [mandatory] [meta, assembly, species, chache_version]

    main:
    versions = Channel.empty()

    vcf.view()
    ENSEMBLVEP_DOWNLOAD ( ensemblvep_info )

    //versions = versions.mix( GATK4_HAPLOTYPECALLER .out.versions )

    emit:

    versions  // channel: [ versions.yml ]
}