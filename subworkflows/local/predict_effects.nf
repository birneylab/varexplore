workflow PREDICT_EFFECTS {
    take:
    vcf // channel: [mandatory] [meta, vcf, tbi]

    main:
    versions = Channel.empty()

    vcf.view()

    //versions = versions.mix( GATK4_HAPLOTYPECALLER .out.versions )

    emit:

    versions  // channel: [ versions.yml ]
}