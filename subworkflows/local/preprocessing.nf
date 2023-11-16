workflow PREPROCESSING {
    take:
    vcf                    // value: [mandatory] vcf_file

    main:
    versions = Channel.empty()

    emit:
    versions              // channel: [ versions.yml ]
}
