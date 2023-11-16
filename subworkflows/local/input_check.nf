//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK   } from '../../modules/local/samplesheet_check'
include { VARIANT_TABLE_CHECK } from '../../modules/local/variant_table_check'

workflow INPUT_CHECK {
    take:
    input         // file: /path/to/samplesheet.csv
    variant_table // file: /path/to/variant_table.csv

    main:
    SAMPLESHEET_CHECK ( input )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_cram_channel(it) }
        .set { reads }

    VARIANT_TABLE_CHECK ( variant_table )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_variant_channel(it) }
        .set { vars }

    emit:
    reads                                     // channel: [ meta, cram, crai ]
    vars                                      // channel: [ meta, chr, start, end ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// funtion to interpolate ${projectDir} inside the samplesheet for testing profiles
// https://stackoverflow.com/questions/37379101/how-to-convert-string-to-gstring-and-replace-placeholder-in-groovy
def interpolate_string(String str) {

    def res = Eval.me('projectDir', projectDir, '"' + str + '"') as String

    return res

}

// Function to get list of [ meta, [ cram, crai ] ]
def create_cram_channel(LinkedHashMap row) {

    // interpolate and read files
    cram = file(interpolate_string(row.cram))
    crai = file(interpolate_string(row.crai))

    // create meta map
    def meta    = [:]
    meta.id     = row.sample
    meta.group  = row.group

    // add path(s) of the cram and crai file(s) to the meta map
    def cram_meta = []
    if (!cram.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Cram file does not exist!\n${row.cram}"
    }
    if (!crai.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Crai file does not exist!\n${row.crai}"
    }
    cram_meta = [ meta, cram, crai ]
    return cram_meta
}

// Function to get list of [ meta, [ chr, start, end ] ]
def create_variant_channel(LinkedHashMap row) {

    chr          = row.chr
    region_start = row.region_start
    region_end   = row.region_end

    // create meta map
    def meta    = [:]
    meta.id     = row.var_id

    var_meta = [ meta, chr, region_start, region_end ]
    return var_meta
}
