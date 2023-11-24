process FORMAT_VEP_OUTPUT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(vep)

    output:
    tuple val(meta), path("*.mut.gz") , emit: mut
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    clean_lines <- sub("## .*", "", readLines("${vep}"))
    clean_lines <- clean_lines[clean_lines != ""]
    clean_lines[1] <- sub("^#", "", clean_lines[1])
    df_in <- read.table(text = clean_lines, sep = "\\t", header = TRUE)
    df_out <- data.frame(
        chr = sub("(.*):(.*)-(.*)", "\\\\1", df_in[["Location"]]),
        start = sub("(.*):(.*)-(.*)", "\\\\2", df_in[["Location"]]),
        end = sub("(.*):(.*)-(.*)", "\\\\3", df_in[["Location"]]),
        patient = df_in[["IND"]],
        type = df_in[["Consequence"]],
        ref = df_in[["REF_ALLELE"]],
        alt = df_in[["Allele"]],
        gene = df_in[["Gene"]],
        biotype = df_in[["BIOTYPE"]],
        symbol = df_in[["SYMBOL"]],
        swissprot = df_in[["SWISSPROT"]],
        trembl = df_in[["TREMBL"]],
        uniparc = df_in[["UNIPARC"]],
        exon = df_in[["EXON"]],
        position_cdna = df_in[["Position.in.cDNA"]],
        position_cds = df_in[["Position.in.CDS"]],
        position_protein = df_in[["Position.in.protein"]],
        aa_change = df_in[["Amino.acid.change"]],
        codon_change = df_in[["Codon.change"]],
        sample_zygosity = df_in[["ZYG"]]
    )
    write.table(df_out, gzfile("${prefix}.mut.gz"), sep = "\\t", row.names = FALSE, quote = FALSE)

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            "END_VERSIONS\\n",
            sep = "\\n"
        )
    )
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mut.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
