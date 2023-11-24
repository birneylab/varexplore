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
    df_in[["Location"]] <- ifelse(
        grepl("^.*:.*-.*\$", df_in[["Location"]]),
        df_in[["Location"]],
        sub("(.*):(.*)", "\\\\1:\\\\2-\\\\2", df_in[["Location"]])
    )
    df_in <- df_in[df_in[["ZYG"]] != "HOMREF",]
    df_out <- data.frame(
        chr = sub("(.*):(.*)-(.*)", "\\\\1", df_in[["Location"]]),
        start = as.numeric(sub("(.*):(.*)-(.*)", "\\\\2", df_in[["Location"]])),
        end = as.numeric(sub("(.*):(.*)-(.*)", "\\\\3", df_in[["Location"]])),
        sample = df_in[["IND"]],
        impact = df_in[["IMPACT"]],
        consequence = df_in[["Consequence"]],
        ref = df_in[["REF_ALLELE"]],
        alt = df_in[["Allele"]],
        gene = df_in[["Gene"]],
        position_cdna = df_in[["cDNA_position"]],
        position_cds = df_in[["CDS_position"]],
        position_protein = df_in[["Protein_position"]],
        aa_change = df_in[["Amino_acids"]],
        codon_change = df_in[["Codons"]],
        sample_zygosity = df_in[["ZYG"]]
    )
    df_out[["end"]] <- ifelse(
        df_out[["end"]] > df_out[["start"]],
        df_out[["end"]],
    ifelse(
        df_out[["end"]] == df_out[["start"]],
        df_out[["end"]] + 1,
        NA
    ))
    stopifnot(all(df_out[["end"]] > df_out[["start"]]))
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
