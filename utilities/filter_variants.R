#!/usr/bin/env Rscript

library("getopt")

spec <- matrix(
  c(
    "help", "h", 0, "logical", "Display help and exit", # nolint: commas_linter, line_length_linter.
    "vcf", "v", 1, "character", "VCF file", # nolint: commas_linter.
    "mut", "m", 1, "character", "Compressed .mut.gz file with ENSEMBL VEP predictions for viewing in IGV", # nolint: commas_linter, line_length_linter.
    "hom1", "a", 1, "character", "File with sample names in the VCF file that should be homozigous for the same allele (no matter which allele, one per line)", # nolint: line_length_linter.
    "hom2", "b", 1, "character", "List of sample names in the VCF file that should be homozigous for the same allele (which must be different from the allele carried by the samples in --hom1, one per line)", # nolint: line_length_linter.
    "het", "c", 1, "character", "List of sample names in the VCF file that should be heterozygous", # nolint: line_length_linter.
    "drop", "d", 1, "character", "List of sample names to omit (one per line)", # nolint: line_length_linter.
    "outname", "o", 1, "character", "Output basename"
  ),
  byrow = TRUE,
  ncol = 5
)

opt <- getopt(spec)

if (
  !is.null(opt[["help"]]) || any(
    sapply(
      c("vcf", "hom1"),
      function(flag) {
        is.null(opt[[flag]])
      }
    )
  ) || all(
    sapply(
      c("hom2", "het"),
      function(flag) {
        is.null(opt[[flag]])
      }
    )
  )
) {
  cat(getopt(spec, usage = TRUE, command = "filter_variants.R"))
  q(status = 1)
}

library("vcfR")
library("tidyverse")


if (!is.null(opt[["outname"]])) {
  outname <- opt[["outname"]]
} else {
  outname <- "filtered_variants"
}

vcf <- read.vcfR(opt[["vcf"]], verbose = FALSE)
df <- extract_gt_tidy(
  vcf,
  format_fields = "GT",
  dot_is_NA = TRUE,
  alleles = FALSE,
  gt_column_prepend = NULL,
  verbose = TRUE
) %>%
  reframe(
    key = Key,
    sample = Indiv,
    gt = str_replace(GT, "\\|", "/"),
    allele_1 = str_remove(gt, "/.*"),
    allele_2 = str_remove(gt, ".*/"),
    allele_1 = ifelse(allele_1 == ".", NA, allele_1),
    allele_2 = ifelse(allele_2 == ".", NA, allele_2)
  )

sample_spec <- read_csv(opt[["hom1"]], col_names = "sample") %>%
  mutate(type = "hom1")

if (!is.null(opt[["hom2"]])) {
  sample_spec <- rbind(
    sample_spec,
    read_csv(opt[["hom2"]], col_names = "sample") %>% mutate(type = "hom2")
  )
}

if (!is.null(opt[["het"]])) {
  sample_spec <- rbind(
    sample_spec,
    read_csv(opt[["het"]], col_names = "sample") %>% mutate(type = "het")
  )
}

sample_names_all <- colnames(vcf@gt)
if (!is.null(opt[["drop"]])) {
  sample_names_to_drop <- read_csv(
    opt[["drop"]],
    col_names = "sample"
  )[["sample"]]
  samples_to_keep <- which(!(sample_names_all %in% sample_names_to_drop))
} else {
  samples_to_keep <- seq_along(sample_names_all)
}

stopifnot(
  length(sample_spec[["sample"]]) == length(unique(sample_spec[["sample"]]))
)

df <- full_join(df, sample_spec, by = "sample") %>%
  filter(!is.na(type))


if (!is.null(opt[["het"]])) {
  key_het_ok <- df %>%
    filter(type == "het") %>%
    group_by(key) %>%
    filter(
      all(
        length(unique(na.omit(c(allele_1, allele_2))) == 2) &
          ((allele_1 != allele_2) | is.na(allele_1) | is.na(allele_2))
      )
    ) %>%
    ungroup() %>%
    reframe(key) %>%
    distinct()

  df <- df %>% inner_join(key_het_ok)
}

if (!is.null(opt[["hom1"]])) {
  key_hom1_ok <- df %>%
    filter(type == "hom1") %>%
    group_by(key) %>%
    filter(
      length(unique(na.omit(c(allele_1, allele_2)))) == 1
    ) %>%
    mutate(allele_hom = unique(na.omit(c(allele_1, allele_2)))) %>%
    ungroup() %>%
    reframe(key, type, allele_hom) %>%
    distinct()

  df_hom <- df %>% inner_join(key_hom1_ok)
  df <- df %>% inner_join(key_hom1_ok %>% reframe(key) %>% distinct())
}

if (!is.null(opt[["hom2"]]) && !is.null(opt[["hom1"]])) {
  key_hom2_ok <- df %>%
    filter(type == "hom2") %>%
    group_by(key) %>%
    filter(
      length(unique(na.omit(c(allele_1, allele_2)))) == 1
    ) %>%
    mutate(allele_hom = unique(na.omit(c(allele_1, allele_2)))) %>%
    ungroup() %>%
    reframe(key, type, allele_hom) %>%
    distinct()

  df_hom <- df %>%
    inner_join(key_hom2_ok) %>%
    bind_rows(df_hom)

  df <- df %>% inner_join(key_hom2_ok %>% reframe(key) %>% distinct())

  key_hom_diff_ok <- df_hom %>%
    reframe(
      type,
      key,
      allele_hom,
    ) %>%
    pivot_wider(
      id_cols = key,
      names_from = type,
      names_prefix = "allele_",
      values_from = allele_hom
    ) %>%
    filter(
      all(
        !is.na(allele_hom1) |
          !is.na(allele_hom2) |
          allele_hom1 != allele_hom2
      )
    ) %>%
    reframe(key)

  df <- df %>% inner_join(key_hom_diff_ok)
}

if (!is.null(opt[["hom2"]]) && is.null(opt[["hom1"]])) {
  stop("Specifying --hom2 and not --hom1 is not allowed")
}

if (
  !is.null(opt[["het"]]) ||
    (!is.null(opt[["hom1"]]) && !is.null(opt[["hom2"]]))
) {
  df <- df %>%
    group_by(key) %>%
    filter(length(unique(na.omit(c(allele_1, allele_2)))) == 2) %>%
    ungroup()
} else {
  df <- df %>%
    group_by(key) %>%
    filter(length(unique(na.omit(c(allele_1, allele_2)))) == 1) %>%
    ungroup()
}

vars_to_keep <- df %>%
  pull(key) %>%
  unique()
vcf_filtered <- vcf[vars_to_keep, samples_to_keep]
write.vcf(vcf_filtered, sprintf("%s.vcf.gz", outname))

if (is.null(opt[["mut"]])) quit(status = 0)

vars_to_keep_ids <- tibble(
  # ids are chr_pos
  id = rownames(extract.gt(vcf_filtered, IDtoRowNames = TRUE))
)
sample_names_to_keep <- tibble(
  sample = colnames(extract.gt(vcf_filtered))
)
mut <- read_tsv(opt[["mut"]]) %>%
  inner_join(vars_to_keep_ids) %>%
  inner_join(sample_names_to_keep)

write_tsv(mut, sprintf("%s.mut.gz", outname))
