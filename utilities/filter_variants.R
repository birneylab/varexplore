#!/usr/bin/env Rscript

library("getopt")

spec <- matrix(
  c(
    "help", "h", 0, "logical"   , "Display help and exit", # nolint: commas_linter, line_length_linter.
    "vcf" , "v", 1, "character", "VCF file", # nolint: commas_linter.
    "mut" , "m", 1, "character", "Compressed .mut.gz file with ENSEMBL VEP predictions for viewing in IGV", # nolint: commas_linter, line_length_linter.
    "hom1", "a", 1, "character", "File with sample names in the VCF file that should be homozigous for the same allele (no matter which allele, one per line)", # nolint: line_length_linter.
    "hom2", "b", 1, "character", "List of sample names in the VCF file that should be homozigous for the same allele (which must be different from the allele carried by the samples in --hom1, one per line)", # nolint: line_length_linter.
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
      c("vcf", "hom1", "hom2"),
      function(flag){
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


if (!is.null(opt[["outname"]])){
  outname <- opt[["outname"]]
} else {
  outname <- "filtered_variants"
}

vcf <- read.vcfR(opt[["vcf"]], verbose = FALSE)
df <- extract_gt_tidy(
  vcf,
  format_fields =  "GT",
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
    allele_2 = str_remove(gt, ".*/")
  )

sample_spec <- rbind(
  read_csv(opt[["hom1"]], col_names = "sample") %>% mutate(type = "hom1"),
  read_csv(opt[["hom2"]], col_names = "sample") %>% mutate(type = "hom2")
)

sample_names_all <- colnames(vcf@gt)
if (!is.null(opt[["drop"]])) {
  sample_names_to_drop <- read_csv(
    opt[["drop"]], col_names = "sample"
  )[["sample"]]
  samples_to_keep <- which(!(sample_names_all %in% sample_names_to_drop))
} else {
  samples_to_keep <- seq_along(sample_names_all)
}

stopifnot(
  length(sample_spec[["sample"]]) == length(unique(sample_spec[["sample"]]))
)

df <- full_join(df, sample_spec, by = "sample") %>%
  filter(!is.na(type)) %>%
  group_by(type, key) %>%
  summarise(
    n_alleles = length(unique(c(allele_1, allele_2))),
    is_hom = all(allele_1 == allele_2),
    allele_if_single = ifelse(n_alleles == 1, unique(c(allele_1, allele_2)), NA)
  ) %>%
  filter(is_hom & n_alleles == 1) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = key,
    names_from = type,
    names_prefix = "allele_",
    values_from = allele_if_single
  ) %>%
  filter(
    !is.na(allele_hom1) & !is.na(allele_hom2) & allele_hom1 != allele_hom2
  ) %>%
  arrange(key)

vars_to_keep <- df %>% pull(unique(key))
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