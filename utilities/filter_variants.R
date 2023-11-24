#!/usr/bin/env Rscript

library("getopt")

spec <- matrix(
  c(
    "help", "h", 0, "logical"   , "Display help and exit", # nolint: commas_linter, line_length_linter.
    "vcf" , "v", 1, "character", "VCF file", # nolint: commas_linter.
    "mut" , "m", 1, "character", "Compressed .mut.gz file with ENSEMBL VEP predictions for viewing in IGV", # nolint: commas_linter, line_length_linter.
    "hom1", "a", 1, "character", "File with sample names in the VCF file that should be homozigous for the same allele (no matter which allele, one per line)", # nolint: line_length_linter.
    "hom2", "b", 1, "character", "List of sample names in the VCF file that should be homozigous for the same allele (which must be different from the allele carried by the samples in --hom1, one per line)", # nolint: line_length_linter.
    "drop", "d", 1, "character", "List of sample names to omit (one per line)" # nolint: line_length_linter.
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

vcf <- read.vcfR(opt[["vcf"]], verbose = FALSE)
df <- extract_gt_tidy(
  vcf,
  format_fields =  "GT",
  dot_is_NA = TRUE,
  alleles = FALSE,
  gt_column_prepend = NULL,
  verbose = TRUE
) %>%
  summarise(
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

stopifnot(
  length(sample_spec[["sample"]]) == length(unique(sample_spec[["sample"]]))
)

df <- full_join(df, sample_spec, by = "sample") %>%
  filter()


df %>% filter(!is.na(type))