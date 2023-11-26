# ![birneylab/varexplore](docs/images/birneylab-varexplore_name_light.png#gh-light-mode-only) ![birneylab/flexlmm](docs/images/birneylab-varexplore_name_dark.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**birneylab/varexplore** is a bioinformatics pipeline that clusters next-generation sequencing data from a set of samples according to the allele that they carry at a given marker, and then calls genetic variants and predicts variant effects on the aggregated meta-samples.
This is useful when dealing with imputed variants from low-coverage sequencing, such as the ones produced by the [birneylab/stitchimpute](https://github.com/birneylab/stitchimpute) pipeline.
Many imputation tools impute only SNPs, and in a GWAS contest the causal variant at locus may be something (i.e., indels) that is in linkage disequilibrium with an imputed marker, but not included in the set of imputed markers themselves.
In such cases, grouping samples according to the allele at a marker of interest (the SNP with the strongest association signal) allows for calling additional variants in the region that are in linkage disequilibrium with the marker.
This would allow to discover for example frameshifts and other mutations with large effects on protein function that are in strong linkage with a marker of of interest.

**Disclaimer**: this pipeline uses the nf-core template but it is not part of nf-core itself.

![birneylab/varexplore_metro_map](docs/images/birneylab_varexplore_drawing.png)

1. Group samples according to the genotypes at the selected markers ([`bcftools`](https://samtools.github.io/bcftools/bcftools.html))
1. Extract a user-defined region from the sequencing files around the variant of interest ([`samtools`](http://www.htslib.org/doc/samtools.html))
1. Merge the sequencing files from samples with the same genotype at the variant of interest ([`samtools`](http://www.htslib.org/doc/samtools.html))
1. Edit the SM tags of the sequencing files so that they are seen as a single sample ([`samtools`](http://www.htslib.org/doc/samtools.html))
1. Joint call variants in the meta-samples for each vairant of interest ([`gatk4`](https://gatk.broadinstitute.org/hc/en-us))
1. Predict variant effects ([`ENSEMBL VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html))

## Integration with [birneylab/stitchimpute](https://github.com/birneylab/stitchimpute)

In order to use a vcf file obtained from the **birneylab/stitchimpute** pipeline, add an ID for the variants that you want to be used for separating the samples using `bcftools`:

```bash
bcftools annotate \
   --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' \
   -O z \
   -o stitchimpute.with_ids.vcf.gz \
   stitchimpute.vcf.gz
```

Then, in `--variant_table` refer to the variants of interest using the same naming scheme. For example, a variant on chromosome 1 position 123 with reference allele A and alternate allele T will have ID:

```
1_123_A_T
```

## Visualising the results

The files `<OUTDIR>/variants/variant_<CHR>_<POS>_<REF>_<ALT>/variant_<CHR>_<POS>_<REF>_<ALT>.mut.gz` containing variant consequences can be directly loaded in [IGV](https://igv.org/) for visualisation. They variants will be coloured according to the ENSEMBL VEP impact ("MODIFIER", "LOW", "MODERATE", "HIGH").

## Filtering variants
The script `utilities/filter_variants.R` can be used as a stand-alone tool to filter variants and ENSEMBL VEP predictions produced by this pipeline according to their allele distribution in different samples. This is useful for example when certain samples should be excluded, or when only variants that have the same genotype in a group of samples and a different genotype in another group of samples should be retained.

It should be run in the folder `<OUTDIR>/variants/variant_<CHR>_<POS>_<REF>_<ALT>` for the grouping variant of interest.
Assuming that the script has been downloaded to the current directory as `filter_variants.R` one can run it as follows:

```bash
filter_variants.R \
   --vcf variant_<CHR>_<POS>_<REF>_<ALT>.vcf.gz \
   --mut variant_<CHR>_<POS>_<REF>_<ALT>.mut.gz \
   --hom1 hom1_samples.txt \
   --hom2 hom2_samples.txt \
   --drop samples_to_drop.txt
```

The files `hom1_samples.txt`, `hom2_samples.txt`, and `samples_to_drop.txt` contain one sample name per line, matching the sample names in the vcf file given in input to the script (those can be obtained running `bcftools query -l variant_<CHR>_<POS>_<REF>_<ALT>.vcf.gz`).
The samples in the file given to the `--drop` argument are excluded from the output.

Variants are retained only if they satisfy all of the following conditions:
- fully homozygous and with the same homozygous genotype in all the `--hom1` samples
- fully homozygous and with the same homozygous genotype in all the `--hom2` samples
- with the genotypes in `--hom1` and `--hom2` samples being not identical

The filtered set of variants is written to a new VCF file `variants_filtered.vcf.gz`. The filtered set of variant consequences is written to a new TSV file `variants_filtered.mut.gz`. The basename of the output can be customised using the `--outname` flag.

## Birneylab-specific information

For ease of use, the ideal settings for stitch for medaka samples have been specified in a profile called `medaka`.
This can be activated with the flag `-profile medaka`.
Always use this profile when working with medaka samples.

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:
```
sample,cram,crai,group
s1,s1.cram,s1.cram.crai,g1
s2,s2.cram,s2.cram.crai,g1
s3,s3.cram,s3.cram.crai,g2
s4,s4.cram,s4.cram.crai,g2
```

Each row represents a sample with associated sequencing files (`cram`) and indexes (`crai`). The `group` column indicates whether the samples should be sub-divided in groups when forming the meta-samples. Simply set all the samples to the same group if you don't want to use this feature. `sample` should match the sample names in the vcf file provided, to assign genotypes to the correct samples.

Prepare also a file containing the variants of interest and the region boundaries that looks as follows:

`variant_table.csv`:
```
var_id,chr,region_start,region_end
1_100_G_A,1,50,150
1_200_C_T,1,150,250
```

Each row represents a variant that should be used to group the samples according to the allele that they carry for such variant.
`var_id` should match the ID field of a variant in the vcf file provided. `chr` should name the chromosome the variant resides in. `region_start` and `region_end` define the region of interest in the meta-samples where new variants in linkage disequilibrium with the variant of interest should be called.

You can run the pipeline using:

```bash
nextflow run birneylab/varexplore \
   -profile <docker/singularity/.../institute> \
   --vcf input.vcf.gz \
   --input sample_table.csv \
   --variant_table variant_table.csv \ 
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

> For more details and further functionality, please refer to the [usage documentation](docs/usage.md) and the [parameter documentation](docs/parameters.md).

> **Warning**:
> It is highly recommended to use the docker or singularity profile. Some processes do not have a working conda configuration.

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Credits

> birneylab/varexplore was originally written by Saul Pierotti.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
