# Incorporating TADs into GWAS Analysis - TAD_Pathways

Gregory P. Way and Casey S. Greene 2016

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.163950.svg)](https://doi.org/10.5281/zenodo.163950)

## Summary

The repository contains methods for manipulating, observing, and visualizing
topologically associating domains (TADs) in the context of SNPs, genes, and 
repeat elements for human (hg19) and mouse (mm9) genomes.

The repository also proposes methods and tools for the incorporation of TAD
domains into the prioritization of GWAS signals through the investigation of
publicly available GWAS data. We introduce **TAD pathways** as a method to
identify the likely causal genes from GWAS independent of distance to sentinel
SNP.

A preprint of our method is available [here on bioRxiv](http://dx.doi.org/10.1101/087718)

## Contact

For all questions and bug reporting please file a
[GitHub issue](https://github.com/greenelab/tad_pathways/issues)

For all other questions contact Casey Greene at csgreene@mail.med.upenn.edu or
Struan Grant at grants@email.chop.edu

## Usage

There are two ways to implement a TAD_Pathways analysis:

1. Disease/Trait Specific - Uses GWAS identified SNPs
2. Custom - Uses custom SNP list

### Disease/Trait Specific

Curates the GWAS catalog and TAD boundaries to visualize TADs and generate
TAD based gene lists. This will also perform a TAD pathways analysis for
Bone Mineral Density GWAS. This will reproduce the analysis and figures used
in the paper.

```sh
# Using python dependencies
conda env create --quiet --force --file environment.yml
source activate tad_pathways

bash scripts/run_pipeline.sh
```

This will download data, perform analyses, and output several genomic figures.
The command will also output TAD based genes for 299 different GWAS traits.
Our TAD_Pathways method can be applied directly using these gene lists.

### Custom

TAD_Pathways is customizable and allows a user to prespecify any SNP list of
interest to test TAD based pathway associations. To perform a custom analysis
create a comma separated file where the first row of each column names the list
of snps below in subsequent rows.

E.g.: `custom_example.csv`

| Group 1 | Group 2 |
| ------- | ------- |
| rs12345 | rs67891 |
| rs19876 | rs54321 |

Then, perform the following steps:

```bash
# Extract locations for SNP list
Rscript --vanilla scripts/tad_util/build_snp_list.R \
        --snp_file "custom_example.csv" \
        --output_file "mapped_results.tsv"

# Build TAD based genelists for each group
python scripts/build_custom_TAD_genelist.py \
       --snp_data_file "mapped_results.tsv" \
       --output_file "custom_tad_genelist.tsv"

# The output file is then ready for the manual "TAD_Pathways" steps below
```


## TAD_Pathways

As a case study to demonstrate the utility of a TAD based approach,
input the TAD based gene list for the Bone Mineral Density (1,297 genes) into a
pathway analysis:

Next, run a
[WebGestalt](http://bioinfo.vanderbilt.edu/webgestalt/ "Pathway Analysis")
pathway analysis on the gene list.

### WebGestalt Parameters

| Parameter | Input |
| --------- | ----- |
| Select gene ID type | *hsapiens__gene_symbol* |
| Enrichment Analysis | *GO Analysis* |
| GO Slim Classification | Yes |
| Reference Set | *hsapiens__genome* |
| Statistical Method | *Hypergeometric* |
| Multiple Test Adjustment | *BH* |
| Significance Level | *Top10* |
| Minimum Number of Genes for a Category | *4*

Note - The output of `scripts/run_pipeline.sh` in *data/TAD_based_genes/* for
all traits is ready for TAD Pathway Analysis.

After performing the WebGestalt analysis, click `Export TSV Only` and save the
file in `data/gestalt/<TRAIT>_gestalt.tsv` where `<TRAIT>` is "BMD" for the
example.

## GWAS/eQTL Integration

### Data Access

* GWAS Catalog (2016-02-25)
* eQTL (2016-05-09)
[eQTL Browser](http://www.ncbi.nlm.nih.gov/projects/gap/eqtl/index.cgi "eQTL")

### Nearest BMD gene GWAS reports

* [Richards _et al._ 2008 Lancet](http://doi.org/10.1016/S0140-6736(08)60599-1)
* [Rivadeneira _et al._ 2009 Nature Genetics](http://doi.org/10.1038/ng.446)
* [Estrada _et al._ 2012 Nature Genetics](http://doi.org/10.1038/ng.2249)
* [Styrkarsdottir _et al._ 2013 Nature](http://doi.org/10.1038/nature12124)

### eQTL Browser Parameters

* Analysis ID (All)
* Association Test Significance Filters (p-value 1 x 10^-1)
* Phenotype Traits  (*Bone Mineral Density*)

## Dependencies

All analyses were performed in the Anaconda python distribution (3.5.1). For
specific package versions please refer to `environment.yml`. R version 3.3.0 was
used for visualization. For more specific environment dependencies refer to our
accompanying docker file at `docker/Dockerfile`

