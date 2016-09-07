# Incorporating TADs into GWAS Analysis - TAD Pathways

## Summary

The repository contains methods for manipulating, observing, and visualizing
topologically associating domains (TADs) in the context of 1000 Genomes Phase
III data, hg19 Gencode genes, and repeat elements. The repository also
proposes methods and tools for the incorporation of TAD domains into the
prioritization of GWAS signals through the investigation of publicly available
GWAS data. We introduce TAD pathways as a method to identify the likely causal
genes from GWAS independent of distance to sentinel SNP.

## Contact

For all questions and bug reporting please file
[GitHub issues](https://github.com/greenelab/tad_pathways/issues)

## Usage

Curate GWAS catalog and TAD boundaries to generate TAD based gene lists

```sh
./ANALYSIS.sh
```

This will download data, perform analyses, and output several figures.

## Tad Pathways

The above command outputs TAD based genes for signal found in 299 different GWAS
traits. As a case study to demonstrate the utility of a TAD based approach,
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

Note - The output of `ANALYSIS.sh` in *data/TAD_based_genes/* for all traits is
ready for TAD Pathway Analysis.

After performing the WebGestalt analysis, click `Export TSV Only` and save the
file in `data/gestalt/<TRAIT>_gestalt.tsv` where `<TRAIT>` is "BMD" for the example.

## GWAS/eQTL Integration

### Data Access  (see `download_data.sh` for more details)

* GWAS Catalog (2016-02-25)
* eQTL (2016-05-09)
[eQTL Browser](http://www.ncbi.nlm.nih.gov/projects/gap/eqtl/index.cgi "eQTL")

### Nearest gene GWAS reports

* [Bone Mineral Density](http://www.ncbi.nlm.nih.gov/pubmed/22504420 "BMD")

### eQTL Browser Parameters

* Analysis ID (All)
* Association Test Significance Filters (p-value 1 x 10^-1)
* Phenotype Traits  (*Bone Mineral Density*)

## Dependencies

All analyses were performed in the Anaconda python distribution (3.5.1). For
specific package versions please refer to `environment.yml`. R version 3.3.0 was
used for visualization. For more specific environment dependencies refer to our
accompanying docker file at `docker/Dockerfile`

