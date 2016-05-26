############################################
# Incorporating TADs into GWAS Analysis - TAD Pathways
(C) 2016 Gregory Way
############################################

#######################
# SUMMARY
#######################
The repository contains methods for manipulating, observing, and visualizing
topologically associating domains (TADs) in the context of 1000 Genomes Phase
III data, and hg19 Gencode genes. The repository also proposes methods and tools
for the incorporation of TAD domains into the prioritization of GWAS signals
through the investigation of publicly available GWAS data. We introduce TAD
pathways as a method to identify the likely causal genes from GWAS independent
of distance to sentinel SNP.

#######################
# CONTACT
#######################
For all questions and bug reporting:
GregWay@mail.med.upenn.edu

#######################
# USAGE
#######################
Curate GWAS catalog and TAD boundaries to generate TAD based gene lists 
~~~~~~~~~~~~~~~~~~
./ANALYSIS.sh
~~~~~~~~~~~~~~~~~~
Downloads data, performs analyses, and outputs several figures.

#######################
# TAD PATHWAY
#######################
The above command outputs TAD based genes for signal found in 299 different GWAS
traits. As a case study to demonstrate the utility of a TAD based approach,
input the TAD based gene lists for the following two diseases into pathway
analysis:

* Bone Mineral Density (1,297 genes for TAD pathway)
* Inflammatory Bowel Disease (2,270 genes)

Run a [WebGestalt](http://bioinfo.vanderbilt.edu/webgestalt/ "Pathway Analysis")
pathway analysis on the gene lists for the above traits.

# WebGestalt Parameters
* Select gene ID type *hsapiens__gene_symbol*
* Enrichment Analysis *GO Analysis*
* GO Slim Classification *Yes*
* Reference Set *hsapiens__genome*
* Statistical Method *Hypergeometric*
* Multiple Test Adjustment *BH*
* Significance Level *Top10*
* Minimum Number of Genes for a Category *4*

Note - The output of ANALYSIS.sh in *data/TAD_based_genes* for all traits is
ready for TAD Pathway Analysis.

#######################
# GWAS/eQTL INTEGRATION
#######################
# Data Access  (see download_data.sh for more details)
* GWAS Catalog (2016-02-25)
* eQTL (2016-05-09) (http://www.ncbi.nlm.nih.gov/projects/gap/eqtl/index.cgi)

# Nearest gene GWAS reports
* Bone Mineral Density (http://www.ncbi.nlm.nih.gov/pubmed/22504420)
* Inflammatory Bowel Disease (http://www.ncbi.nlm.nih.gov/pubmed/?term=26192919)

# eQTL Browser Parameters
* Analysis ID (All)
* Association Test Significance Filters (p-value 1 x 10^-1)
* Phenotype Traits  (*Bone mineral density*, *Inflammatory bowel disease*)

#######################
# DEPENDENCIES
#######################
# Python 2.7.6
* pandas (0.17.1)
* vcf (0.6.7)
* scipy (0.14.0)
* scikit-learn (0.16.1)
* pyliftover (0.3)

# R 3.3.0
* readr (0.2.2)
* VennDiagram (1.6.17)

# Linux (Ubuntu 14.04)
* bcftools (1.3)
* vcftools (0.1.11)

