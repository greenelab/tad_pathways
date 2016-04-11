############################################
# Incorporating TADs into GWAS Analysis
(C) 2016 Gregory Way
############################################

#############
# SUMMARY
#############
The repository contains methods for manipulating, observing, and visualizing
topologically associating domains (TADs) in the context of 1000 Genomes Phase
III data, and hg19 Gencode genes. The repository also proposes methods and tools
for the incorporation of TAD domains into the prioritization of GWAS signals
through the investigation of publicly available GWAS data.

#############
# CONTACT
#############
For all questions and bug reporting:
GregWay@mail.med.upenn.edu

#############
# USAGE
#############
The analysis can be performed in full through
a single command:
~~~~~~~~~~~~~~~~~~
./ANALYSIS.sh
~~~~~~~~~~~~~~~~~~
which downloads data, performs analyses, and outputs all figures.

#######################
# DEPENDENCIES
#######################
# Python 2.7
* pandas (0.17.1)
* vcf (0.6.7)
* scipy (0.14.0)
* scikit-learn (0.16.1)

# R 3.2.3
* readr (0.2.2)

# Linux (Ubuntu 14.04)
* plink (1.07)
* vcftools (0.1.11)
* Hubber (Sleipnir Library)(1.0)
