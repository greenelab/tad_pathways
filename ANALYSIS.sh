#!/bin/bash

##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TAD Pathways - Analytical code for TAD based analysis and visualization
# (C) 2016 Gregory Way
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################

##############################
# PART 1: Download Data
##############################
# Download 1000G Phase III data, hg19 Gencode genes, NHGRI-EBI GWAS Catalog,
# and hESC TAD domain boundaries
./download_data.sh

##############################
# PART 2: Process Data
##############################
# Extract common SNPs from 1000 Genomes data
python generate_common-snps.py

# Generate index files (maps to TAD identifiers to enable fast lookup)
# 1000G SNP / genes / repeat elements
python generate_index_files.py

# Convert GWAS Catalog to hg19
python bin/convert_GWAS_catalog_hg19.py

# Extract data from NHGRI-EBI GWAS catalog
Rscript NHGRI-EBI_GWAS_summary.R

##############################
# PART 3: Visualize SNPs and Genes in TADs
##############################
# Output histograms and line graphs of SNP/Gene/Repeat locations in TADs
python visualize_TAD_locations.py

##############################
# PART 4: TAD Pathway Analysis
##############################
# Extract TAD based genes
python build_TAD_genelists.py

# Input the following trait specific genes into WebGestalt
# http://bioinfo.vanderbilt.edu/webgestalt/
# 1) Bone Mineral Density
# 2) Inflammatory Bowel Disease

# Obtain the pathway with the lowest p value enrichment
