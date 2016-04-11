#!/bin/bash

##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exploring TADs - Analytical code for TAD based analysis and visualization
# (C) 2016 Gregory Way
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################

##############################
# PART 1: DOWNLOAD DATA
##############################
# Download 1000G Phase III data, hg19 Gencode genes, NHGRI-EBI GWAS Catalog,
# and hESC TAD domain boundaries
./download_data.sh

# Extract common SNPs from 1000 Genomes data
python generate_common-snps.py

##############################
# PART 2: PROCESS DATA
##############################
# Generate index files (maps to TAD identifiers to enable fast lookup)
# (A) Generate 1000G SNP based index file
python generate_snpidx.1000G.py

# (B) Generate gene based index file
python generate_geneidx.py

# (C) Convert GWAS Catalog to hg19
python convert_GWAS_catalog_hg19.py

##############################
# PART 3: Extract data from NHGRI-EBI GWAS catalog
##############################
# (A) Summarize the catalog by subsetting specific journals
Rscript NHGRI-EBI_GWAS_summary.R

# (B) LD correction -- Manual Step
# (https://www.broadinstitute.org/mpg/snap/ldsearchpw.php)
# Chromosome specific SNP files are saved in '/data/SNAP/' and must be 
# fed manually into SNAP with the following parameters:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SNP data set: 1000 Genomes Pilot 1
# r2 threshold: No limit
# Population panel: CEU
# Distance limit: 500
# Download to: File
# Include each query snp as a proxy for itself: yes
# Suppress warning messages in output: no
# Filter by Array: No
# Output Columns: D', Genome Coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These files should also be saved in '/data/SNAP/' and renamed to
# '/data/SNAP/SNAP_Results_chr{NUM}.txt'

# (C) Extract only independent SNPs (r2 < 0.2)
Rscript NHGRI-EBI_SNAP_summary.R

##############################
# PART 4: Visualize SNPs and Genes in TADs
##############################
# Output histograms and line graphs of SNP/Gene locations in TADs
python visualize_SNPS_1000G.py
python visualize_genes.py
python visualize_independent_snps.py

##############################
# PART 5: Network Prioritization of Trait Genes using TADs
##############################
# Extract TAD based genes
python build_TAD_genelists.py

# Network prioritization using these gene lists
python build_TAD_networks.py
