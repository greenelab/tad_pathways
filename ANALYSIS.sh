#!/bin/bash

##############################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exploring TADs - Analytical code for TAD based analysis and visualization
# Written by: Greg Way
#
# Way, G., Greene, C., Grant, S.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################

##############################
# PART 1: DOWNLOAD DATA
##############################
# (A) 1000G project SNPs ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ (1000 Genomes SNPs)
sh download/download_1000G.sh   # Download the data
python generate_common-snps.py  # Extract common SNPs from the data

# (B) Get gene coordinates ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
sh download/download_gencodegenes.sh

# (C) NHGRI-EBI GWAS Catalog  https://www.ebi.ac.uk/gwas/docs/downloads 
# The catalog is already downloaded and is placed in data/gwas_catalog_v1.0.tsv
# SNP rsID mapped to build 144

# (D) hESC Topologically Associated Domains http://compbio.med.harvard.edu/modencode/webpage/hic/hESC_domains_hg19.bed
sh download/download_TADs.sh

##############################
# Part 2: GENERATE INDEX FILES
##############################
# Each index file will map to TAD identifiers and enable fast lookup
# (A) Generate 1000G SNP based index files
python generate_snpidx.1000G.py

# (B) Generate gene based index files
python generate_geneidx.py

##############################
# Part 3: Visualize SNP and Genes in TADs
##############################
# This section will output histograms and line graphs of SNP/Gene locations in TADs
python visualize_SNPS_1000G.py
python visualize_genes.py  # Also outputs a chi square test for gene start enrichment near boundaries

##############################
# Part 4: Extract data from NHGRI-EBI GWAS catalog
##############################
# (A) Summarize the catalog by subsetting specific journals
R --no-save < generate_NHGRI-EBI_GWAS_summary.R

# (B) Linkage disequilibrium correction (https://www.broadinstitute.org/mpg/snap/ldsearchpw.php)
# Following this analysis, chromosome specific SNP files are saved in '/data/SNAP/' and must be 
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
# These files should also be saved in '/data/SNAP/' and renamed to '/data/SNAP/SNAP_Results_chr{NUM}.txt'

# (C) Extract only independent SNPs (r2 < 0.2)
R --no-save < generate_NHGRI-EBI_SNAP_summary.R

# (D) Visualize independent signals in TADs
python TAD_independent_snps.py
