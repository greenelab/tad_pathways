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

# Output boxplot of GC distribution across TADs
python gc_content_distribution.py

##############################
# PART 4: TAD Pathway Analysis
##############################
# Extract TAD based genes
python build_TAD_genelists.py

# ----------------------------
# Manual Step - WebGestalt Analysis
# ----------------------------
# Input the following trait specific genes into WebGestalt
# http://bioinfo.vanderbilt.edu/webgestalt/
# 1) Bone Mineral Density (BMD)
# 2) Type 2 Diabetes (T2D)

# For WebGestalt Parameters see README
# Click 'Export TSV Only' and save as 'data/gestalt/<TRAIT>_gestalt.tsv'
# Where <TRAIT> is either 'BMD' or 'T2D'

# Process WebGestalt Parameters for BMD and T2D
python parse_gestalt.py -t BMD
python parse_gestalt.py -t T2D

##############################
# PART 5: Download eQTL data for trait of interest
##############################
# ----------------------------
# Manual Step - See README for instructions on how to download eQTL data
# ----------------------------

##############################
# PART 6: Integrative GWAS/eQTL/TAD Analysis for BMD and T2D
##############################
# Output evidence tables
python construct_evidence.py -t BMD -e data/eqtl/eqtl_BMD_genelist.tsv \
-g data/gwas_based_genes/Bone_mineral_density_hg19_SNPs_GWAS_genelists.tsv \
-p 'skeletal system development'

python construct_evidence.py -t T2D -e data/eqtl/eqtl_T2D_genelist.tsv \
-g data/gwas_based_genes/Type_2_diabetes_hg19_SNPs_GWAS_genelists.tsv \
-p 'insulin secretion'

# Summarize GWAS/eQTL/TAD evidence
python assign_evidence_to_TADs.py -e tad_pathway/BMD_gene_evidence.csv \
-t data/gwas_TAD_location/Bone_mineral_density_hg19_SNPs.tsv \
-o tad_pathway/BMD_evidence_summary_v2.tsv

python assign_evidence_to_TADs.py -e tad_pathway/T2D_gene_evidence.csv \
-t data/gwas_TAD_location/Type_2_diabetes_hg19_SNPs.tsv \
-o tad_pathway/T2D_evidence_summary.tsv

# Output venn diagrams and gene lists for both traits
R --no-save --args tad_pathway/BMD_gene_evidence.csv \
BMD < integrative_summary.R

R --no-save --args tad_pathway/T2D_gene_evidence.csv \
T2D < integrative_summary.R

