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
# Download human and mouse genome data
python scripts/download/download_hg.py
python scripts/download/download_mm.py

##############################
# PART 2: Process Data
##############################
# Extract common SNPs from 1000 Genomes and mouse genome 
python scripts/generate_common-snps.py --genome 'hg'
python scripts/generate_common-snps.py --genome 'mm'

# Generate index files mapping 1000G SNPs, genes, and repeat elements
python scripts/generate_index_files.py --TAD-Boundary 'hESC'
python scripts/generate_index_files.py --TAD-Boundary 'IMR90'
python scripts/generate_index_files.py --TAD-Boundary 'mESC'
python scripts/generate_index_files.py --TAD-Boundary 'cortex'

##############################
# PART 3: Visualize SNPs and Genes in TADs
##############################
# Output histograms and line graphs of SNP/Gene/Repeat locations in TADs
# and gc content distribution across human and mouse tads
sh scripts/tad_util/viz/visualize.sh

##############################
# PART 4: TAD Pathway Analysis
##############################
# Convert human GWAS catalog to hg19
python scripts/convert_GWAS_catalog_hg19.py

# Extract data from NHGRI-EBI GWAS catalog
Rscript scripts/NHGRI-EBI_GWAS_summary.R

# Extract TAD based genes
python scripts/build_TAD_genelists.py --TAD-Boundary 'hESC'

# ----------------------------
# Manual Step - WebGestalt Analysis
# ----------------------------
# Input the Bone Mineral Density (BMD) specific genes into WebGestalt
# http://bioinfo.vanderbilt.edu/webgestalt/

# For WebGestalt Parameters see README
# Click 'Export TSV Only' and save as 'data/gestalt/BMD_gestalt.tsv'

# Process WebGestalt Parameters for BMD and T2D
python scripts/parse_gestalt.py --trait BMD

##############################
# PART 5: Download eQTL data for trait of interest
##############################
# ----------------------------
# Manual Step - See README for instructions on how to download eQTL data
# ----------------------------

##############################
# PART 6: Integrative GWAS/eQTL/TAD Analysis
##############################
# Output evidence tables
python scripts/construct_evidence.py \
--trait BMD \
--eqtl data/eqtl/eqtl_BMD_genelist.tsv \
--gwas data/gwas_based_genes/Bone_mineral_density_hg19_SNPs_GWAS_genelists.tsv \
--pathway 'skeletal system development'

# Summarize GWAS/eQTL/TAD evidence
python scripts/assign_evidence_to_TADs.py \
--evidence-fh tad_pathway/BMD_gene_evidence.csv \
--tad-snp-path data/gwas_TAD_location/Bone_mineral_density_hg19_SNPs.tsv \
--output-fh tad_pathway/BMD_evidence_summary.tsv

# Output venn diagrams and gene lists for both traits
R --no-save --args tad_pathway/BMD_gene_evidence.csv \
BMD < integrative_summary.R

