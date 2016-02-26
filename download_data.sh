#!/bin/bash

# (C) 2016 Gregory Way

# Description: 
# Downloads the following data:
# 1) 1000 Genomes Phase III
# 2) hg19 Gencode Genes
# 3) NHGRI-EBI GWAS Catalog
# 4) hESC TAD domain boundaries

# Usage:
# This script is run by 'ANALYSIS.sh' but can also be run directly by calling
# from the parent directory './download_data.sh'. Only run once.

# Output:
# 1) 24 distinct raw data files separated by chromosome ('data/raw1000G/')
# 2) GTF file of all hg19 Gencode genes ('data/')
# 3) Catalog of all significant GWAS findings ('data/')
# 4) BED file of hg19 hESC TAD boundaries (Dixon et al. 2012; 'data/')

KG_download_folder='data/raw1000G/'
download_folder='data/'

######################
# 1000 Genomes Phase III
######################
base_url='ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr'
end_url='.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# Download somatic chromosomes
for chrom in $(seq 22)
do
  wget '--directory-prefix='$KG_download_folder $base_url$chrom$end_url
done

# Download sex chromosomes
wget '--directory-prefix='$KG_download_folder $base_url'X.phase3_shapeit2_'\
'mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
wget '--directory-prefix='$KG_download_folder $base_url'Y.phase3_integrated'\
'_v1b.20130502.genotypes.vcf.gz'

######################
# Gencode
######################
base_url='ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/'\
'gencode.v19.annotation.gtf.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip 'data/gencode.v19.annotation.gtf.gz'

######################
# GWAS Catalog
######################
base_url='https://bitbucket.org/gwaygenomics/download/raw/'\
'2e0bd4b7462581f6cf68d69aa51e288d1fa8943a/gwas/'\
'gwas_catalog_v1.0.1-downloaded_2016-02-25.tsv'

wget '--directory-prefix='$download_folder $base_url

######################
# hESC TAD Boundaries
######################
base_url='http://compbio.med.harvard.edu/modencode/webpage/hic/'\
'hESC_domains_hg19.bed'

wget '--directory-prefix='$download_folder $base_url

