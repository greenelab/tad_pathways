#!/bin/bash

# Description: Script will download 1000 Genomes Phase III data from the dedicated FTP site below. Data are separated
# by chromosome and the script parses this information to download 24 distinct raw data files.

# Usage: Run from command line from the parent directory 'download/download_1000G.sh'

# Output: 24 distinct raw data files separated by chromosome

download_folder='data/raw1000G/'
base_url='ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr'
end_url='.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# Download somatic chromosomes
for chrom in $(seq 22)
do
  wget '--directory-prefix='$download_folder $base_url$chrom$end_url
done

# Download sex chromosomes
wget '--directory-prefix='$download_folder $base_url'X.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
wget '--directory-prefix='$download_folder $base_url'Y.phase3_integrated_v1b.20130502.genotypes.vcf.gz'
