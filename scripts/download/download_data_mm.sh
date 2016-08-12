#!/bin/bash

# (C) 2016 Gregory Way

# Description: 
# Downloads the following mouse genome data:
# 1) Mouse Genome Project SNPs (v3)
# 2) hg19 Gencode Genes
# 3) RepeatMasker mm9 repeat elements
# 4) TAD domain boundaries (mESC and Cortex)
# 5) mm9 sequence

# Usage:
# This script is run once by 'ANALYSIS.sh'

# Output:
# 1) 24 distinct raw data files separated by chromosome ('data/raw_snps/')
# 2) GTF file of all mm9 Gencode genes ('data/')
# 3) Fasta file of all imputed genomic repeat elements ('data/')
# 4) BED file (mm9 mESC and cortex) TAD boundaries (Dixon et al. 2012; 'data/')
# 5) FASTA files for each mm9 chromosome

download_folder='data/mm/'
######################
# 1) Mouse Genome Project v5
######################
base_url='ftp://ftp-mouse.sanger.ac.uk/REL-1211-SNPs_Indels/'\
'mgp.v2.snps.annot.reformat.vcf.gz'

# Download full VCF file
wget '--directory-prefix='$download_folder $base_url

######################
# 2) Gencode Mouse
######################
base_url='ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/'\
'gencode.vM9.annotation.gtf.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip $download_folder'gencode.v19.annotation.gtf.gz'

######################
# 3) Repeat Elements
######################
base_url='www.repeatmasker.org//genomes/mm9/RepeatMasker-rm328-db20090604/'\
'mm9.fa.out.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip $download_folder'mm9.fa.out.gz'
sed 's/ \+/\t/g' $download_folder'mm9.fa.out' > $download_folder'mm9.fa.out.tsv'

######################
# 4) TAD Boundaries
######################
base_url='http://chromosome.sdsc.edu/mouse/hi-c/'
mESC='mESC.domain.tar.gz'
cortex='cortex.domain.tar.gz'

wget '--directory-prefix='$download_folder $base_url$mESC
wget '--directory-prefix='$download_folder $base_url$cortex
tar -xf $download_folder$mESC -C $download_folder
tar -xf $download_folder$cortex -C $download_folder
mv 'data/mm/mESC/HindIII_combined/total.HindIII.combined.domain' \
'data/mESC_domains.mm9.bed'
mv 'data/mm/cortex/combined/total.combined.domain' 'data/cortex_domains_mm9.bed'

######################
# 5) hg19 FASTA files
######################
mm9_download_folder='data/mm/mm9_fasta/'
base_url='http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr'
end_url='.fa.gz'

# Download somatic chromosomes
for chrom in $(seq 22)
do
  wget '--directory-prefix='$mm9_download_folder $base_url$chrom$end_url
  gunzip $mm9_download_folder'chr'$chrom$end_url
done

# Download sex chromosomes
wget '--directory-prefix='$mm9_download_folder $base_url'X.fa.gz'
gunzip $mm9_download_folder'chrX.fa.gz'
wget '--directory-prefix='$mm9_download_folder $base_url'Y.fa.gz'
gunzip $mm9_download_folder'chrY.fa.gz'

