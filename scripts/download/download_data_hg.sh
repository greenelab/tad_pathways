#!/bin/bash

# (C) 2016 Gregory Way

# Description: 
# Downloads the following human genome data:
# 1) 1000 Genomes Phase III
# 2) hg19 Gencode Genes
# 3) RepeatMasker hg19 repeat elements
# 4) NHGRI-EBI GWAS Catalog
# 5) UCSC hg38/hg19 LiftOver Chain File
# 6) TAD domain boundaries (hESC and IMR90)
# 7) hg19 sequence

# Usage:
# This script is run once by 'ANALYSIS.sh'

# Output:
# 1) 24 distinct raw data files separated by chromosome ('data/hg/raw1000G/')
# 2) GTF file of all hg19 Gencode genes ('data/hg/')
# 3) Fasta file of all imputed genomic repeat elements ('data/hg/')
# 4) Catalog of all significant GWAS findings ('data/hg/')
# 5) UCSC coordinates mapping file ('data/hg/')
# 6) BED file (hg19 hESC/IMR90) TAD boundaries (Dixon et al. 2012; 'data/hg/')
# 7) FASTA files for each hg19 chromosome

download_folder='data/hg/'
######################
# 1) 1000 Genomes Phase III
######################
KG_download_folder='data/hg/raw1000G/'
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
# 2) Gencode
######################
base_url='ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/'\
'gencode.v19.annotation.gtf.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip $download_folder'gencode.v19.annotation.gtf.gz'

######################
# 3) Repeat Elements
######################
base_url='www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/'\
'hg19.fa.out.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip $download_folder'hg19.fa.out.gz'
sed 's/ \+/\t/g' $download_folder'hg19.fa.out' > $download_folder'hg19.fa.out.tsv'

######################
# 4) GWAS Catalog
######################
base_url='https://bitbucket.org/gwaygenomics/download/raw/'\
'2e0bd4b7462581f6cf68d69aa51e288d1fa8943a/gwas/'\
'gwas_catalog_v1.0.1-downloaded_2016-02-25.tsv'

wget '--directory-prefix='$download_folder $base_url

######################
# 5) UCSC LiftOver Chain
######################
url='http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/'\
'hg38ToHg19.over.chain.gz'

wget '--directory-prefix='$download_folder $url

######################
# 6) hESC TAD Boundaries
######################
base_url='http://compbio.med.harvard.edu/modencode/webpage/hic/'
hESC='hESC_domains_hg19.bed'
IMR90='IMR90_domains_hg19.bed'

wget '--directory-prefix='$download_folder $base_url$hESC
wget '--directory-prefix='$download_folder $base_url$IMR90

######################
# 7) hg19 FASTA files
######################
hg19_download_folder='data/hg/hg19_fasta/'
base_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr'
end_url='.fa.gz'

# Download somatic chromosomes
for chrom in $(seq 22)
do
  wget '--directory-prefix='$hg19_download_folder $base_url$chrom$end_url
  gunzip $hg19_download_folder'chr'$chrom$end_url
done

# Download sex chromosomes
wget '--directory-prefix='$hg19_download_folder $base_url'X.fa.gz'
gunzip $hg19_download_folder'chrX.fa.gz'
wget '--directory-prefix='$hg19_download_folder $base_url'Y.fa.gz'
gunzip $hg19_download_folder'chrY.fa.gz'

