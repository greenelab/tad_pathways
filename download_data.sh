#!/bin/bash

# (C) 2016 Gregory Way

# Description: 
# Downloads the following data:
# 1) 1000 Genomes Phase III
# 2) hg19 Gencode Genes
# 3) RepeatMasker hg19 repeat elements
# 4) NHGRI-EBI GWAS Catalog
# 5) UCSC hg38/hg19 LiftOver Chain File
# 6) hESC TAD domain boundaries

# Usage:
# This script is run once by 'ANALYSIS.sh'

# Output:
# 1) 24 distinct raw data files separated by chromosome ('data/raw1000G/')
# 2) GTF file of all hg19 Gencode genes ('data/')
# 3) Fasta file of all imputed genomic repeat elements ('data/')
# 4) Catalog of all significant GWAS findings ('data/')
# 5) UCSC coordinates mapping file ('data/')
# 6) BED file of hg19 hESC TAD boundaries (Dixon et al. 2012; 'data/')

######################
# 1) 1000 Genomes Phase III
######################
KG_download_folder='data/raw1000G/'
download_folder='data/'
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

gunzip 'data/hg19.fa.out.gz'

######################
# 3) Repeat Elements
######################
base_url='www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/'\
'hg19.fa.out.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip 'data/gencode.v19.annotation.gtf.gz'
sed 's/ \+/\t/g' data/hg19.fa.out > data/hg19.fa.out.tsv

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
base_url='http://compbio.med.harvard.edu/modencode/webpage/hic/'\
'hESC_domains_hg19.bed'

wget '--directory-prefix='$download_folder $base_url

