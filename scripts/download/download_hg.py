# Download Human Data

# The script will retreive human data used in the TAD pathways pipeline

# 1. 1000 Genomes Phase III SNPs
# 2. hg19 Gencode Genes
# 3. RepeatMasker hg19 Repeat Elements
# 4. The NHGRI-EBI GWAS Catalog
# 5. The UCSC hg38 to hg19 LiftOver Chain File
# 6. TAD Domain Boundaries (hESC and IMR90)
# 7. The Full hg19 Sequence

import os
from download_util import download_file

# Make new folder if it doesn't exist already
download_folder = 'data/hg/'
if not (os.path.exists(download_folder)):
    os.makedirs(download_folder)

# 1. 1000 Genomes Phase III SNPs
KG_folder = 'data/hg/raw1000G/'
base_url = 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/'

# Create download directory
if not os.path.exists(KG_folder):
    os.makedirs(KG_folder)

# Download somatic chromosomes
for chrom in range(1, 23):
    fh = 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.' \
         'genotypes.vcf.gz'.format(chrom)
    download_file(base_url, fh, KG_folder)

# Download sex chromosomes
sex_chrm = ['ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.'
            'genotypes.vcf.gz',
            'ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz']
for chrom in sex_chrm:
    download_file(base_url, chrom, KG_folder)

# 2. hg19 Gencode Genes
base_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/'
fh = 'gencode.v19.annotation.gtf.gz'
download_file(base_url, fh, download_folder)

# 3. RepeatMasker hg19 Repeat Elements
base_url = 'http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-' \
           'db20140131/'
fh = 'hg19.fa.out.gz'
download_file(base_url, fh, download_folder)

# 4. The NHGRI-EBI GWAS Catalog
# GWAS Catalog downloaded on February 25th, 2016
base_url = 'https://bitbucket.org/gwaygenomics/download/raw/'\
           '2e0bd4b7462581f6cf68d69aa51e288d1fa8943a/gwas/'
fh = 'gwas_catalog_v1.0.1-downloaded_2016-02-25.tsv'
download_file(base_url, fh, download_folder)

# 5. TAD Domain Boundaries (hESC and IMR90)
base_url = 'http://compbio.med.harvard.edu/modencode/webpage/hic/'
hESC = 'hESC_domains_hg19.bed'
IMR90 = 'IMR90_domains_hg19.bed'
download_file(base_url, hESC, download_folder)
download_file(base_url, IMR90, download_folder)

# 6. The Full hg19 Sequence
hg19_download_folder = 'data/hg/hg19_fasta/'
base_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/'

# Create download directory
if not os.path.exists(hg19_download_folder):
    os.makedirs(hg19_download_folder)

# Download chromosome sequences
for chrom in list(range(1, 23)) + ['X', 'Y']:
    fh = 'chr{}.fa.gz'.format(chrom)
    download_file(base_url, fh, hg19_download_folder)
