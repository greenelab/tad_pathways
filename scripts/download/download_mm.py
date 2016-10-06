# Download Mouse Data

# The script will retreive mouse data used in the TAD pathways pipeline

# 1. Mouse Genome Project SNPs (v2)
# 2. mm9 Gencode Genes
# 3. RepeatMasker mm9 repeat elements
# 4. TAD Domain Boundaries (mESC and Cortex)
# 5. mm9 sequence

import os
from download_util import download_file, process_repeats

# Make new folder if it doesn't exist already
genome = 'mm'
download_folder = 'data/mm/'
if not (os.path.exists(download_folder)):
    os.makedirs(download_folder)

# 1. Mouse Genome Project SNPs v2
base_url = 'ftp://ftp-mouse.sanger.ac.uk/REL-1211-SNPs_Indels/'
filename = 'mgp.v2.snps.annot.reformat.vcf.gz'
download_file(base_url, filename, download_folder)

# 2. Gencode Mouse
base_url = 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M9/'
filename = 'gencode.vM9.annotation.gtf.gz'
download_file(base_url, filename, download_folder)

# 3. RepeatMasker mm9 Repeat Elements
base_url = 'http://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328' \
           '-db20090604/'
filename = 'mm9.fa.out.gz'
download_file(base_url, filename, download_folder)
process_repeats(download_folder, filename, genome)

# 4. TAD Domain Boundaries
base_url = 'http://chromosome.sdsc.edu/mouse/hi-c/'
mESC = 'mESC.domain.tar.gz'
cortex = 'cortex.domain.tar.gz'
download_file(base_url, mESC, download_folder)
download_file(base_url, cortex, download_folder)
# TODO - process these domain files (extract and take combined)

# 5. mm9 Full Sequence
mm9_download_folder = 'data/mm/mm9_fasta/'
base_url = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/'

# Create download directory
if not os.path.exists(mm9_download_folder):
    os.makedirs(mm9_download_folder)

# Download chromosome sequences
for chrom in list(range(1, 20)) + ['X', 'Y']:
    filename = 'chr{}.fa.gz'.format(chrom)
    download_file(base_url, filename, mm9_download_folder)
