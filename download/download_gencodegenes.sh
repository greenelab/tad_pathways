#!/bin/bash

# Description: Script will download hg19 gencode genes from the dedicated FTP site below.

# Usage: Run from command line from the parent directory 'download/download_gencodegenes.sh'

# Output: an unzipped gtf file of all gencode genes mapped to hg19 

download_folder='data/'
base_url='ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'

wget '--directory-prefix='$download_folder $base_url

gunzip 'data/gencode.v19.annotation.gtf.gz'
