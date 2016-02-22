#!/bin/bash

# Description: Script will download the hESC domains identified by Dixon et al. 2012 but converted to hg19
# by the UCSC liftover tool (as processed through the site below)

# Usage: Run from command line from the parent directory 'download/download_1000G.sh'

# Output: A BED file of all hg19 hESC TAD boundaries 

download_folder='data/'
base_url='http://compbio.med.harvard.edu/modencode/webpage/hic/hESC_domains_hg19.bed'

wget '--directory-prefix='$download_folder $base_url
