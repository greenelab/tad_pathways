"""
generate_snpidx_1000G.py
Author: Greg Way

Description:
Map GWAS SNP information to TADs based on genomic locations. We are considering
the GRCh37/hg19 genome build for both SNP positions/TAD locations respectively.

Usage:
Run from command line 'python generate_snpidx_1000G.py'

Output:
The script will output a pickle file that will store a python dictionary to
enable fast lookups of SNP information as they fall in TADs
"""

import sys
sys.path.insert(1, 'bin/')
from util import initialize_TAD_dictionary
import pandas as pd
import pickle


####################################
# Load Constants
####################################
TAD_LOC = 'data/hESC_domains_hg19.bed'
REF_SNPS = 'data/common-snps.tsv'

####################################
# Initialize a dictionary with the structure:
# 1st key: chromosome, 2nd key: TAD ID:start-stop
####################################
TADdict = initialize_TAD_dictionary(TAD_LOC)

####################################
# Once dict is loaded, make each value a pandas dataframe for easy lookup
####################################
COL_ID = ['POSITION', 'CHROMOSOME', 'MAF', 'RS']

# Look over each key in the TAD dictionary
for key in TADdict.iterkeys():
    if "Boundary" not in key:
        for lockey, locval in TADdict[key].iteritems():
            if len(locval) > 0:
                TADdict[key][lockey] = pd.DataFrame(locval, columns=COL_ID)

TADdict['Boundary'] = pd.DataFrame(TADdict['Boundary'], columns=COL_ID)

####################################
# Save dictionary in picklefile
####################################
pickle.dump(TADdict, open('output/SNPindex_1000G_hg19.p', 'wb'))
