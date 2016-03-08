"""
(C) 2016 Gregory Way
generate_snpidx_1000G.py

Description:
Map GWAS SNP information to TADs based on genomic locations. We are considering
the GRCh37/hg19 genome build for both SNP positions/TAD locations respectively.

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python generate_snpidx_1000G.py'

Output:
The script will output a pickle file that will store a python dictionary to
enable fast lookups of SNP information as they fall in TADs
"""

import sys
sys.path.insert(1, 'bin/')
from util import initialize_TAD_dictionary
from util import parse_TAD_name
import pandas as pd
import pickle
import csv

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
# Populate TAD dictionary with SNP data
####################################
with open(REF_SNPS) as snp_fh:
    snpfile = csv.reader(snp_fh, delimiter='\t')
    next(snpfile)  # Skip header
    for snprow in snpfile:
        snprow.pop()
        snp_chrom = snprow[0]
        snp_pos = int(snprow[1])
        
        # Subset lookup to only the given chromosome
        match_idx = 0
        for tad in TADdict[str(snp_chrom)]:
            tadinfo = parse_TAD_name(tad)
            tad_start = tadinfo[1]
            tad_end = tadinfo[2]

            if tad_start <= snp_pos < tad_end:
                TADdict[str(snp_chrom)][tad].append(snprow)
                match_idx = 0  # reset match counter
                continue
            match_idx += 1

        # If there are no matches, the SNP must be in a boundary
        if match_idx == len(TADdict[str(snp_chrom)]):
            TADdict['Boundary'].append(snprow)

####################################
# Once dict is loaded, make each value a pandas dataframe for easy lookup
####################################
COL_ID = ['CHROMOSOME', 'POSITION', 'MAF', 'RS']

# Look over each key in the TAD dictionary
for key in TADdict.iterkeys():
    if "Boundary" not in key:
        if len(TADdict[key]) > 0:
            for lockey, locval in TADdict[key].iteritems():
                if len(locval) > 0:
                    TADdict[key][lockey] = pd.DataFrame(locval, columns=COL_ID)

TADdict['Boundary'] = pd.DataFrame(TADdict['Boundary'], columns=COL_ID)

####################################
# Save dictionary in picklefile
####################################
pickle.dump(TADdict, open('output/SNPindex_1000G_hg19.p', 'wb'))
