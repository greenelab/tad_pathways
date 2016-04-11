"""
(C) 2016 Gregory Way
generate_geneidx.py

Description:
Observe SNP distributions and locations across TADs (hg19)

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python visualize_SNPS_1000G.py'

Output:
Two matplotlib plots (SNP distribution and histogram across TADs)
"""

import sys
sys.path.insert(1, 'bin/')
from util import parse_TAD_name
from util import parse_SNP_position
from util import ID_TAD_bins
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare

####################################
# Load Constants
####################################
NUM_BINS = 50

####################################
# Load Data
####################################
KG_TADs = pd.read_pickle('output/SNPindex_1000G_hg19.p')

####################################
# Observe location of SNPs in TADs
####################################
# Initialize a dictionary that will store location information for gene types
TADdict = range(NUM_BINS)
TADdict = dict.fromkeys(TADdict)
for key in TADdict.iterkeys():
    TADdict[key] = 0

# Initialize a chromosome specific dictionary holding bin specific info
Chromosome_Specific_TAD = dict.fromkeys(KG_TADs.keys())
for chromkey in Chromosome_Specific_TAD.keys():
    # Create a list of zeros of len(NUM_BINS)
    Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Also observe distribution of amount of SNPs in each TAD
howmanySNPsinTAD = []

# Loop through each TAD and SNP to fill the histogram
SNP_side = [0, 0]
for key in KG_TADs.iterkeys():
    if key not in ['NoTAD', 'Boundary', 'X', 'Y']:
        for tadkey, tadval in KG_TADs[key].iteritems():
            howmanySNPsinTAD.append(len(tadval))
            if len(tadval) > 0:
                tadname = parse_TAD_name(tadkey)
                for snp in range(len(tadval)):
                    snp_loc = parse_SNP_position(tadval.ix[snp, :])
                    bin_assign = ID_TAD_bins(tadname, NUM_BINS, snp_loc)
                    TADdict[bin_assign] += 1

                    # Count the number of SNPs in left/right
                    if bin_assign < 25:
                        SNP_side[0] += 1
                    else:
                        SNP_side[1] += 1

                    # Add a count for chromsome specificity
                    Chromosome_Specific_TAD[key][bin_assign] += 1

# Convert to DataFrame
SNPLocations = pd.DataFrame.from_dict(TADdict, orient='index')
SNPLoc_chrom = pd.DataFrame.from_dict(Chromosome_Specific_TAD, orient='index')
SNPLoc_chrom = SNPLoc_chrom.transpose()

####################
# Perform ChiSquare test of SNPs to 'right' of TAD
####################
TAD_SNP_sig = chisquare(SNP_side)

with open('output/TAD_SNP_rightness_Chisquare.txt', 'w') as chisq_fh:
    chisq_fh.write('SNPs in the left vs. right of TAD\n')
    chisq_fh.write('left,right\n')
    chisq_fh.write(','.join('%s' % x for x in SNP_side) + '\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_SNP_sig))

####################################
# Visualization
####################################
# SNP Distribution
SNPLocations.plot()
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('SNP location inside TADs')
plt.grid(True)
plt.savefig('figures/SNPlocations_TADs_1000G_hg19.png')
plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
for chrom in Chromosome_Specific_TAD.keys():
    if chrom not in ['NoTAD', 'Boundary', 'X', 'Y']:
        nameplot = 'figures/chrom/SNPdistribution_chrom_' + chrom + '_' + str(NUM_BINS) + 'bins.png'
        chrom_loc = SNPLoc_chrom[chrom]
        chrom_loc.plot()
        plt.xlabel('Bins (Normalized TAD length)')
        plt.ylabel('Frequency')
        plt.title('SNP location inside TADs\nChromosome ' + str(chrom))
        plt.grid(True)
        plt.savefig(nameplot)
        plt.close()

# Histogram
n, bins, patches = plt.hist(howmanySNPsinTAD, NUM_BINS, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('Number of SNPs')
plt.ylabel('Frequency')
plt.title('Distribution of number of SNPs in TADs')
plt.axis([min(howmanySNPsinTAD), max(howmanySNPsinTAD), 0, 0.0005])
plt.grid(True)
addtext = 'mean = ' + str(round(np.mean(howmanySNPsinTAD), 2)) + '\n' + \
          'std = ' + str(round(np.std(howmanySNPsinTAD), 2))
plt.text(max(howmanySNPsinTAD), .0003, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig('figures/SNPcounts_TADs_1000G_hg19.png')
plt.close()
