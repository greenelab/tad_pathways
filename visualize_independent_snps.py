"""
(C) 2016 Gregory Way
TAD_independent_snps.py

Description:
Visualize the independent SNPs in the NHGRI-EBI GWAS catalog

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python visualize_independent_snps.py'

Output:
1) Several matplotlib plots describing the independent SNPs in TADs
2) Coordinates for TADs with the most SNPs (independent and all)
"""

import sys
sys.path.insert(1, 'bin/')
from util import parse_TAD_name
from util import parse_SNP_position
from util import ID_TAD_bins
import csv
import pandas as pd
import matplotlib.pyplot as plt
import operator

####################################
# Load Constants
####################################
NUM_BINS = 50

####################################
# Load picklefile
####################################
KG_TADs = pd.read_pickle('output/SNPindex_1000G_hg19.p')

tad_id_dict = {}
for key in KG_TADs.iterkeys():
    if key not in ['NoTAD', 'Boundary', 'X', 'Y']:
        for tadkey in KG_TADs[key].iterkeys():
            tad_id = parse_TAD_name(tadkey)[0]
            tad_id_dict[tad_id] = 0

rs_set = set()
cool_snps = []
with open('output/all_independent_snps_across_traits.csv', 'r') as snp_fh:
    snp_data = csv.reader(snp_fh)
    for snp in snp_data:
        rs = snp[0]
        chromosome = snp[1]
        KG_TAD_Sub = KG_TADs[chromosome]
        for KG_key in KG_TAD_Sub.keys():
            if len(KG_TAD_Sub[KG_key]) > 0:
                if rs in KG_TAD_Sub[KG_key]['RS'].tolist():
                    rs_set.add(rs)
                    KG_key_id = parse_TAD_name(KG_key)[0]
                    tad_id_dict[KG_key_id] += 1
                    if KG_key_id == 1617:
                        cool_snps.append(rs)

######################################
# Generate Plots
######################################
plt.plot(tad_id_dict.keys(), tad_id_dict.values(), 'ro')
plt.xlabel('TAD ID')
plt.ylabel('Number of Independent SNP signals')
plt.title('Independent SNP signals across TADS')
plt.axis([0, len(tad_id_dict), 0, max(tad_id_dict.values()) + 5])
plt.grid(True)
plt.savefig('figures/independent_SNPS_across_TADS_hg19.png')
plt.close()

# Get Maximum Value TAD
max_tad_id = max(tad_id_dict.iteritems(), key=operator.itemgetter(1))[0]

for key in KG_TADs.iterkeys():
    if key not in ['NoTAD', 'Boundary', 'X', 'Y']:
        for tadkey in KG_TADs[key].iterkeys():
            tad_info = parse_TAD_name(tadkey)
            tad_id = tad_info[0]
            if tad_id == max_tad_id:
                chrom = 'chr' + key
                location = str(tad_info[1]) + '-' + str(tad_info[2])

# Write the maximum independent signals to file
with open('output/max_independent_tad.txt', 'w') as max_fh:
    max_fh.write(chrom + ':' + location)

######################################
# All SNPs (dependent + Independent)
######################################
tad_id_dict_all = {}
for key in KG_TADs.iterkeys():
    if key not in ['NoTAD', 'Boundary', 'X', 'Y']:
        for tadkey in KG_TADs[key].iterkeys():
            tad_id = parse_TAD_name(tadkey)[0]
            tad_id_dict_all[tad_id] = 0

with open('output/all_unique_trait_snps_chrID.csv', 'r') as snp_fh:
    snp_data = csv.reader(snp_fh)
    for snp in snp_data:
        if snp[1] not in ['23', 'NA']:
            rs = snp[0]
            chromosome = snp[1]
            KG_TAD_Sub = KG_TADs[chromosome]
            for KG_key in KG_TAD_Sub.keys():
                if len(KG_TAD_Sub[KG_key]) > 0:
                    if rs in KG_TAD_Sub[KG_key]['RS'].tolist():
                        KG_key_id = parse_TAD_name(KG_key)[0]
                        tad_id_dict_all[KG_key_id] += 1

# Generate plots
plt.plot(tad_id_dict_all.keys(), tad_id_dict_all.values(), 'ro')
plt.xlabel('TAD ID')
plt.ylabel('Number of SNP signals (All SNP Signals)')
plt.title('ALL SNP signals across TADS')
plt.axis([0, len(tad_id_dict_all), 0, max(tad_id_dict_all.values()) + 5])
plt.grid(True)
plt.savefig('figures/all_signal_SNPS_across_TADS_hg19.png')
plt.close()

# Get Maximum Value TAD
max_tad_id_all = max(tad_id_dict_all.iteritems(), key=operator.itemgetter(1))[0]

for key in KG_TADs.iterkeys():
    if key not in ['NoTAD', 'Boundary', 'X', 'Y']:
        for tadkey in KG_TADs[key].iterkeys():
            tad_info = parse_TAD_name(tadkey)
            tad_id = tad_info[0]
            if tad_id == max_tad_id_all:
                chrom = 'chr' + key
                location = str(tad_info[1]) + '-' + str(tad_info[2])

# Write the maximum independent signals to file
with open('output/max_all_snps_tad.txt', 'w') as max_fh:
    max_fh.write(chrom + ':' + location)

# Plot independent vs. dependent
plt.plot(tad_id_dict_all.values(), tad_id_dict.values(), 'ro')
plt.xlabel('Number of SNP signals (All SNP Signals)')
plt.ylabel('Number of independent SNP signals')
plt.title('All Signals vs. Independent Signals')
plt.grid(True)
plt.axis([0, max(tad_id_dict_all.values()) + 5, 0,
         max(tad_id_dict.values()) + 5])
plt.savefig('figures/all_vs_independent_SNPS_across_TADS_hg19.png')
plt.close()

###################################
# Visualize distribution of independent SNPs
###################################
# Observe location of independent SNPs in TADs
TADdict = range(NUM_BINS)
TADdict = dict.fromkeys(TADdict)
for key in TADdict.iterkeys():
    # Create a list of zeros of len(NUM_BINS)
    TADdict[key] = 0

# Loop through each TAD and SNP to fill the histogram
for key in KG_TADs.iterkeys():
    if key not in ['NoTAD', 'Boundary', 'X', 'Y']:
        for tadkey, tadval in KG_TADs[key].iteritems():
            if len(tadval) > 0:
                tadname = parse_TAD_name(tadkey)
                for snp in range(len(tadval)):
                    if tadval.ix[snp, :]['RS'] in rs_set:
                        snp_loc = parse_SNP_position(tadval.ix[snp, :])
                        bin_assign = ID_TAD_bins(tadname, NUM_BINS, snp_loc)
                        TADdict[bin_assign] += 1

SNPLocations = pd.DataFrame.from_dict(TADdict, orient='index')

# Plot locations of SNPs
SNPLocations.plot()
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('Independent SNP location inside TADs')
plt.grid(True)
plt.savefig('figures/Independent_SNPlocations_1000Genomes_hg19.png')
plt.close()
