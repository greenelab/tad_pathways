"""
(C) 2016 Gregory Way
visualize_TAD_locations.py

Description:
Observe SNP/gene/repeat distributions and locations across TADs (hg19)

Usage:
Is called by 'ANALYSIS.sh'

Output:
1) SNP distribution and histogram across TADs
2) Gene distribution
    A) Several matplotlib plots in the 'figures/' folder
    B) Several gene location plots across TADs in 'figures/gene/'
    C) Chisquare analysis results for gene start sites near boundaries
3) Repeat distribution
    A) Several matplotlib plots in the 'figures/' folder
    B) Several repeat location plots across TADs in 'figures/repeats/'
"""

import sys
sys.path.insert(1, 'bin/')
from util import parse_TAD_name
from util import parse_SNP_position
from util import parse_gene_info
from util import parse_repeat_info
from util import ID_TAD_bins
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.figure.max_open_warning = 0
from scipy.stats import chisquare

####################################
# Load Constants
####################################
NUM_BINS = 50
# Gene types described here: http://www.gencodegenes.org/gencode_biotypes.html
gene_types = ['protein_coding', 'pseudogene', 'miRNA', 'snRNA', 'snoRNA',
              'lincRNA', 'rRNA', 'all', 'misc_RNA', 'sense_intronic',
              'processed_transcript', 'sense_overlapping', 'IG_pseudogene',
              'TR_pseudogene', 'TR_gene', 'polymorphic_pseudogene',
              'antisense', 'IG_gene', '3prime_overlapping_ncrna']
TR_genes = ['TR_V_gene', 'TR_J_gene', 'TR_D_gene', 'TR_C_gene']
IG_genes = ['IG_V_gene', 'IG_C_gene', 'IG_C_gene', 'IG_J_gene', 'IG_D_gene']
TR_pseud = ['TR_V_pseudogene', 'TR_J_pseudogene']
IG_pseud = ['IG_V_pseudogene', 'IG_C_pseudogene', 'IG_J_pseudogene']

####################################
# Load Data
####################################
KG_TADs = pd.read_pickle('index/SNPindex_1000G_hg19.p')
TADdictgenes = pd.read_pickle('index/Gene_Assignments_In_TADs_hg19.p')
TADdictrepeats = pd.read_pickle('index/Repeats_In_TADs_hg19.p')

####################################
# PART 1 - SNPs ####################
####################################

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
    if key in ['NoTAD', 'Boundary', 'X', 'Y']:
        continue
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

with open('results/TAD_SNP_rightness_Chisquare.txt', 'w') as chisq_fh:
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
plt.grid(False)
plt.savefig('figures/SNPlocations_TADs_1000G_hg19.png')
plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['NoTAD', 'Boundary', 'X', 'Y']:
        continue
    nameplot = 'figures/chrom/SNPdistribution_chrom_' + chrom + '_' + \
               str(NUM_BINS) + 'bins.png'
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

####################################
# PART 2 - Genes ###################
####################################

####################################
# Observe location of genes in TADs
####################################
# Initialize a dictionary that will store location information for gene types
TADgenepos = dict.fromkeys(gene_types)
for genekey in gene_types:
    # Create a list of zeros of len(NUM_BINS)
    TADgenepos[genekey] = [0] * (NUM_BINS)

# Initialize a chromosome specific dictionary holding bin specific info
Chromosome_Specific_TAD = dict.fromkeys(TADdictgenes.keys())
for chromkey in Chromosome_Specific_TAD.keys():
    if chromkey in ['Boundary']:
        continue
    # Create a list of zeros of len(NUM_BINS)
    Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Determine significance of 'leftness' of TADs
TAD_near_boundary = [0, 0]
TAD_nb_protein = [0, 0]
TAD_nb_rRNA = [0, 0]
left_bin = NUM_BINS / 2

# Also, we're interested in the genes near TAD boundaries
gene_TAD_boundary = []
# Obtain the gene position information
for key, value in TADdictgenes.iteritems():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey, tadval in TADdictgenes[key].iteritems():
        tadkey = parse_TAD_name(tadkey)
        for gene in range(len(tadval)):
            # gene_info stores [gene_type, chromosome, start_pos]
            gene_info = parse_gene_info(tadval.ix[gene, :])
            gene_name = tadval.ix[gene, :]['gene']
            binID = ID_TAD_bins(tadkey, NUM_BINS, gene_info[2],
                                IDtype='gene')
            if binID == -1:
                continue

            # Add the bin number to the appropriate dictionary
            gene_type = gene_info[0]
            if gene_type in TR_genes:
                TADgenepos['TR_gene'][binID] += 1
            elif gene_type in IG_genes:
                TADgenepos['IG_gene'][binID] += 1
            elif gene_type in TR_pseud:
                TADgenepos['TR_pseudogene'][binID] += 1
            elif gene_type in IG_pseud:
                TADgenepos['IG_pseudogene'][binID] += 1
            else:
                TADgenepos[gene_type][binID] += 1

            # Add the bin number of all genes
            TADgenepos['all'][binID] += 1

            # Determine which half of bin gene belongs in and store
            if binID in [0, 1, 48, 49]:
                if gene_type == 'protein_coding':
                    TAD_nb_protein[0] += 1
                elif gene_type == 'rRNA':
                    TAD_nb_rRNA[0] += 1
                TAD_near_boundary[0] += 1
                gene_TAD_boundary.append([gene_name, binID])

            else:
                if gene_type == 'protein_coding':
                    TAD_nb_protein[1] += 1
                elif gene_type == 'rRNA':
                    TAD_nb_rRNA[1] += 1
                TAD_near_boundary[1] += 1

            # Add a count for chromsome specificity
            Chromosome_Specific_TAD[key][binID] += 1

# Convert to pandas dataframe
geneLocations = pd.DataFrame(TADgenepos)
chromLocations = pd.DataFrame(Chromosome_Specific_TAD)
gene_TAD_boundary = pd.DataFrame(gene_TAD_boundary, columns=['gene', 'bin'])

# Write gene TAD boundary to csv
gene_TAD_boundary.to_csv('tables/genes_at_TAD_boundaries.tsv', sep='\t')

####################
# Perform ChiSquare test of TAD near boundary
####################
# Expected values for each bin
ex_all = sum(TAD_near_boundary) / NUM_BINS
ex_pro = sum(TAD_nb_protein) / NUM_BINS
ex_rRNA = sum(TAD_nb_rRNA) / NUM_BINS

TAD_nb_sig = chisquare(TAD_near_boundary, f_exp=[ex_all * 4, ex_all * 46])
TAD_nb_sig_pro = chisquare(TAD_nb_protein, f_exp=[ex_pro * 4, ex_pro * 46])
TAD_nb_sig_rRNA = chisquare(TAD_nb_rRNA, f_exp=[ex_rRNA * 4, ex_rRNA * 46])

with open('results/TADBoundary_Chisquare.txt', 'w') as chisq_fh:
    chisq_fh.write('Near boundaries of all genes\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_near_boundary) + '\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_nb_sig) + '\n')
    chisq_fh.write('Near boundaries of protein coding genes\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_nb_protein) + '\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_nb_sig_pro) + '\n')
    chisq_fh.write('Near boundaries of rRNA genes\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_nb_rRNA) + '\n')
    chisq_fh.write(','.join('%s' % x for x in TAD_nb_sig_rRNA) + '\n')

#####################################
# Plot line graphs for each gene type
#####################################
for g in gene_types:
    nameplot = 'figures/gene/genelocation_' + g + '_' + str(NUM_BINS) + '.png'
    g_loc = geneLocations[g]
    fig, ax = plt.subplots(1)
    g_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Gene Location across TADs\n' + g)
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
# Separate Graphs
for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['Boundary']:
        continue
    nameplot = 'figures/chrom/GeneDistribution_chrom_' + chrom + '_' + \
               str(NUM_BINS) + 'bins.png'
    chrom_loc = chromLocations[chrom]
    fig, ax = plt.subplots(1)
    chrom_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Gene location (Start Site) across TADs\nChromosome ' +
              str(chrom))
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

# Combined Graph
fig, ax = plt.subplots(1)
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('Gene location (Start Site) across TADs\nby Chromosome')
plt.grid(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['Boundary']:
        continue
    chrom_loc = chromLocations[chrom]
    chrom_loc.plot()

plt.savefig('figures/chrom/GeneDistribution_chromosomes.png')
plt.close()

####################################
# Obtain the number of genes in each TAD
####################################
howmanygenesinTAD = []
TADlength = []
for key, value in TADdictgenes.iteritems():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey, tadval in TADdictgenes[key].iteritems():
        howmanygenesinTAD.append(len(tadval))
        tadname = parse_TAD_name(tadkey)
        TADlength.append((tadname[2] - tadname[1])/1000)

####################################
# Output Histograms
####################################
# Plot number of genes in tads
n, bins, patches = plt.hist(howmanygenesinTAD, NUM_BINS, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('Number of genes')
plt.ylabel('Frequency')
plt.title('Distribution of number of genes in TADs')
plt.axis([min(howmanygenesinTAD), max(howmanygenesinTAD), 0, 0.06])
plt.grid(False)
addtext = 'mean = ' + str(round(np.mean(howmanygenesinTAD), 2)) + '\n' + \
          'std = ' + str(round(np.std(howmanygenesinTAD), 2))
plt.text(130, .055, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig('figures/TAD_Gene_histogram_hg19.png')
plt.close()

# Plot how long TADs are
n, bins, patches = plt.hist(TADlength, NUM_BINS * 2, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('TAD length (kb)')
plt.ylabel('Frequency')
plt.title('How long are TADs?')
plt.axis([min(TADlength), max(TADlength), 0, 0.0023])
plt.grid(False)
addtext = 'mean = ' + str(round(np.mean(TADlength), 2)) + '\n' + 'std = ' + \
          str(round(np.std(TADlength), 2))
plt.text(4000, .0021, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig('figures/TAD_length_histogram_hg19.png')

####################################
# PART 3 - Repeats #################
####################################

####################################
# Process repeat picklefile
####################################
all_repeat_types = []
for key in TADdictrepeats.keys():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey in TADdictrepeats[key].keys():
        repeat_elements = TADdictrepeats[key][tadkey]
        repeat_elements = repeat_elements['repeat'].tolist()
        all_repeat_types = all_repeat_types + repeat_elements

all_repeat_types = set(all_repeat_types)
repeat_types = []
for repeat in all_repeat_types:
    if '?' not in repeat:
        repeat_types.append(repeat)

repeat_types.append('all')

####################################
# Observe location of repeats in TADs
####################################
# Initialize a dictionary that will store location information for repeats
TADrepeatpos = dict.fromkeys(repeat_types)
for repeatkey in repeat_types:
    # Create a list of zeros of len(NUM_BINS)
    TADrepeatpos[repeatkey] = [0] * (NUM_BINS)

# Initialize a chromosome specific dictionary holding bin specific info
Chromosome_Specific_TAD = dict.fromkeys(TADdictrepeats.keys())
for chromkey in Chromosome_Specific_TAD.keys():
    if chromkey not in ['Boundary']:
        # Create a list of zeros of len(NUM_BINS)
        Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Obtain the repeat position information
for key, value in TADdictrepeats.iteritems():
    if key in ['Y', 'Boundary']:
        continue
    for tadkey, tadval in TADdictrepeats[key].iteritems():
        tadkey = parse_TAD_name(tadkey)
        for repeat in range(len(tadval)):
            # repeat_info stores [gene_type, chromosome, start_pos]
            repeat_info = parse_repeat_info(tadval.ix[repeat, :])
            repeat_name = repeat_info[0]

            if '?' not in repeat_name:
                binID = ID_TAD_bins(tadkey, NUM_BINS, repeat_info[2],
                                    IDtype='gene')
                if binID == -1:
                    continue

                # if the start is exactly on a TAD boundary, put it in the
                # previous bin
                if binID == 50:
                    binID = 49

                # Add to bin of specific repeat
                TADrepeatpos[repeat_name][binID] += 1

                # Add the bin number of all repeats
                TADrepeatpos['all'][binID] += 1

                # Add a count for chromsome specificity
                Chromosome_Specific_TAD[key][binID] += 1

# Convert to pandas dataframe
repeatLocations = pd.DataFrame(TADrepeatpos)
chromLocations = pd.DataFrame(Chromosome_Specific_TAD)

#####################################
# Plot line graphs for each gene type
#####################################
for r in repeat_types:
    fn = r.replace('/', '_')
    nameplot = 'figures/repeats/repeat_' + fn + '_' + str(NUM_BINS) + '.png'
    g_loc = repeatLocations[r]
    fig, ax = plt.subplots(1)
    g_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Repeat Location across TADs\n' + r)
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
# Separate Graphs
for chrom in Chromosome_Specific_TAD.keys():
    if chrom in ['Y', 'Boundary']:
        continue
    nameplot = 'figures/chrom/RepeatDistribution_chrom_' + chrom + '_' + \
               str(NUM_BINS) + 'bins.png'
    chrom_loc = chromLocations[chrom]
    fig, ax = plt.subplots(1)
    chrom_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Repeat location across TADs\nChromosome ' +
              str(chrom))
    plt.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    plt.savefig(nameplot)
    plt.close()

# Combined Graph
fig, ax = plt.subplots(1)
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('Frequency')
plt.title('Repeat location (Start Site) across TADs\nby Chromosome')
plt.grid(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

for chrom in Chromosome_Specific_TAD.keys():
    if chrom not in ['Y', 'Boundary']:
        chrom_loc = chromLocations[chrom]
        chrom_loc.plot()

plt.savefig('figures/chrom/RepeatDistribution_chromosomes.png')
plt.close()

####################################
# Obtain the number of Repeats in each TAD
####################################
howmanyrepeatsinTAD = []
tad_names = []
for key, value in TADdictrepeats.iteritems():
    if key not in ['Y', 'Boundary']:
        for tadkey, tadval in TADdictrepeats[key].iteritems():
            howmanyrepeatsinTAD.append(len(tadval))
            tad_names.append(tadkey)

howmanyrepeatsinTAD = np.array(howmanyrepeatsinTAD)

####################################
# Output Histograms
####################################
# Plot number of genes in tads
n, bins, patches = plt.hist(howmanyrepeatsinTAD, NUM_BINS, normed=1,
                            facecolor='green', alpha=0.75)
plt.xlabel('Number of Repeat Elements')
plt.ylabel('Frequency')
plt.title('Distribution of number of Repeat Elements in TADs')
plt.axis([min(howmanyrepeatsinTAD), max(howmanyrepeatsinTAD), 0, 0.001])
plt.grid(False)
addtext = 'mean = ' + str(round(np.mean(howmanyrepeatsinTAD), 2)) + '\n' + \
          'std = ' + str(round(np.std(howmanyrepeatsinTAD), 2))
plt.text(5000, .00055, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig('figures/TAD_Repeat_histogram_hg19.png')
plt.close()
