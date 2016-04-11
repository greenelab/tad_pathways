"""
(C) 2016 Gregory Way
visualize_genes.py

Description:
Observe gene distributions and locations across TADs (hg19)

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python visualize_genes.py'

Output:
1) Several matplotlib plots in the 'figures/' folder
2) Several gene location plots across TADs in 'figures/gene/'
3) Chisquare analysis results for gene start sites near boundaries
"""

import sys
sys.path.insert(1, 'bin/')
from util import parse_TAD_name
from util import ID_TAD_bins
from util import parse_gene_info
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare

####################################
# Load picklefile
####################################
TADdictgenes = pd.read_pickle('output/Gene_Assignments_In_TADs_hg19.p')
print TADdictgenes['5']['1054:76404244-76804244']
####################################
# Load Constants
####################################
NUM_BINS = 50
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
    if chromkey not in ['Boundary']:
        # Create a list of zeros of len(NUM_BINS)
        Chromosome_Specific_TAD[chromkey] = [0] * (NUM_BINS)

# Determine significance of 'leftness' of TADs
TAD_near_boundary = [0, 0]
TAD_nb_protein = [0, 0]
TAD_nb_rRNA = [0, 0]
left_bin = NUM_BINS / 2

# Obtain the gene position information
for key, value in TADdictgenes.iteritems():
    if key not in ['Boundary']:
        for tadkey, tadval in TADdictgenes[key].iteritems():
            tadkey = parse_TAD_name(tadkey)
            for gene in range(len(tadval)):
                # gene_info stores [gene_type, chromosome, start_pos]
                gene_info = parse_gene_info(tadval.ix[gene, :])
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

with open('output/TADBoundary_Chisquare.txt', 'w') as chisq_fh:
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
    g_loc.plot()
    plt.xlabel('Bins (Normalized TAD length)')
    plt.ylabel('Frequency')
    plt.title('Gene location inside TADs\n' + g)
    plt.grid(True)
    plt.savefig(nameplot)
    plt.close()

#####################################
# Plot Chromosome Specific Line Graphs
#####################################
# Separate Graphs
for chrom in Chromosome_Specific_TAD.keys():
    if chrom not in ['Boundary']:
        nameplot = 'figures/chrom/GeneDistribution_chrom_' + chrom + '_' + str(NUM_BINS) + 'bins.png'
        chrom_loc = chromLocations[chrom]
        chrom_loc.plot()
        plt.xlabel('Bins (Normalized TAD length)')
        plt.ylabel('Frequency')
        plt.title('Gene location (Start Site) inside TADs\nChromosome ' + str(chrom))
        plt.grid(True)
        plt.savefig(nameplot)
        plt.close()

# Combined Graph
for chrom in Chromosome_Specific_TAD.keys():
    if chrom not in ['Boundary']:
        chrom_loc = chromLocations[chrom]
        chrom_loc.plot()
        plt.xlabel('Bins (Normalized TAD length)')
        plt.ylabel('Frequency')
        plt.title('Gene location (Start Site) inside TADs\nby Chromosome')
        plt.grid(True)

plt.savefig('figures/chrom/GeneDistribution_chromosomes.png')
plt.close()

####################################
# Obtain the number of genes in each TAD
####################################
howmanygenesinTAD = []
TADlength = []
for key, value in TADdictgenes.iteritems():
    if key not in ['Boundary']:
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
plt.grid(True)
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
plt.grid(True)
addtext = 'mean = ' + str(round(np.mean(TADlength), 2)) + '\n' + 'std = ' + \
          str(round(np.std(TADlength), 2))
plt.text(4000, .0021, addtext, verticalalignment='top',
         horizontalalignment='right')
plt.savefig('figures/TAD_length_histogram_hg19.png')
