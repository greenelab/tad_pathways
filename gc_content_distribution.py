"""
(C) 2016 Gregory Way
gc_content_distribution.py

Description:
Observe GC Content distributions and locations across TADs (hg19)

Usage:
Is called by 'ANALYSIS.sh'

Output:
GC Content distribution and histogram across TADs
"""

import sys
sys.path.insert(1, 'bin/')
from util import parse_TAD_name
import random
import pandas as pd
import matplotlib.pyplot as plt
plt.figure.max_open_warning = 0
from Bio import SeqIO

random.seed(123)
####################################
# Define Functions
####################################


def load_fasta(chrom, fasta_loc):
    """
    Retrieve fasta file

    Arguments:
    :param chrom: the chromosome of interest (format 'chr#')
    :param fasta_loc: the location that stores the hg19 fasta files

    Output:
    fasta file for the given chromosome
    """

    chrom_fa = fasta_loc + chrom + '.fa'
    record = SeqIO.read(open(chrom_fa), 'fasta')

    nucleotides = str(record.seq)
    # FASTA file has upper and lowercase letters
    # lower case = repetative elements (ignore for now)
    nucleotides = nucleotides.upper()

    return nucleotides


def parse_fasta_tad(tadname, nuc):
    """
    Retrieve nucleotide sequence of tad

    Arguments:
    :param tadname: of the format tadid:tadstart-tadend
    :param nuc: the nucleotides for the same chromosome as the TAD

    Output:
    a list of nucleotides in the TAD
    """

    tadname = parse_TAD_name(tadname)
    start = tadname[1]
    end = tadname[2]

    tadsequence = nuc[start:end]

    return tadsequence


def determine_gc_content(seq, start, end):
    """
    Determine the gc content for a given sequence given bin coordinates

    Arguments:
    :param seq: a nucleotide sequence
    :param start: where to subset the sequence
    :param end: where to subset the sequence

    Output:
    A count of GC content within the specific coordinates
    """

    if len(seq) < end:
        end = len(seq)
    subset_seq = seq[start:end]

    A = subset_seq.count('A')
    C = subset_seq.count('C')
    G = subset_seq.count('G')
    T = subset_seq.count('T')
    N = subset_seq.count('N')

    known_length = A + C + G + T

    if known_length + N == N:
        GC = 0.5
    else:
        GC = (G + C) / float(known_length)

    return float(GC)


def split_TAD_bins(tadlength, num_bins):
    """
    Return a list of coordinates to partition the TAD

    Arguments:
    :param tadlength: how long the TAD is
    :param num_bins: how many bins to split

    Output:
    a list of tuples with the locations for starting and ending bins
    """

    import random

    avgbin = tadlength / num_bins
    remainder = tadlength % num_bins
    if remainder > 0:
        randadd = random.sample(range(0, num_bins), remainder)
    else:
        randadd = []

    return_list = []
    current_idx = 0
    for binID in range(0, num_bins):
        next_idx = current_idx + avgbin
        if binID in randadd:
            return_list.append([current_idx + 1, next_idx + 1])
            current_idx = next_idx + 1
        else:
            return_list.append([current_idx + 1, next_idx])
            current_idx = next_idx
    return return_list


####################################
# Load Constants
####################################
NUM_BINS = 50
FASTA_LOC = 'data/hg19_fasta/'

####################################
# Load Data
####################################
TADdictgenes = pd.read_pickle('index/Gene_Assignments_In_TADs_hg19.p')

####################################
# Analysis
####################################
gc_distribution = []
for chrom in TADdictgenes.keys():
    if chrom in ['Boundary', 'Y']:
        continue
    lookup_chrom = 'chr' + str(chrom)
    # Get the TAD sequence:
    chr_sequence = load_fasta(lookup_chrom, FASTA_LOC)
    for TAD in TADdictgenes[chrom].keys():
        tad_sequence = parse_fasta_tad(TAD, chr_sequence)
        tad_length = len(tad_sequence)

        # Get the TAD bins
        tad_bins = split_TAD_bins(tad_length, NUM_BINS)

        # Now, loop over the TAD bins and extract GC Content
        tad_gc = []
        coord_idx = 0

        for coord in tad_bins:
            start = tad_bins[coord_idx][0]
            end = tad_bins[coord_idx][1]
            gc = determine_gc_content(tad_sequence, start, end)
            tad_gc.append(gc)
            coord_idx += 1

        gc_distribution.append(tad_gc)

####################################
# Plot
####################################
gc_distribution = pd.DataFrame(gc_distribution)

fig, ax = plt.subplots()
gc_distribution.boxplot(column=range(0, 50), notch=True, positions=range(0, 50),
                        return_type='dict')
ax.set_xticklabels(0)
ax.tick_params(axis=u'both', which=u'both', length=0)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#fig.set_size_inches(15, 10)
plt.grid(False)
plt.xlabel('Bins (Normalized TAD length)')
plt.ylabel('GC Percentage')
plt.title('GC Content Across TADs')
plt.savefig('figures/gc_distribution.png')
