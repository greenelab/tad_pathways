"""
2016 Gregory Way
scripts/visualize_gc_and_divergence.py

Description:
Observe GC Content distributions and locations across TADs

Usage:
Is called by 'scripts/visualize.sh' which is run inside of
'scripts/run_pipeline.sh':

        python gc_content_distribution.py --TAD-Boundary 'hESC'

Output:
GC Content distribution and histogram across TADs as a .png file
"""

import random
import argparse

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

from tad_util.util import assign_bin, load_tad

plt.figure.max_open_warning = 0
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("paper", rc={"font.size": 20, "axes.titlesize": 20,
                             "axes.labelsize": 20, "xtick.labelsize": 12,
                             "ytick.labelsize": 12})
random.seed(123)

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--TAD-Boundary', help='boundary cell type. The'
                    'options can be "hESC", "IMR90", "mESC", or "cortex"')
args = parser.parse_args()

# Define Constants
NUM_BINS = 50
xlab = [''] * NUM_BINS
for x in range(0, 50, 10):
    xlab[x] = x
TAD_CELL = args.TAD_Boundary

# Generate file names depending on input
if TAD_CELL in ['hESC', 'IMR90']:
    GENOME = 'hg19'
    BASE_DIR = 'data/hg/'
elif TAD_CELL in ['mESC', 'cortex']:
    GENOME = 'mm9'
    BASE_DIR = 'data/mm/'
else:
    raise ValueError('Please input: "hESC", "IMR90", "mESC", or "cortex"')

BASE_FH = GENOME + '_' + TAD_CELL
REPEAT_INDEX = 'index/REPEATS_index_' + BASE_FH + '.tsv.bz2'
GC_FIG_FH = 'figures/' + GENOME + '/gc_distribution_' + BASE_FH + '.png'
DIV_FIG_FH = 'figures/' + GENOME + '/repeat_divergence_' + BASE_FH + '.png'
ALU_FIG_FH = 'figures/' + GENOME + '/alu_divergence_' + BASE_FH + '.png'
TAD_LOC = BASE_DIR + TAD_CELL + '_domains_' + GENOME + '.bed'

if GENOME == 'hg19':
    FASTA_LOC = 'data/hg/hg19_fasta/'
elif GENOME == 'mm9':
    FASTA_LOC = 'data/mm/mm9_fasta/'


def load_fasta(chrom, fasta_loc):
    """
    Retrieve fasta file

    Arguments:
    :param chrom: the chromosome of interest (format 'chr#')
    :param fasta_loc: the location that stores the hg19 fasta files

    Output:
    fasta file for the given chromosome
    """

    chrom_fa = fasta_loc + 'chr' + chrom + '.fa'
    record = SeqIO.read(open(chrom_fa), 'fasta')

    nucleotides = str(record.seq)
    # FASTA file has upper and lowercase letters
    # lower case = repetative elements
    nucleotides = nucleotides.upper()

    return nucleotides


def split_TAD_bins(tadlength, num_bins):
    """
    Return a list of coordinates to partition the TAD

    Arguments:
    :param tadlength: how long the TAD is
    :param num_bins: how many bins to split

    Output:
    a list of tuples with the locations for starting and ending bins
    """

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
            return_list.append((current_idx + 1, next_idx + 1))
            current_idx = next_idx + 1
        else:
            return_list.append((current_idx + 1, next_idx))
            current_idx = next_idx
    return return_list


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

    import collections
    
    start = int(start)
    end = int(end)

    subset_seq = seq[start:end]

    c = collections.Counter(subset_seq)

    # Known length will be zero if entire sequence is not known or `N`
    known_length = sum(c[base] for base in 'ACTG')
    GC = (c['G'] + c['C']) / known_length if known_length else 0.5

    return GC


def get_gc_content(tad, seq, bins):
    """
    Determine the gc content of all TADs across TAD bins

    Arguments:
    :param tad: a row in a TAD boundary DataFrame
    :param seq: the SeqIO.seq object for the chromosome fasta
    :param bins: int, number of bins to distribute gc content of TAD sequence

    Output:
    A pandas series of gc content of all TADs across bins
    """

    tad_start = int(tad['start'])
    tad_end = int(tad['end'])

    tad_length = tad_end - tad_start

    # Get the TAD bins
    tad_bins = split_TAD_bins(tad_length, bins)

    tad_sequence = seq[tad_start:tad_end]

    # Now, loop over the TAD bins and extract GC Content
    tad_gc = []

    for coord in tad_bins:
        start, end = coord

        gc = determine_gc_content(tad_sequence, start, end)
        tad_gc.append(gc)

    return pd.Series(tad_gc)

# Load Data
tad_df = load_tad(TAD_LOC)
repeat_df = pd.read_table(REPEAT_INDEX, index_col=0)
repeat_df = repeat_df.ix[~pd.isnull(repeat_df['TAD_id'])]
bin_r = repeat_df.apply(lambda x: assign_bin(x, bins=NUM_BINS, ID='repeat'),
                        axis=1)
repeat_df = repeat_df.assign(tad_bin=bin_r)
repeat_df = repeat_df[repeat_df['tad_bin'] != -1]
alu_df = repeat_df[repeat_df['repeat'] == 'SINE/Alu']

# Plot divergence
p = sns.boxplot(x=repeat_df['tad_bin'], y=repeat_df['div'],
                color=sns.xkcd_rgb['medium green'])
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Repeat Divergence', xlabel='TAD Bins')
p.set_title('')
plt.tight_layout()
plt.savefig(DIV_FIG_FH)
plt.close()

# ALU divergence
p = sns.boxplot(x=alu_df['tad_bin'], y=alu_df['div'],
                color=sns.xkcd_rgb['medium green'])
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='ALU Repeat Divergence', xlabel='TAD Bins')
p.set_title('')
plt.tight_layout()
plt.savefig(ALU_FIG_FH)
plt.close()

# Plot GC content
gc_content_df = pd.DataFrame()
for chrom in tad_df['chromosome'].unique():
    tad_sub = tad_df[tad_df['chromosome'] == chrom]
    fasta = load_fasta(str(chrom), FASTA_LOC)
    gc_content = tad_sub.apply(lambda x: get_gc_content(x, seq=fasta,
                                                        bins=NUM_BINS), axis=1)
    gc_content_df = gc_content_df.append(gc_content, ignore_index=True)

p = sns.boxplot(data=gc_content_df, color=sns.xkcd_rgb['medium green'])
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='GC Content', xlabel='TAD Bins')
p.set_title('')
plt.tight_layout()
plt.savefig(GC_FIG_FH)
plt.close()
