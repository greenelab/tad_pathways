"""
2016 Gregory Way
scripts/viz/generate_TAD_locations.py

Description:
Summarizes the location of genomic elements across TADs

Usage:
Is called by 'scripts/viz/visualize.sh' which is run inside of
'scripts/run_pipeline.sh'. This particular script will output the location
of genomic elements in a given input TAD

      python scripts/viz/generate_TAD_locations.py --TAD-Boundary 'hESC'

Output:
Several .png plots in "figures/GENOME/" and chisquare analyses of the
"rightness" of SNPs in TADs and protein coding genes near boundaries.
"""

import argparse

import pandas as pd
from scipy.stats import chisquare
import matplotlib.pyplot as plt
import seaborn as sns

from tad_util.util import assign_bin

plt.figure.max_open_warning = 0
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("paper", rc={"font.size": 20, "axes.titlesize": 20,
                             "axes.labelsize": 20, "xtick.labelsize": 12,
                             "ytick.labelsize": 12})

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--TAD-Boundary', help='boundary cell type. The'
                    'options can be "hESC", "IMR90", "mESC", or "cortex"')
args = parser.parse_args()

# Load Constants
NUM_BINS = 50
TAD_CELL = args.TAD_Boundary
xlab = [''] * NUM_BINS
for x in range(0, 50, 10):
    xlab[x] = x

if TAD_CELL in ['hESC', 'IMR90']:
    GENOME = 'hg19'
elif TAD_CELL in ['mESC', 'cortex']:
    GENOME = 'mm9'
else:
    raise ValueError('Please input: "hESC", "IMR90", "mESC", or "cortex"')

# Input files
BASE_FH = GENOME + '_' + TAD_CELL

SNP_INDEX = 'index/SNP_index_' + BASE_FH + '.tsv.bz2'
GENE_INDEX = 'index/GENE_index_' + BASE_FH + '.tsv.bz2'
REPEAT_INDEX = 'index/REPEATS_index_' + BASE_FH + '.tsv.bz2'

# Output files
FIG_BASE = 'figures/' + GENOME + '/'

SNP_HIST = FIG_BASE + 'snp_hist_' + BASE_FH + '.png'
SNP_DIST = FIG_BASE + 'snp_tad_distribution_' + BASE_FH + '.png'
SNP_CHROM = FIG_BASE + 'chrom/snp_tad_distribution_chrom_' + BASE_FH
SNP_CHI = 'results/tad_snp_rightness_chi_' + BASE_FH + '.txt'

GENE_HIST = FIG_BASE + 'gene_hist_' + BASE_FH + '.png'
GENE_CHROM = FIG_BASE + 'chrom/gene_tad_distribution_chrom_' + BASE_FH
GENE_TYPE = FIG_BASE + 'gene_type_' + BASE_FH + '_'
GENE_CHI = 'results/tad_gene_boundary_chi_' + BASE_FH + '.txt'

TAD_LEN_HIST = FIG_BASE + TAD_CELL + '_tad_length_histogram.png'

REPEAT_HIST = FIG_BASE + 'repeat_hist_' + BASE_FH + '.png'
REPEAT_TYPE = FIG_BASE + 'repeat_type_' + BASE_FH + '_'
REPEAT_DIST = FIG_BASE + 'repeat_type_all_distribution_' + BASE_FH + '.png'

# Load Data
gene_types_df = pd.read_table('tables/gene_classification.tsv')
snp_df = pd.read_table(SNP_INDEX, index_col=0)
gene_df = pd.read_table(GENE_INDEX, index_col=0)
repeat_df = pd.read_table(REPEAT_INDEX, index_col=0)

#########################
# PART 1 - SNPs
#########################
# Process SNP dataframe
snp_df = snp_df[snp_df['TAD_id'] != 'Boundary']
bin_s = snp_df.apply(lambda x: assign_bin(x, bins=NUM_BINS, ID='SNP'), axis=1)
snp_df = snp_df.assign(tad_bin=bin_s)

# Histogram of number of SNPs per TAD
snp_hist = snp_df.groupby('TAD_id').tad_bin.count()
p = sns.distplot(snp_hist, kde=False, color=sns.xkcd_rgb['medium green'])
sns.despine()
p.set(ylabel='', xlabel='Number of SNPs')
p.set_title('Histogram of number of SNPs in TADs')
plt.tight_layout()
plt.savefig(SNP_HIST)
plt.close()

# Distribution of SNPs across TADs
summary_snp = snp_df['tad_bin'].value_counts(sort=False)
p = sns.pointplot(x=summary_snp.index, y=summary_snp,
                  color=sns.xkcd_rgb["medium green"], scale=0.5)
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Number of SNPs', xlabel='TAD Bins')
p.set_title('Distribution of SNPs across TADs')
plt.tight_layout()
plt.savefig(SNP_DIST)
plt.close()

# Chromosome-specific distribution
snp_chrom = snp_df.groupby('chromosome').tad_bin.value_counts(sort=False).\
  unstack(level=0)
for chrom, chrom_df in snp_chrom.iteritems():
    fh = SNP_CHROM + '_' + str(chrom) + '.png'
    p = sns.pointplot(x=chrom_df.index, y=chrom_df,
                      color=sns.xkcd_rgb["medium green"], scale=0.5)
    sns.despine()
    p.set(xticklabels=xlab)
    p.set(ylabel='Number of SNPs', xlabel='TAD Bins')
    p.set_title('SNP Distribution in Chromosome ' + str(chrom))
    plt.tight_layout()
    plt.savefig(fh)
    plt.close()

# SNPs appear to be more concentrated on the right side of TADs
snp_side = [snp_df[snp_df['tad_bin'] < 25].shape[0],
            snp_df[snp_df['tad_bin'] >= 25].shape[0]]
tad_snp_sig = chisquare(snp_side)

with open(SNP_CHI, 'w') as chisq_fh:
    chisq_fh.write('SNPs in the left vs. right of ' + TAD_CELL + ' TAD\n')
    chisq_fh.write('left,right\n')
    chisq_fh.write(','.join('%s' % x for x in snp_side) + '\n')
    chisq_fh.write(','.join('%s' % x for x in tad_snp_sig))

#########################
# PART 2 - Genes
#########################
# Process genes
gene_df = gene_df[gene_df['TAD_id'] != 'Boundary']
bin_assign_gene = gene_df.apply(lambda x: assign_bin(x, bins=NUM_BINS,
                                                     ID='gene'), axis=1)
gene_df = gene_df.assign(tad_bin=bin_assign_gene)
gene_df = gene_df[gene_df['tad_bin'] != -1]

# Histogram of number of Genes per TAD
gene_hist = gene_df.groupby('TAD_id').tad_bin.count()
p = sns.distplot(gene_hist, kde=False, color=sns.xkcd_rgb['medium green'])
sns.despine()
p.set(ylabel='', xlabel='Number of Genes')
p.set_title('Histogram of number of Genes in TADs')
plt.tight_layout()
plt.savefig(GENE_HIST)
plt.close()

# Chromosome specific distribution of genes across TADs
gene_chrom = gene_df.groupby('chromosome').tad_bin.value_counts(sort=False).\
  unstack(level=0)
for chrom, chrom_df in gene_chrom.iteritems():
    fh = GENE_CHROM + '_' + str(chrom) + '.png'
    p = sns.pointplot(x=chrom_df.index, y=chrom_df,
                      color=sns.xkcd_rgb["medium green"], scale=0.5)
    sns.despine()
    p.set(xticklabels=xlab)
    p.set(ylabel='Number of Genes', xlabel='TAD Bins')
    p.set_title('Gene Distribution in Chromosome ' + str(chrom))
    plt.tight_layout()
    plt.savefig(fh)
    plt.close()

# Gene-type specific distribution across TADs
gene_types_df = gene_types_df[gene_types_df[GENOME] == 1]
summary_gene_classes = []
for idx, gene in gene_types_df.iterrows():
    gene_class = gene['gene_class']
    gene_type = gene['gene_type']

    if gene_class in ['tr_gene', 'ig_gene', 'tr_pseud', 'ig_pseud']:
        gene_type = gene_types_df[gene_types_df['gene_class'] ==
                                  gene_class]['gene_type']
        gene_sub_df = gene_df[gene_df['gene_type'].isin(gene_type)]
        plot_title = gene_class

        if gene_class in summary_gene_classes:
            continue
        else:
            summary_gene_classes.append(gene_class)

    elif gene_class == 'std' and gene_type != 'all':
        gene_sub_df = gene_df[gene_df['gene_type'] == gene_type]
        plot_title = gene_type
    elif gene_type == 'all':
        gene_sub_df = gene_df
        plot_title = 'Distribution of Genes across TADs'

    fh = GENE_TYPE + plot_title + '.png'
    sum_gene_sub = gene_sub_df['tad_bin'].value_counts(sort=False).sort_index()

    p = sns.pointplot(x=sum_gene_sub.index, y=sum_gene_sub,
                      color=sns.xkcd_rgb["medium green"], scale=0.5)
    sns.despine()
    p.set(xticklabels=xlab)
    p.set(ylabel='Number of Genes', xlabel='TAD Bins')
    p.set_title(plot_title)
    plt.tight_layout()
    plt.savefig(fh)
    plt.close()

# Chisquare of genes on TAD boundaries
protein_coding = gene_df[gene_df['gene_type'] == 'protein_coding']
bin_list = list(range(NUM_BINS))[0:2] + list(range(NUM_BINS))[-2:]
boundary_df = protein_coding[protein_coding['tad_bin'].isin(bin_list)]
num_genes_b = boundary_df.shape[0]
num_genes_c = protein_coding.shape[0] - num_genes_b
chi_test = [num_genes_b, num_genes_c]
exp = protein_coding.shape[0] / NUM_BINS
bound_chi = chisquare(chi_test, f_exp=[exp * len(bin_list),
                                       exp * (NUM_BINS - len(bin_list))])

with open(GENE_CHI, 'w') as chisq_fh:
    chisq_fh.write('Genes at boundaries vs. center of ' + TAD_CELL + ' TAD\n')
    chisq_fh.write('bound,center\n')
    chisq_fh.write(','.join('%s' % x for x in chi_test) + '\n')
    chisq_fh.write(','.join('%s' % x for x in bound_chi))

# Distribution of TAD length
tad_df = gene_df[['TAD_id', 'TAD_start', 'TAD_end']].drop_duplicates()
tad_hist = (tad_df['TAD_end'] - tad_df['TAD_start']) / 1000

p = sns.distplot(tad_hist, kde=False, color=sns.xkcd_rgb['medium green'])
p.set(ylabel='', xlabel='TAD Length (kb)')
p.set_title('Histogram of number of Genes in TADs')
sns.despine()
plt.tight_layout()
plt.savefig(TAD_LEN_HIST)
plt.close()

#########################
# PART 3 - Repeats
#########################
# Process Repeats
repeat_df = repeat_df.fillna('Boundary')
repeat_df = repeat_df[repeat_df['TAD_id'] != 'Boundary']
bin_assign_repeat = repeat_df.apply(lambda x: assign_bin(x, bins=NUM_BINS,
                                                         ID='repeat'), axis=1)
repeat_df = repeat_df.assign(tad_bin=bin_assign_repeat)
repeat_df = repeat_df[repeat_df['tad_bin'] != -1]

# Histogram of number of repeats per TAD
repeat_hist = repeat_df.groupby('TAD_id').tad_bin.count()
p = sns.distplot(repeat_hist, kde=False, color=sns.xkcd_rgb['medium green'])
sns.despine()
p.set(ylabel='', xlabel='Number of Repeats')
p.set_title('Histogram of number of Repeats in TADs')
plt.tight_layout()
plt.savefig(REPEAT_HIST)
plt.close()

# Distribution of different classes of repeats across TADs
for repeat_type in repeat_df['repeat'].unique():
    if '?' not in repeat_type:
        repeat_fh = repeat_type.replace('/', '_')
        fh = REPEAT_TYPE + repeat_fh + '.png'
        rep_sub_df = repeat_df[repeat_df['repeat'] == repeat_type]
        sum_rep = rep_sub_df['tad_bin'].value_counts(sort=False).sort_index()
        p = sns.pointplot(x=sum_rep.index, y=sum_rep,
                          color=sns.xkcd_rgb["medium green"], scale=0.5)
        sns.despine()
        p.set(xticklabels=xlab)
        p.set(ylabel='Number of Repeats', xlabel='TAD Bins')
        p.set_title(repeat_type + ' Distribution')
        plt.tight_layout()
        plt.savefig(fh)
        plt.close()

# Distribution of all repeats
sum_repeat = repeat_df['tad_bin'].value_counts(sort=False).sort_index()
p = sns.pointplot(x=sum_repeat.index, y=sum_repeat,
                  color=sns.xkcd_rgb["medium green"], scale=0.5)
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Number of Repeats', xlabel='TAD Bins')
p.set_title('All Repeats Distribution')
plt.tight_layout()
plt.savefig(REPEAT_DIST)
plt.close()
