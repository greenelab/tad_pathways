"""
2016 Gregory Way
scripts/visualize_genomic_elements.py

Description:
Summarizes the location of genomic elements across TADs

Usage:
Is called by 'scripts/visualize.sh' which is run inside of
'scripts/run_pipeline.sh'. This particular script will output the location
of genomic elements in a given input TAD

      python scripts/visualize_gwas_distribution.py

Output:
Two .pdf plots in "figures/" describing the GWAS signal distribution across
hESC and IMR90 TADs
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tad_util.util import load_tad, assign_bin


def curate_gwas_elements(tad_df, input_df):
    """
    Loop through TAD boundaries to assign GWAS snps to TADs

    Arguments:
    :param tad_df: pandas dataframe, TAD boundary locations
    :param input_df: pandas dataframe, the gwas locations

    Output:
    Large summary dataframes of the TAD assignments for GWAS signals
    """

    big_out_df = pd.DataFrame()

    for tad in tad_df.itertuples():
        tad_id, chrom, start, end = list(tad)
        start, end = int(start), int(end)

        chm_sub_df = input_df[input_df['chromosome'] == chrom]

        elem_sub_df = chm_sub_df[(chm_sub_df['position'] >= start) &
                                 (chm_sub_df['position'] < end)]

        # Assign TAD information
        elem_sub_df['TAD_id'] = tad_id
        elem_sub_df['TAD_start'] = start
        elem_sub_df['TAD_end'] = end

        big_out_df = big_out_df.append(elem_sub_df, ignore_index=True)

    return big_out_df

# Set plotting defaults
plt.figure.max_open_warning = 0
sns.set_style("whitegrid")
sns.set_style("ticks")
sns.set_context("paper", rc={"font.size": 20, "axes.titlesize": 20,
                             "axes.labelsize": 20, "xtick.labelsize": 12,
                             "ytick.labelsize": 12})

# Load Constants
num_bins = 50
xlab = [''] * num_bins
for x in range(0, 50, 10):
    xlab[x] = x

# Load and process data
gwas_df = pd.read_table('data/gwas_catalog_hg19.tsv')

hesc_df = load_tad('data/hg/hESC_domains_hg19.bed')
imr_df = load_tad('data/hg/IMR90_domains_hg19.bed')

gwas_df = gwas_df.dropna(subset=['CHR_POS', 'CHR_ID'])
gwas_df = gwas_df[gwas_df['CHR_ID'] != 'Not Mapped']
gwas_df = gwas_df.assign(chromosome=gwas_df.CHR_ID.str[3:].astype(str))
gwas_df = gwas_df.assign(position=gwas_df.CHR_POS.astype(int))

# Map content
hesc_tad_df = curate_gwas_elements(hesc_df, gwas_df)
imr_tad_df = curate_gwas_elements(imr_df, gwas_df)

# Assign bins to gwas signals
hesc_bins = hesc_tad_df.apply(lambda x: assign_bin(x, bins=num_bins, ID='SNP'),
                              axis=1)
imr_bins = imr_tad_df.apply(lambda x: assign_bin(x, bins=num_bins, ID='SNP'),
                            axis=1)
hesc_tad_df = hesc_tad_df.assign(tad_bin=hesc_bins)
imr_tad_df = hesc_tad_df.assign(tad_bin=imr_bins)

# Plot and output distributions
summary_hesc = hesc_tad_df['tad_bin'].value_counts(sort=False)
summary_imr = imr_tad_df['tad_bin'].value_counts(sort=False)

p = sns.pointplot(x=summary_hesc.index, y=summary_hesc,
                  color=sns.xkcd_rgb["medium green"], scale=0.5)
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Number of GWAS signals', xlabel='TAD Bins')
p.set_title('Distribution of SNPs across TADs')
plt.tight_layout()
plt.savefig('figures/hesc_gwas_distribution.pdf')
plt.close()

p = sns.pointplot(x=summary_imr.index, y=summary_imr,
                  color=sns.xkcd_rgb["medium green"], scale=0.5)
sns.despine()
p.set(xticklabels=xlab)
p.set(ylabel='Number of GWAS signals', xlabel='TAD Bins')
p.set_title('Distribution of SNPs across TADs')
plt.tight_layout()
plt.savefig('figures/imr90_gwas_distribution.pdf')
plt.close()
