"""
(C) 2016 Gregory Way
generate_index_files.py

Description:
1) Map GWAS SNP information to TADs based on genomic locations. We are
considering the GRCh37/hg19 genome build for both SNP positions/TAD locations
respectively.

2) Map reference genes to TADs based on genomic locations. For this iteration we
are considering the hg19 assembly.

3) Map repeatmasker elements to TADs based on genomic locations (hg19)

Usage:
Is called by 'ANALYSIS.sh'

Output:
The script will output a pickle file that will store a python dictionary to
enable fast lookups of SNP/gene/repeat information as they fall in TADs
"""

import sys
sys.path.insert(1, 'bin/')
from util import initialize_TAD_dictionary
from util import parse_TAD_name
from util import parse_gene_gtf
import pandas as pd
import pickle
import csv

####################################
# Load Constants
####################################
TAD_LOC = 'data/hESC_domains_hg19.bed'
REF_SNPS = 'data/common-snps.tsv'
REF_GENE_GTF = 'data/gencode.v19.annotation.gtf'
TAD_LOC = 'data/hESC_domains_hg19.bed'
REPEAT_FH = 'data/hg19.fa.out.tsv'

####################################
# PART 1 - SNPs ####################
####################################

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
pickle.dump(TADdict, open('index/SNPindex_1000G_hg19.p', 'wb'))

####################################
# PART 2 - Genes ###################
####################################

####################################
# Initialize a dictionary with the structure:
# 1st key: chromosome, 2nd key: TAD ID + Location
####################################
TADdict = initialize_TAD_dictionary(TAD_LOC)

####################################
# Load and process data into the dict
####################################
spanned_genes = []
only_spanned_gene_names = []
with open(REF_GENE_GTF) as refGene_fh:
    refGene = csv.reader(refGene_fh, delimiter='\t')
    next(refGene)  # Skip over header information
    next(refGene)
    next(refGene)
    next(refGene)
    next(refGene)

    for row in refGene:

        # We are only interested in genes here
        if row[2] == 'gene':
            g_info = parse_gene_gtf(row)

            # Extract gene info and assign
            gene_name = g_info[0]
            gene_chrom = g_info[2]
            gene_start = g_info[3]
            gene_end = g_info[4]

            # Currently, there are no Y or M chromosome TADs
            if gene_chrom == 'Y' or gene_chrom == 'M':
                continue

            # Loop over all the TADs in the given chromosome
            # to find which TAD the gene is found in
            for key in TADdict[gene_chrom].iterkeys():
                t_info = parse_TAD_name(key)  # Extract TAD info
                tad_id = t_info[0]
                tad_start = t_info[1]
                tad_end = t_info[2]

                # If the gene is entirely in the TAD, it belongs in
                if gene_start >= tad_start and gene_end <= tad_end:
                    TADdict[gene_chrom][key].append(g_info)
                    continue

                # If the gene is partially in, still put it in
                if (tad_start <= gene_start <= tad_end or
                        tad_start <= gene_end <= tad_end):
                    TADdict[gene_chrom][key].append(g_info)

                    if gene_name in only_spanned_gene_names:
                        continue
                    else:
                        spanned_genes.append(g_info)

                    # Keep track of spanning genes but only store once
                    only_spanned_gene_names.append(gene_name)
                    TADdict['Boundary'].append(g_info)
                    continue

####################################
# Make each value a pandas dataframe for easy lookup
####################################
for chm in TADdict.iterkeys():
    if chm not in ['Boundary', 'Y']:
        for key, value in TADdict[chm].iteritems():
            # Some TADs do not have any genes in them
            if len(value) > 0:
                TADdict[chm][key] = pd.DataFrame(value,
                                                 columns=('gene', 'type',
                                                          'chrom', 'start',
                                                          'end', 'strand'))
    else:
        TADdict[chm] = pd.DataFrame(TADdict[chm], columns=('gene', 'type',
                                                           'chrom', 'start',
                                                           'end', 'strand'))

# Write the spanned genes as a text file
with open('tables/SpannedGenesAcrossTADs.tsv', 'w') as span_fh:
    csv_file = csv.writer(span_fh, delimiter='\t')
    csv_file.writerow(['gene', 'type', 'chrom', 'start', 'end', 'strand'])
    for gene in spanned_genes:
        csv_file.writerow(gene)

####################################
# Save dictionary in picklefile
####################################
pickle.dump(TADdict, open('index/Gene_Assignments_In_TADs_hg19.p', 'wb'))

####################################
# PART 3 - Repeat Elements #########
####################################

####################################
# Load Data - Messy Repeat Data
####################################
RepeatElements = pd.read_csv(REPEAT_FH, delimiter='\t',
                             skiprows=[0, 1, 2], header=None,
                             names=['chrom', 'begin', 'end', 'repeat'],
                             index_col=False, usecols=[5, 6, 7, 11])

####################################
# Initialize a dictionary with the structure:
# 1st key: chromosome, 2nd key: TAD ID + Location
####################################
TADdict = initialize_TAD_dictionary(TAD_LOC)

####################################
# Load and process data into the dict
####################################
for idx, row in RepeatElements.iterrows():

    # Extract repeat info and assign
    chrom = row['chrom'][3:]
    start = str(row['begin'])
    end = str(row['end'])

    if chrom not in TADdict.keys():
        continue

    if chrom in ['Y', 'Boundary']:
        continue

    if (')' not in start) and (')' not in end):
        # Loop over all the TADs in the given chromosome
        # to find which TAD the gene is found in
        for key in TADdict[chrom].iterkeys():
            t_info = parse_TAD_name(key)  # Extract TAD info
            tad_start = t_info[1]
            tad_end = t_info[2]

            # If the gene is entirely in the TAD, it belongs in
            if int(start) >= tad_start and int(end) <= tad_end:
                TADdict[chrom][key].append([chrom, int(start),
                                            int(end), row['repeat']])
                continue

            # If the gene is partially in, still put it in
            if (tad_start <= int(start) <= tad_end or
                    tad_start <= int(end) <= tad_end):
                TADdict[chrom][key].append([chrom, int(start),
                                            int(end), row['repeat']])
                continue

####################################
# Make each value a pandas dataframe for easy lookup
####################################
for chm in TADdict.iterkeys():
    if chm not in ['Y', 'Boundary']:
        for key, value in TADdict[chm].iteritems():
            # Some TADs do not have any genes in them
            if len(value) > 0:
                TADdict[chm][key] = pd.DataFrame(value,
                                                 columns=('chrom', 'begin',
                                                          'end', 'repeat'))
####################################
# Save dictionary in picklefile
####################################
pickle.dump(TADdict, open('index/Repeats_In_TADs_hg19.p', 'wb'))
