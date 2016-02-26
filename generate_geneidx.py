"""
(C) 2016 Gregory Way
generate_geneidx.py

Description:
Map reference genes to TADs based on genomic locations. For this iteration we
are considering the hg19 assembly.

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python generate_geneidx.py'

Output:
The script will output a pickle file that will store a python dictionary to
enable fast lookups of genes in TADs. The script will also output a csv file
of all genes that span TADs.
"""

import sys
sys.path.insert(1, 'bin/')
from util import initialize_TAD_dictionary
from util import parse_gene_gtf
from util import parse_TAD_name
import csv
import pandas as pd
import pickle

####################################
# Load Constants
####################################
REF_GENE_GTF = 'data/gencode.v19.annotation.gtf'
TAD_LOC = 'data/hESC_domains_hg19.bed'

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
    if chm not in ['Boundary']:
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
with open('output/SpannedGenesAcrossTADs.csv', 'w') as span_fh:
    csv_file = csv.writer(span_fh)
    csv_file.writerow(['gene', 'type', 'chrom', 'start', 'end', 'strand'])
    for gene in spanned_genes:
        csv_file.writerow(gene)

####################################
# Save dictionary in picklefile
####################################
pickle.dump(TADdict, open('output/Gene_Assignments_In_TADs_hg19.p', 'wb'))
