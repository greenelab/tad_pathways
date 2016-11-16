"""
2016 Gregory Way
scripts/construct_evidence.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line

        'python scripts/construct_evidence.py -t <trait> -e <eqtl> -g <gwas>'

Output:
a .csv evidence file (two columns: [gene, type of evidence])
"""

import os
import argparse
import pandas as pd

# Load Command Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--trait', help='symbol for trait data')
parser.add_argument("-e", "--eqtl", help="location of eqtl data")
parser.add_argument("-g", "--genelist", help="location of gwas genelist")
parser.add_argument("-p", "--pathway", help="pathway of interest")
args = parser.parse_args()

# Load Constants
trait = args.trait
trait_file = os.path.join('data', 'gestalt',
                          '{}_complete_gestalt.tsv'.format(trait))
eqtl_file = args.eqtl
gwas_file = args.genelist
pathway = args.pathway
output_file = os.path.join('tad_pathway', '{}_gene_evidence.csv'.format(trait))

# Load Data
pathway_genes = pd.read_csv(trait_file, delimiter='\t')
pathway_genes = pathway_genes.loc[pathway_genes['go_name'] == pathway, :]
pathway_genes = pathway_genes['symbol'].tolist()

eqtl_genes = pd.read_csv(eqtl_file, delimiter='\t')
eqtl_genes = list(set(eqtl_genes['Gene'].tolist()))

gwas_genes_df = pd.read_table(gwas_file)
gwas_genes = gwas_genes_df.MAPPED_GENE

# Process GWAS genes
full_gwas_genes = []
for gene in gwas_genes:
    gene = gene.replace(' - ', ', ').split(', ')
    full_gwas_genes += gene

# Process all genes
all_genes = set(full_gwas_genes + eqtl_genes + pathway_genes)
all_assignments = []

# Logic to determine genic evidence
for gene in all_genes:
    if gene in full_gwas_genes:
        assignment = 'gwas'
        if gene in eqtl_genes:
            assignment = '{}_eqtl'.format(assignment)
        if gene in pathway_genes:
            assignment = '{}_tad'.format(assignment)
        if assignment == 'gwas_eqtl_tad':
            assignment = 'all'
    else:
        if gene in eqtl_genes:
            assignment = 'eqtl'
            if gene in pathway_genes:
                assignment = '{}_tad'.format(assignment)
        else:
            assignment = 'tad'

    all_assignments.append(assignment)

evidence = pd.DataFrame([all_genes, all_assignments]).T
evidence.columns = ['gene', 'evidence']
evidence.to_csv(output_file, delimiter=',', index=False)
