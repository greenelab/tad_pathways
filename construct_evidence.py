"""
(C) 2016 Gregory Way
construct_evidence.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line 'python construct_evidence.py -t <trait> -e <eqtl> -g <gwas>'

Output:
a .csv evidence file (two columns: [gene, type of evidence])
"""

from optparse import OptionParser
import pandas as pd
import pickle

####################################
# Load Command Arguments
####################################
parser = OptionParser()  # Load command line options
parser.add_option("-t", "--trait", dest="trait",
                  help="symbol for trait data", type="string")
parser.add_option("-e", "--eQTL-location", dest="eqtl",
                  help="location of eqtl data", type="string")
parser.add_option("-g", "--GWAS-location", dest="gwas",
                  help="location of gwas data", type="string")
parser.add_option("-p", "--pathway", dest="pathway",
                  help="pathway of interest", type="string")
(options, args) = parser.parse_args()

####################################
# Load Constants
####################################
TRAIT = options.trait
TRAIT_PICKLE = 'data/gestalt/' + TRAIT + '_pathways.p'
EQTL_FH = options.eqtl
GWAS_FH = options.gwas
PATHWAY = options.pathway
OUTPUT_FH = 'tad_pathway/' + TRAIT + '_gene_evidence.csv'

####################################
# Load Data
####################################
pathway_genes = pd.read_pickle(TRAIT_PICKLE)
pathway_genes = pathway_genes[PATHWAY][8]['gene'].tolist()

eqtl_genes = pd.read_csv(EQTL_FH, delimiter='\t')
eqtl_genes = list(set(eqtl_genes['Gene'].tolist()))

gwas_genes = pd.read_csv(GWAS_FH, header=None)
gwas_genes = gwas_genes[0].tolist()

####################################
# Process Genes
####################################
all_genes = set(gwas_genes + eqtl_genes + pathway_genes)
all_assignments = []

# Logic to determine genic evidence
for gene in all_genes:
    if gene in gwas_genes:
        assignment = 'gwas'
        if gene in eqtl_genes:
            assignment = assignment + '_eqtl'
        if gene in pathway_genes:
            assignment = assignment + '_tad'
        if assignment == 'gwas_eqtl_tad':
            assignment = 'all'

    if gene in eqtl_genes:
        if gene not in gwas_genes:
            assignment = 'eqtl'
            if gene in pathway_genes:
                assignment = assignment + '_tad'

    if gene in pathway_genes:
        if gene not in gwas_genes:
            if gene not in eqtl_genes:
                assignment = 'tad'

    all_assignments.append(assignment)

evidence = pd.DataFrame([all_genes, all_assignments]).T
evidence.columns = ['gene', 'evidence']
evidence.to_csv(OUTPUT_FH, delimiter=',', index=False)
