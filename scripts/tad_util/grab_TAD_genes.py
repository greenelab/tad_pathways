"""
(C) 2016 Gregory Way
grab_TAD_genes.py

Description:
Uses intemediate files and gene assignment in TADs index file to curate
a final list of all trait-specific genes that fall in signal TADs

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python build_TAD_genes.py'

Output:
Trait specific .tsv files of one column each indicating all the genes that fall
in signal TADs
"""

from optparse import OptionParser
import pandas as pd
import math

####################################
# Load Command Arguments
####################################
parser = OptionParser()  # Load command line options
parser.add_option("-f", "--TAD-GWAS-File", dest="TG_file",
                  help="Location of TAD/GWAS data", type="string")
parser.add_option("-o", "--output-TAD-file", dest="output",
                  help="Name of the full genelist output file", type="string")
parser.add_option("-g", "--output-GWAS-file", dest="output_gwas",
                  help="Name of the GWAS genelist output file", type="string")
(options, args) = parser.parse_args()

####################################
# Load Constants
####################################
TAD_GWAS_LOC = options.TG_file
OUTPUT_FH = options.output
OUTPUT_GWAS_FH = options.output_gwas

####################################
# Define Functions
####################################


def buildTADkey(gwas_snp):
    """
    gwas_snp - a single row subset of the TAD-GWAS input file

    output - The lookup info in the TAD gene dictionary
    i.e. [Chromosome, TAD_ID:TAD_Start-TAD_End]
    """
    chrom = gwas_snp['chrom'][3:]
    start = int(gwas_snp['TADStart'])
    end = int(gwas_snp['TADEnd'])
    tad_num = int(gwas_snp['TADidx'])
    output = str(tad_num) + ':' + str(start) + '-' + str(end)
    return [str(chrom), output]


def appendtogenelist(newgenes, genelist):
    """
    newgenes - a list of genes to append to genelist
    genelist - the current list of genes

    output - genelist with appended newgenes
    """
    for gene in newgenes:
        genelist.append(gene)
    return genelist

####################################
# Load picklefile (dictionary)
####################################
TADdictgenes = pd.read_pickle('index/Gene_Assignments_In_TADs_hg19.p')

####################################
# Load data
####################################
TADGWAS = pd.read_csv(TAD_GWAS_LOC, sep='\t')

####################################
# Get all TAD based genes
####################################
# GWAS genes are all genes nearest to each significant SNP
gwas_genes = []
# GWAS signal are the set of genes that fall in enriched TADs
GWASsignal = set()

for tadix in range(len(TADGWAS)):
    # Subset the pandas File
    tad_sub = TADGWAS.ix[tadix, :]

    # Get the gene identified in the GWAS
    gene = str(tad_sub['gene']).split(',')
    if type(gene) is list:
        gwas_genes = appendtogenelist(gene, gwas_genes)
    else:
        gwas_genes.append(gene)

    # If the gene does fall in a TAD
    if not math.isnan(tad_sub['TADidx']):

        # Build the key to lookup TAD in dict and lookup
        chrom, TADkey = buildTADkey(tad_sub)
        FullTADinfo = TADdictgenes[chrom][TADkey]

        # if the TAD is not empty
        if len(FullTADinfo) > 0:

            # Get the genes in the TAD and add to the set
            GWAS_genes = FullTADinfo.gene[0:].tolist()
            for g in GWAS_genes:
                GWASsignal.add(g)

# Save the TAD based genes to a file
full_gene_list = pd.Series(list(GWASsignal))
with open(OUTPUT_FH, 'wt') as w_file:
    full_gene_list.to_csv(w_file, index=False, header=False, sep='\t')

# Save the nearest genes to a file
gwas_gene_list = pd.Series(list(set(gwas_genes)))
with open(OUTPUT_GWAS_FH, 'wt') as w_file:
    gwas_gene_list.to_csv(w_file, index=False, header=False, sep='\t')
