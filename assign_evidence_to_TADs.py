"""
(C) 2016 Gregory Way
assign_evidence_to_TADs.py

Description:
Takes in genes and evidence support and assigns each gene to the TAD

Usage:
python assign_evidence_to_TADs.py -e path/to/evidence.csv -t path/to/TADassign

Output:
Trait specific .tsv files of one column each indicating all the genes that fall
in signal TADs
"""
import sys
sys.path.insert(1, 'bin/')
from optparse import OptionParser
import pandas as pd
import math

####################################
# Load Command Arguments
####################################
parser = OptionParser()  # Load command line options
parser.add_option("-e", "--evidence-fh", dest="evidence_fh",
                  help="Location of evidence file", type="string")
parser.add_option("-t", "--path-to-tad", dest="tad",
                  help="path to TAD/GWAS data", type="string")
parser.add_option("-o", "--output-fh", dest="out_fh",
                  help="location to write results", type="string")
(options, args) = parser.parse_args()

####################################
# Load Constants
####################################
EVIDENCE_FH = options.evidence_fh
TAD_GWAS_FH = options.tad
OUT_FH = options.out_fh

####################################
# Load picklefile (dictionary)
####################################
TADdictgenes = pd.read_pickle('index/Gene_Assignments_In_TADs_hg19.p')

####################################
# Load data
####################################
TADGWAS = pd.read_csv(TAD_GWAS_FH, sep='\t')
EVIDENCE = pd.read_csv(EVIDENCE_FH)

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


def parse_ev_key(tadkey):
    """
    tadkey - the key for the evidence dictionary

    output - Information stored in the evidence key
    i.e. [ID, Chromosome, Start, End, UCSC lookup]
    """
    chrom = tadkey.split(':')
    ID = chrom[1]
    start = chrom[2].split('-')
    end = start[1]
    start = start[0]
    chrom = chrom[0]
    ucsc = chrom + ':' + start + '-' + end
    return [ID, chrom, start, end, ucsc]

####################################
# Investigate each significant TADs
####################################
evidence_dict = {}
for tadix in range(len(TADGWAS)):
    # Subset the pandas File
    tad_sub = TADGWAS.ix[tadix, :]

    # Get the gene identified in the GWAS
    gene = str(tad_sub['gene']).split(',')
    gene = [x.strip(' ') for x in gene]

    # If the gene does fall in a TAD
    if not math.isnan(tad_sub['TADidx']):

        # Build the key to lookup TAD in dict and lookup
        chrom, TADkey = buildTADkey(tad_sub)
        FullTADinfo = TADdictgenes[chrom][TADkey]

        # Build evidence key
        e_key = 'chr' + str(chrom) + ':' + TADkey

        if e_key not in evidence_dict.keys():
            evidence_dict[e_key] = []

        # if the TAD is not empty
        if len(FullTADinfo) > 0:

            # Get the genes in the TAD
            TAD_genes = FullTADinfo.gene[0:].tolist()

            # Subset the evidence data frame to genes in TAD
            evidence_sub = EVIDENCE.ix[EVIDENCE['gene'].isin(TAD_genes), :]

            # Loop over each of the rows
            TAD_gene_evidence = []
            for ev in range(0, evidence_sub.shape[0]):
                gene = evidence_sub['gene'].tolist()[ev]
                ev_type = evidence_sub['evidence'].tolist()[ev]
                TAD_gene_evidence.append([gene, ev_type])

            # Add this gene and evidence to the dictionary
            evidence_dict[e_key].append(TAD_gene_evidence)

####################################
# Write out results to file
####################################
with open(OUT_FH, 'w') as out_fh:
    out_fh.write('ID\tCHROMOSOME\tSTART\tEND\tUCSC\tGENE\tEVIDENCE\n')
    for tadkey in evidence_dict.keys():
        info = parse_ev_key(tadkey)
        ID = info[0]
        chrom = info[1]
        start = info[2]
        end = info[3]
        ucsc = info[4]
        loop_idx = 0
        for evidence in evidence_dict[tadkey][0]:
            gene = evidence[0]
            e = evidence[1]
            if loop_idx == 0:
                out_fh.write('\t'.join([ID, chrom, start, end, ucsc,
                                        gene, e]) + '\n')
            else:
                out_fh.write('\t\t\t\t\t' + gene + '\t' + e + '\n')
            loop_idx += 1
