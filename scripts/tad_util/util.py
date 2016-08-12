"""
util.py
Author: Greg Way

Description:
Location for useful methods for interacting with TADs/SNPs/Genes

Usage:
Import into python scripts - do not envoke directly

Output:
See specific functions
"""

##############################
# Define functions
##############################


def initialize_TAD_dictionary(TAD_LOC):
    """
    Parameters:
    TAD_LOC: the location of the TAD boundary coordinate file

    Output:
    an empty dictionary with TAD information as keys. Key structure
    is 'TAD_ID:TADstart-TADend'.
    """

    import csv

    TADdict = {}
    TAD_idx = 0  # Will become the TAD based ID
    with open(TAD_LOC) as tadfile:
        tadreader = csv.reader(tadfile, delimiter='\t')
        for row in tadreader:
            # Build the name for the 1st key
            chrom = row[0][3:]

            # Initialize a dictionary with chromosome key
            if chrom not in TADdict:
                TADdict[str(chrom)] = {}

            # Build the name for the 2nd key
            TADStart = str(row[1])
            TADEnd = str(row[2])
            TADkey = str(TAD_idx) + ':' + TADStart + '-' + TADEnd
            TADdict[str(chrom)][TADkey] = []
            TAD_idx += 1

    # Add a place to store SNPs that fall in TAD boundaries
    TADdict['Boundary'] = []

    # Add Y chromosome
    TADdict['Y'] = []

    return TADdict


def parse_TAD_name(tad_name):
    """
    Parameters:
    tad_name - a single dictionary key with format 'TAD_ID:TADstart-TADend'

    Output:
    a list with the parsed TAD information
    """
    tad_name = str(tad_name)
    tad_name_ID_pos = tad_name.split(':')
    tad_position = tad_name_ID_pos[1].split('-')
    return [tad_name_ID_pos[0], int(tad_position[0]), int(tad_position[1])]


def parse_gene_gtf(gene_info):
    """
    Parameters:
    gene_info - a row from the input gtf file. File format is described here:
    "useast.ensemble.org/info/website/upload/gff.html?redirect=no"

    Output:
    Parsed gene information from the gtf file with the format:
    [gene_name, gene_type, chromosome, gene_start, gene_end, strand]
    """
    chrom = gene_info[0][3:]
    l_start = gene_info[3]
    l_end = gene_info[4]
    strand = gene_info[6]
    gene_info = str(gene_info[8])  # The 8th value has ';' separated attributes
    attrb = gene_info.split(';')
    gene_name = attrb[4].split(' ')[2].strip('"')
    gene_type = attrb[2].split(' ')[2].strip('"')
    return (gene_name, gene_type, str(chrom), int(l_start), int(l_end), strand)


def parse_SNP_position(snp_info):
    '''
    Parameters:
    snp_info - a row within a single pandas DataFrame in the pickle dictionary

    Output:
    The position of the SNP
    '''
    return int(snp_info['POSITION'])


def ID_TAD_bins(tadkey, bins, start_pos, IDtype='SNP'):
    '''
    Parameters:
    tadkey: the parsed output of parse_TAD_name()
    bins: the number of bins to distribute each gene into
    start_pos: starting position of element
    IDtype: either 'SNP' (default) or 'Gene'

    Output:
    the appropriate bin number for the given input tadkey and SNP.
    If IDtype='gene', the output will be -1 if the start site is not within TAD
    (i.e. - the gene is overlapping boundary regions)
    '''
    # get information from the tadkey
    tadstart = int(tadkey[1])
    tadend = int(tadkey[2])

    # Determine the bin the gene belongs in
    diff = float(start_pos - int(tadstart))
    if IDtype == 'SNP':
        bin_assign = int(diff/(tadend - tadstart) * bins)
    elif IDtype == 'gene':
        if diff < 0:
            bin_assign = -1
        elif tadend < start_pos:
            bin_assign = -1
        else:
            bin_assign = int(diff/(tadend - tadstart) * bins)

    # Return the bin assignment
    return bin_assign


def parse_gene_info(gene_info):
    '''
    Parameters:
    gene_info - a row within a single pandas DataFrame in the pickle dictionary

    Output:
    gene information in the form [type, chromosome, start_site]
    '''
    gtype = gene_info['type']
    chrom = gene_info['chrom']
    strnd = gene_info['strand']
    if strnd == '+':
        start = int(gene_info['start'])
    else:
        start = int(gene_info['end'])
    output = [gtype, chrom, start]
    return output


def parse_repeat_info(repeat_info):
    '''
    Parameters:
    repeat_info - a row in repeat pandas DataFrame

    Output:
    repeat information [type, chromosome, start]
    '''
    repeat_type = repeat_info['repeat']
    chrom = repeat_info['chrom']
    start = repeat_info['begin']
    return [repeat_type, chrom, start]
