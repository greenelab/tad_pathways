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


def load_tad(TAD_LOC):
    """
    Parameters:
    TAD_LOC: the location of the TAD boundary coordinate file

    Output:
    an empty dictionary with TAD information as keys. Key structure
    is 'TAD_ID:TADstart-TADend'.
    """

    import pandas as pd

    snp_header = ['chromosome', 'start', 'end']
    TAD_df = pd.read_table(TAD_LOC, names=snp_header)
    TAD_df['chromosome'] = TAD_df['chromosome'].map(lambda x: x[3:])

    return TAD_df


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
    gene_info - a row from the input gtf file.

    Output:
    Parsed gene information from the gtf file [gene_type, gene_name]
    """

    import pandas as pd

    info = gene_info['info'].split(';')
    build_info = {}
    for attrb in range(1, (len(info) - 1)):
        attrb_s = info[attrb].split(' ')
        attrb_name = attrb_s[1].strip('"')
        new_attrb = attrb_s[2].strip('"')
        build_info[attrb_name] = new_attrb
    return_list = [build_info['gene_type'], build_info['gene_name']]

    return pd.Series(return_list)


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
