"""
(C) 2016 Gregory Way
parse_gestalt.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line 'python parse_gestalt.py -t <TRAIT>'

    -t <TRAIT> is predefined by the manual WebGestalt step (see README)
    e.g. <TRAIT> is BMD for Bone Mineral Density

Output:
pickle files with pathway specific dictionaries and summary of pathway p values
"""

from optparse import OptionParser
import pandas as pd
import StringIO

###################
# Define Functions
###################


def read_gestalt(fh):
    """
    Read in the gestalt data and parse pathways

    Arguments:
    :param fh: the location of the file to read

    Output:
    1) Pandas Dataframe of entire TRAIT pathways
    2) Pathway summary text file
    """

    # Read in the entire file as string
    with open(fh, 'r') as trait_fh:
        gestalt = trait_fh.read()

    # Remove trailing whitespace from file
    gestalt = gestalt.rstrip()

    # Triple new lines demarcate pathways, skip header info
    stanzas = gestalt.split('\n\n\n')[1:]

    # Initialize empty list to append pathway DataFrames
    stanzas_df = list()

    # Loop through each pathway and combine into large DataFrame
    pathway_names = []
    pathway_adjp = []
    pathway_numgenes = []
    for pathway in stanzas:
        pathway_file = StringIO.StringIO(pathway)

        # Parse pathway and enrichment information
        go_info = next(pathway_file).rstrip('\n').split('\t')
        enrichment = next(pathway_file).rstrip('\n').split(';')
        go_class = go_info[0]
        go_name = go_info[1]
        go_id = go_info[2]
        adj_p = float(enrichment[5][5:])

        # Read in the io file
        pathway_header = ['symbol', 'NA', 'symbol2', 'name', 'entrez_id',
                          'ensemble_id']
        pathway_df = pd.read_table(pathway_file, skip_blank_lines=True,
                                   names=pathway_header)

        # Add GO information to the pathway DataFrame
        pathway_df['go_class'] = go_class
        pathway_df['go_name'] = go_name
        pathway_df['go_id'] = go_id
        pathway_df['adjP'] = adj_p
        stanzas_df.append(pathway_df)

        # Add to general information
        pathway_names.append(go_name)
        pathway_adjp.append(adj_p)
        pathway_numgenes.append(pathway_df.shape[0])

    return_df = pd.concat(stanzas_df)
    gen_info = pd.DataFrame([pathway_names, pathway_adjp, pathway_numgenes]).T
    gen_info.columns = ['name', 'adjP', 'num_genes']

    return return_df, gen_info

# Run the command if called by name
if __name__ == '__main__':
    ####################################
    # Load Command Arguments
    ####################################
    parser = OptionParser()  # Load command line options
    parser.add_option("-t", "--trait", dest="trait",
                      help="symbol for trait data", type="string")
    (options, args) = parser.parse_args()

    ####################################
    # Load Constants
    ####################################
    TRAIT = options.trait
    TRAIT_FH = 'data/gestalt/' + TRAIT + '_gestalt.tsv'
    OUT_TRAIT_FH = 'data/gestalt/' + TRAIT + '_complete_gestalt.tsv'
    OUT_P_FH = 'data/gestalt/' + TRAIT + '_pvalues.tsv'

    ####################################
    # Load and Save Data
    ####################################
    gestalt_data, gestalt_p = read_gestalt(TRAIT_FH)

    # Write files
    gestalt_data = gestalt_data.sort_values(by='adjP')
    gestalt_data.to_csv(OUT_TRAIT_FH, sep='\t')
    gestalt_p = gestalt_p.sort_values(by='adjP')
    gestalt_p.to_csv(OUT_P_FH, sep='\t')
