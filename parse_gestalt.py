"""
(C) 2016 Gregory Way
parse_gestalt.py

Description:
Take as input the .tsv results from the WebGestalt analysis and outputs trait
specific TAD pathway genes

Usage:
Command line 'python parse_gestalt.py -t <TRAIT>'

Output:
pickle files with pathway specific dictionaries and summary of pathway p values
"""

from optparse import OptionParser
import pandas as pd
import pickle

###################
# Define Functions
###################


def read_gestalt(fh):
    """
    Read in the gestalt data and parse pathways

    Arguments:
    :param fh: the location of the file to read

    Output:
    Two dictionaries:
    1) Pathway names as keys, holding enrichment info
    2) Pathway names as keys, holding p value info
    """

    pathway_dict = {}
    pathway_pvalue = {}
    pathway_number = 0
    gene_DataFrame = pd.Series()
    pathway_name = ''
    with open(fh, 'r') as trait_fh:
        # Skip 10 header lines
        for x in xrange(10):
            next(trait_fh)

        for line in trait_fh:
            line = line.strip('\n').split('\t')

            # Test the line type and perform different operations accordingly
            # Skip the line if the first element is blank
            if len(line[0]) == 0:
                continue

            # This line has the GO category and pathway name and ID
            if len(line[0].split(' ')) == 2:

                # Add the pathway to the dictionary
                if pathway_number > 0:
                    pathway_dict[pathway_name].append(gene_DataFrame.T)
                    pathway_pvalue[pathway_name].append(pathway_number)
                    pathway_number = 0

                pathway_cat = line[0]
                pathway_name = line[1]
                pathway_id = line[2]
                pathway_dict[pathway_name] = [pathway_cat, pathway_id]
                pathway_pvalue[pathway_name] = []
                continue

            # This line has specific attributes
            if len(line[0].split(';')) > 1:
                attributes = line[0].split(';')
                for attr in attributes:
                    pathway_dict[pathway_name].append(attr)
                    if attr.split('=')[0] == 'adjP':
                        adjp = float(attr.split('=')[1])
                        pathway_pvalue[pathway_name].append(adjp)
                continue

            # This line has gene info
            if line[1] == 'NA':
                gene = pd.Series(line, index=['gene', 'NA', 'symbol', 'descr',
                                              'entrez', 'ensemble'])
                if pathway_number == 0:
                    gene_DataFrame = pd.Series()

                gene_DataFrame = pd.concat([gene_DataFrame, gene], axis=1)
                pathway_number += 1

        # Add the final pathway outside the for loop
        pathway_dict[pathway_name].append(gene_DataFrame.T)
        pathway_pvalue[pathway_name].append(pathway_number)

    return pathway_dict, pathway_pvalue

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
    OUT_FH = 'data/gestalt/' + TRAIT + '_pathways.p'
    OUT_P_FH = 'data/gestalt/' + TRAIT + '_pvalues.tsv'

    ####################################
    # Load and Save Data
    ####################################
    gestalt_data, gestalt_p = read_gestalt(TRAIT_FH)

    # Save pickle file for pathway lookup
    pickle.dump(gestalt_data, open(OUT_FH, 'wb'))

    # Process p value data and output to tsv
    gestalt_p = pd.DataFrame.from_dict(gestalt_p, orient='index')
    gestalt_p.columns = ['adjP', 'num_genes']
    gestalt_p = gestalt_p.sort_values(by='adjP')
    gestalt_p.to_csv(OUT_P_FH, sep='\t')
