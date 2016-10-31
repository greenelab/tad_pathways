"""
2016 Gregory Way
TAD Pathways
scripts/build_TAD_genelists.py

Description:
Calls two scripts from bin/ to build genelists of all genes in TADs harboring
significant SNPs

Usage:
Is called by 'scripts/run_pipeline.sh'

Output:
Trait specific .tsv files of one column each indicating all the genes that fall
in signal TADs
"""

import os
from subprocess import call
import re


def check_file(path):
    """
    Make directory if path does not exist.
    """
    if not os.path.exists(path):
        os.makedirs(path)

gwas_results = os.path.join('data', 'gwas_catalog')
tad_gwas_loc = os.path.join('data', 'gwas_TAD_location')
tad_genes_loc = os.path.join('data', 'TAD_based_genes')

gwas_files = os.listdir(gwas_results)

identify_script = os.path.join('scripts', 'tad_util', 'Identify_TAD_signal.py')
grab_genes_script = os.path.join('scripts', 'tad_util', 'grab_TAD_genes.py')
tad_data = os.path.join('data', 'hg', 'hESC_domains_hg19.bed')

check_file(gwas_results)
check_file(tad_gwas_loc)
check_file(tad_genes_loc)

# Output GWAS info with corresponding SNP signal TAD
for gwas in gwas_files:
    file_replace = re.sub(r'\b.tsv\b', '_SNPs.tsv', gwas)
    output_file = os.path.join(tad_gwas_loc, file_replace)
    print('Building SNP list for: ' + gwas)
    call(['python', identify_script,
          '--tad_data', tad_data,
          '--gwas_data', os.path.join(gwas_results, gwas),
          '--output_file', output_file])

# Output subset TAD files that indicate all genes that fall in signal TADs
# for the given input genes
tad_gwas_files = os.listdir(tad_gwas_loc)
for tad in tad_gwas_files:
    out_tad_replace = re.sub(r'\b.tsv\b', '_TAD_genelists.tsv', tad)
    out_gwas_replace = re.sub(r'\b.tsv\b', '_GWAS_genelists.tsv', tad)
    output_file = os.path.join(tad_genes_loc, out_tad_replace)
    print('Grabbing TAD based genes for: ' + tad)
    call(['python', grab_genes_script,
          '--tad_gwas_file', os.path.join(tad_gwas_loc, tad),
          '--output_tad_file', output_file])
