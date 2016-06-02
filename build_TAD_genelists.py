"""
(C) 2016 Gregory Way
build_TAD_genelists.py

Description:
Calls two scripts from bin/ to build genelists of all genes in TADs harboring
significant SNPs

Usage:
Is called by 'ANALYSIS.sh' but can also be run through:
Command line 'python build_TAD_genelists.py'

Output:
Trait specific .txt files of one column each indicating all the genes that fall
in signal TADs
"""

import os
from subprocess import call
import re

GWAS_RESULTS = 'data/gwas_catalog/'
TAD_GWAS_LOC = 'data/gwas_TAD_location/'
TAD_GENES_LOC = 'data/TAD_based_genes/'
GWAS_LOC = 'data/gwas_based_genes/'

GWAS_FILES = os.listdir(GWAS_RESULTS)

for gwas in GWAS_FILES:
    output_fh = TAD_GWAS_LOC + re.sub(r'\b.tsv\b', '_SNPs.tsv', gwas)
    print output_fh
    call(['python', 'bin/Identify_TAD_signal.py', '-t',
          'data/hESC_domains_hg19.bed', '-g', GWAS_RESULTS + gwas,
          '-o', output_fh])

TAD_GWAS_FILES = os.listdir(TAD_GWAS_LOC)

for tad in TAD_GWAS_FILES:
    output_fh = TAD_GENES_LOC + re.sub(r'\b.tsv\b', '_TAD_genelists.tsv', tad)
    output_gwas_fh = GWAS_LOC + re.sub(r'\b.tsv\b', '_GWAS_genelists.tsv', tad)
    print output_fh
    call(['python', 'bin/grab_TAD_genes.py', '-f', TAD_GWAS_LOC + tad,
         '-o', output_fh, '-g', output_gwas_fh])
