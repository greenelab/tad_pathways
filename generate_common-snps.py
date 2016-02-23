"""
generate_common-snps.py
Author: Greg Way
Modified from https://github.com/dhimmel/kg/blob/master/common-snps.ipynb

Description:
Extract SNPs from VCF files, subset into new SNP specific VCF files
stratified by chromosome and output a csv of all 1000 Genomes Phase 3
common SNPs

Usage:
Command line 'python generate_common-snps.py' (must have 1000 genomes VCF
files in data/raw1000G)

Output:
Will keep common SNPs in the 1000 Genomes Project Phase 3 data and write out
new subset vcf.gz files in 'data/subsetvcf'
"""
import os
import subprocess
import vcf

##################
# Define Constants
##################
DOWNLOAD_DIR = 'data/raw1000G/'
SUBSET_DIR = 'data/subsetvcf/'

# Maximum major allele frequency to allow
MAJOR_AF_MAX = 0.95

# Compute vcf file names
VCF_FILES = [
    'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'
    '.vcf.gz'.format(x)
    for x in range(1, 23)]
VCF_FILES += [
    'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz',
    'ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz']

##################
# Define Functions
##################


def get_majar_af(aafs):
    """
    Returns the major allele frequency from
    the list of alternate allele frequencies.
    """
    allele_freqs = aafs + [1 - sum(aafs)]
    return max(allele_freqs)

##################
# Create subset vcf files using bcftools
##################
for vcf_file in VCF_FILES:
    subprocess.call(['bcftools', 'view', '--exclude-types',
                     'indels,mnps,other', '--apply-filters',
                     'PASS', '--max-af',
                     str(MAJOR_AF_MAX) + ':major', '--drop-genotypes',
                     '--output-type', 'z', '--output-file',
                     SUBSET_DIR + vcf_file, DOWNLOAD_DIR + VCF_FILES])

##################
# Write out common SNPs
##################
with open('data/common-snps.tsv', 'w') as writer_fh:
    writer_fh.write('\t'.join(['chromosome', 'position',
                    'major_allele_frequency', 'rsid', '\n']))
    for vcf_file in VCF_FILES:
        path = os.path.join(SUBSET_DIR, vcf_file)

        for r in vcf.Reader(filename=path):
            # Quality control
            if r.FILTER:
                print(r.FILTER)
                continue

            # Exclude non-SNPs
            if not r.is_snp:
                print(r.var_type)
                continue

            # Major allele frequency check
            major_af = get_majar_af(r.INFO['AF'])
            if major_af > MAJOR_AF_MAX:
                print(major_af)
                continue

            # Write SNP to tsv
            row = [str(r.CHROM), str(r.POS), str(round(major_af, 6)),
                   str(r.ID), '\n']
            writer_fh.write('\t'.join(row))
