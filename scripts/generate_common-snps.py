"""
2016 Gregory Way - TAD Pathways
scripts/generate_common-snps.py
Modified from https://github.com/dhimmel/kg/blob/master/common-snps.ipynb

Description:
Extract SNPs from VCF files, subset into new SNP specific VCF files
stratified by chromosome and output a csv of all 1000 Genomes Phase 3
common SNPs or common mouse genome SNPs

Usage:
Is called by 'scripts/run_pipeline.sh' but can also be run through command line

           python generate_common-snps.py --genome <GENOME>

Where <GENOME> can be either 'hg' or 'mm' for human and mouse respectively

Output:
Will keep common SNPs in the 1000 Genomes Project Phase 3 data or Mouse Genomes
Project data and write out new subset vcf.gz files
"""

import os
import subprocess
import argparse
import csv

import vcf

# 1. Define Constants
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', help='Can be either hg or mm')
args = parser.parse_args()
GENOME = args.genome

# Maximum major allele frequency to allow
MAJOR_AF_MAX = 0.95
RAW_DIR = 'data/' + GENOME

if GENOME == 'hg':
    RAW_DIR = RAW_DIR + '/raw1000G/'
    VCF_FH = [
        'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes'
        '.vcf.gz'.format(x)
        for x in range(1, 23)]
    VCF_FH += [
        'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes' +
        '.vcf.gz', 'ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz']

elif GENOME == 'mm':
    VCF_FH = ['mgp.v2.snps.annot.reformat.vcf.gz']

else:
    raise ValueError('Please input either "hg" or "mm" for --genome')

# Build output directories
SUBSET_DIR = 'data/' + GENOME + '/subsetvcf/'
OUTPUT_FH = 'data/' + GENOME + '/' + GENOME + '_common-snps.tsv'

# Check if the subset directory exists
if not os.path.exists(SUBSET_DIR):
    os.makedirs(SUBSET_DIR)

# 2. Define Functions


def get_major_af(aafs, genome):
    """
    Returns the major allele frequency from the alternate allele frequencies.
    VCF format specifications: https://samtools.github.io/hts-specs/VCFv4.1.pdf

    Arguments:
    :param aafs: r.INFO in the vcf file Reader
    :param genome: the genome (either hg or mm)

    Output:
    The major allele frequency for the variant
    """
    if genome == 'hg':
        aafs = aafs['AF']
        allele_freqs = aafs + [1 - sum(aafs)]
    elif genome == 'mm':
        AC = aafs['AC'][0]
        AN = aafs['AN']
        AF = [float(AC) / float(AN)]
        allele_freqs = [1 - AF[0]] + AF
    return max(allele_freqs)

# 3. Create subset vcf files using bcftools
for vcf_file in VCF_FH:
    VCF_OUTPUT = os.path.join(SUBSET_DIR, vcf_file)

    # Only run bcftools if the file doesn't already exists
    if not os.path.exists(VCF_OUTPUT):
        RAW_FH = os.path.join(RAW_DIR, vcf_file)
        subprocess.call(['bcftools',
                         'view',
                         '--exclude-types', 'indels,mnps,other',
                         '--apply-filters', 'PASS',
                         '--max-af', str(MAJOR_AF_MAX) + ':major',
                         '--drop-genotypes',
                         '--output-type', 'z',
                         '--output-file', VCF_OUTPUT,
                         RAW_FH])

# 4. Write out common SNPs
with open(OUTPUT_FH, 'w') as writer_fh:
    snp_writer = csv.writer(writer_fh, delimiter='\t')
    snp_writer.writerow(['chromosome', 'position', 'major_allele_frequency',
                         'rsid'])
    for vcf_file in VCF_FH:
        path = os.path.join(SUBSET_DIR, vcf_file)

        for r in vcf.Reader(filename=path):
            # Quality control
            if r.FILTER:
                continue

            # Exclude non-SNPs
            if not r.is_snp:
                continue

            # Major allele frequency check
            major_af = get_major_af(r.INFO, genome=GENOME)
            if major_af > MAJOR_AF_MAX:
                continue

            # Exclude SNPs without rsid
            if r.ID is None:
                continue

            # Write SNP to tsv
            snp_writer.writerow([str(r.CHROM), str(r.POS),
                                 str(round(major_af, 6)), str(r.ID)])
