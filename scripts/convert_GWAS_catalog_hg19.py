"""
2016 Gregory Way
scripts/convert_GWAS_catalog_hg19.py

Description:
The GWAS catalog is downloaded in hg38 coordinates. Use pyliftover to convert
coordinates to hg19.

Usage:
Called by scripts/run_pipeline.sh

Output:
The GWAS Catalog but with hg19 coordinates
"""

import csv
from pyliftover import LiftOver

gwas_file = 'data/hg/gwas_catalog_v1.0.1-downloaded_2016-02-25.tsv'
gwas_output_file = 'data/hg/gwas_catalog_hg19.tsv'


def convert_coord(chromosome, position, liftover):
    """
    Convert SNPs across genome builds based on chromosome and position

    Arguments:
    chromosome - a string indicating the chromosome number
    position - the physical location of the SNP in the chromosome
    liftover - a pyliftover object mapping across genomic builds

    Output:
    lifted-over genomic elements
    """

    if position == '':
        return 'NA', 'NA'

    coord = int(position) + 1

    # Get chromosome coordinates
    chrom = 'chr' + str(chromosome)

    # Convert from hg38 to hg19
    converted = liftover.convert_coordinate(chrom, coord)

    # If chromosome name is incorrect
    if converted is None:
        return 'Check_Chromosome', 'Check_Chromosome'

    # If locus is deleted in conversion
    if len(converted) == 0:
        return 'Locus_Deleted', 'Locus_Deleted'

    # Convert to updated genome build
    new_chrom, new_coord, _, _ = converted[0]

    return new_chrom, new_coord

# Convert GWAS catalog to hg19
lo = LiftOver('data/hg/hg38ToHg19.over.chain.gz')

with open(gwas_file) as gwas_fh, open(gwas_output_file, 'w') as gwas_out_fh:
    reader = csv.DictReader(gwas_fh, delimiter='\t')
    writer = csv.DictWriter(gwas_out_fh, delimiter='\t',
                            fieldnames=reader.fieldnames)
    writer.writeheader()

    for row in reader:

        # X chromosomes are coded as '23' in the EBI-NHGRI catalog
        if row['CHR_ID'] == '23':
            row['CHR_ID'] = 'X'

        row['CHR_ID'], row['CHR_POS'] = convert_coord(row['CHR_ID'],
                                                      row['CHR_POS'], lo)
        writer.writerow(row)
