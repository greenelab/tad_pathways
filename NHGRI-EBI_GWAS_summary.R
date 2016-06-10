# (C) 2016 Gregory Way
# NHGRI-EBI_GWAS_summary.R

# Description: 
# Process the entire NHGRI-EBI GWAS catalog and parse out disease/trait specific
# information for subsets of GWAS findings (only those found in replication
# required journals)

# Usage:
# This script is run by 'ANALYSIS.sh' but can also be run directly by calling
# from the parent directory 'R --no-save < NHGRI-EBI_GWAS_summary.R'

# Output:
# Several chromosome specific tables and a large summary table of all findings
# and a histogram of number of SNPs per trait/disease

# Load data
gwas_catalog <- readr::read_tsv("data/gwas_catalog_hg19.tsv")

# Filter data to only replication required journals
repl_journals <- c("N Engl J Med", "Science", "Nature", "Nat Genet", "Lancet")
gwas_filter <- gwas_catalog[gwas_catalog$JOURNAL %in% repl_journals, ]

# Summarize the data and write to file
gwas_traitinfo = list()
traits <- unique(gwas_filter$`DISEASE/TRAIT`)
num_snps <- c()
for (trait in traits) {
  gwas_subset <- gwas_filter[gwas_filter$`DISEASE/TRAIT` == trait, ]
  num_snps <- c(num_snps, nrow(gwas_subset))
  cat(trait, ' ', nrow(gwas_subset), '\n')
  gwas_traitinfo[[trait]] <- gwas_subset[ , c('SNPS', 'CHR_ID', 'CHR_POS', 
                                              'DISEASE/TRAIT', 
                                              'REPORTED GENE(S)', 'MAPPED_GENE',
                                              'PUBMEDID')]
  trait_name <- gsub(' ', '_', trait)
  trait_name <- gsub('[/:(),.-]', '_', trait_name)
  filename = paste0('data/gwas_catalog/', trait_name, '.tsv')
  print(filename)
  readr::write_tsv(gwas_traitinfo[[trait]], filename)
}

num_snp_order <- order(table(gwas_filter$`DISEASE/TRAIT`), decreasing = T)
num_snps_trait <- table(gwas_filter$`DISEASE/TRAIT`)[num_snp_order]

# Output histogram of number of SNPs per trait/disease
png('figures/num_snps_in_trait.png', height = 800, width = 1000)
hist(num_snps_trait, breaks = 100, xlab = 'Number of SNPs', 
     main = 'NHGRI-EBI GWAS\nNumber of Significant SNPs per Disease/Trait',
     cex = 5)
dev.off()

# Save the total SNP information to a file separated by chromosome
gwas_use_data <- gwas_filter[, c('SNPS', 'CHR_ID', 'CHR_POS', 'DISEASE/TRAIT',
                                 'REPORTED GENE(S)', 'MAPPED_GENE', 'PUBMEDID')]

gwas_use_data <- gwas_use_data[gwas_use_data$SNPS != '', ]
gwas_use_data_chrom <- list()
for (chrom in 1:length(unique(gwas_use_data$CHR_ID))){
  gwas_use_data_chrom[[chrom]] <- gwas_use_data[gwas_use_data$CHR_ID == chrom, ]
  write.table(gwas_use_data_chrom[[chrom]], 
              paste0('data/SNAP/all_unique_snps_chr', chrom, '.tsv'), 
              col.names = F, row.names = F, quote = F, sep = '\t')
}

# Save to file with chromosome information
gwas_snps_unique_ch <- cbind(gwas_use_data$SNPS, gwas_use_data$CHR_ID)
write.table(gwas_snps_unique_ch, 'data/SNAP/all_unique_trait_snps_chrID.csv', 
            sep = ',', col.names = F, row.names = F, quote = F)

# Chromosome specific files are input into:
# https://www.broadinstitute.org/mpg/snap/ldsearchpw.php
# Parameters: 
# SNP data set: 1000 Genomes Pilot 1
# r2 threshold: No limit
# Population panel: CEU
# Distance limit: 500
# Download to: File
# Include each query snp as a proxy for itself: yes
# Suppress warning messages in output: no
# Filter by Array: No
# Output Columns: D', Genome Coordinates
