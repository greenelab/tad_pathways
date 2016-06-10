# (C) 2016 Gregory Way
# NHGRI-EBI_GWAS_summary.R

# Description: 
# Process the entire NHGRI-EBI GWAS catalog and parse out disease/trait specific
# information for subsets of GWAS findings (only those found in replication
# required journals)

# Usage:
# This script is run by 'ANALYSIS.sh'

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
  gwas_traitinfo[[trait]] <- gwas_subset[ , c('SNPS', 'CHR_ID', 'CHR_POS', 
                                              'DISEASE/TRAIT', 
                                              'REPORTED GENE(S)', 'MAPPED_GENE',
                                              'PUBMEDID')]
  trait_name <- gsub(' ', '_', trait)
  trait_name <- gsub('[/:(),.-]', '_', trait_name)
  filename <- paste0('data/gwas_catalog/', trait_name, '_hg19.tsv')
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