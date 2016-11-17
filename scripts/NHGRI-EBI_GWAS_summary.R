# 2016 Gregory Way
# TAD Pathways
# scripts/NHGRI-EBI_GWAS_summary.R

# Description: 
# Process the entire NHGRI-EBI GWAS catalog and parse out disease/trait specific
# information for subsets of GWAS findings (only those found in replication
# required journals)

# Usage:
# This script is run by "scripts/run_pipeline.sh"

# Output:
# Several chromosome specific tables and a large summary table of all findings
# and a histogram of number of SNPs per trait/disease

library("checkpoint")
checkpoint("2016-02-25")

# Load data
gwas_catalog <- readr::read_tsv("data/hg/gwas_catalog_hg19.tsv")

# Filter data to only replication required journals
repl_journals <- c("N Engl J Med", "Science", "Nature", "Nat Genet", "Lancet")
gwas_filter <- gwas_catalog[gwas_catalog$JOURNAL %in% repl_journals, ]

dir.create(file.path("data", "gwas_catalog"))

# Summarize the data and write to file
traits <- unique(gwas_filter$`DISEASE/TRAIT`)
num_snps <- c()
for (trait in traits) {
  gwas_subset <- gwas_filter[gwas_filter$`DISEASE/TRAIT` == trait, ]
  num_snps <- c(num_snps, nrow(gwas_subset))
  gwas_traitinfo <- gwas_subset[ , c("SNPS", "CHR_ID", "CHR_POS", 
                                              "DISEASE/TRAIT", 
                                              "REPORTED GENE(S)", "MAPPED_GENE",
                                              "PUBMEDID")]
  
  # Remove SNPs that are not mapped
  gwas_traitinfo <- gwas_traitinfo[complete.cases(gwas_traitinfo), ]
  
  # Generate file names
  if (nrow(gwas_traitinfo) > 0) {
    trait_name <- gsub(" ", "_", trait)
    trait_name <- gsub("[/:(),.-]", "_", trait_name)
    filename <- file.path("data", "gwas_catalog", paste0(trait_name,
                                                         "_hg19.tsv"))
    
    readr::write_tsv(gwas_traitinfo, filename)
  }
}

num_snp_order <- order(table(gwas_filter$`DISEASE/TRAIT`), decreasing = T)
num_snps_trait <- table(gwas_filter$`DISEASE/TRAIT`)[num_snp_order]

# Output histogram of number of SNPs per trait/disease
png("figures/num_snps_in_trait.png", height = 800, width = 1000)
hist(num_snps_trait, breaks = 100, xlab = "Number of SNPs", 
     main = "NHGRI-EBI GWAS\nNumber of Significant SNPs per Disease/Trait",
     cex = 5)
dev.off()