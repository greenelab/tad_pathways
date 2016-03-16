# (C) 2016 Gregory Way
# NHGRI-EBI_SNAP_summary.R

# Description: 
# Take as input the chromosome specific SNAP files, which includes linkage
# disequilibrium estimates for each disease associated SNP, and output only
# independent SNPs (in LD < 0.2)

# Usage:
# This script is run by 'ANALYSIS.sh' but can also be run directly by calling
# from the parent directory 'R --no-save < NHGRI-EBI_SNAP_summary.R'

# Output:
# A table of independent SNPs and corresponding genomic locations
# A scatterplot of number of independent SNPs per chromosome

# Load all SNAP results -- See 'NHGRI_EBI_GWAS_summary.R'
snap_file_dir <- c('data/SNAP/')
snap_files <- list.files(snap_file_dir)[grepl('SNAP_Results', 
                                              list.files(snap_file_dir))]

# Load function
find_independent_snps <- function(snap_file, chromosome) {
  
  # Arguments:
  #    snap_file - results of broad institute SNAP analysis for particular chrom
  #    chromosome - chromosome corresponding to the snap_file
  # Output:
  #    dataframe of independent SNP signals across the chrom and the chrm number
  
  unique_test_snps <- unique(c(snap_file$SNP, snap_file$Proxy))
  independent_snps <- c()
  redundant_snps <- c()
  for (snp in unique_test_snps) {
    test_sub <- snap_file[snap_file$SNP == snp | snap_file$Proxy == snp, ]
    for (relation in 1:nrow(test_sub)) {
      test_snp <- test_sub$SNP[relation]
      test_prox <- test_sub$Proxy[relation]
      if (test_snp == test_prox) {
        next
      }
      test_LD <- test_sub$RSquared[relation]
      if (test_LD < 0.2) {
        if (length(grep(test_snp, redundant_snps)) == 0) {
          independent_snps <- unique(c(independent_snps, test_snp))
        } 
        if (length(grep(test_prox, redundant_snps)) == 0) {
          independent_snps <- unique(c(independent_snps, test_prox))
        }
        
      } else {
        if (length(grep(test_snp, redundant_snps)) >= 1) {
          if (length(grep(test_prox, redundant_snps)) == 0) {
            independent_snps <- unique(c(independent_snps, test_snp))
          } 
          if (length(grep(test_snp, independent_snps)) >= 1) {
            independent_snps <- setdiff(independent_snps, test_snp)
          }
        } else if (length(grep(test_prox, redundant_snps)) >= 1) {
          if (length(grep(test_snp, redundant_snps)) == 0) {
            independent_snps <- unique(c(independent_snps, test_prox))
          } 
          if (length(grep(test_prox, independent_snps)) >= 1) {
            independent_snps <- setdiff(independent_snps, test_prox)
          }
        } else if (length(grep(test_snp, independent_snps)) >= 1) {
          redundant_snps <- unique(c(redundant_snps, test_prox))
          if (length(grep(test_prox, independent_snps)) >= 1) {
            independent_snps <- setdiff(independent_snps, test_prox)
          }
        } else if (length(grep(test_prox, independent_snps)) >= 1) {
          redundant_snps <- unique(c(redundant_snps, test_snp))
          if (length(grep(test_snp, independent_snps)) >= 1) {
            independent_snps <- setdiff(independent_snps, test_snp)
          }
        }
      }
    }
  }
  chrom <- rep(chromosome, length(independent_snps))
  return(cbind(independent_snps, chrom))
}

# Analysis and write output to file
all_independent_snps <- c()
for (snap_fh in snap_files) {
  snap <- readr::read_tsv(paste(snap_file_dir, snap_fh, sep = ""))
  chrom <- unlist(strsplit(unlist(strsplit(snap_fh, '_'))[3], '[.]'))[1]
  chrom <- substring(chrom, 4, nchar(chrom))
  all_independent_snps <- rbind(all_independent_snps, 
                                find_independent_snps(snap, chrom))
}

write.table(all_independent_snps, 
            'output/all_independent_snps_across_traits.csv', sep = ',', 
            col.names = T, row.names = F)

summary_snps <- as.data.frame(table(all_independent_snps[,2]))
png("figures/independent_SNP_by_chromosome.png", height = 500, width = 600)
plot(as.numeric(paste(summary_snps[,1])), summary_snps[,2], pch = 16,
     xlab = 'Chromosome', 
     ylab = 'Number of independent SNPs', 
     main = paste0("Summary of Independent SNPs across each chromosome \n",
            "(299 Traits/Diseases in Replication-Required Journals)"))
dev.off()
