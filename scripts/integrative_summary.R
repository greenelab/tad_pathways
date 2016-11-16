# 2016 Gregory Way
# scripts/Integrative_summary.R

# Description: 
# Take as input the genes identified by the TAD pathway analysis, the eQTL
# analysis, and the nearest gene GWAS to determine evidence overlaps

# Usage:
# The script is run by scripts/run_pipeline with evidence file and trait as
# command arguments

# Output:
# Venn diagrams of disease specific evidence overlaps

library("checkpoint")
checkpoint("2016-02-25")

# Load in command arguments
args <- commandArgs(trailingOnly = T)

# Parse command argumentsx
evidence <- args[1]
trait <- args[2]
venn_output_file <- file.path('figures', paste0('venn_', trait, '.tiff'))

# Read in Data
evidence_genes <- readr::read_csv(evidence)

# Prepare for VennDiagram input
tad_genes <- c()
eqtl_genes <- c()
gwas_genes <- c()
for (gene in 1:nrow(evidence_genes)) {
  assignment <- evidence_genes$evidence[gene]
  
  if (assignment == 'tad') {
    tad_genes <- c(tad_genes, gene)
  } else if (assignment == 'gwas') {
    gwas_genes <- c(gwas_genes, gene)
  } else if (assignment == 'eqtl') {
    eqtl_genes <- c(eqtl_genes, gene)
  } else if (assignment == 'gwas_tad') {
    gwas_genes <- c(gwas_genes, gene)
    tad_genes <- c(tad_genes, gene)
  } else if (assignment == 'gwas_eqtl') {
    gwas_genes <- c(gwas_genes, gene)
    eqtl_genes <- c(eqtl_genes, gene)
  } else if (assignment == 'eqtl_tad') {
    tad_genes <- c(tad_genes, gene)
    eqtl_genes <- c(eqtl_genes, gene)
  } else if (assignment == 'all') {
    tad_genes <- c(tad_genes, gene)
    eqtl_genes <- c(eqtl_genes, gene)
    gwas_genes <- c(gwas_genes, gene)
  }
  
}

venn_list <- list('eQTL' = eqtl_genes, 'GWAS' = gwas_genes,
                  'TAD Pathway' = tad_genes)

# Output Venn Diagram
VennDiagram::venn.diagram(x = venn_list,
                          filename = VENN_FH,
                          fill = c("red", "blue", "yellow"),
                          height = 1500,
                          width = 1500,
                          euler.d = F,
                          scaled = F)
