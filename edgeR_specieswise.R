
#!/usr/bin/env Rscript

library(edgeR)
library(tximport)
library(readr)
library(dplyr)
library(pheatmap)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
metadata_file <- args[1]

# Read metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Read all quant.sf files with tximport
files <- setNames(metadata$path, metadata$sample_id)
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Store results in list
results_list <- list()

# Loop per species
for (sp in unique(metadata$species)) {
  cat("Processing species:", sp, "\n")
  sp_meta <- metadata[metadata$species == sp, ]
  sp_samples <- sp_meta$sample_id
  sp_counts <- txi$counts[, sp_samples]

  # Prepare DGEList
  group <- factor(sp_meta$condition)
  dge <- DGEList(counts = sp_counts)
  dge <- calcNormFactors(dge, method = "TMMwsp")

  # Design matrix: individual + condition
  sp_meta$individual <- factor(sp_meta$individual)
  design <- model.matrix(~ individual + condition, data = sp_meta)

  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)

  # Test for conditionSummer vs conditionSpring
  qlf <- glmQLFTest(fit, coef = "conditionSummer")

  # Save results
  result_file <- paste0("edgeR_results_", gsub(" ", "_", sp), ".csv")
  topTags(qlf, n = Inf) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    write.csv(file = result_file, row.names = FALSE)

  cat("Saved:", result_file, "\n")
}
