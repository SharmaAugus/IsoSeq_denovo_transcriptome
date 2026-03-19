#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(DESeq2)
  library(edgeR)
  library(EnhancedVolcano)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

# Read in the count matrix and metadata
countData <- read.table("RSEM_counts_matrix.tsv", header = TRUE, row.names = 1, sep = "\t")
colData <- read.table("metadata.tsv", header = TRUE, row.names = 1, sep = "\t")

# Ensure column order matches between count matrix and metadata
countData <- countData[, rownames(colData)]

# Use top 500 features by variance for visualization
ntop <- 500

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ group)
dds <- DESeq(dds)

# List of all desired comparisons
comparisons <- list(
  c("group", "Sp_W", "Su_W"),
  c("group", "Sp_W", "Sp_C"),
  c("group", "Su_W", "Su_C"),
  c("group", "Sp_C", "Su_C"),
  c("group", "Sp_C", "Su_W")
)

# Function to perform DESeq2 for each comparison and save results
for (comp in comparisons) {
  res <- results(dds, contrast = comp)
  resOrdered <- res[order(res$padj), ]
  comp_name <- paste(comp[2], "vs", comp[3], sep = "_")
  outFile <- paste0("DESeq2_", comp_name, "_results.tsv")
  write.table(as.data.frame(resOrdered), file = outFile, sep = "\t", quote = FALSE)
  
  # Count DEGs
  degs <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)
  cat(comp_name, ": ", nrow(degs), " significant DEGs\n")
  
  # Volcano plot
  pdf(paste0("Volcano_", comp_name, ".pdf"))
  EnhancedVolcano(resOrdered,
                  lab = rownames(resOrdered),
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = comp_name,
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  pointSize = 2.0,
                  labSize = 3.0)
  dev.off()
}

# edgeR analysis
group <- factor(colData$group)
dge <- DGEList(counts = countData, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep,,keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~0 + group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

# Create all pairwise contrasts
contrast_list <- makeContrasts(
  Sp_W_vs_Su_W = groupSp_W - groupSu_W,
  Sp_W_vs_Sp_C = groupSp_W - groupSp_C,
  Su_W_vs_Su_C = groupSu_W - groupSu_C,
  Sp_C_vs_Su_C = groupSp_C - groupSu_C,
  Sp_C_vs_Su_W = groupSp_C - groupSu_W,
  levels = design
)

for (contrast_name in colnames(contrast_list)) {
  lrt <- glmLRT(fit, contrast = contrast_list[, contrast_name])
  res <- topTags(lrt, n = Inf)$table
  outFile <- paste0("edgeR_", contrast_name, "_results.tsv")
  write.table(res, file = outFile, sep = "\t", quote = FALSE)
  
  # Count DEGs
  degs <- subset(res, FDR < 0.05 & abs(logFC) > 1)
  cat(contrast_name, ": ", nrow(degs), " significant DEGs\n")
  
  # Volcano plot
  pdf(paste0("Volcano_", contrast_name, ".pdf"))
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'logFC',
                  y = 'FDR',
                  title = contrast_name,
                  pCutoff = 0.05,
                  FCcutoff = 1.0,
                  pointSize = 2.0,
                  labSize = 3.0)
  dev.off()
}

cat("All DESeq2 and edgeR analyses completed successfully.\n")
