expr <- read.table("combined_expression_TPM_matrix.txt",
                   header=TRUE,
                   row.names=1,
                   sep="\t")

dim(expr)

### Log transformation
expr_log <- log2(expr + 1)

### Filter Low Expression Genes
keep <- rowSums(expr > 1) >= 5
expr_filtered <- expr_log[keep, ]

### Confirm Matrix Orientation
dim(expr_filtered)
76317 × 18

###Variance Filtering
library(matrixStats)
vars <- rowVars(expr_filtered)

> class(expr_filtered)
[1] "data.frame"
> expr_filtered <- as.matrix(expr_filtered)
> class(expr_filtered)
[1] "matrix" "array" 

library(matrixStats)

### Calculate Variance
vars <- rowVars(expr_filtered)
summary(vars)


### keep top 25,000 transcripts (select most variable transcripts)
top <- order(vars, decreasing = TRUE)[1:25000]

expr_final <- expr_filtered[top, ]
### Check size:
dim(expr_final)

### Prepare Trait Metadata
sampleNames <- rownames(datExpr)
sampleNames

season <- ifelse(grepl("V$", sampleNames), 1, 0)
data.frame(sampleNames, season)

### Extract Individual IDs. Remove the last letter (P or V).
individual <- sub("[VP]$", "", sampleNames)
data.frame(sampleNames, individual)

### Build Trait Metadata Table
traitData <- data.frame(
  Season = season,
  Individual = individual
  )
rownames(traitData) <- sampleNames
traitData

### Convert Individual to Factor
traitData$Individual <- as.factor(traitData$Individual)
# Now R treats individuals as categorical variables.

### Sample Clustering (First WGCNA Diagnostic). Just if any samples are outliers before network construction.
sampleTree <- hclust(dist(datExpr), method = "average")

plot(sampleTree,
  main="Sample clustering to detect outliers",
  xlab="Samples",
  sub="")
traitColors <- data.frame(
  Season = numbers2colors(traitData$Season,
  colors = c("skyblue","orange"))
  )

plotDendroAndColors(sampleTree,
  traitColors,
  groupLabels="Season")

# Check 
dim(datExpr)

### Choose Soft Threshold Power. This determines how the gene correlation network is constructed.
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr,
       powerVector = powers, verbose = 5)
# For plots and choosing the best power

dev.off()

par(mfrow=c(1,2))
plot(sft$fitIndices[,1],
        -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit (R^2)",
         type="n")

text(sft$fitIndices[,1],
   -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,
    col="red")
abline(h=0.8,col="blue")

plot(sft$fitIndices[,1],
    sft$fitIndices[,5],
    xlab="Soft Threshold",
    ylab="Mean Connectivity",
    type="n")

text(sft$fitIndices[,1],
      sft$fitIndices[,5],
      labels=powers,
       col="red")


### blockwiseModules calculation 
net <- blockwiseModules(
  datExpr,
  power = 6,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "Flower_TOM",
  verbose = 3)

# Here selected power is 6
# This computes the adjacency matrix.
# But WGCNA usually goes directly to TOM (Topological Overlap Matrix) because it is more biologically meaningful.
#Parameter                Meaning
# power              soft threshold (6)
# minModuleSize      minimum genes in module
# mergeCutHeight     merge similar modules
#TOMType               network type

#Check detected modules:
table(net$colors)

### Convert Module Labels to Colors WGCNA modules are usually shown as colors.
moduleColors <- labels2colors(net$colors)

### Plot Gene Dendrogram with Modules
plotDendroAndColors(
net$dendrograms[[1]],
moduleColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE,
hang = 0.03)

table(moduleColors)
moduleColors

### Extract Module Eigengenes
MEs <- net$MEs
dim(MEs)

### Correlate Modules With Traits. 
### 1) Our trait table: Season and Individual

moduleTraitCor <- cor(MEs, traitData$Season, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples=18)

### Create Module–Trait Heatmap
par(mar = c(6,8,3,3))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Season",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  zlim = c(-1,1)
)

### Use Proper Module Colors
#You already have:
moduleColors <- labels2colors(net$colors)

# Now align them with MEs:
MEs_col <- moduleEigengenes(datExpr, moduleColors)$eigengenes

# Correlation + text
moduleTraitCor <- cor(MEs_col, traitData$Season, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 18)
textMatrix <- paste(signif(moduleTraitCor,2), "\n(",
           signif(moduleTraitPvalue,1), ")", sep="")

# Remove Text Labels 
yLabels = moduleColorNames
yLabels = FALSE

#Add Module Color Bar 
moduleColorNames <- substring(names(MEs_col), 3)

#labeledHeatmap does NOT directly support side color bars but pheatmap
library(pheatmap)
# Row annotation (module colors)
annotation_row <- data.frame(Module = moduleColorNames)
rownames(annotation_row) <- rownames(moduleTraitCor)

#Force textMatrix to be a matrix with same dimensions:
textMatrix <- matrix(textMatrix, 
            nrow = nrow(moduleTraitCor),
            ncol = ncol(moduleTraitCor))

# Fix Annotation Colors 
ann_colors <- list(Module = setNames(moduleColorNames, moduleColorNames))

library(pheatmap)
pheatmap(
  moduleTraitCor,
  color = colorRampPalette(c("blue","white","red"))(50),
  display_numbers = textMatrix,
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize = 7,
  fontsize_number = 5,
  border_color = NA
  )

### Show Only Significant Values
textMatrixSig <- textMatrix
textMatrixSig[moduleTraitPvalue > 0.05] <- ""

> library(pheatmap)
> # Ensure textMatrix is proper matrix
  > textMatrix <- matrix(textMatrix,
                         +                      nrow = nrow(moduleTraitCor),
                         +                      ncol = ncol(moduleTraitCor))
# Row annotation
annotation_row <- data.frame(Module = moduleColorNames)
rownames(annotation_row) <- rownames(moduleTraitCor)
# Correct color mapping
 ann_colors <- list(
 Module = setNames(moduleColorNames, moduleColorNames)
  )
# Plot
pheatmap(
  moduleTraitCor,
  color = colorRampPalette(c("blue","white","red"))(50),
  display_numbers = textMatrixSig,
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize = 7,
  fontsize_number = 6,
  border_color = NA
 )

### 2) Model Seasonal Effect While Controlling Individual. (Simple & Powerful): Residual Approach
#Remove individual effect first.
# Convert Individual to factor
traitData$Individual <- as.factor(traitData$Individual)
# Remove individual effect from each module eigengene
MEs_resid <- MEs_col
for(i in 1:ncol(MEs_col)){
   fit <- lm(MEs_col[,i] ~ traitData$Individual)
   MEs_resid[,i] <- residuals(fit)
   }

### Now Correlate With Season
moduleTraitCor <- cor(MEs_resid, traitData$Season, use="p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 18)
# Recreate Heatmap (Same as Before)
textMatrix <- paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,1), ")", sep="")
pheatmap(
  moduleTraitCor,
  color = colorRampPalette(c("blue","white","red"))(50),
  display_numbers = textMatrixSig,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize = 7
  )


######3) Direct Paired Difference
# Instead of regression: Split samples
spring <- datExpr[traitData$Season == 0, ]
summer <- datExpr[traitData$Season == 1, ]

# Ensure same order
spring <- spring[order(traitData$Individual[traitData$Season == 0]), ]
summer <- summer[order(traitData$Individual[traitData$Season == 1]), ]
 
# Compute difference
diffExpr <- summer - spring
MEs_diff <- moduleEigengenes(diffExpr, moduleColors)$eigengenes

modulePvalues <- apply(MEs_diff, 2, function(x) {
 t.test(x)$p.value
 })

moduleMeans <- colMeans(MEs_diff)

moduleTraitCor <- matrix(moduleMeans, ncol = 1)

rownames(moduleTraitCor) <- colnames(MEs_diff)
colnames(moduleTraitCor) <- "SeasonEffect"

textMatrix <- paste(signif(moduleMeans,2), "\n(",
                signif(modulePvalues,1), ")", sep="")

# Convert to matrix
textMatrix <- matrix(textMatrix,
             nrow = nrow(moduleTraitCor),
             ncol = ncol(moduleTraitCor))
# Only significant 
textMatrixSig <- textMatrix
textMatrixSig[moduleTraitPvalue > 0.05] <- ""

library(pheatmap)

> # Row annotation
annotation_row <- data.frame(Module = moduleColorNames)
rownames(annotation_row) <- rownames(moduleTraitCor)

# Correct color mapping
  ann_colors <- list(
  Module = setNames(moduleColorNames, moduleColorNames)
  )
pheatmap(
  moduleTraitCor,
  color = colorRampPalette(c("blue","white","red"))(50),
  display_numbers = textMatrixSig,
  number_color = "black",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize = 7,
  fontsize_number = 4,
  border_color = NA
  )

# strongest seasonal modules
orderModules <- order(abs(moduleMeans), decreasing = TRUE)

topModules <- data.frame(
  Module = names(moduleMeans)[orderModules],
  Effect = moduleMeans[orderModules],
  Pvalue = modulePvalues[orderModules]
)
head(topModules, 10)

######### Selected method 2 For further analysis
#  recompute these using MEs_resid:
moduleTraitCor <- cor(MEs_resid, traitData$Season)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 18)

## Identify Top Modules
topModules <- data.frame(
  Module = rownames(moduleTraitCor)[orderModules],
  Correlation = moduleTraitCor[orderModules,1],
  Pvalue = moduleTraitPvalue[orderModules,1]
  )
 
head(topModules, 10)

# Top Module is Meivory
module <- "ivory"
# Get all genes in this module
genes_ivory <- colnames(datExpr)[moduleColors == module]
length(genes_ivory)

head(genes_ivory)
length(genes_ivory)

## Calculate Module Membership (kME). This identifies hub genes (most important genes in module)
#Calculate correlation with module eigengene
kME <- cor(datExpr, MEs_resid, use = "p")
# Extract kME for Meivory
kME_ivory <- kME[, "Meivory"]


#### Identify TOP Hub Genes. (Hub Genes most strongly connected within module → potential regulators)
hubGenes_ivory <- names(sort(abs(kME_ivory), decreasing = TRUE))[1:30] 
hubGenes_ivory

##Combine Everything into a Table
ivory_table <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  kME = kME_ivory
)

ivory_table <- ivory_table[moduleColors == "ivory", ]
ivory_table <- ivory_table[order(-abs(ivory_table$kME)), ]
 
head(ivory_table, 20)

## Add Gene Significance (GS). This connects genes to season
GS <- cor(datExpr, traitData$Season, use = "p") 
ivory_table$GS <- GS[moduleColors == "ivory"]

ivory_table <- ivory_table[order(-abs(ivory_table$GS)), ]
topHubGenes <- ivory_table[
  abs(ivory_table$kME) > 0.8 & abs(ivory_table$GS) > 0.6,
]

topHubGenes

plot(abs(ivory_table$kME),
     abs(ivory_table$GS),
     xlab = "Module Membership (kME)",
     ylab = "Gene Significance (GS)",
     main = "MEivory: Hub Gene Selection")


###### Annotation
# We have ids like MOREFL_tx_... in Annotation folder while here it is like transcript/…, so first match those ids.
ivory_table$Gene <- gsub("transcript/", "MORE24FL_tx_", ivory_table$Gene)

#Load Annotation File in R
annot <- read.delim("/scratch/MORE24/Annotation/Prueba2-ok/MORE24-FLOR-Mar-Trinotate.tsv",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)
head(annot)
colnames(annot)

### merge top Module (ivory_table) with annotation file
annot_ivory <- merge(ivory_table,
                     annot,
                     by.x = "Gene",
                     by.y = "X.gene_id",
                     all.x = TRUE)
sum(!is.na(annot_ivory$sprot_Top_BLASTX_hit))
nrow(annot_ivory)
### 204 but in ivory table it is 139 means IDs are NOT matching correctly OR duplicates exist in annotation
sum(duplicated(annot$X.gene_id))
### YES → duplicates exist → causing row expansion
### Fix Annotation Table, remove duplicates
annot_unique <- annot[!duplicated(annot$X.gene_id), ]

### Now merge again
annot_ivory <- merge(ivory_table,
                     annot_unique,
                     by.x = "Gene",
                     by.y = "X.gene_id",
                     all.x = TRUE)
nrow(annot_ivory)
### Check again annotation coverage
sum(!is.na(annot_ivory$sprot_Top_BLASTX_hit))

head(annot_ivory[, c("Gene","kME","GS","sprot_Top_BLASTX_hit","Pfam")])

### TOP HUB GENES
topHubAnnotated <- annot_ivory[order(-abs(annot_ivory$kME)), ][1:20, ]

topHubAnnotated[, c("Gene","kME","GS","sprot_Top_BLASTX_hit","Pfam")]

### FINAL HUB GENES
finalHubGenes <- annot_ivory[
  abs(annot_ivory$kME) > 0.85 & abs(annot_ivory$GS) > 0.7,
]

finalHubGenes[, c("Gene","sprot_Top_BLASTX_hit","Pfam")]
head(topHubAnnotated[, c("sprot_Top_BLASTX_hit","Pfam")], 10)

write.csv(annot_ivory,
          "MEivory_full_annotated.csv",
          row.names = FALSE)
write.csv(topHubAnnotated,
          "MEivory_topHubGenes.csv",
          row.names = FALSE)
write.csv(finalHubGenes,
          "MEivory_finalHubGenes.csv",
          row.names = FALSE)

finalHubGenes <- annot_ivory[
  abs(annot_ivory$kME) > 0.85 & abs(annot_ivory$GS) > 0.6,
]

### Plot Network in R
library(WGCNA)

### Export Network to Cytoscape
library(WGCNA)
module <- "ivory"
inModule <- moduleColors == module
### Compute adjacency
adjacency_matrix <- adjacency(datExpr[, inModule], power = 6)
# TOM similarity
TOM <- TOMsimilarity(adjacency_matrix)
# Dissimilarity
dissTOM <- 1 - TOM
# Clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
# Plot
plot(geneTree, main = "Gene clustering (MEivory)")
TOM <- TOMsimilarityFromExpr(datExpr[, inModule], power = 6)
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering (MEivory)")


#datExpr coloumn has old name ids so annotation gene name ids will match to them.
#From annotation file 
mapping <- annot_ivory[, c("Gene", "sprot_Top_BLASTX_hit")]

# Extract readable gene names
mapping$GeneName <- sub("\\^.*", "", mapping$sprot_Top_BLASTX_hit)
head(mapping$GeneName)

# Match with datExpr IDs
mapping$Gene <- gsub("MORE24FL_tx_", "transcript/", mapping$Gene)
# Replace column names
match_idx <- match(colnames(datExpr), mapping$Gene)
new_names <- mapping$GeneName[match_idx]

Prepare annotation fallback to original if NA
new_names[is.na(new_names)] <- colnames(datExpr)[is.na(new_names)]
colnames(datExpr) <- new_names

## Hub Gene Expression Plot
# Prepare annotation
annotation_col <- data.frame(
  Season = traitData$Season,
  Individual = traitData$Individual
)
rownames(annotation_col) <- rownames(datExpr)
# plot
library(pheatmap)

hub <- finalHubGenes$Gene
hub <- gsub("MORE24FL_tx_", "transcript/", hub)

# convert again to NEW names
hub_names <- mapping$GeneName[match(hub, mapping$Gene)]

mat <- datExpr[, hub_names]

pheatmap(
  mat,
  scale = "row",
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue","white","red"))(50),
  fontsize_row = 8,
  main = "Hub gene expression (MEivory)"
)
hub_names
sum(hub %in% mapping$Gene)

# Missing annotation (".") and HS704_ARATH" appears twice.All hub genes matched,But names still messy
# Remove bad gene names
valid <- hub_names != "." & !is.na(hub_names)

hub_names <- hub_names[valid]
hub <- hub[valid]
# Remove duplicates
hub_names_unique <- unique(hub_names)

#Extract matrix
mat <- datExpr[, hub_names_unique]
dim(mat)

library(pheatmap)

# Order samples properly
annotation_col <- data.frame(
  Season = traitData$Season,
  Individual = traitData$Individual
)
rownames(annotation_col) <- rownames(datExpr)

ord <- order(annotation_col$Season)
mat <- mat[ord, ]
annotation_col <- annotation_col[ord, ]

pheatmap(
  mat,
  scale = "row",
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue","white","red"))(50),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "MEivory Hub Genes (Heat Stress Module)"
)

# there is error and followinh things checked and changed
rowSds <- apply(mat, 1, sd)
summary(rowSds)
dim(mat_filtered)
sum(is.na(mat))

all(rownames(mat) == rownames(annotation_col))

str(mat)

mat <- as.matrix(mat)
mode(mat) <- "numeric"
range(mat)

apply(mat, 2, sd)
mat_filtered <- mat[, apply(mat, 2, sd) > 1e-6]
# this worked mat_t <- t(mat). Rows = genes → scaled (correct) and Columns = samples → annotated

pheatmap(
  mat_t,
  scale = "row",
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue","white","red"))(50),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "MEivory Hub Genes (Heat Stress Module)"
)

# Cytoscape (export using gene names)
nodeNames <- colnames(datExpr)[inModule]

exportNetworkToCytoscape(
  TOM,
  edgeFile = "Cytoscape_edges_ivory.txt",
  nodeFile = "Cytoscape_nodes_ivory.txt",
  weighted = TRUE,
  threshold = 0.1,
  nodeNames = nodeNames,
  nodeAttr = moduleColors[inModule]
)
# Error came, I used nodeNames <- colnames(datExpr)[inModule], this gives ALL genes in module (139) 
# TOM <- TOMsimilarityFromExpr(datExpr[, topGenes], power = 6), was calculated for top genes (25)
# For clean hub-centered network, not all 139 genes

topGenes <- ivory_table$Gene[abs(ivory_table$kME) > 0.8]

# Ensure they exist in datExpr
topGenes <- topGenes[topGenes %in% colnames(datExpr)]
# Subset expression
expr_sub <- datExpr[, topGenes]
# Compute TOM
library(WGCNA)

# TOM <- TOMsimilarityFromExpr(expr_sub, power = 6) not worked with mismatched ids
TOM <- TOMsimilarityFromExpr(datExpr[, inModule], power = 6)

exportNetworkToCytoscape(
  TOM,
  edgeFile = "Cytoscape_edges_ivory.txt",
  nodeFile = "Cytoscape_nodes_ivory.txt",
  weighted = TRUE,
  threshold = 0.3,
  nodeNames = colnames(datExpr)[inModule],
  nodeAttr = moduleColors[inModule]
)
# FIX duplicate edges

edges <- read.delim("Cytoscape_edges_ivory.txt")

edges_unique <- edges[!duplicated(
  t(apply(edges[,1:2], 1, sort))
), ]

write.table(edges_unique,
            "Cytoscape_edges_ivory_clean.txt",
            sep="\t", row.names=FALSE, quote=FALSE)


# Get genes present in edges
edge_unique <- unique(c(edges$fromNode, edges$toNode))

# Filter node table
nodes <- read.delim("Cytoscape_nodes_ivory.txt")
nodes_clean <- nodes[nodes$nodeName %in% edge_genes, ]

# Save new node file
write.table(nodes_clean,
            "Cytoscape_nodes_ivory_clean.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)


### GO ENRICHMENT
# Get genes present in edges
edge_genes <- unique(c(edges$fromNode, edges$toNode))

# Filter node table
nodes_clean <- nodes[nodes$nodeName %in% edge_genes, ]

# Save new node file
write.table(nodes_clean,
            "Cytoscape_nodes_ivory_clean.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)ICHMENT
# module genes
ivory_genes <- ivory_table$Gene

# load GO file
go <- read.delim("/scratch/MORE24/Annotation/Prueba2-ok/go_annotations.txt",
                 header = FALSE,
                 stringsAsFactors = FALSE)
colnames(go) <- c("Gene","GO")
# Filter
go_ivory <- go[go$Gene %in% ivory_genes, ]
# Quick enrichment
table(go_ivory$GO)

install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.At.tair.db")  # optional Arabidopsis reference

ivory_genes <- ivory_table$Gene
# Extract GO annotation
go <- read.delim("/scratch/MORE24/Annotation/Prueba2-ok/go_annotations.txt",
                 header = FALSE,
                 stringsAsFactors = FALSE)

colnames(go) <- c("Gene", "GO")
# filer for ivory Module
go_ivory <- go[go$Gene %in% ivory_genes, ]

# Rub enrichment
library(clusterProfiler)

ego <- enricher(
  gene = go_ivory$Gene,
  TERM2GENE = go,
  pvalueCutoff = 0.05
)
barplot(ego, showCategory = 20)

dotplot(ego, showCategory = 20)
head(go)
# we have seen, Gene   GO:0008150^GO:0003674^GO:0005575 , split GO terms
library(tidyr)
library(dplyr)

# Rename columns properly
colnames(go) <- c("Gene", "GO")

# Split multiple GO terms
go_long <- go %>%
  separate_rows(GO, sep = "[,;^ ]+")

# clean data
go_long <- go_long[go_long$GO != "" & !is.na(go_long$GO), ]
# Prepare TERM2GENE
TERM2GENE <- go_long[, c("GO", "Gene")]

library(clusterProfiler)

ego <- enricher(
  gene = ivory_genes,
  TERM2GENE = TERM2GENE,
  pvalueCutoff = 0.05
)
ego 
barplot(ego, showCategory = 20)
dotplot(ego, showCategory = 20)

all_genes <- unique(go_long$Gene)

ego <- enricher(
  gene = ivory_genes,
  TERM2GENE = TERM2GENE,
  universe = all_genes,
  pvalueCutoff = 0.05
)
head(go)
head(go_long)




