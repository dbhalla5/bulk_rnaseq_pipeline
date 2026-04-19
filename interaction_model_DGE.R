### R version 4.4.3

###### bulk RNA-seq downstream analysis #############

### Load required libraries
library(tximport)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(magrittr)
library(ggrepel)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

### Set working directory
setwd("/path/to/RNAseq_analysis")

### Load sample information
sample_info <- read.table("sample_sheet.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$SampleName
sample_info$type <- factor(sample_info$type, levels = c("DMSO", "Treatment_1", "Treatment_2", "Treatment_1_2"))

### Load quantification files
files <- file.path("Salmon_output", sample_info$SampleName, "quant.sf")
names(files) <- sample_info$SampleName
tx2gene <- read_tsv("Salmon_output/tx2gene.tsv", col_names = FALSE)

### Import transcript counts
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

### Build DESeq2 dataset
#dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~type)

# Convert counts to integers by rounding
txi_counts_int <- round(txi$counts)

# Now create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = txi_counts_int,
                              colData = sample_info,
                              design = ~type)



### Filter genes with low counts
keep <- rowSums(counts(dds)) > 5
dds <- dds[keep,]

### Run DESeq2
dds <- DESeq(dds)

### Variance Stabilizing Transformation
vsd <- vst(dds, blind = FALSE)
vst_counts <- assay(vsd)

### Define colors for treatments
treatment_colors <- setNames(c("#E5E5E5", "#EEA9B8", "#FFA500", "#9ACD32"),
                             levels(sample_info$type))

### Log2-transformed raw counts for quick QC plots
logcounts <- log2(counts(dds, normalized = TRUE) + 1)

### QC Plot: Boxplot of log2 counts
boxplot(logcounts,
        xlab = "", ylab = "Log2(Counts)",
        las = 2, col = treatment_colors[sample_info$type],
        main = "Log2(Normalized Counts)")
abline(h = median(logcounts), col = "blue")

### QC Plot: Boxplot of VST counts
boxplot(vst_counts,
        xlab = "", ylab = "VST Counts",
        las = 2, col = treatment_colors[sample_info$type],
        main = "VST-Transformed Counts")
abline(h = median(vst_counts), col = "blue")

### QC Plot: SD vs Mean (VST)
plot(rowMeans(vst_counts), rowSds(vst_counts),
     main = "VST Counts: SD vs Mean",
     xlab = "Mean Expression", ylab = "Standard Deviation")

### PCA Plot
pca_data <- plotPCA(vsd, intgroup = "type", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = type, label = rownames(pca_data))) +
  geom_point(size = 3) +
  geom_text_repel() +
  scale_color_manual(values = treatment_colors) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of Samples")

### Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Distance Heatmap")


################# 1. For Treatment_1 #############################

res_Treatment_1 <- results(dds, contrast = c("type", "Treatment_1", "DMSO"))  # Example: Treatment_1 vs DMSO
res_Treatment_1 <- res_Treatment_1[order(res_Treatment_1$padj), ]  # Sort by adjusted p-value
summary(res_Treatment_1)


# Convert DESeqResults to a dataframe
res_Treatment_1_df <- as.data.frame(res_Treatment_1)

# Ensure there are no NA values in padj
res_Treatment_1_df <- dplyr::filter(res_Treatment_1_df, !is.na(padj))

# Compute -log10(padj)
res_Treatment_1_df$log10padj <- -log10(res_Treatment_1_df$padj)

# Volcano Plot
ggplot(res_Treatment_1_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

plotMA(res_Treatment_1, main = "Treatment_1 MA Plot", ylim = c(-5, 5))


top_genes <- rownames(head(res_Treatment_1[order(res_Treatment_1$padj), ], 50))  # Top 50 genes
pheatmap(assay(vsd)[top_genes, ], 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = sample_info, 
         main = "Treatment_1 Heatmap of Top DEGs")


# Convert gene names to Entrez IDs
gene_list_Treatment_1 <- res_Treatment_1$log2FoldChange
names(gene_list_Treatment_1) <- rownames(res_Treatment_1)
gene_list_Treatment_1 <- sort(gene_list_Treatment_1, decreasing = TRUE)

# Run GSEA
gsea_results_Treatment_1 <- gseGO(geneList = gene_list_Treatment_1, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
dotplot(gsea_results_Treatment_1)

############################################################


################# 2. For Treatment_2 #############################

res_Treatment_2 <- results(dds, contrast = c("type", "Treatment_2", "DMSO"))  # Example: Treatment_2 vs DMSO
res_Treatment_2 <- res_Treatment_2[order(res_Treatment_2$padj), ]  # Sort by adjusted p-value
summary(res_Treatment_2)


# Convert DESeqResults to a dataframe
res_Treatment_2_df <- as.data.frame(res_Treatment_2)

# Ensure there are no NA values in padj
res_Treatment_2_df <- dplyr::filter(res_Treatment_2_df, !is.na(padj))

# Compute -log10(padj)
res_Treatment_2_df$log10padj <- -log10(res_Treatment_2_df$padj)

# Volcano Plot
ggplot(res_Treatment_2_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

plotMA(res_Treatment_2, main = "Treatment_2 MA Plot", ylim = c(-5, 5))


top_genes <- rownames(head(res_Treatment_2[order(res_Treatment_2$padj), ], 50))  # Top 50 genes
pheatmap(assay(vsd)[top_genes, ], 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = sample_info, 
         main = "Treatment_2 Heatmap of Top DEGs")


# Convert gene names to Entrez IDs
gene_list_Treatment_2 <- res_Treatment_2$log2FoldChange
names(gene_list_Treatment_2) <- rownames(res_Treatment_2)
gene_list_Treatment_2 <- sort(gene_list_Treatment_2, decreasing = TRUE)

# Run GSEA
gsea_results_Treatment_2 <- gseGO(geneList = gene_list_Treatment_2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
dotplot(gsea_results_Treatment_2)

############################################################


################# 3. For Treatment_1_2 #############################

res_Treatment_1_2 <- results(dds, contrast = c("type", "Treatment_1_2", "DMSO"))  # Example: Treatment_1_2 vs DMSO
res_Treatment_1_2 <- res_Treatment_1_2[order(res_Treatment_1_2$padj), ]  # Sort by adjusted p-value
summary(res_Treatment_1_2)


# Convert DESeqResults to a dataframe
res_Treatment_1_2_df <- as.data.frame(res_Treatment_1_2)

# Ensure there are no NA values in padj
res_Treatment_1_2_df <- dplyr::filter(res_Treatment_1_2_df, !is.na(padj))

# Compute -log10(padj)
res_Treatment_1_2_df$log10padj <- -log10(res_Treatment_1_2_df$padj)

# Volcano Plot
ggplot(res_Treatment_1_2_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

plotMA(res_Treatment_1_2, main = "Treatment_1_2 MA Plot", ylim = c(-5, 5))


top_genes <- rownames(head(res_Treatment_1_2[order(res_Treatment_1_2$padj), ], 50))  # Top 50 genes
pheatmap(assay(vsd)[top_genes, ], 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = sample_info, 
         main = "Treatment_1_2 Heatmap of Top DEGs")


# Convert gene names to Entrez IDs
gene_list_Treatment_1_2 <- res_Treatment_1_2$log2FoldChange
names(gene_list_Treatment_1_2) <- rownames(res_Treatment_1_2)
gene_list_Treatment_1_2 <- sort(gene_list_Treatment_1_2, decreasing = TRUE)

# Run GSEA
gsea_results_Treatment_1_2 <- gseGO(geneList = gene_list_Treatment_1_2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")
dotplot(gsea_results_Treatment_1_2)

############################################################

################# 4. Interaction model of Treatment_1_2 #############################


# Extract treatment group from column names
dds$treatment <- factor(sub("_.*", "", colnames(dds)),
                        levels = c("DMSO", "Treatment_1", "Treatment_2", "Treatment_1_2"))

# Define presence of each drug
dds$drugA <- factor(ifelse(dds$treatment %in% c("Treatment_1", "Treatment_1_2"), "yes", "no"), levels = c("no", "yes"))
dds$drugB <- factor(ifelse(dds$treatment %in% c("Treatment_2",  "Treatment_1_2"), "yes", "no"), levels = c("no", "yes"))

# Set up interaction design
design(dds) <- ~ drugA + drugB + drugA:drugB

# Run DESeq
dds <- DESeq(dds)


vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pca_data <- plotPCA(vsd, intgroup = "type", returnData = TRUE)
# OUTPUT OF ABOVE:pca_data <- plotPCA(vsd, intgroup = "type", returnData = TRUE)
# using ntop=500 top features by variance

sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists), 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Distance Heatmap")


# Extract interaction term (synergistic effect of Treatment_1_2)
res_interaction <- results(dds, name = "drugAyes.drugByes")
res_interaction <- res_interaction[order(res_interaction$padj),] # Sort by adjusted p-value

summary(res_interaction)

# OUTPUT OF THE COMMAND 'summary(res_interaction)'

  # out of 16583 with nonzero total read count
  # adjusted p-value < 0.1
  # LFC > 0 (up)       : 70, 0.42%
  # LFC < 0 (down)     : 110, 0.66%
  # outliers [1]       : 17, 0.1%
  # low counts [2]     : 2570, 15%
  # (mean count < 5)
  # [1] see 'cooksCutoff' argument of ?results
  # [2] see 'independentFiltering' argument of ?results



# Convert DESeqResults to a dataframe
res_interaction_df <- as.data.frame(res_interaction)

# Ensure there are no NA values in padj
res_interaction_df <- dplyr::filter(res_interaction_df, !is.na(padj))

# Compute -log10(padj)
res_interaction_df$log10padj <- -log10(res_interaction_df$padj)

# write the dataframe to csv file.
write.csv(as.data.frame(res_interaction_df), "Treatment_1_2_df_interaction_model.csv")

# Volcano Plot
ggplot(res_interaction_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()

# MA plot
plotMA(res_interaction, main = "interaction model MA Plot", ylim = c(-5, 5))

# heatmap drawn from top 50 genes
top_genes <- rownames(head(res_interaction[order(res_interaction$padj), ], 50))  # Top 50 genes
pheatmap(assay(vsd)[top_genes, ], 
         scale = "row", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = sample_info, 
         main = "interaction model Heatmap of Top DEGs")


# Convert gene names to Entrez IDs
gene_list_interaction <- res_interaction_df$log2FoldChange
names(gene_list_interaction) <- rownames(res_interaction_df)
gene_list_interaction <- sort(gene_list_interaction, decreasing = TRUE)


# Run GSEA
gsea_results_interaction <- gseGO(
  geneList = gene_list_interaction,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2  # less stringent
)


# Run GSEA
# gsea_results_interaction <- gseGO(geneList = gene_list_interaction, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")


# OUTPUT OF THE COMMAND ABOVE 'gsea_results_interaction <- gseGO(geneList = gene_list_interaction, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH")'
# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
# preparing geneSet collections...
# GSEA analysis...
# no term enriched under specific pvalueCutoff...
# Warning message:
#   In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                There are ties in the preranked stats (7.32% of the list).
#                              The order of those tied genes will be arbitrary, which may produce unexpected results.
# 

# Plot only if enriched
if (nrow(gsea_results_interaction@result) > 0) {
  dotplot(gsea_results_interaction)
} else {
  message("No enriched GO terms found for interaction model (even with relaxed cutoff).")
}

# The following is only to be used in 
# ----- overrepresentation analysis
# enrichGO
# enrichKEGG
# enrichPathway
# ----- These methods ask: "Are known pathways overrepresented among my significant DEGs?" ----

res_interaction_sig <- res_interaction[which(res_interaction$padj < 0.05), ]

write.csv(as.data.frame(res_interaction_sig), "Interaction_model_DEGs_padj_005.csv")

# 1. Overrepresentation Analysis (ORA) on significant DEGs
sig_genes <- rownames(res_interaction_sig)

# KEGG
kegg_results <- enrichKEGG(
  gene         = bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)
  
  # OUTPUT OF THE ABOVE COMMAND 'kegg_results...
  # Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
  # Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
  # 'Treatment_2ect()' returned 1:1 mapping between keys and columns
  # Warning message:
  #   In bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) :
  #   5.15% of input gene IDs are fail to map...

# Reactome
# reactome_results <- enrichPathway(
#   gene         = bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
#   organism     = "human",
#   pvalueCutoff = 0.05
# )

# OUTPUT OF THE ABOVE COMMAND 'reactome_results...'

  # 'Treatment_2ect()' returned 1:1 mapping between keys and columns
  # Warning message:
  #   In bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) :
  #   6.11% of input gene IDs are fail to map..



# 2. Visualization
dotplot(kegg_results, showCategory = 50)
# dotplot(reactome_results, showCategory = 20)

############

geneList <- res_interaction_sig$log2FoldChange
names(geneList) <- rownames(res_interaction_sig)

# Run GO analysis
ego <- enrichGO(gene = names(geneList),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
barplot(ego, showCategory = 20)

write.csv(as.data.frame(ego), "GO_BP_enrichGO_interaction.csv")

## kegg pathway analysis for synergistic effect (only relevant for Treatment_1_2)

# Convert gene IDs
geneList <- res_interaction_sig$log2FoldChange
names(geneList) <- rownames(res_interaction_sig)

# Run GO analysis

ego <- enrichKEGG(
  gene         = bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)


barplot(ego, showCategory = 20)

write.csv(as.data.frame(ego), "GO_BP_enrichKEGG_interaction.csv")







