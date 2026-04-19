######  bulk RNA-seq downstream analysis #############

#Loading required libraries:

library(AnnotationDbi)
library(biomaRt)
library(tximport)
library(tidyverse)
library(ggfortify)
library(DESeq2)
library(sva)
library(RColorBrewer)
library(patchwork)
library(ggdendro)
library(corrplot)
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(ggplot2)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(ReportingTools)
library(magrittr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)  # Use correct organism database


setwd("/path/to/bulk_RNA_seq_analysis")

D3_sampleinfo <- read.table("sample_sheet.txt", header = T, stringsAsFactors = T)
rownames(D3_sampleinfo) <- D3_sampleinfo$SampleName
D3_sampleinfo$type %<>% relevel("DMSO")

files <- str_c("Salmon_output/", D3_sampleinfo$SampleName, "/quant.sf")
files <- set_names(files, D3_sampleinfo$SampleName)
tx2gene <- read_tsv("Salmon_output/tx2gene.tsv", col_names = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

D3_rawCounts <- round(txi$counts, 0) 
# check dimension of count matrix
dim(D3_rawCounts)
all(rownames(D3_sampleinfo) == colnames(D3_rawCounts))

#Keep genes with a minimum of 5 fragments (arbitrary) across all samples:
keep <- rowSums(D3_rawCounts) > 5
table(keep, useNA="always")
D3_filtCounts <- D3_rawCounts[keep,]
dim(D3_filtCounts)

# Get log2 counts
logcounts <- log2(D3_filtCounts + 1)
summary(logcounts[,1]) # summary for first column 
summary(logcounts) # summary for each column
# make a colour vector

# Ensure 'type' is character
D3_sampleinfo$type <- as.character(D3_sampleinfo$type)

# Define color mapping
TreatmentCols <- setNames(c("#E5E5E5", "#EEA9B8", "#FFA500", "#9ACD32"), 
                          c("DMSO", "Treatment_1", "Treatment_2", "Treatment_1_2"))[D3_sampleinfo$type]

# Debug to ensure correct color mapping
print(unique(D3_sampleinfo$type))  # Should be exactly "DMSO", "Treatment_1", "Treatment_2", "Treatment_1_2"
print(unique(TreatmentCols))  # Should not have NA

# Generate boxplot
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=TreatmentCols,
        main="Log2(Counts)")

# Let's add a blue horizontal line that corresponds to the median
abline(h=apply(as.data.frame(apply(logcounts, 2, median)), 2, median), col="blue")

vst_counts <- vst(D3_rawCounts)
# Check distributions of samples using boxplots
boxplot(vst_counts,
        xlab="",
        ylab="VST counts",
        las=2,
        col=TreatmentCols)
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(vst_counts), col="blue")
# VST counts standard deviation (sd) vs mean expression
plot(rowMeans(vst_counts), rowSds(vst_counts),
     main='VST counts: sd vs mean')

rlogcounts <- rlog(D3_filtCounts)
# run PCA
pcDat <- prcomp(t(rlogcounts)) # plot PCA

#################

# Read sample information
sample_info <- read.table("sample_sheet.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(sample_info) <- sample_info$SampleName
sample_info$type <- as.factor(sample_info$type)
sample_info$type <- relevel(sample_info$type, ref = "DMSO")  # Set DMSO as control

# Define file paths for Salmon outputs
files <- file.path("Salmon_output", sample_info$SampleName, "quant.sf")
names(files) <- sample_info$SampleName

# Load transcript-to-gene mapping
tx2gene <- read_tsv("Salmon_output/tx2gene.tsv", col_names = FALSE)

# Import transcript-level quantification
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = sample_info, design = ~type)
dds <- DESeq(dds)  # Run DESeq2 normalization and differential expression


vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pca_data <- plotPCA(vsd, intgroup = "type", returnData = TRUE)
ggplot(pca_data, aes(PC1, PC2, color = type)) +
  geom_point(size = 4) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()


sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists), 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample Distance Heatmap")

################# For Treatment_1 #############################

res_Treatment_1 <- results(dds, contrast = c("type", "Treatment_1", "DMSO"))  # Example: Treatment_1 vs DMSO
res_Treatment_1 <- res_Treatment_1[order(res_Treatment_1$padj), ]  # Sort by adjusted p-value
summary(res_Treatment_1)


# Convert DESeqResults to a dataframe
res_Treatment_1_df <- as.data.frame(res_Treatment_1)

# Ensure there are no NA values in padj
res_Treatment_1_df <- dplyr::filter(res_Treatment_1_df, !is.na(padj))

# Compute -log10(padj)
res_Treatment_1_df$log10padj <- -log10(res_Treatment_1_df$padj)

# write the dataframe to csv file.
write.csv(as.data.frame(res_Treatment_1_df), "Treatment_1_df_addditive_model.csv")

# ------- get mapping to Ensembl gene ID, biotype ----------- #
# Connect to Ensembl (human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract gene symbols from DE results
gene_symbols_Treatment_1 <- rownames(res_Treatment_1_df)

# Get annotations
annot_Treatment_1 <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "gene_biotype"),
               filters = "hgnc_symbol",
               values = gene_symbols_Treatment_1,
               mart = mart)

# Add gene symbol column to DE result for joining
res_Treatment_1_df$hgnc_symbol <- rownames(res_Treatment_1_df)

# Merge
res_annotated_Treatment_1 <- left_join(res_Treatment_1_df, annot_Treatment_1, by = "hgnc_symbol")

# ------------------------------------------------------------#

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


################# For Treatment_2 #############################

res_Treatment_2 <- results(dds, contrast = c("type", "Treatment_2", "DMSO"))  # Example: Treatment_2 vs DMSO
res_Treatment_2 <- res_Treatment_2[order(res_Treatment_2$padj), ]  # Sort by adjusted p-value
summary(res_Treatment_2)


# Convert DESeqResults to a dataframe
res_Treatment_2_df <- as.data.frame(res_Treatment_2)

# Ensure there are no NA values in padj
res_Treatment_2_df <- dplyr::filter(res_Treatment_2_df, !is.na(padj))

# Compute -log10(padj)
res_Treatment_2_df$log10padj <- -log10(res_Treatment_2_df$padj)

# write the dataframe to csv file.
write.csv(as.data.frame(res_Treatment_2_df), "Treatment_2_df_addditive_model.csv")

# ------- get mapping to Ensembl gene ID, biotype ----------- #

# Connect to Ensembl (human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract gene symbols from DE results
gene_symbols_Treatment_2 <- rownames(res_Treatment_2_df)

# Get annotations
annot_Treatment_2 <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "gene_biotype"),
                    filters = "hgnc_symbol",
                    values = gene_symbols_Treatment_2,
                    mart = mart)

# Add gene symbol column to DE result for joining
res_Treatment_2_df$hgnc_symbol <- rownames(res_Treatment_2_df)

# Merge
res_annotated_Treatment_2 <- left_join(res_Treatment_2_df, annot_Treatment_2, by = "hgnc_symbol")

# ------------------------------------------------------------#


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


################# For Treatment_1_2 #############################

res_Treatment_1_2 <- results(dds, contrast = c("type", "Treatment_1_2", "DMSO"))  # Example: Treatment_1_2 vs DMSO
res_Treatment_1_2 <- res_Treatment_1_2[order(res_Treatment_1_2$padj), ]  # Sort by adjusted p-value
summary(res_Treatment_1_2)


# Convert DESeqResults to a dataframe
res_Treatment_1_2_df <- as.data.frame(res_Treatment_1_2)

# Ensure there are no NA values in padj
res_Treatment_1_2_df <- dplyr::filter(res_Treatment_1_2_df, !is.na(padj))

# Compute -log10(padj)
res_Treatment_1_2_df$log10padj <- -log10(res_Treatment_1_2_df$padj)

# write the dataframe to csv file.
write.csv(as.data.frame(res_Treatment_1_2_df), "Treatment_1_2_df_addditive_model.csv")

# ------- get mapping to Ensembl gene ID, biotype ----------- #
# Connect to Ensembl (human)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Extract gene symbols from DE results
gene_symbols_Treatment_1_2 <- rownames(res_Treatment_1_2_df)

# Get annotations
annot_Treatment_1_2 <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "gene_biotype"),
                    filters = "hgnc_symbol",
                    values = gene_symbols_Treatment_1_2,
                    mart = mart)

# Add gene symbol column to DE result for joining
res_Treatment_1_2_df$hgnc_symbol <- rownames(res_Treatment_1_2_df)

# Merge
res_annotated_Treatment_1_2 <- left_join(res_Treatment_1_2_df, annot_Treatment_1_2, by = "hgnc_symbol")

# ------------------------------------------------------------#




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

########   PCA plot 

# Perform PCA on variance-stabilized data
vsd <- vst(dds, blind = TRUE)  # Variance-stabilizing transformation
pca_data <- plotPCA(vsd, intgroup = "type", returnData = TRUE)

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 4) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme_minimal()

######   Heatmap for sample clustering 

# Treatment_2ect top variable genes for clustering
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

# Plot heatmaps
# 1. color scale represents raw variance-stabilized counts (VST units).
pheatmap(assay(vsd)[topVarGenes,], 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = as.data.frame(colData(dds)),
         main = "each cell in heatmap = variance-stabilized expression value"
         )

#2. 
pheatmap(assay(vsd)[topVarGenes,], 
         scale = "row",
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = as.data.frame(colData(dds)),
         main = "Relative expression per gene (z-scored)"
         )


############   Gene Expression Tables (FPKM)

# Convert raw counts to FPKM
dds <- estimateSizeFactors(dds)
fpkm_values <- fpkm(dds)  # Get FPKM values

# Save as table
write.csv(as.data.frame(fpkm_values), "FPKM_gene_expression.csv")


#############  Differentially Expressed Gene Lists (LFC & FDR)

# Run DESeq2 differential analysis
dds <- DESeq(dds)

# Get results for each comparison
res_Treatment_1 <- results(dds, contrast = c("type", "Treatment_1", "DMSO"))
res_Treatment_2  <- results(dds, contrast = c("type", "Treatment_2", "DMSO"))
res_Treatment_1_2 <- results(dds, contrast = c("type", "Treatment_1_2", "DMSO"))

# Filter significant DEGs (FDR < 0.05)
res_Treatment_1_sig <- res_Treatment_1[which(res_Treatment_1$padj < 0.05), ]
res_Treatment_2_sig <- res_Treatment_2[which(res_Treatment_2$padj < 0.05), ]
res_Treatment_1_2_sig <- res_Treatment_1_2[which(res_Treatment_1_2$padj < 0.05), ]

# Save DEGs
write.csv(res_Treatment_1_sig, "DEGs_Treatment_1_vs_DMSO.csv")
write.csv(res_Treatment_2_sig, "DEGs_Treatment_2_vs_DMSO.csv")
write.csv(res_Treatment_1_2_sig, "DEGs_Treatment_1_2_vs_DMSO.csv")


#############    Volcano Plots for Each Comparison 

####  Treatment_1.  #####

# Convert results to data frame
res_Treatment_1_df <- as.data.frame(res_Treatment_1) %>% filter(!is.na(padj))
res_Treatment_1_df$log10padj <- -log10(res_Treatment_1_df$padj)

# Volcano Plot
ggplot(res_Treatment_1_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Treatment_1 Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()


#######  Treatment_2 #######

# Convert results to data frame
res_Treatment_2_df <- as.data.frame(res_Treatment_2) %>% filter(!is.na(padj))
res_Treatment_2_df$log10padj <- -log10(res_Treatment_2_df$padj)

# Volcano Plot
ggplot(res_Treatment_2_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Treatment_2 Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()


########   Treatment_1_2 #####

# Convert results to data frame
res_Treatment_1_2_df <- as.data.frame(res_Treatment_1_2) %>% filter(!is.na(padj))
res_Treatment_1_2_df$log10padj <- -log10(res_Treatment_1_2_df$padj)

# Volcano Plot
ggplot(res_Treatment_1_2_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Treatment_1_2 Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
  theme_minimal()


##############################################################
#########  Pathway Analysis (GO, GSEA, KEGG)
######    Gene Ontology (GO)

### Treatment_1 #######

# Convert gene IDs
geneList <- res_Treatment_1_sig$log2FoldChange
names(geneList) <- rownames(res_Treatment_1_sig)

# Run GO analysis
ego <- enrichGO(gene = names(geneList),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
barplot(ego, showCategory = 20)

write.csv(as.data.frame(ego), "GO_BP_enrichGO_Treatment_1.csv")

  #### GSEA

# Create named numeric vector: log2FoldChange values named by gene symbols
geneList <- res_Treatment_1_df$log2FoldChange
names(geneList) <- rownames(res_Treatment_1_df)

# Sort in decreasing order (most upregulated genes first)
geneList <- sort(geneList, decreasing = TRUE)

gsea_res <- gseGO(geneList = geneList,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
dotplot(gsea_res)

write.csv(as.data.frame(gsea_res), "GO_BP_gseGO_Treatment_1.csv")

  #### KEGG pathway analysis


# Convert Gene Symbols to Entrez IDs
gene_df <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Ensure geneList has Entrez IDs as names
geneList_entrez <- geneList[names(geneList) %in% gene_df$SYMBOL]  # Filter mapped genes
names(geneList_entrez) <- gene_df$ENTREZID[match(names(geneList_entrez), gene_df$SYMBOL)]  # Rename with Entrez IDs

# Run KEGG enrichment
kk <- enrichKEGG(gene = names(geneList_entrez), 
                 organism = "hsa",  # Use "mmu" for mouse
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.05)

# Check results
if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
  print("No significant KEGG pathways found. Try increasing pvalueCutoff (e.g., 0.1).")
} else {
  dotplot(kk)
}

write.csv(as.data.frame(kk), "GO_BP_enrichKEGG_Treatment_1.csv")


#### Treatment_2 ###############

# Convert gene IDs
geneList <- res_Treatment_2_sig$log2FoldChange
names(geneList) <- rownames(res_Treatment_2_sig)

# Run GO analysis
ego <- enrichGO(gene = names(geneList),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
barplot(ego, showCategory = 20)

write.csv(as.data.frame(ego), "GO_BP_enrichGO_Treatment_2.csv")

#### GSEA

# Create named numeric vector: log2FoldChange values named by gene symbols
geneList <- res_Treatment_2_df$log2FoldChange
names(geneList) <- rownames(res_Treatment_2_df)

# Sort in decreasing order (most upregulated genes first)
geneList <- sort(geneList, decreasing = TRUE)

gsea_res <- gseGO(geneList = geneList,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
dotplot(gsea_res)

write.csv(as.data.frame(gsea_res), "GO_BP_gseGO_Treatment_2.csv")

#### KEGG pathway analysis


# Convert Gene Symbols to Entrez IDs
gene_df <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Ensure geneList has Entrez IDs as names
geneList_entrez <- geneList[names(geneList) %in% gene_df$SYMBOL]  # Filter mapped genes
names(geneList_entrez) <- gene_df$ENTREZID[match(names(geneList_entrez), gene_df$SYMBOL)]  # Rename with Entrez IDs

# Run KEGG enrichment
kk <- enrichKEGG(gene = names(geneList_entrez), 
                 organism = "hsa",  # Use "mmu" for mouse
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.05)

# Check results
if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
  print("No significant KEGG pathways found. Try increasing pvalueCutoff (e.g., 0.1).")
} else {
  dotplot(kk)
}

write.csv(as.data.frame(kk), "GO_BP_enrichKEGG_Treatment_2.csv")

#### Treatment_1_2 ##############################

# Convert gene IDs
geneList <- res_Treatment_1_2_sig$log2FoldChange
names(geneList) <- rownames(res_Treatment_1_2_sig)

# Run GO analysis
ego <- enrichGO(gene = names(geneList),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
barplot(ego, showCategory = 20)

write.csv(as.data.frame(ego), "GO_BP_enrichGO_Treatment_1_2.csv")

#### GSEA

# Create named numeric vector: log2FoldChange values named by gene symbols
geneList <- res_Treatment_1_2_df$log2FoldChange
names(geneList) <- rownames(res_Treatment_1_2_df)

# Sort in decreasing order (most upregulated genes first)
geneList <- sort(geneList, decreasing = TRUE)

gsea_res <- gseGO(geneList = geneList,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
dotplot(gsea_res)

write.csv(as.data.frame(gsea_res), "GSEA_BP_gseGO_Treatment_1_2.csv")

#### KEGG pathway analysis


# Convert Gene Symbols to Entrez IDs
gene_df <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Ensure geneList has Entrez IDs as names
geneList_entrez <- geneList[names(geneList) %in% gene_df$SYMBOL]  # Filter mapped genes
names(geneList_entrez) <- gene_df$ENTREZID[match(names(geneList_entrez), gene_df$SYMBOL)]  # Rename with Entrez IDs

# Run KEGG enrichment
kk <- enrichKEGG(gene = names(geneList_entrez), 
                 organism = "hsa",  # Use "mmu" for mouse
                 pAdjustMethod = "BH", 
                 pvalueCutoff = 0.05)

# Check results
if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
  print("No significant KEGG pathways found. Try increasing pvalueCutoff (e.g., 0.1).")
} else {
  dotplot(kk)
}

write.csv(as.data.frame(kk), "GO_BP_enrichKEGG_Treatment_1_2.csv")


############ ###################### #####################

# Summary
# Method	Use res_interaction_sig?	Notes
# gseGO()	 No	Uses ranked gene list, not a threshold
# gseKEGG()	 No	Same as above
# enrichGO()	 Yes	Requires a gene set of significant genes
# enrichKEGG()	 Yes	Same
# enrichPathway()	Yes	Same


