# ==============================================================================
# GSE100026: RAPGEF1 Analysis in CML (QC + Partial Correlation)
# Repository Publication Script
# ==============================================================================

# Required packages installation guide (uncomment if needed):
# install.packages(c("svglite", "ggplot2", "reshape2", "pheatmap", "limma", "ggrepel", "dplyr", "ggpubr", "clusterProfiler", "org.Hs.eg.db"))

library(svglite) 
library(ggplot2)
library(reshape2)
library(pheatmap)
library(limma)
library(ggrepel)
library(dplyr)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

# ==============================================================================
# 1. PATH CONFIGURATION & DATA LOADING
# ==============================================================================

# Universal paths setup: Place the data file in your current R working directory
file_path  <- "GSE100026_expressed_gene_RPKM.txt.gz"
image_path <- "./Export"
dir.create(image_path, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(file_path)) {
  stop(paste("Input file not found at:", getwd(), "/", file_path, "\nPlease place the data file in your working directory."))
}

raw_data <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE)
rownames(raw_data) <- make.unique(as.character(raw_data[,1]))
rpkm_data <- raw_data[,-1]

sample_columns <- c("CML-BC11", "CML-BC12", "CML-BC13", "CML-BC14", "CML-BC15", 
                    "CML-CP10", "CML-CP6", "CML-CP7", "CML-CP8", "CML-CP9", 
                    "Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4", "Ctrl5")

numeric_matrix <- as.matrix(rpkm_data[, sample_columns])

metadata <- data.frame(
  SampleID = colnames(numeric_matrix),
  Condition = factor(rep(c("Blast_Crisis", "Chronic_Phase", "Control"), each = 5),
                     levels = c("Control", "Chronic_Phase", "Blast_Crisis"))
)
rownames(metadata) <- metadata$SampleID

# Filter and Log2 transformation
expressed_genes <- rowSums(numeric_matrix > 0.5) >= 3
filtered_matrix <- numeric_matrix[expressed_genes, ]
log_matrix <- log2(filtered_matrix + 1)

target_gene <- "RAPGEF1"
if(!target_gene %in% rownames(log_matrix)) stop(paste("Target gene", target_gene, "not found in matrix."))

# ==============================================================================
# 2. QUALITY CONTROL (QC) METRICS
# ==============================================================================

# QC 1: Total Library Size (Sum of Expression)
df_lib <- data.frame(
  Sample = metadata$SampleID,
  Condition = factor(metadata$Condition, labels = c("Healthy Control", "CML - CP", "CML - BC")),
  Total_Expression = colSums(numeric_matrix)
)

p_lib <- ggplot(df_lib, aes(x = Sample, y = Total_Expression, fill = Condition)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("#27ae60", "#2980b9", "#c0392b")) +
  labs(title = "QC: Total Sample Expression", x = "Samples", y = "Total RPKM Sum") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

ggsave(filename = file.path(image_path, "QC_1_Library_Size.svg"), plot = p_lib, width = 8, height = 5)

# QC 2: Expression Distribution (Boxplots)
df_melt <- melt(log_matrix)
colnames(df_melt) <- c("Gene", "Sample", "Log2RPKM")
df_melt$Condition <- metadata$Condition[match(df_melt$Sample, metadata$SampleID)]
levels(df_melt$Condition) <- c("Healthy Control", "CML - CP", "CML - BC")

p_dist <- ggplot(df_melt, aes(x = Sample, y = Log2RPKM, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.6) +
  scale_fill_manual(values = c("#27ae60", "#2980b9", "#c0392b")) +
  labs(title = "QC: Gene Expression Distribution", x = "Samples", y = "Log2(RPKM + 1)") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16))

ggsave(filename = file.path(image_path, "QC_2_Expression_Distribution.svg"), plot = p_dist, width = 9, height = 5)

# QC 3: Sample-to-Sample Correlation Heatmap
sample_cor <- cor(log_matrix, method = "pearson")
df_ann_qc <- data.frame(Condition = factor(metadata$Condition, labels = c("Healthy Control", "CML - CP", "CML - BC")))
rownames(df_ann_qc) <- metadata$SampleID
ann_colors_qc <- list(Condition = c("Healthy Control" = "#27ae60", "CML - CP" = "#2980b9", "CML - BC" = "#c0392b"))

svglite(filename = file.path(image_path, "QC_3_Sample_Correlation.svg"), width = 9, height = 8)
pheatmap(sample_cor, 
         annotation_col = df_ann_qc, 
         annotation_colors = ann_colors_qc,
         display_numbers = TRUE, number_color = "black", 
         fontsize = 12, fontsize_number = 9, angle_col = 45,
         color = colorRampPalette(c("white", "#e74c3c"))(100),
         main = "QC: Sample-to-Sample Pearson Correlation")
dev.off()

# ==============================================================================
# 3. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ==============================================================================

pca_res <- prcomp(t(log_matrix), scale. = FALSE) 
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- factor(metadata$Condition, 
                       levels = c("Control", "Chronic_Phase", "Blast_Crisis"),
                       labels = c("Healthy Control", "CML - Chronic Phase", "CML - Blast Crisis"))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  geom_point(size = 5, shape = 21, stroke = 1.2, color = "white") +
  scale_color_manual(values = c("#27ae60", "#2980b9", "#c0392b")) +
  scale_fill_manual(values = c("#27ae60", "#2980b9", "#c0392b")) +
  labs(
    title = "Principal Component Analysis (PCA)",
    subtitle = "Dataset: GSE100026 | RNA-seq Transcriptomic Profile",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)"),
    color = "Clinical Status", fill = "Clinical Status"
  ) +
  theme_bw(base_size = 15) +
  theme(text = element_text(family = "sans"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey30"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2))

ggsave(filename = file.path(image_path, "PCA_GSE100026_Publication.svg"), plot = p_pca, width = 9, height = 7)

# ==============================================================================
# 4. DIFFERENTIAL EXPRESSION & VOLCANO/MA PLOTS
# ==============================================================================

design <- model.matrix(~0 + metadata$Condition)
colnames(design) <- levels(metadata$Condition) 

fit <- lmFit(log_matrix, design)
contrast_matrix <- makeContrasts(
  CP_vs_Ctrl = Chronic_Phase - Control,
  BC_vs_Ctrl = Blast_Crisis - Control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE) 

stats_CP <- topTable(fit2, coef = "CP_vs_Ctrl", number = Inf)
stats_BC <- topTable(fit2, coef = "BC_vs_Ctrl", number = Inf)

# Volcano Plot Function
plot_volcano <- function(stats_df, plot_title, hl_gene) {
  stats_df$Status <- "Not Significant"
  stats_df$Status[stats_df$logFC > 1 & stats_df$adj.P.Val < 0.05] <- "Up-regulated"
  stats_df$Status[stats_df$logFC < -1 & stats_df$adj.P.Val < 0.05] <- "Down-regulated"
  
  max_x <- max(abs(stats_df$logFC), na.rm = TRUE)
  gene_hl <- stats_df[rownames(stats_df) == hl_gene, ]
  
  ggplot(stats_df, aes(x = logFC, y = -log10(adj.P.Val), color = Status)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_point(data = gene_hl, aes(x = logFC, y = -log10(adj.P.Val)), 
               color = "black", fill = "yellow", size = 5, shape = 21, stroke = 1.5) +
    geom_text_repel(data = gene_hl, aes(label = hl_gene), 
                    color = "black", fontface = "bold.italic", size = 6, box.padding = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_x_continuous(limits = c(-max_x, max_x)) +
    scale_color_manual(values = c("Down-regulated" = "#3498db", "Not Significant" = "grey80", "Up-regulated" = "#e74c3c")) +
    labs(title = plot_title, x = expression(bold(Log[2]~"Fold Change")), y = expression(bold("-Log"[10]~"Adjusted P-value (FDR)"))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          legend.position = "bottom", legend.text = element_text(size = 12))
}

p_cp <- plot_volcano(stats_CP, "DE: Chronic Phase vs Control", target_gene)
p_bc <- plot_volcano(stats_BC, "DE: Blast Crisis vs Control", target_gene)

ggsave(filename = file.path(image_path, "Volcano_CP_vs_Ctrl.svg"), plot = p_cp, width = 8, height = 8)
ggsave(filename = file.path(image_path, "Volcano_BC_vs_Ctrl.svg"), plot = p_bc, width = 8, height = 8)

# QC 4: MA Plot Function
plot_ma <- function(stats_df, plot_title, hl_gene) {
  stats_df$Status <- "Not Significant"
  stats_df$Status[stats_df$logFC > 1 & stats_df$adj.P.Val < 0.05] <- "Up-regulated"
  stats_df$Status[stats_df$logFC < -1 & stats_df$adj.P.Val < 0.05] <- "Down-regulated"
  
  gene_hl <- stats_df[rownames(stats_df) == hl_gene, ]
  
  ggplot(stats_df, aes(x = AveExpr, y = logFC, color = Status)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_point(data = gene_hl, aes(x = AveExpr, y = logFC), 
               color = "black", fill = "yellow", size = 5, shape = 21, stroke = 1.5) +
    geom_text_repel(data = gene_hl, aes(label = hl_gene), 
                    color = "black", fontface = "bold.italic", size = 6, box.padding = 1) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    scale_color_manual(values = c("Down-regulated" = "#3498db", "Not Significant" = "grey80", "Up-regulated" = "#e74c3c")) +
    labs(title = plot_title, x = "Average Log2 Expression", y = "Log2 Fold Change") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          legend.position = "bottom", legend.text = element_text(size = 12))
}

ma_cp <- plot_ma(stats_CP, "QC MA Plot: Chronic Phase vs Control", target_gene)
ma_bc <- plot_ma(stats_BC, "QC MA Plot: Blast Crisis vs Control", target_gene)

ggsave(filename = file.path(image_path, "QC_4_MA_CP_vs_Ctrl.svg"), plot = ma_cp, width = 8, height = 7)
ggsave(filename = file.path(image_path, "QC_4_MA_BC_vs_Ctrl.svg"), plot = ma_bc, width = 8, height = 7)

# ==============================================================================
# 5. TARGET GENE EXPRESSION BOXPLOT
# ==============================================================================

df_rapgef1 <- data.frame(
  SampleID = metadata$SampleID,
  Condition = factor(metadata$Condition, labels = c("Healthy Control", "CML - CP", "CML - BC")),
  Expression = as.numeric(log_matrix[target_gene, ])
)

# Statistical metrics parsing from Limma results
stat.test <- data.frame(
  group1 = "Healthy Control",
  group2 = c("CML - CP", "CML - BC"),
  p.adj = c(stats_CP[target_gene, "adj.P.Val"], stats_BC[target_gene, "adj.P.Val"])
)
stat.test$p.adj.signif <- symnum(stat.test$p.adj, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                 symbols = c("****", "***", "**", "*", "ns"))
y_max <- max(df_rapgef1$Expression)
stat.test$y.position <- c(y_max + 0.4, y_max + 0.9)

p_rapgef1 <- ggplot(df_rapgef1, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6, linewidth = 0.6) + 
  geom_jitter(width = 0.15, size = 3.5, shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = c("#27ae60", "#2980b9", "#c0392b")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) + 
  labs(
    title = bquote(italic(.(target_gene)) ~ "Expression Profile"),
    y = expression(bold(Log[2]~"(RPKM + 1)")), x = ""
  ) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", size = 8, bracket.size = 1) + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        legend.position = "none", 
        axis.text.x = element_text(face = "bold", color = "black", angle = 45, hjust = 1, size = 14))

ggsave(filename = file.path(image_path, paste0(target_gene, "_Boxplot.svg")), plot = p_rapgef1, width = 5, height = 6)

# ==============================================================================
# 6. PARTIAL CO-EXPRESSION ANALYSIS (Residual Matrix)
# ==============================================================================

# Step A: Regress out the effect of clinical condition
design_cond <- model.matrix(~ Condition, data = metadata)
fit_cond <- lmFit(log_matrix, design_cond)

# Step B: Extract residual variation across all genes
res_matrix <- residuals(fit_cond, log_matrix)
target_res <- as.numeric(res_matrix[target_gene, ])

# Step C: Spearman correlation on residuals
cor_vals <- suppressWarnings(cor(t(res_matrix), target_res, method = "spearman"))
p_vals <- apply(res_matrix, 1, function(x) {
  suppressWarnings(cor.test(x, target_res, method = "spearman")$p.value)
})

stats_partial <- data.frame(
  Gene = rownames(res_matrix),
  Partial_Correlation = as.numeric(cor_vals),
  P_Value = as.numeric(p_vals)
)

stats_partial$Adj_P_Value <- p.adjust(stats_partial$P_Value, method = "fdr")
stats_partial <- na.omit(stats_partial)

# Filtering for highly significant interactions
sig_partial <- stats_partial %>% 
  filter(Adj_P_Value < 0.05, abs(Partial_Correlation) > 0.45) %>% 
  pull(Gene)

if(length(sig_partial) >= 2) {
  ordered_samples <- metadata %>% arrange(Condition) %>% pull(SampleID)
  
  ordered_matrix_raw <- log_matrix[sig_partial, ordered_samples]
  ordered_matrix_res <- res_matrix[sig_partial, ordered_samples]
  
  df_ann <- data.frame(Condition = metadata$Condition[match(ordered_samples, metadata$SampleID)])
  rownames(df_ann) <- ordered_samples
  levels(df_ann$Condition) <- c("Healthy Control", "CML - CP", "CML - BC")
  
  ann_colors <- list(Condition = c("Healthy Control" = "#27ae60", "CML - CP" = "#2980b9", "CML - BC" = "#c0392b"))
  
  # Heatmap 1: Raw Expression Matrix
  col_dist <- dist(t(ordered_matrix_raw))
  col_hclust <- hclust(col_dist, method = "ward.D2")
  sample_weights <- ifelse(df_ann$Condition == "Healthy Control", -1000, 100)
  col_dend <- as.hclust(reorder(as.dendrogram(col_hclust), wts = sample_weights))
  
  svglite(filename = file.path(image_path, "Heatmap_1_Raw_Signature.svg"), width = 8, height = 9)
  pheatmap(ordered_matrix_raw, scale = "row", clustering_method = "ward.D2",
           cluster_cols = col_dend,  
           annotation_col = df_ann, annotation_colors = ann_colors, 
           show_rownames = FALSE, 
           fontsize = 12, fontsize_col = 11, angle_col = 45,
           color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
           main = "RAPGEF1 Network: Raw Expression (Disease Effect Present)")
  dev.off()
  
  # Heatmap 2: Residual Expression Matrix
  svglite(filename = file.path(image_path, "Heatmap_2_Partial_Signature.svg"), width = 8, height = 9)
  pheatmap(ordered_matrix_res, scale = "row", clustering_method = "ward.D2",
           cluster_cols = TRUE, 
           annotation_col = df_ann, annotation_colors = ann_colors, 
           show_rownames = FALSE, 
           fontsize = 12, fontsize_col = 11, angle_col = 45,
           color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
           main = "RAPGEF1 Network: Partial Co-expression (Confounder Removed)")
  dev.off()
  
} else {
  message("Warning: Insufficient genes passed filtering thresholds for co-expression profiling.")
}

# ==============================================================================
# 7. GENE ONTOLOGY (GO) ENRICHMENT ANALYSIS
# ==============================================================================

pos_genes <- stats_partial %>% filter(Partial_Correlation > 0, Adj_P_Value < 0.05) %>% pull(Gene) 

if(length(pos_genes) >= 10) {
  genes_entrez <- bitr(pos_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  go_res <- enrichGO(gene = genes_entrez, OrgDb = org.Hs.eg.db, ont = "BP", 
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  
  if(!is.null(go_res) && nrow(go_res) > 0) {
    p_dot <- enrichplot::dotplot(go_res, showCategory = 15) +
      labs(title = "GO Enrichment: Pure RAPGEF1 Correlated Genes") + 
      theme_bw(base_size = 14) +
      theme(axis.text.y = element_text(size = 12, face = "bold"), 
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
    
    ggsave(file.path(image_path, "Figure_GO_Dotplot_Partial.svg"), plot = p_dot, width = 10, height = 8)
  }
}

# ==============================================================================
# 8. SUB-PATHWAY GENE EXTRACTION
# ==============================================================================

go_df <- as.data.frame(go_res)
target_pathway <- "myeloid cell differentiation"
myeloid_data <- go_df %>% filter(Description == target_pathway)

if(nrow(myeloid_data) > 0) {
  myeloid_genes_string <- myeloid_data$geneID[1]
  myeloid_genes_list <- unlist(strsplit(myeloid_genes_string, "/"))
  
  message("=====================================================")
  message(paste("Genes associated with:", target_pathway))
  message("=====================================================")
  print(myeloid_genes_list)
  
} else {
  message("The specified pathway was not found in the significant GO results.")
}