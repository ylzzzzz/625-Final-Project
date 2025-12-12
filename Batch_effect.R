mydata <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata.csv")
gene_data <- mydata[, 7:ncol(mydata)]
clinical_data <- mydata[, 1:6]
batch_labels <- clinical_data$batch_id
disease_labels <- clinical_data$Alzheimers_Disease
length(disease_labels)

# library(sva)
# 
# combat_corrected <- ComBat(dat = as.matrix(t(genedata_clean)),
#                            batch = batch_labels,
#                            mod = disease_labels,
#                            par.prior = TRUE, 
#                            prior.plots = FALSE)
# data_batch_corrected <- t(combat_corrected)

library(limma)

batch_corrected <- removeBatchEffect(
  t(genedata_clean),
  batch = clinical_data$batch_id,
  design = model.matrix(~ clinical_data$Alzheimers_Disease)
)
data_batch_corrected <- t(batch_corrected)

library(uwot)

set.seed(123)
umap_result <- umap(genedata_clean, n_neighbors = 20, min_dist = 0.3, metric = "euclidean")
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$batch <- factor(batch_labels)
umap_df$disease <- factor(disease_labels, levels = c(0,1), labels = c("Healthy", "Diseased"))

library(ggplot2)

ggplot(umap_df, aes(UMAP1, UMAP2, color = batch)) +
  geom_point(size = 2, alpha = 0.8)  +
  theme_minimal() +
  labs(title = "UMAP projection of gene expression",
       color = "Disease Status") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(umap_df, aes(UMAP1, UMAP2, color = disease)) +
  geom_point(size = 2, alpha = 0.8)  +
  theme_minimal() +
  labs(title = "UMAP projection of gene expression",
       color = "Disease Status") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------------after batch effect correction
library(uwot)

set.seed(123)
umap_result <- umap(data_batch_corrected, n_neighbors = 20, min_dist = 0.3, metric = "euclidean")
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$batch <- factor(batch_labels)
umap_df$disease <- factor(disease_labels, levels = c(0,1), labels = c("Healthy", "Diseased"))

library(ggplot2)

ggplot(umap_df, aes(UMAP1, UMAP2, color = batch)) +
  geom_point(size = 2, alpha = 0.8)  +
  theme_minimal() +
  labs(title = "UMAP projection of gene expression",
       color = "Disease Status") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(umap_df, aes(UMAP1, UMAP2, color = disease)) +
  geom_point(size = 2, alpha = 0.8)  +
  theme_minimal() +
  labs(title = "UMAP projection of gene expression",
       color = "Disease Status") +
  theme(plot.title = element_text(hjust = 0.5))


#---------------------------DEA------------
design <- model.matrix(~ Alzheimers_Disease, data = clinical_data)
fit <- lmFit(batch_corrected, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "Alzheimers_Disease", adjust = "BH", number = Inf)
head(results)

sig_genes <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(sig_genes) 

library(EnhancedVolcano)
EnhancedVolcano(results,
                lab = rownames(results),
                x = "logFC",
                y = "adj.P.Val",
                title = "AD vs Control",
                pCutoff = 0.05,
                FCcutoff = 1)
library(pheatmap)
library(RColorBrewer)

# 提取表达矩阵
topGenes <- rownames(head(results, 10))
mat <- data_batch_corrected[, topGenes]
mat_scaled <- t(scale(t(mat)))          # 按基因Z-score标准化
mat_scaled[is.na(mat_scaled)] <- 0

# 重新定义注释表
annotation_col <- data.frame(
  Disease = factor(clinical_data$Alzheimers_Disease,
                   levels = c(0, 1),
                   labels = c("Control", "AD"))
)
rownames(annotation_col) <- clinical_data$Sample_ID
annotation_col <- annotation_col[colnames(mat_scaled), , drop = FALSE]

# 手动设定颜色映射
ann_colors <- list(
  Disease = c(Control = "#4daf4a", AD = "#e41a1c")
)

# 绘制热图
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,      
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Top 10 DEGs (Batch-corrected)"
)




