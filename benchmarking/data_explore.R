mydata <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata_final.csv")
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
  t(gene_data),
  batch = clinical_data$batch_id,
  design = model.matrix(~ clinical_data$Alzheimers_Disease)
)
data_batch_corrected <- t(batch_corrected)

library(uwot)

set.seed(123)
umap_result <- umap(gene_data, n_neighbors = 20, min_dist = 0.3, metric = "euclidean")
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$batch <- factor(batch_labels)
umap_df$disease <- factor(disease_labels, levels = c(0,1), labels = c("Healthy", "Diseased"))

library(ggplot2)

ggplot(umap_df, aes(UMAP1, UMAP2, color = batch)) +
  geom_point(size = 2, alpha = 0.8)  +
  theme_minimal() +
  labs(title = "UMAP projection of gene expression",
       color = "Batch Status") +
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
       color = "Batch Status") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(umap_df, aes(UMAP1, UMAP2, color = disease)) +
  geom_point(size = 2, alpha = 0.8)  +
  theme_minimal() +
  labs(title = "UMAP projection of gene expression",
       color = "Disease Status") +
  theme(plot.title = element_text(hjust = 0.5))



