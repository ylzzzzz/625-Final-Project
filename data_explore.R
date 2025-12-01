mydata <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata.csv")
gene_data <- mydata[, 7:ncol(mydata)]
clinical_data <- mydata[, 1:6]
batch_labels <- clinical_data$batch_id
disease_labels <- clinical_data$Alzheimers_Disease
length(disease_labels)

# deal with NA ------------------------------------------------------------
clean_NA <- function(X, max_na = 0.5, method = c("median", "mean")) {
  
  method <- match.arg(method)
  
  if (!is.matrix(X)) {
    X <- as.data.frame(X)
  }
  
  na_prop <- colMeans(is.na(X))
  too_many_na <- na_prop > max_na
  if (any(too_many_na)) {
    cat(sprintf("Removed %d variables with > %.0f%% missing values:\n", 
                sum(too_many_na), 100 * max_na))
    print(names(X)[too_many_na])
    X <- X[, !too_many_na, drop = FALSE]
  }
  
  
  cat("Imputing remaining NA values using", method, "...\n")
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[[j]])) {
      if (method == "median") {
        fill_value <- median(X[[j]], na.rm = TRUE)
      } else if (method == "mean") {
        fill_value <- mean(X[[j]], na.rm = TRUE)
      }
      X[[j]][is.na(X[[j]])] <- fill_value
    }
  }
  
  X <- as.matrix(sapply(X, as.numeric))
  
  if (anyNA(X)) {
    stop("There are still NA values after imputation!")
  } else {
    cat("all missing values handled successfully.\n")
    cat(sprintf("Final dimensions: %d samples x %d features\n", nrow(X), ncol(X)))
  }
  
  return(X)
}

genedata_clean <- clean_NA(gene_data, max_na = 0.5, method = "median")

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
