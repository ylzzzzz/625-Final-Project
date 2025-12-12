library(limma)
library(EnhancedVolcano)

mydata = read.csv("C:/Users/trelo/Downloads/Alzheimers_Disease_cleandata.csv")
gene_data = mydata[, 7:ncol(mydata)]
clinical_data = mydata[, 1:6]
batch_labels = clinical_data$batch_id
disease_labels = clinical_data$Alzheimers_Disease
length(disease_labels)

# deal with NA ------------------------------------------------------------
clean_NA = function(X, max_na = 0.5, method = c("median", "mean")) {
  
  method = match.arg(method)
  
  if (!is.matrix(X)) {
    X = as.data.frame(X)
  }
  
  na_prop = colMeans(is.na(X))
  too_many_na = na_prop > max_na
  if (any(too_many_na)) {
    cat(sprintf("Removed %d variables with > %.0f%% missing values:\n", 
                sum(too_many_na), 100 * max_na))
    print(names(X)[too_many_na])
    X = X[, !too_many_na, drop = FALSE]
  }
  
  
  cat("Imputing remaining NA values using", method, "...\n")
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[[j]])) {
      if (method == "median") {
        fill_value = median(X[[j]], na.rm = TRUE)
      } else if (method == "mean") {
        fill_value = mean(X[[j]], na.rm = TRUE)
      }
      X[[j]][is.na(X[[j]])] = fill_value
    }
  }
  
  X = as.matrix(sapply(X, as.numeric))
  
  if (anyNA(X)) {
    stop("There are still NA values after imputation!")
  } else {
    cat("all missing values handled successfully.\n")
    cat(sprintf("Final dimensions: %d samples x %d features\n", nrow(X), ncol(X)))
  }
  
  return(X)
}

genedata_clean = clean_NA(gene_data, max_na = 0.5, method = "median")



#Loop for each DEA
# Identify unique cohorts
unique_batches = unique(clinical_data$batch_id)

# Create a list to store the results for each cohort
cohort_results_list = list()

# Iterate through each cohort
for (cohort in unique_batches) {
  cat(paste0("\nProcessing Cohort (Batch): ", cohort, "...\n"))
  
  # Get indices for this specific batch
  idx = which(clinical_data$batch_id == cohort)
  
  # Subset clinical and gene data
  sub_clinical = clinical_data[idx, ]
  # Note: genedata_clean is Samples x Genes, so we slice rows
  sub_gene = genedata_clean[idx, ] 
  
  # We need both Healthy (0) and Diseased (1) samples to run a comparison
  if (length(unique(sub_clinical$Alzheimers_Disease)) < 2) {
    cat(paste("  [!] Skipping cohort", cohort, "- missing one of the study groups (needs both Case and Control).\n"))
    next
  }
  
  # Run Limma (DEA)
  # Create design matrix for THIS cohort only
  # We do NOT include batch_id in the design because we are inside a single batch
  design_sub = model.matrix(~ Alzheimers_Disease, data = sub_clinical)
  
  # Transpose sub_gene because limma expects Genes as Rows, Samples as Columns
  fit_sub = lmFit(t(sub_gene), design_sub)
  fit_sub = eBayes(fit_sub)
  
  # Get top table
  res_sub = topTable(fit_sub, coef = "Alzheimers_Disease", adjust = "BH", number = Inf)
  
  # --- D. Store and Print Summary ---
  cohort_results_list[[as.character(cohort)]] = res_sub
  
  sig_count = sum(res_sub$adj.P.Val < 0.05 & abs(res_sub$logFC) > 1)
  cat(paste("  -> Finished. Significant genes found:", sig_count, "\n"))
  
  # Generate Volcano Plot per Cohort
  # This will create a plot in your 'Plots' pane for each loop iteration
  #print(
  #  EnhancedVolcano(res_sub,
  #                  lab = rownames(res_sub),
  #                  x = "logFC",
  #                  y = "adj.P.Val",
  #                  title = paste0("Cohort ", cohort, ": AD vs Control"),
  #                  subtitle = paste(sig_count, "sig genes"),
  #                  pCutoff = 0.05,
  #                  FCcutoff = 1,
  #                  legendPosition = 'right',
  #                  caption = "")
  #)
}


################################################################################
# Results
# Check which cohorts were analyzed
names(cohort_results_list)

# View the top results for a specific batch (e.g., Batch "1")
# Replace "1" with whatever your actual batch IDs are
head(cohort_results_list[["1"]])

# Example: Write all results to separate CSV files
for (cohort_name in names(cohort_results_list)) {
  write.csv(cohort_results_list[[cohort_name]], 
            file = paste0("DEA_Results_Cohort_", cohort_name, ".csv"))
}

################################################################################
# Comparison of common h=genes
# Extract significant gene names from each cohort
sig_genes_list = lapply(cohort_results_list, function(df) {
  rownames(df)[df$adj.P.Val < 0.05 & abs(df$logFC) > 1]
})

# Find genes common to all analyzed cohorts
common_genes = Reduce(intersect, sig_genes_list)

cat("\nGenes significant in ALL cohorts:\n")
print(common_genes)

################################################################################
library(EnhancedVolcano)
library(cowplot) # Install if needed: install.packages("cowplot")

# 1. Filter: Keep only cohorts with significant genes
# We define significance as adj.P.Val < 0.05 and abs(logFC) > 1
plots_list = list()

for (cohort_name in names(cohort_results_list)) {
  
  res = cohort_results_list[[cohort_name]]
  
  # Check if there are ANY genes meeting the criteria
  sig_genes = subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
  
  if (nrow(sig_genes) > 0) {
    
    # Create the plot object but DO NOT print it yet. Save it to a list.
    p = EnhancedVolcano(res,
                         lab = rownames(res),
                         x = "logFC",
                         y = "adj.P.Val",
                         title = paste0("Cohort ", cohort_name),
                         subtitle = paste(nrow(sig_genes), "sig genes"),
                         pCutoff = 0.05,
                         FCcutoff = 1,
                         legendPosition = 'none', # Hide legend to save space in the grid
                         caption = "")
    
    plots_list[[cohort_name]] = p
  } else {
    cat(paste("Skipping Cohort", cohort_name, "- 0 significant genes found.\n"))
  }
}

# 2. Plot All Together
if (length(plots_list) > 0) {
  cat(paste("Plotting", length(plots_list), "cohorts with significant findings...\n"))
  
  # Arrange in a grid
  # ncol = 2 means 2 plots per row. Adjust as needed.
  final_grid = plot_grid(plotlist = plots_list, ncol = 2, labels = "AUTO")
  
  # Print the grid
  print(final_grid)
  
} else {
  cat("No cohorts had significant genes based on current cutoffs.")
}

################################################################################
library(pheatmap)

# 1. Identify all unique significant genes across the filtered cohorts
# (Using the list we filtered above)
all_sig_genes = c()
cohorts_with_sig = names(plots_list) # Uses the names from the previous step

for (coh in cohorts_with_sig) {
  df = cohort_results_list[[coh]]
  sigs = rownames(df)[df$adj.P.Val < 0.05 & abs(df$logFC) > 1]
  all_sig_genes = c(all_sig_genes, sigs)
}
unique_sig_genes = unique(all_sig_genes)

# 2. Build a matrix of LogFC values
# Rows = Genes, Columns = Cohorts
logfc_matrix = matrix(NA, nrow = length(unique_sig_genes), ncol = length(cohorts_with_sig))
rownames(logfc_matrix) = unique_sig_genes
colnames(logfc_matrix) = cohorts_with_sig

for (coh in cohorts_with_sig) {
  # Extract LogFC for these genes from the full results of that cohort
  full_res = cohort_results_list[[coh]]
  
  # Match gene names
  matches = full_res[unique_sig_genes, "logFC"]
  logfc_matrix[, coh] = matches
}

# 3. Plot Heatmap of Log Fold Changes
# Red = Upregulated in AD, Blue = Downregulated in AD
pheatmap(logfc_matrix,
         main = "LogFC Comparison Across Cohorts (Sig. Genes Only)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_cols = FALSE, # Keep cohort order (or TRUE to group similar cohorts)
         cluster_rows = TRUE,
         show_rownames = FALSE, # Hide row names if too many genes
         na_col = "grey") # Grey if a gene was missing in a specific cohort

