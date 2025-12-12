library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

mydata = read.csv("./data/Alzheimers_Disease_final.csv")
gene_data = mydata[, 7:ncol(mydata)]
clinical_data = mydata[, 1:6]
batch_labels = clinical_data$batch_id
disease_labels = clinical_data$Alzheimers_Disease
length(disease_labels)

#Loop for each DEA
# Identify unique cohorts
unique_batches = unique(clinical_data$batch_id)
#unique_batches = unique_batches[]

# Create a list to store the results for each cohort
cohort_results_list = list()

# Iterate through each cohort
for (cohort in unique_batches) {
  cat(paste0("\nProcessing Cohort (Batch): ", cohort, "...\n"))
  
  # 1. Subset Data
  idx = which(clinical_data$batch_id == cohort)
  sub_clinical = clinical_data[idx, ]
  sub_gene = gene_data[idx, ] # Changed from genedata_clean to gene_data
  
  # Check for both groups
  if (length(unique(sub_clinical$Alzheimers_Disease)) < 2) {
    cat(paste("  [!] Skipping cohort", cohort, "- needs both Case and Control.\n"))
    next
  }
  
  # 2. Run Limma (DEA)
  design_sub = model.matrix(~ Alzheimers_Disease, data = sub_clinical)
  fit_sub = lmFit(t(sub_gene), design_sub)
  fit_sub = eBayes(fit_sub)
  res_sub = topTable(fit_sub, coef = "Alzheimers_Disease", adjust = "BH", number = Inf)
  
  # Store results
  cohort_results_list[[as.character(cohort)]] = res_sub
  sig_count = sum(res_sub$adj.P.Val < 0.05 & abs(res_sub$logFC) > 1)
  cat(paste("  -> Finished. Significant genes found:", sig_count, "\n"))
  
  # 3. Generate Volcano Plot
  volcano_plot <- EnhancedVolcano(res_sub,
                                  lab = rownames(res_sub),
                                  x = "logFC",
                                  y = "adj.P.Val",
                                  title = paste0("Cohort ", cohort, ": AD vs Control"),
                                  subtitle = paste(sig_count, "sig genes"),
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  legendPosition = 'right',
                                  caption = "")
  
  ggsave(filename = paste0("./figures/Cohort_", cohort, "_AD_vs_Control_volcano.png"), 
         plot = volcano_plot, width = 8, height = 6, dpi = 300)
  
  # --- 4. Heatmap Section (FIXED) ---
  
  # A. Select Top Genes
  top_n <- 20
  top_genes <- rownames(head(res_sub, top_n))
  
  # B. Prepare Matrix [Genes x Samples]
  # We transpose because pheatmap expects Genes as Rows
  mat <- t(sub_gene[, top_genes])
  
  # C. Filter Zero Variance Rows (Fixes warning/scaling issues)
  # Calculate variance for each gene (row)
  gene_vars <- apply(mat, 1, var)
  # Keep only genes with variance > 0
  mat <- mat[gene_vars > 0, , drop = FALSE]
  
  # D. Scale (Z-score normalization)
  # Scale works on columns, so we transpose twice: t(scale(t(mat)))
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0 # Safety replace NAs
  
  # E. Prepare Annotation
  # Force input to character to ensure factor levels match '0' and '1' safely
  clean_disease <- as.character(sub_clinical$Alzheimers_Disease)
  
  annotation_col <- data.frame(
    Disease = factor(clean_disease,
                     levels = c("0", "1"), 
                     labels = c("Control", "AD"))
  )
  
  # CRITICAL: Match rownames of annotation to column names of matrix
  rownames(annotation_col) <- colnames(mat_scaled)
  
  # F. Define Colors
  # The list name "Disease" MUST match the column name in annotation_col "Disease"
  ann_colors <- list(
    Disease = c(Control = "#4daf4a", AD = "#e41a1c")
  )
  
  # G. Generate Heatmap
  heatmap_filename <- paste0("./figures/Cohort_", cohort, "_Top", top_n, "_Heatmap.png")
  
  tryCatch({
    pheatmap(
      mat_scaled,
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = FALSE, 
      main = paste0("Cohort ", cohort, ": Top ", top_n, " DEGs"),
      filename = heatmap_filename,
      width = 8,
      height = 6
    )
    Sys.sleep(5)
  }, error = function(e) {
    cat(paste("  [!] Heatmap failed for cohort", cohort, ":", e$message, "\n"))
  })
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



heatmap_filename <- paste0("./figures/Cohort_", cohort, "_Top", top_n, "_Heatmap.png")

tryCatch({
  pheatmap(logfc_matrix,
           main = "LogFC Comparison Across Cohorts (Sig. Genes Only)",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cluster_cols = FALSE, # Keep cohort order (or TRUE to group similar cohorts)
           cluster_rows = TRUE,
           show_rownames = TRUE, # Hide row names if too many genes
           na_col = "grey") # Grey if a gene was missing in a specific cohort
}, error = function(e) {
  cat(paste("  [!] Heatmap failed for cohort", cohort, ":", e$message, "\n"))
})
