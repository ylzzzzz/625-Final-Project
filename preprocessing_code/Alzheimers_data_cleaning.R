cran_packages = c("jsonlite", "dplyr", "tidyr", "data.table")
bioc_packages = c("sva")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
missing_cran = cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)
missing_bioc = bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc)

library(jsonlite)
library(dplyr)
library(tidyr)
library(data.table)
library(sva)

#######################################################################################

# Load and apply gene_synonym.json
gene_map_path = "gene_synonym.json" 
gene_map_list = fromJSON(gene_map_path)

gene_map_dt = data.table(
  alias = names(gene_map_list),
  canonical = unlist(gene_map_list),
  stringsAsFactors = FALSE
)
setkey(gene_map_dt, alias)
all_known_genes = unique(c(gene_map_dt$alias, gene_map_dt$canonical))

# Define clinical columns
clinical_columns = c(
  "Alzheimers_Disease", 
  "Age", 
  "Gender"
)



base_path = "data/Alzheimers_Disease"

# Find all .csv files directly in this folder
combined_files = list.files(base_path, pattern = "\\.csv$", full.names = FALSE)

# Filter out any sub-folder CSVs if they exist (seperate clinical_data and gene_data files)
combined_files = combined_files[!grepl("clinical_data/|gene_data/", combined_files)]

if (length(combined_files) == 0) {
  stop(paste("No .csv files found directly in:", base_path))
}

# Get the Batch IDs from the filenames
batch_ids = sub("\\.csv$", "", combined_files)

message(sprintf("Found %d combined datasets to process: %s", 
                length(batch_ids), paste(batch_ids, collapse = ", ")))

load_csv = function(batch_id, base_path, clinical_list) {
  message(sprintf("Processing: %s", batch_id))
  # Define file path
  file_path = file.path(base_path, paste0(batch_id, ".csv"))
  
  if (!file.exists(file_path)) {
    warning(sprintf("  SKIPPING: File not found: %s", file_path))
    return(NULL)
  }
  # Load Combined Data (Samples x (Features + Genes))
  data = as.data.frame(fread(file_path))
  
  # Set row names (Sample IDs) from the first column
  if (any(duplicated(data[[1]]))) {
    warning(sprintf("  SKIPPING %s: Duplicate sample IDs found in first column.", batch_id))
    return(NULL)
  }
  rownames(data) = data[[1]]
  data = data[, -1, drop = FALSE]
  
  # Separate clinical and gene columns
  all_cols = colnames(data)
  
  # Identify clinical columns: any column name in our master list
  clinical_cols = all_cols[all_cols %in% clinical_list]
  
  # Identify gene columns: everything *else*
  gene_cols = all_cols[!all_cols %in% clinical_list]
  
  if (length(gene_cols) == 0) {
    warning(sprintf("  SKIPPING %s: No columns were identified as genes. (Check clinical_columns list)", batch_id))
    return(NULL)
  }
  if (length(clinical_cols) == 0) {
    warning(sprintf("  SKIPPING %s: No clinical columns found. (Check clinical_columns list)", batch_id))
    return(NULL)
  }
  
  # Create the two data objects
  clinical_data = data[, clinical_cols, drop = FALSE]
  gene_data = data[, gene_cols, drop = FALSE]
  
  clinical_data$batch_id = batch_id
  
  # transpose gene data
  expr_data = t(gene_data)
  
  # Ensure all gene expression data is numeric
  expr_data = apply(expr_data, 2, as.numeric)
  rownames(expr_data) = colnames(gene_data) # Transpose swaps row/col names
  
  message(sprintf("  Loaded %s: %d samples, %d genes, %d clinical features.", 
                  batch_id, ncol(expr_data), nrow(expr_data), ncol(clinical_data) - 1))
  
  # Return both as a list
  return(list(
    expr_data = expr_data,
    clinical_data = clinical_data
  ))
}

# Load all datasets
message("--- Starting data load loop ---")

# We now pass clinical_columns to the function
all_data_list = lapply(batch_ids, load_csv, 
                        base_path = base_path, 
                        clinical_list = clinical_columns)

# Filter out any NULLs (from folders that failed to load)
all_data_list = all_data_list[!sapply(all_data_list, is.null)]

if (length(all_data_list) == 0) {
  stop("Failed to load any datasets. Check file paths and clinical_columns list.")
}
message(sprintf("Successfully loaded %d datasets.", length(all_data_list)))


# Feature Alignment
message("Starting feature alignment...")

# This is the same alignment function as before
align_and_aggregate_dt = function(expr_matrix, gene_map_dt) {
  expr_dt = as.data.table(expr_matrix, keep.rownames = "alias")
  expr_dt = gene_map_dt[expr_dt, on = "alias"] # Left join
  expr_dt[is.na(canonical), gene_final := alias]
  expr_dt[!is.na(canonical), gene_final := canonical]
  expr_dt[, `:=`(alias = NULL, canonical = NULL)]
  expr_aggregated = expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_final]
  
  expr_df = as.data.frame(expr_aggregated)
  rownames(expr_df) = expr_df$gene_final
  expr_df$gene_final = NULL
  return(as.matrix(expr_df))
}

# Apply the alignment function to all datasets
aligned_expr_list = lapply(all_data_list, function(data) {
  batch_name = data$clinical_data$batch_id[1]
  message(sprintf("  Aligning batch: %s", batch_name)) 
  align_and_aggregate_dt(data$expr_data, gene_map_dt)
})

# Find the UNION of all genes
message("Finding the union (all unique genes) across all datasets...")
all_genes = Reduce(union, lapply(aligned_expr_list, rownames))
message(sprintf("Found %d unique genes in total across all datasets.", length(all_genes)))

# Re-index all matrices to the new master gene list
# This will fill missing genes with 'NA'
message("Re-indexing matrices... (This may take a moment)")
final_expr_list = lapply(aligned_expr_list, function(expr_matrix) {
  
  sample_names = colnames(expr_matrix)
  
  # Create a new, empty matrix filled with NAs
  new_matrix = matrix(NA, 
                       nrow = length(all_genes), 
                       ncol = length(sample_names),
                       dimnames = list(all_genes, sample_names))
  
  # Find the genes that *this* matrix actually has
  common_genes = intersect(rownames(expr_matrix), all_genes)
  
  # Fill in the data we have
  new_matrix[common_genes, ] = expr_matrix[common_genes, ]
  
  return(new_matrix)
})

message("Feature alignment and union complete.")


# Combine data
message("Combining data...")

# Combine all expression matrices by column (cbind)
# This works now because all matrices have the exact same rows (all_genes)
combined_expression_matrix = do.call(cbind, final_expr_list)
message(sprintf("Combined Expression Matrix (with NAs) dimensions: %s", 
                paste(dim(combined_expression_matrix), collapse = " x ")))

# Combine all clinical data frames by row (bind_rows)
clinical_data_list = lapply(all_data_list, `[[`, "clinical_data")
combined_clinical_data = bind_rows(clinical_data_list)
message(sprintf("Combined Clinical Data dimensions: %s", 
                paste(dim(combined_clinical_data), collapse = " x ")))

# Re-order expression matrix columns to match clinical data rows
sample_order = rownames(combined_clinical_data)
combined_expression_matrix = combined_expression_matrix[, sample_order]

if (ncol(combined_expression_matrix) != nrow(combined_clinical_data)) {
  stop("Sample mismatch after combining. This should not happen.")
}
message("Sample counts match. Data is not yet ready for ComBat.")



# Standardize biological variable 
message("Standardizing biological variable...")

messy_column = "Alzheimers_Disease" 

if (!messy_column %in% names(combined_clinical_data)) {
  warning(sprintf("Column '%s' not found. Using the first clinical column as a guess.", messy_column))
  messy_column = names(combined_clinical_data)[1] 
  message(sprintf("Using fallback column: %s", messy_column))
}

# Create a new 'bio_group' column
combined_clinical_data = combined_clinical_data %>%
  mutate(
    bio_group = case_when(
      grepl("0|control|normal|healthy", .data[[messy_column]], ignore.case = TRUE) ~ "Control",
      grepl("1|alzheimer|ad|dementia", .data[[messy_column]], ignore.case = TRUE) ~ "AD",
      TRUE ~ NA_character_
    )
  )

message("Groups after standardization (pre-filtering):")
print(table(combined_clinical_data$bio_group, useNA = "ifany"))

# Filter out NA samples from the clinical data
# We are NOT filtering the expression matrix here.
is_na_rows = is.na(combined_clinical_data$bio_group)
n_removed = sum(is_na_rows)

if (n_removed > 0) {
  message(sprintf("Note: Removing %d samples that had 'NA' in bio_group.", n_removed))
  combined_clinical_data = combined_clinical_data[!is_na_rows, ]
  combined_expression_matrix = combined_expression_matrix[, !is_na_rows]
}


message("Creating final merged dataframe... (ComBat was NOT run)")

# For modeling, transpose the expression data so samples are rows
# This matrix will be (Samples x Genes) and will contain NAs
final_expression_df = t(combined_expression_matrix)
final_expression_df = as.data.frame(final_expression_df)

message(sprintf("Final expression data frame dimensions: %s",
                paste(dim(final_expression_df), collapse = " x ")))

# Merge with clinical data (by rownames)
final_modeling_df = merge(combined_clinical_data, 
                           final_expression_df,
                           by = 0) # 'by = 0' merges by rownames

# Tidy up the 'Row.names' column created by merge
final_modeling_df = final_modeling_df %>%
  rename(Sample_ID = Row.names)

message(sprintf("Successfully created final modeling data frame with %d samples and %d total features.",
                nrow(final_modeling_df), ncol(final_modeling_df)))

message("--- SCRIPT COMPLETE ---")
message("This data contains ALL genes, including those with NA values.")
message("Showing first 6 rows and 10 columns:")
print(head(final_modeling_df[, 1:10]))



