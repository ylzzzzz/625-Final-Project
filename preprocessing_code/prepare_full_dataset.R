# full dataset compiler

# -----------------------------------------------------------------------------
# SETUP & PACKAGES
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

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

# !!! IMPORTANT !!!
# You must list ALL potential disease status column names here.
# If a folder has a column "Acute_Myeloid_Leukemia" and it is NOT in this list,
# the script will try to treat it as a gene and fail.
clinical_columns = c(
  "Alzheimers_Disease", 
  "Acute_Myeloid_Leukemia", # Example: Add your other folder column names here
  "Adrenocortical_Cancer",  # Example
  "Allergies",              # Example
  "Age", 
  "Gender"
)

# -----------------------------------------------------------------------------
# DATA LOADING FUNCTIONS
# -----------------------------------------------------------------------------

# Function to load a single CSV batch
load_csv = function(batch_id, current_folder_path, clinical_list) {
  # Define file path
  file_path = file.path(current_folder_path, paste0(batch_id, ".csv"))
  
  if (!file.exists(file_path)) {
    warning(sprintf("  SKIPPING: File not found: %s", file_path))
    return(NULL)
  }
  
  # Load Combined Data (Samples x (Features + Genes))
  # suppressor warnings for column type guessing during initial fread
  data = as.data.frame(fread(file_path, showProgress = FALSE))
  
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
  
  # Add batch_id and also track which Disease Folder this came from
  clinical_data$batch_id = batch_id
  clinical_data$source_folder = basename(current_folder_path)
  
  # transpose gene data
  expr_data = t(gene_data)
  
  # Ensure all gene expression data is numeric
  expr_data = apply(expr_data, 2, as.numeric)
  rownames(expr_data) = colnames(gene_data) 
  
  message(sprintf("  Loaded %s: %d samples, %d genes.", 
                  batch_id, ncol(expr_data), nrow(expr_data)))
  
  return(list(
    expr_data = expr_data,
    clinical_data = clinical_data
  ))
}

# Function to process an entire folder (e.g., "data/Alzheimers_Disease")
process_folder = function(folder_path, clinical_list) {
  message(sprintf("Entering folder: %s", folder_path))
  
  # Find all .csv files directly in this folder
  combined_files = list.files(folder_path, pattern = "\\.csv$", full.names = FALSE)
  
  # Filter out sub-folder CSVs (metadata files)
  combined_files = combined_files[!grepl("clinical_data/|gene_data/", combined_files)]
  
  if (length(combined_files) == 0) {
    warning(sprintf("No valid .csv files found in %s", folder_path))
    return(NULL)
  }
  
  # Get Batch IDs
  batch_ids = sub("\\.csv$", "", combined_files)
  
  # Apply load_csv to every file in THIS folder
  folder_data = lapply(batch_ids, load_csv, 
                       current_folder_path = folder_path, 
                       clinical_list = clinical_list)
  
  return(folder_data)
}

# -----------------------------------------------------------------------------
# MAIN EXECUTION LOOP
# -----------------------------------------------------------------------------

message("--- Starting Global Data Load ---")

# 1. Define Root Path
root_data_path = "data"

# 2. Find all sub-directories in 'data'
# recursive=FALSE ensures we only get the immediate disease folders
all_disease_folders = list.dirs(root_data_path, full.names = TRUE, recursive = FALSE)

if (length(all_disease_folders) == 0) {
  stop("No subdirectories found in 'data'.")
}

message(sprintf("Found %d disease folders: %s", 
                length(all_disease_folders), 
                paste(basename(all_disease_folders), collapse = ", ")))

# 3. Loop through every folder and load data
# This returns a Nested List: List of Folders -> List of Batches
nested_data_list = lapply(all_disease_folders, process_folder, clinical_list = clinical_columns)

# 4. Flatten the list
# We want one big list of batches, regardless of which folder they came from
all_data_list = unlist(nested_data_list, recursive = FALSE)

# Filter out NULLs (failed loads)
all_data_list = all_data_list[!sapply(all_data_list, is.null)]

if (length(all_data_list) == 0) {
  stop("Failed to load any datasets. Check file paths and clinical_columns list.")
}

message(sprintf("Successfully loaded %d total datasets across all folders.", length(all_data_list)))


# -----------------------------------------------------------------------------
# FEATURE ALIGNMENT
# -----------------------------------------------------------------------------
message("Starting feature alignment...")

align_and_aggregate_dt = function(expr_matrix, gene_map_dt) {
  expr_dt = as.data.table(expr_matrix, keep.rownames = "alias")
  expr_dt = gene_map_dt[expr_dt, on = "alias"] 
  expr_dt[is.na(canonical), gene_final := alias]
  expr_dt[!is.na(canonical), gene_final := canonical]
  expr_dt[, `:=`(alias = NULL, canonical = NULL)]
  expr_aggregated = expr_dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_final]
  
  expr_df = as.data.frame(expr_aggregated)
  rownames(expr_df) = expr_df$gene_final
  expr_df$gene_final = NULL
  return(as.matrix(expr_df))
}

aligned_expr_list = lapply(all_data_list, function(data) {
  # Silent alignment to reduce console spam
  align_and_aggregate_dt(data$expr_data, gene_map_dt)
})

message("Finding the union of all genes...")
all_genes = Reduce(union, lapply(aligned_expr_list, rownames))
message(sprintf("Found %d unique genes in total.", length(all_genes)))

message("Re-indexing matrices... (This may take a moment)")
final_expr_list = lapply(aligned_expr_list, function(expr_matrix) {
  sample_names = colnames(expr_matrix)
  new_matrix = matrix(NA, 
                      nrow = length(all_genes), 
                      ncol = length(sample_names),
                      dimnames = list(all_genes, sample_names))
  common_genes = intersect(rownames(expr_matrix), all_genes)
  new_matrix[common_genes, ] = expr_matrix[common_genes, ]
  return(new_matrix)
})


# -----------------------------------------------------------------------------
# COMBINE DATA
# -----------------------------------------------------------------------------
message("Combining data...")

combined_expression_matrix = do.call(cbind, final_expr_list)

clinical_data_list = lapply(all_data_list, `[[`, "clinical_data")
combined_clinical_data = bind_rows(clinical_data_list)

# Re-order expression matrix columns to match clinical data rows
sample_order = rownames(combined_clinical_data)
combined_expression_matrix = combined_expression_matrix[, sample_order]

message("Data combined.")

# -----------------------------------------------------------------------------
# STANDARDIZATION
# -----------------------------------------------------------------------------
message("Standardizing biological variable...")

# Because we have multiple diseases now, we need to find the target column dynamically
# We look for the first column that matches one of our known disease names
# (excluding Age, Gender, batch_id, source_folder)
possible_disease_cols = setdiff(names(combined_clinical_data), 
                                c("Age", "Gender", "batch_id", "source_folder"))

# We assume the first valid disease column found for a row is the one we want
# This creates a coalesced column "Disease_Status_Raw"
combined_clinical_data = combined_clinical_data %>%
  mutate(Disease_Status_Raw = coalesce(!!!syms(possible_disease_cols)))

message("Raw status labels found:")
print(table(combined_clinical_data$Disease_Status_Raw))

# Standardize to Control vs Case
combined_clinical_data = combined_clinical_data %>%
  mutate(
    bio_group = case_when(
      grepl("0|control|normal|healthy", Disease_Status_Raw, ignore.case = TRUE) ~ "Control",
      grepl("1|alzheimer|ad|dementia|case|tumor|cancer|leukemia|positive", Disease_Status_Raw, ignore.case = TRUE) ~ "Case",
      TRUE ~ NA_character_
    )
  )

message("Groups after standardization:")
print(table(combined_clinical_data$bio_group, useNA = "ifany"))

is_na_rows = is.na(combined_clinical_data$bio_group)
if (sum(is_na_rows) > 0) {
  message(sprintf("Removing %d samples with unknown status.", sum(is_na_rows)))
  combined_clinical_data = combined_clinical_data[!is_na_rows, ]
  combined_expression_matrix = combined_expression_matrix[, !is_na_rows]
}

# -----------------------------------------------------------------------------
# FINAL OUTPUT
# -----------------------------------------------------------------------------

final_expression_df = as.data.frame(t(combined_expression_matrix))

final_modeling_df = merge(combined_clinical_data, 
                          final_expression_df,
                          by = 0) 

final_modeling_df = final_modeling_df %>%
  rename(Sample_ID = Row.names)

message("--- SCRIPT COMPLETE ---")
print(head(final_modeling_df[, 1:10]))