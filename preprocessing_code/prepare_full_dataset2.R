# MEMORY SAFE DATASET COMPILER

# -----------------------------------------------------------------------------
# SETUP
# -----------------------------------------------------------------------------
cran_packages = c("jsonlite", "dplyr", "tidyr", "data.table")
bioc_packages = c("sva")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
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

# Load gene map
gene_map_path = "gene_synonym.json" 
gene_map_list = fromJSON(gene_map_path)
gene_map_dt = data.table(
  alias = names(gene_map_list),
  canonical = unlist(gene_map_list),
  stringsAsFactors = FALSE
)
setkey(gene_map_dt, alias)

# Clinical columns to extract
clinical_columns = c(
  "Alzheimers_Disease", 
  "Acute_Myeloid_Leukemia", 
  "Adrenocortical_Cancer",  
  "Allergies",              
  "Age", 
  "Gender"
)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS (FIXED)
# -----------------------------------------------------------------------------

process_single_file = function(file_info, clinical_list, map_dt) {
  
  full_path = file_info$path
  batch_id = file_info$batch_id
  folder = file_info$folder
  
  message(sprintf("Processing: %s", batch_id))
  
  # 1. Fast Load
  # We read as data.frame to ensure safe handling before conversion
  data = as.data.frame(fread(full_path, showProgress = FALSE))
  
  # 2. Check duplicates
  if (any(duplicated(data[[1]]))) {
    warning(sprintf("  SKIPPING %s: Duplicate samples.", batch_id))
    return(NULL)
  }
  rownames(data) = data[[1]]
  data = data[, -1, drop = FALSE]
  
  # 3. Separate Columns
  all_cols = colnames(data)
  clinical_cols = all_cols[all_cols %in% clinical_list]
  gene_cols = all_cols[!all_cols %in% clinical_list]
  
  if (length(gene_cols) == 0 || length(clinical_cols) == 0) return(NULL)
  
  # 4. Extract Clinical Data
  clin_df = data[, clinical_cols, drop = FALSE]
  clin_df$batch_id = batch_id
  clin_df$source_folder = folder
  clin_df$Sample_ID = rownames(clin_df)
  
  # 5. Process Expression Data (THE FIX)
  # Extract gene data
  gene_df = data[, gene_cols, drop = FALSE]
  
  # Convert to matrix first (this keeps it as character if there are mixed types)
  expr_mat = as.matrix(gene_df)
  
  # Transpose (Genes x Samples)
  expr_mat = t(expr_mat)
  
  # SAFE NUMERIC CONVERSION
  # storage.mode is faster than apply() and forces conversion. 
  # Non-numeric strings will become NA (which is what we want).
  storage.mode(expr_mat) = "numeric"
  
  # 6. Align Features (Using data.table)
  # keep.rownames creates a column named "alias"
  expr_dt = as.data.table(expr_mat, keep.rownames = "alias")
  
  # Clean up raw data to free RAM
  rm(data, gene_df, expr_mat) 
  
  # Join with gene map
  expr_dt = map_dt[expr_dt, on = "alias"] 
  
  # Handle canonical names
  expr_dt[is.na(canonical), gene_final := alias]
  expr_dt[!is.na(canonical), gene_final := canonical]
  
  # Aggregate (Mean)
  # We explicitly specify .SDcols to ensure we only average the sample columns
  sample_cols = names(expr_dt)[3:(ncol(expr_dt)-1)] # Skip alias, canonical. gene_final is last.
  
  # Note: sample_cols are the sample IDs.
  # We need to average columns 3 to (ncol-1). 
  # Actually, let's just exclude the ID columns explicitly.
  cols_to_mean = setdiff(names(expr_dt), c("alias", "canonical", "gene_final"))
  
  expr_aggregated = expr_dt[, lapply(.SD, mean, na.rm = TRUE), 
                            by = gene_final, 
                            .SDcols = cols_to_mean]
  
  # Convert back to matrix
  final_mat = as.matrix(expr_aggregated[, -1])
  rownames(final_mat) = expr_aggregated$gene_final
  
  # Force Garbage Collection
  rm(expr_dt, expr_aggregated)
  gc()
  
  return(list(clinical = clin_df, expr = final_mat))
}

# -----------------------------------------------------------------------------
# EXECUTION
# -----------------------------------------------------------------------------

message("--- Starting Memory-Optimized Load ---")

# 1. Build File List first (Do not load yet)
root_data_path = "data"
all_disease_folders = list.dirs(root_data_path, full.names = TRUE, recursive = FALSE)
file_manifest = list()

for (folder in all_disease_folders) {
  files = list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  files = files[!grepl("clinical_data/|gene_data/", files)]
  if(length(files) > 0) {
    batch_ids = sub("\\.csv$", "", basename(files))
    for(i in seq_along(files)) {
      file_manifest[[length(file_manifest)+1]] = list(
        path = files[i], 
        batch_id = batch_ids[i], 
        folder = basename(folder)
      )
    }
  }
}

message(sprintf("Found %d files to process.", length(file_manifest)))

# 2. Loop and Process (One by One)
processed_list = list()

for (i in seq_along(file_manifest)) {
  res = process_single_file(file_manifest[[i]], clinical_columns, gene_map_dt)
  if (!is.null(res)) {
    processed_list[[length(processed_list)+1]] = res
  }
  # CRITICAL: Clean up memory after every iteration
  gc() 
}

if (length(processed_list) == 0) stop("No data loaded.")

# -----------------------------------------------------------------------------
# COMBINE RESULTS
# -----------------------------------------------------------------------------
message("Combining Clinical Data...")
all_clinical = do.call(bind_rows, lapply(processed_list, `[[`, "clinical"))
rownames(all_clinical) = all_clinical$Sample_ID

message("Unifying Gene Lists...")
all_matrices = lapply(processed_list, `[[`, "expr")
all_genes = Reduce(union, lapply(all_matrices, rownames))
message(sprintf("Total Unique Genes: %d", length(all_genes)))

# -----------------------------------------------------------------------------
# MATRIX MERGE
# -----------------------------------------------------------------------------
message("Merging Expression Matrices...")

total_samples = nrow(all_clinical)
final_matrix = matrix(NA, nrow = length(all_genes), ncol = total_samples)
rownames(final_matrix) = all_genes
colnames(final_matrix) = rownames(all_clinical)

for (mat in all_matrices) {
  sample_names = colnames(mat)
  common_genes = intersect(rownames(mat), all_genes)
  col_indices = match(sample_names, colnames(final_matrix))
  final_matrix[common_genes, col_indices] = mat[common_genes, ]
  rm(mat) 
}

rm(all_matrices, processed_list)
gc()

# -----------------------------------------------------------------------------
# STANDARDIZATION
# -----------------------------------------------------------------------------
message("Standardizing...")

possible_disease_cols = setdiff(names(all_clinical), 
                                c("Age", "Gender", "batch_id", "source_folder", "Sample_ID"))

all_clinical = all_clinical %>%
  mutate(Disease_Status_Raw = coalesce(!!!syms(possible_disease_cols))) %>%
  mutate(
    bio_group = case_when(
      grepl("0|control|normal|healthy", Disease_Status_Raw, ignore.case = TRUE) ~ "Control",
      grepl("1|alzheimer|ad|dementia|case|tumor|cancer|leukemia|positive", Disease_Status_Raw, ignore.case = TRUE) ~ "Case",
      TRUE ~ NA_character_
    )
  )

valid_samples = !is.na(all_clinical$bio_group)
all_clinical = all_clinical[valid_samples, ]
final_matrix = final_matrix[, valid_samples]

# -----------------------------------------------------------------------------
# MEMORY-SURGICAL SAVE (Writes 1 CSV)
# -----------------------------------------------------------------------------
message("--- Starting Surgical Merge & Write ---")

# 1. Transpose the Matrix
message("Transposing matrix... (Hold on)")
t_expr = t(final_matrix)

rm(final_matrix)
gc() 

# 2. Convert to data.table
message("Converting to data.table...")
final_dt = as.data.table(t_expr)

rm(t_expr)
gc()

if (!all(rownames(all_clinical) == rownames(all_clinical))) {
  stop("Row mismatch! Do not proceed.")
}

message("Injecting clinical data...")
final_dt[, Sample_ID := rownames(all_clinical)]
final_dt[, Bio_Group := all_clinical$bio_group]
final_dt[, Age := all_clinical$Age]
final_dt[, Gender := all_clinical$Gender]
# You can uncomment below if you want Batch ID in the final CSV
# final_dt[, Batch := all_clinical$batch_id] 

setcolorder(final_dt, c("Sample_ID", "Bio_Group", "Age", "Gender"))

message("Writing massive CSV to disk...")
fwrite(final_dt, "Final_Combined_Dataset.csv")

message("--- SUCCESS! Data saved to Final_Combined_Dataset.csv ---")