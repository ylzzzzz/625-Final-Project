# ==============================================================================
# CONFIGURATION & SETUP
# ==============================================================================
cran_packages = c("jsonlite", "dplyr", "tidyr", "data.table")
bioc_packages = c("sva")

# Efficient package loading
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Parameters
MIN_SAMPLE_SIZE = 51           # Only merge if sample size > 50
GENE_MAP_PATH   = "gene_synonym.json"
BASE_PATH       = "data/Alzheimers_Disease"
USE_7ZIP        = TRUE         # Set to TRUE to compress output with 7zip

# Define clinical columns (Keep this list accurate!)
CLINICAL_COLS = c("Alzheimers_Disease", "Age", "Gender", "batch_id", "Sample_ID", "bio_group")

# ==============================================================================
# DATA LOADING & PRE-PROCESSING
# ==============================================================================

# 1. Load Gene Map
gene_map_list = fromJSON(GENE_MAP_PATH)
gene_map_dt = data.table(
  alias = names(gene_map_list),
  canonical = unlist(gene_map_list),
  stringsAsFactors = FALSE
)
setkey(gene_map_dt, alias)

# 2. Identify Files
combined_files = list.files(BASE_PATH, pattern = "\\.csv$", full.names = FALSE)
combined_files = combined_files[!grepl("/", combined_files)] # Exclude subfolders

if (length(combined_files) == 0) stop(paste("No .csv files found in:", BASE_PATH))
batch_ids = sub("\\.csv$", "", combined_files)

# 3. Load Function
load_csv = function(batch_id, base_path, clinical_list) {
  file_path = file.path(base_path, paste0(batch_id, ".csv"))
  
  if (!file.exists(file_path)) return(NULL)
  
  dt = fread(file_path)
  
  if (any(duplicated(dt[[1]]))) {
    warning(sprintf("SKIPPING %s: Duplicate sample IDs.", batch_id))
    return(NULL)
  }
  
  # Set rownames manually
  sample_ids = dt[[1]]
  data_df = as.data.frame(dt[, -1, with = FALSE])
  rownames(data_df) = sample_ids
  
  # Identify columns
  all_cols = colnames(data_df)
  clinical_found = intersect(all_cols, clinical_list)
  gene_cols = setdiff(all_cols, clinical_list)
  
  if (length(gene_cols) == 0 || length(clinical_found) == 0) {
    warning(sprintf("SKIPPING %s: Missing gene or clinical data.", batch_id))
    return(NULL)
  }
  
  # --- Filter: Sample Size > 50 ---
  if (nrow(data_df) <= 50) {
    message(sprintf("SKIPPING %s: Sample size %d is <= 50.", batch_id, nrow(data_df)))
    return(NULL)
  }
  
  # Split Data
  clinical_data = data_df[, clinical_found, drop = FALSE]
  clinical_data$batch_id = batch_id
  
  # Transpose Gene Data (Genes x Samples) for alignment
  expr_data = t(data_df[, gene_cols, drop = FALSE])
  
  message(sprintf("Loaded %s: %d samples.", batch_id, ncol(expr_data)))
  return(list(expr_data = expr_data, clinical_data = clinical_data))
}

# Load Datasets
message("--- Loading Datasets ---")
all_data_list = lapply(batch_ids, load_csv, base_path = BASE_PATH, clinical_list = CLINICAL_COLS)
all_data_list = all_data_list[!sapply(all_data_list, is.null)] 

if (length(all_data_list) == 0) stop("No datasets remained after filtering.")

# ==============================================================================
# ALIGNMENT & MERGING
# ==============================================================================

align_and_aggregate = function(expr_matrix, map_dt) {
  dt = as.data.table(expr_matrix, keep.rownames = "alias")
  dt[map_dt, on = "alias", gene_final := i.canonical]
  dt[is.na(gene_final), gene_final := alias]
  dt_agg = dt[, lapply(.SD, mean, na.rm = TRUE), by = gene_final, .SDcols = colnames(expr_matrix)]
  
  mat = as.matrix(dt_agg[, -1, with = FALSE])
  rownames(mat) = dt_agg$gene_final
  return(mat)
}

message("--- Aligning Genes ---")
aligned_expr_list = lapply(all_data_list, function(d) {
  align_and_aggregate(d$expr_data, gene_map_dt)
})

# Union of all genes
all_genes = Reduce(union, lapply(aligned_expr_list, rownames))
message(sprintf("Total unique genes aligned: %d", length(all_genes)))

# Re-index to Master Gene List (Fill NAs)
message("--- Re-indexing Matrices ---")
final_expr_list = lapply(aligned_expr_list, function(mat) {
  new_mat = matrix(NA, nrow = length(all_genes), ncol = ncol(mat), 
                   dimnames = list(all_genes, colnames(mat)))
  common = intersect(rownames(mat), all_genes)
  new_mat[common, ] = mat[common, ]
  return(new_mat)
})

# Combine
expr_combined = do.call(cbind, final_expr_list)
clinical_combined = bind_rows(lapply(all_data_list, `[[`, "clinical_data"))
expr_combined = expr_combined[, rownames(clinical_combined)]

# ==============================================================================
# CLINICAL STANDARDIZATION
# ==============================================================================
message("--- Standardizing Clinical Data ---")

target_col = "Alzheimers_Disease"
if (!target_col %in% names(clinical_combined)) target_col = names(clinical_combined)[1]

clinical_combined = clinical_combined %>%
  mutate(bio_group = case_when(
    grepl("0|control|normal|healthy", .data[[target_col]], ignore.case = TRUE) ~ "Control",
    grepl("1|alzheimer|ad|dementia", .data[[target_col]], ignore.case = TRUE) ~ "AD",
    TRUE ~ NA_character_
  ))

# Filter Clinical NAs
keep_rows = !is.na(clinical_combined$bio_group)
clinical_combined = clinical_combined[keep_rows, ]
expr_combined = expr_combined[, keep_rows]

# Merge final
final_df = merge(clinical_combined, t(expr_combined), by = "row.names")
final_df = final_df %>% rename(Sample_ID = Row.names)

# ==============================================================================
# NA IMPUTATION (Written Out Explicitly)
# ==============================================================================

# Dynamically separate gene data
clinical_cols_final = intersect(names(final_df), CLINICAL_COLS)
gene_cols_final = setdiff(names(final_df), CLINICAL_COLS)

gene_data = final_df[, gene_cols_final]

clean_NA = function(X, max_na = 0.5, method = c("median", "mean")) {
  method = match.arg(method)
  
  if (!is.data.frame(X)) X = as.data.frame(X)
  
  # 1. Remove columns with too many NAs
  na_prop = colMeans(is.na(X))
  too_many_na = na_prop > max_na
  
  if (any(too_many_na)) {
    cat(sprintf("Removed %d variables with > %.0f%% missing values.\n", 
                sum(too_many_na), 100 * max_na))
    X = X[, !too_many_na]
  }
  
  # 2. Impute remaining NAs (Explicit Loop)
  cat("Imputing remaining NA values using", method, "...\n")
  
  for (j in seq_len(ncol(X))) {
    # Only process if there are NAs in this specific column
    if (anyNA(X[[j]])) {
      
      # Calculate statistic based on method
      if (method == "median") {
        fill_value = median(X[[j]], na.rm = TRUE)
      } else if (method == "mean") {
        fill_value = mean(X[[j]], na.rm = TRUE)
      }
      
      # Replace NAs
      X[[j]][is.na(X[[j]])] <- fill_value
    }
  }
  
  # Ensure strictly numeric matrix output
  # (suppress warnings in case coercion happens, though data should be numeric)
  X_mat = suppressWarnings(as.matrix(sapply(X, as.numeric)))
  
  if (anyNA(X_mat)) {
    stop("There are still NA values after imputation!")
  } else {
    cat("All missing values handled successfully.\n")
    cat(sprintf("Final dimensions: %d samples x %d features\n", nrow(X_mat), ncol(X_mat)))
  }
  
  return(X_mat)
}

# Run cleaning
genedata_clean = clean_NA(gene_data, max_na = 0.5, method = "median")
# Clean clinical data:
# 1. Subset the actual data (fixing the immediate error)
clinical_data_subset = final_df[, clinical_cols_final, drop = FALSE]

# 2. Separate Numeric vs Categorical clinical data
# (clean_NA is distinctively designed for numeric matrices)
clin_num_vars = sapply(clinical_data_subset, is.numeric)
clin_cat_vars = !clin_num_vars

# 3. Clean only the numeric clinical data (e.g., Age)
if (any(clin_num_vars)) {
  # Subset numeric columns
  clin_num_data = clinical_data_subset[, clin_num_vars, drop = FALSE]
  
  # Run your custom function
  clin_num_clean = clean_NA(clin_num_data, max_na = 0.5, method = "median")
  
  # Convert back to data frame to merge later
  clin_num_clean = as.data.frame(clin_num_clean)
} else {
  clin_num_clean = data.frame(row.names = rownames(clinical_data_subset))
}

# 4. Handle categorical data (Optional: fill with Mode or specific value)
if (any(clin_cat_vars)) {
  clin_cat_data = clinical_data_subset[, clin_cat_vars, drop = FALSE]
  
  # Simple imputation for categorical: Replace NA with "Unknown" or Mode
  # Here we just ensure they aren't lost
  clin_cat_data[is.na(clin_cat_data)] <- "Unknown" 
} else {
  clin_cat_data = data.frame(row.names = rownames(clinical_data_subset))
}

# 5. Recombine
clinical_clean_final = cbind(clin_cat_data, clin_num_clean)

# Verify rows match before binding (Critical for clinical+gene data)
if (!all(rownames(clinical_clean_final) == rownames(genedata_clean))) {
  stop("Row names do not match! Sort them before binding.")
}

# Bind the two cleaned objects directly
final_clean_df = cbind(clinical_clean_final, as.data.frame(genedata_clean))

message(sprintf("Final dimensions: %d samples x %d features", 
                nrow(final_clean_df), ncol(final_clean_df)))

# ==============================================================================
# SAVING & 7ZIP
# ==============================================================================

out_csv = "./data/Alzheimers_Disease_cleandata_final.csv"

# Combine clean genes back with clinical data
final_clean_df = cbind(final_df[, clinical_clean_final], as.data.frame(genedata_clean))

message(sprintf("Writing CSV to %s ...", out_csv))
write.csv(final_clean_df, out_csv, row.names = TRUE)

if (USE_7ZIP) {
  # Check if 7z is installed/in PATH
  seven_zip_path = Sys.which("7z")
  
  if (nzchar(seven_zip_path)) {
    out_archive = paste0(out_csv, ".7z")
    message(sprintf("Compressing to %s ...", out_archive))
    
    # 7z a [archive_name] [file_to_compress]
    # -mx=9 is max compression
    system2("7z", args = c("a", "-mx=9", out_archive, out_csv))
    
    if (file.exists(out_archive)) {
      message("Compression successful.")
    } else {
      warning("7-Zip command ran but archive was not found.")
    }
  } else {
    warning("7-Zip executable ('7z') not found in system PATH. Skipping compression.")
  }
}

message("--- SCRIPT COMPLETE ---")
