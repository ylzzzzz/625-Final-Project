# ==============================================================================
# Project Setup: Dependencies and Environment
# ==============================================================================

# CRAN Packages (Standard R packages)
cran_packages <- c(
  "jsonlite", "dplyr", "tidyr", "data.table", "scales",       # Data Manipulation
  "ggplot2", "pheatmap", "RColorBrewer", "DiagrammeR",        # Visualization
  "randomForest", "caTools", "caret", "pROC", "uwot",         # ML & Stats
  "keras3", "tensorflow"                                      # Deep Learning
)

# Bioconductor Packages (Genomics/Bioinformatics)
bioc_packages <- c(
  "sva", 
  "limma", 
  "EnhancedVolcano"
)

# Install BiocManager (if missing)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("Installing BiocManager...")
  install.packages("BiocManager")
}

# Install Missing CRAN Packages
installed_cran <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed_cran)

if (length(missing_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

# Install Missing Bioconductor Packages
installed_bioc <- rownames(installed.packages())
missing_bioc <- setdiff(bioc_packages, installed_bioc)

if (length(missing_bioc) > 0) {
  message("Installing missing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE) # update=FALSE prevents long prompts
}

# Load All Libraries
all_packages <- c(cran_packages, bioc_packages)

message("Loading libraries...")
sapply(all_packages, require, character.only = TRUE)


# ==============================================================================
# INSTRUCTIONS & BACKEND SETUP
# ==============================================================================

# NOTE: Python Backend Setup
# Since you are using keras3 and tensorflow, you must ensure the Python backend 
# is configured. If this is your first time running this on this machine, 
# uncomment and run the line below once:

# keras3::install_keras() 

# NOTE: Data Dependency
# Ensure you have run the following scripts before proceeding with modeling:
# 1. Alzheimers_data_cleaning.R
# 2. utils.R
#
# Check: exists("final_modeling_df_fully_cleaned") should be TRUE