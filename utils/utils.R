remove_missing_columns <- function(df, threshold = 0.1, verbose = TRUE) {
  # drop columns with percent of total NA entries above threshold
  missing_proportions <- colMeans(is.na(df))
  columns_to_keep_mask <- missing_proportions <= threshold
  columns_dropped <- names(df)[!columns_to_keep_mask]
  df_cleaned <- df[, columns_to_keep_mask, drop = FALSE]
  if (verbose) {
    cat(paste("Columns with >", threshold * 100, "% missing values (Dropped): ", length(columns_dropped), "\n"))
    cat(paste("Columns Kept: ", sum(columns_to_keep_mask), "\n"))
  }
  return(df_cleaned)
}

remove_missing_rows <- function(df, threshold = 0.1, verbose = TRUE) {
  # drop rows with percent of NA entries above threshold
  missing_proportions <- rowMeans(is.na(df))
  rows_to_keep_mask <- missing_proportions <= threshold
  rows_dropped_count <- sum(!rows_to_keep_mask)
  df_cleaned <- df[rows_to_keep_mask, , drop = FALSE]
  if (verbose) {
    cat(paste("Rows with >", threshold * 100, "% missing values (Dropped): ", rows_dropped_count, "\n"))
    cat(paste("Rows Kept: ", sum(rows_to_keep_mask), "\n"))
  }
  return(df_cleaned)
}

final_modeling_df_fully_cleaned <- remove_missing_columns(final_modeling_df, threshold = 0.1)
final_modeling_df_fully_cleaned <- remove_missing_rows(final_modeling_df_fully_cleaned, threshold = 0.001)


# subset of final_modeling_df_fully_cleaned with only rows and columns that contain NA
rows_with_na <- apply(final_modeling_df_fully_cleaned, 1, anyNA)
cols_with_na <- apply(final_modeling_df_fully_cleaned, 2, anyNA)
na_subset <- final_modeling_df_fully_cleaned[rows_with_na, cols_with_na]
View(na_subset)
