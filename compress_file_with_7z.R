compress_file_with_7z = function(df, file_name) {
  
  # 1. Ensure file_name ends in .csv
  if (!grepl("\\.csv$", file_name, ignore.case = TRUE)) {
    csv_path = paste0(file_name, ".csv")
  } else {
    csv_path = file_name
  }
  
  # 2. Write the dataframe to CSV first
  message(sprintf("Writing CSV to: %s", csv_path))
  write.csv(df, csv_path, row.names = FALSE)
  
  # 3. Find 7-Zip Executable
  # First, check if it is in the System PATH
  seven_zip_exe = Sys.which("7z")
  
  # If not found in PATH, check the default Windows install location
  if (!nzchar(seven_zip_exe)) {
    default_path = "C:/Program Files/7-Zip/7z.exe"
    if (file.exists(default_path)) {
      seven_zip_exe = default_path
    }
  }
  
  # 4. Run Compression
  if (nzchar(seven_zip_exe)) {
    # Define archive name
    archive_path = paste0(csv_path, ".7z")
    message(sprintf("Found 7-Zip at: %s", seven_zip_exe))
    message(sprintf("Compressing to: %s ...", archive_path))
    
    # Run 7z command
    # We use shQuote to handle spaces in "Program Files" safely
    system2(command = seven_zip_exe, 
            args = c("a", "-mx=9", shQuote(archive_path), shQuote(csv_path)))
    
    if (file.exists(archive_path)) {
      message("Compression successful.")
    } else {
      warning("7-Zip command ran but archive was not created.")
    }
    
  } else {
    warning("Could not find 7z.exe in PATH or 'C:/Program Files/7-Zip/'. Saved as CSV only.")
  }
}

# --- Usage ---
#compress_file_with_7z(final_clean_df, "./data/Alzheimers_Disease_final")
