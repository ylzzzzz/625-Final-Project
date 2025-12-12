library(data.table)
library(vroom)
library(archive)
library(bench)
library(readr)
library(ggplot2)

# Define the file paths
csv_path <- "./data/Alzheimers_Disease_final.csv"
zip_path <- "./data/Alzheimers_Disease_final.csv.7z"

# Find 7z
seven_zip_exe = Sys.which("7z")

# If not found in PATH, check the default Windows install location
if (!nzchar(seven_zip_exe)) {
  default_path = "C:/Program Files/7-Zip/7z.exe"
  if (file.exists(default_path)) {
    seven_zip_exe = default_path
  }
}

# Run the benchmark
results <- bench::mark(
  "Base R" = read.csv(csv_path),
  
  "data.table" = fread(csv_path),
  
  "vroom" = vroom(csv_path, show_col_types = FALSE),
  
  "7-Zip" = fread(cmd = paste(shQuote(seven_zip_exe), "e -so ", zip_path)),
  
  # Options
  check = FALSE,
  iterations = 5,
  filter_gc = FALSE
)

# Re-order for specific sorting
results$expression <- factor(results$expression, levels = c("data.table", "vroom", "7-Zip", "Base R"))

# Plot results
autoplot(results) +
  labs(title = "Benchmark: CSV Reading Methods", y = "Method") +
  theme_minimal()
