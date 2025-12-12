library(data.table)
library(vroom)
library(archive)
library(bench)

# Find 7z
seven_zip_exe = Sys.which("7z")

# If not found in PATH, check the default Windows install location
if (!nzchar(seven_zip_exe)) {
  default_path = "C:/Program Files/7-Zip/7z.exe"
  if (file.exists(default_path)) {
    seven_zip_exe = default_path
  }
}

# Define the file paths
csv_path <- "./data/Alzheimers_Disease_final.csv"
zip_path <- "./data/Alzheimers_Disease_final.csv.7z"

# Run the benchmark
results <- bench::mark(
  # 1. Base R
  read_csv = read.csv(csv_path),
  
  # 2. data.table
  fread = fread(csv_path),
  
  # 3. vroom
  vroom = vroom(csv_path, show_col_types = FALSE),
  
  # 4. zip read
  fread_zip <- fread(cmd = paste(shQuote(seven_zip_exe), "e -so ", zip_path)),
  
  # Options
  check = FALSE,
  iterations = 5,
  filter_gc = FALSE
)

print(results)
plot(results)
