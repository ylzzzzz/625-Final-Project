library(ggplot2)
library(dplyr)
library(tidyr)
library(bench)
library(data.table)
library(vroom)

# --- [Setup Phase - same as your code] ---
csv_path <- "./data/Alzheimers_Disease_final.csv"
zip_path <- "./data/Alzheimers_Disease_final.csv.7z"

seven_zip_exe = Sys.which("7z")
if (!nzchar(seven_zip_exe)) {
  default_path = "C:/Program Files/7-Zip/7z.exe"
  if (file.exists(default_path)) {
    seven_zip_exe = default_path
  }
}

# --- [Run Benchmark] ---
results <- bench::mark(
  "Base R" = read.csv(csv_path),
  "data.table" = fread(csv_path),
  "vroom" = vroom(csv_path, show_col_types = FALSE),
  "7-Zip" = fread(cmd = paste(shQuote(seven_zip_exe), "e -so ", zip_path)),
  check = FALSE,
  iterations = 1,     # Note: With 1 iteration, you won't get error bars/distributions
  filter_gc = FALSE
)

# --- [Plot 1: The "Plain" Plot (Built-in)] ---
# bench::autoplot knows how to handle the gc column automatically.
# It uses a violin/beeswarm plot by default.
p_plain <- autoplot(results) +
  labs(title = "bench::autoplot() Results")

print(p_plain)

# --- [Plot 2: The Custom ggplot] ---
# 1. Convert to tibble to drop 'bench_mark' class (avoids the "Assigned data" error)
# 2. Select only relevant columns (drops 'gc')
# 3. Unnest the time column
plot_data <- results %>%
  as_tibble() %>%
  select(expression, time) %>%
  unnest(cols = c(time)) %>%
  mutate(expression = factor(expression, levels = c("data.table", "vroom", "7-Zip", "Base R")))

# Create the custom plot
# Note: Since iterations=1, we use geom_col (bar) instead of geom_boxplot
p_custom <- ggplot(plot_data, aes(x = expression, y = as.numeric(time), fill = expression)) +
  geom_col(alpha = 0.7) +
  coord_flip() + # Horizontal bars are often easier to read for benchmarks
  labs(
    y = "Time (seconds)", 
    x = "Method",
    title = "Custom Benchmark: Load Times (1 Iteration)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

print(p_custom)