library(keras3)
library(tensorflow)
library(caTools)
library(caret)
library(DiagrammeR)
library(ggplot2)
library(tidyr)
library(scales)

# INSTRUCTIONS!!!
# install backend
# tensorflow::install_tensorflow() or keras::install_keras()
# Run Alzheimers_data_cleaning.R and utils.R first to ensure 'final_modeling_df_fully_cleaned' dataframe exists in the environment
model_data = final_modeling_df_fully_cleaned

# Get original factor levels
factor_target = as.factor(model_data$Alzheimers_Disease)
original_levels = levels(factor_target)

# Check if we have exactly two levels (for binary classification)
if (length(original_levels) != 2) {
  stop("This script is designed for binary classification (2 levels), but the target has ", length(original_levels), " levels.")
}

# Get the level names
negative_class_name = make.names(original_levels[1]) # e.g., "X0" or "Healthy"
positive_class_name = make.names(original_levels[2]) # e.g., "X1" or "Case"

print(paste("Target variable levels:", negative_class_name, "(will be 0) and", positive_class_name, "(will be 1)"))

# Convert target to numeric 0 (for the first level) and 1 (for the second level)
model_data$Alzheimers_Disease = as.numeric(factor_target) - 1


# Predictor Preprocessing
# Remove columns that are not predictors
cols_to_remove = c("Sample_ID", "batch_id", "bio_group")
model_data = model_data[, !(names(model_data) %in% cols_to_remove)]

print("Data dimensions after removing non-predictor columns:")
print(dim(model_data))


# Data Splitting
# Split Data into Training (80%) and Test Sets (20%)
# We can use the numeric 0/1 target for splitting
split = sample.split(model_data$Alzheimers_Disease, SplitRatio = 0.80)

train_set = subset(model_data, split == TRUE)
test_set = subset(model_data, split == FALSE)

print(paste("Training set size:", nrow(train_set)))
print(paste("Test set size:", nrow(test_set)))

# Create x/y splits
y_train = train_set$Alzheimers_Disease
x_train = train_set[, !(names(train_set) %in% c("Alzheimers_Disease"))]
y_test = test_set$Alzheimers_Disease
x_test = test_set[, !(names(test_set) %in% c("Alzheimers_Disease"))]

print("Created x/y train/test splits.")
print(paste("x_train dimensions:", paste(dim(x_train), collapse = " x ")))


# Preprocessing for Keras/TensorFlow
# Standardization
# Convert data frames to matrices.

# Calculate scaling factors (mean, sd) from the TRAINING data only
x_train_mean = apply(x_train, 2, mean, na.rm = TRUE)
x_train_sd = apply(x_train, 2, sd, na.rm = TRUE)

# Handle columns with zero variance (sd=0) to avoid NaN
x_train_sd[x_train_sd == 0] = 1

# Apply scaling (standardization) to both training and test sets
x_train_scaled = sweep(x_train, 2, x_train_mean, "-")
x_train_scaled = sweep(x_train_scaled, 2, x_train_sd, "/")

x_test_scaled = sweep(x_test, 2, x_train_mean, "-")
x_test_scaled = sweep(x_test_scaled, 2, x_train_sd, "/")

# Handle any potential NaNs that might have resulted
x_train_scaled[is.na(x_train_scaled)] = 0
x_test_scaled[is.na(x_test_scaled)] = 0

print("x_train and x_test have been standardized.")

# Convert to matrix format
x_train_matrix = as.matrix(x_train_scaled)
y_train_matrix = as.matrix(y_train)

x_test_matrix = as.matrix(x_test_scaled)
y_test_matrix = as.matrix(y_test)

# Get number of features for the model's input layer
n_features = ncol(x_train_matrix)
print(paste("Data converted to matrix. Ready for Keras with", n_features, "features."))

# Build model
inputs = layer_input(shape = c(n_features))

outputs = inputs %>%
  layer_dense(
    units = 64, 
    activation = 'relu', 
    kernel_regularizer = regularizer_l1_l2(l1 = 0.001, l2 = 0.001)
  ) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation = 'sigmoid')

model = keras_model(inputs = inputs, outputs = outputs)

print("Keras model defined using Functional API.")
summary(model)

# Compile the Model
model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam(learning_rate = 0.0001),
  metrics = c('accuracy')
)

print("Keras model compiled.")


# Train the Model
print("Starting model training...")
weight_list = list("0" = 1, "1" = 4)

history = model %>% fit(
  x = x_train_matrix,
  y = y_train_matrix,
  class_weight = weight_list,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2,
  view_metrics = FALSE
)

print("Model training complete.")
print(plot(history))


# Predict and Evaluate the Model
pred_probs = model %>% predict(x_test_matrix)

# Convert probabilities to class labels (0 or 1) with 0.5 category cutoff
pred_classes_numeric = ifelse(pred_probs > 0.5, 1, 0)

# Convert numeric predictions (0, 1) and test labels (0, 1) back to factors
pred_factor = factor(
  pred_classes_numeric,
  levels = c(0, 1),
  labels = c(negative_class_name, positive_class_name)
)

y_test_factor = factor(
  y_test, # y_test is already the numeric 0/1 vector
  levels = c(0, 1),
  labels = c(negative_class_name, positive_class_name)
)

# Use the same 'positive' class name
print(paste("Using", positive_class_name, "as the 'positive' class for evaluation."))

# Generate the confusion matrix
conf_matrix = confusionMatrix(
  data = pred_factor,
  reference = y_test_factor,
  positive = positive_class_name
)

print("--- Model Evaluation on Test Set (TensorFlow/Keras) ---")
print(conf_matrix)







################################################################################
# Model Architecture plot

# Nodes for Input, Hidden Layers, Dropouts, and Output
viz_code = "
digraph neural_network {
  
  # Graph direction: Left to Right
  graph [rankdir = LR, splines=polyline]
  
  # Node Styles
  node [shape = box, style = filled, fontname = Helvetica, fontsize=10]
  
  # Define Nodes (Layers) with colors
  # Input Layer
  input [label = 'Input Layer\n(n_features)', fillcolor = '#E1F5FE']
  
  # Hidden Layer 1 (Dense)
  h1 [label = 'Dense Layer 1\n64 Units\nReLU + L1/L2 Reg', fillcolor = '#FFF9C4']
  
  # Dropout 1
  d1 [label = 'Dropout\n40%', fillcolor = '#FFCCBC']
  
  # Hidden Layer 2 (Dense)
  h2 [label = 'Dense Layer 2\n32 Units\nReLU', fillcolor = '#FFF9C4']
  
  # Dropout 2
  d2 [label = 'Dropout\n30%', fillcolor = '#FFCCBC']
  
  # Output Layer
  out [label = 'Output Layer\n1 Unit\nSigmoid', fillcolor = '#C8E6C9']
  
  # Edges (Connections)
  input -> h1
  h1 -> d1
  d1 -> h2
  h2 -> d2
  d2 -> out
}
"
grViz(viz_code)
################################################################################
# Training history plot

# Convert Keras history to a Data Frame
history_df = as.data.frame(history)

# Reshape for ggplot
history_long = history_df %>%
  pivot_longer(
    cols = -epoch, 
    names_to = "metric", 
    values_to = "value"
  ) %>%
  mutate(
    data_type = ifelse(grepl("val_", metric), "Validation", "Training"),
    metric_type = ifelse(grepl("loss", metric), "Loss (Binary Crossentropy)", "Accuracy")
  )

# Plot
ggplot(history_long, aes(x = epoch, y = value, color = data_type)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~metric_type, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Training" = "#1f77b4", "Validation" = "#ff7f0e")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Neural Network Training Performance",
    subtitle = "Monitoring Loss and Accuracy over Epochs",
    x = "Epoch",
    y = "Value",
    color = "Data Split"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
################################################################################
# Confusion matrix heatmap

# Extract the table from existing conf_matrix object
cm_table = as.data.frame(conf_matrix$table)

# Plot Heatmap
ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 8, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#a8dadc", high = "#1d3557") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Confusion Matrix",
    subtitle = paste("Model Accuracy:", round(conf_matrix$overall['Accuracy'], 3)),
    x = "Actual Label (Truth)",
    y = "Predicted Label"
  ) +
  coord_fixed()
################################################################################
# Significant calculations
# Number of Genes (Features) Used
n_genes = ncol(x_train_matrix)

# Training Time
timing = system.time({
  history = model %>% fit(
    x = x_train_matrix,
    y = y_train_matrix,
    class_weight = list("0" = 1, "1" = 4),
    epochs = 50,
    batch_size = 32,
    validation_split = 0.2,
    verbose = 0 # Turn off text log for cleaner timing
  )
})
train_time_seconds = timing[["elapsed"]]

# Prediction Accuracy
train_acc = tail(as.numeric(history$metrics$accuracy), 1)

# Testing Accuracy (From Confusion Matrix)
test_acc = as.numeric(conf_matrix$overall["Accuracy"])

# RMSE (Root Mean Squared Error)
# Training Data
train_probs = predict(model, x_train_matrix)
rmse_train = sqrt(mean((train_probs - y_train_matrix)^2))

# Testing Data
test_probs = predict(model, x_test_matrix)
rmse_test = sqrt(mean((test_probs - y_test_matrix)^2))


# PRINT REPORT
cat("================================================\n")
cat("          MODEL PERFORMANCE REPORT              \n")
cat("================================================\n")
cat(sprintf("1. Genes Selected (Features): %d\n", n_genes))
cat(sprintf("2. Training Time:             %.2f seconds\n", train_time_seconds))
cat("------------------------------------------------\n")
cat("3. Prediction Accuracy:\n")
cat(sprintf("   - Training Set:            %.2f%%\n", train_acc * 100))
cat(sprintf("   - Testing Set:             %.2f%%\n", test_acc * 100))
cat("------------------------------------------------\n")
cat("4. RMSE (Root Mean Squared Error):\n")
cat(sprintf("   - Training Set:            %.4f\n", rmse_train))
cat(sprintf("   - Testing Set:             %.4f\n", rmse_test))
cat("================================================\n")