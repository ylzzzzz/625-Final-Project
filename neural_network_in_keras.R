library(keras)
library(tensorflow)
library(caTools)
library(caret)

# --- Installation Note ---
# install the backend
# tensorflow::install_tensorflow() or keras::install_keras()

# !!!!!!!
# Run Alzheimers_data_cleaning.R and utils.R first to ensure 'final_modeling_df_fully_cleaned' exists in the environment
model_data <- final_modeling_df_fully_cleaned

# Target Variable Preprocessing
# Keras needs a numeric target (0/1) for binary classification.

# Get original factor levels
factor_target <- as.factor(model_data$Alzheimers_Disease)
original_levels <- levels(factor_target)

# Check if we have exactly two levels (for binary classification)
if (length(original_levels) != 2) {
  stop("This script is designed for binary classification (2 levels), but the target has ", length(original_levels), " levels.")
}

# Get the level names
negative_class_name <- make.names(original_levels[1]) # e.g., "X0" or "Healthy"
positive_class_name <- make.names(original_levels[2]) # e.g., "X1" or "Case"

print(paste("Target variable levels:", negative_class_name, "(will be 0) and", positive_class_name, "(will be 1)"))

# Convert target to numeric 0 (for the first level) and 1 (for the second level)
model_data$Alzheimers_Disease <- as.numeric(factor_target) - 1


# Predictor Preprocessing
# Remove columns that are not predictors
cols_to_remove <- c("Sample_ID", "batch_id", "bio_group")
model_data <- model_data[, !(names(model_data) %in% cols_to_remove)]

print("Data dimensions after removing non-predictor columns:")
print(dim(model_data))


# Data Splitting
# Split Data into Training (80%) and Test Sets (20%)
# We can use the numeric 0/1 target for splitting
split <- sample.split(model_data$Alzheimers_Disease, SplitRatio = 0.80)

train_set <- subset(model_data, split == TRUE)
test_set <- subset(model_data, split == FALSE)

print(paste("Training set size:", nrow(train_set)))
print(paste("Test set size:", nrow(test_set)))

# Create x/y splits (target 'Alzheimers_Disease' is now 0/1)
y_train <- train_set$Alzheimers_Disease
x_train <- train_set[, !(names(train_set) %in% c("Alzheimers_Disease"))]
y_test <- test_set$Alzheimers_Disease
x_test <- test_set[, !(names(test_set) %in% c("Alzheimers_Disease"))]

print("Created x/y train/test splits.")
print(paste("x_train dimensions:", paste(dim(x_train), collapse = " x ")))


# Preprocessing for Keras/TensorFlow
# Scale the data (e.g., standardization).
# Convert data frames to matrices.

# Calculate scaling factors (mean, sd) from the TRAINING data only
x_train_mean <- apply(x_train, 2, mean, na.rm = TRUE)
x_train_sd <- apply(x_train, 2, sd, na.rm = TRUE)

# Handle columns with zero variance (sd=0) to avoid NaN
x_train_sd[x_train_sd == 0] <- 1

# Apply scaling (standardization) to both training and test sets
x_train_scaled <- sweep(x_train, 2, x_train_mean, "-")
x_train_scaled <- sweep(x_train_scaled, 2, x_train_sd, "/")

x_test_scaled <- sweep(x_test, 2, x_train_mean, "-")
x_test_scaled <- sweep(x_test_scaled, 2, x_train_sd, "/")

# Handle any potential NaNs that might have resulted (e.g., from all-NA columns)
x_train_scaled[is.na(x_train_scaled)] <- 0
x_test_scaled[is.na(x_test_scaled)] <- 0

print("x_train and x_test have been standardized.")

# Convert to matrix format
x_train_matrix <- as.matrix(x_train_scaled)
y_train_matrix <- as.matrix(y_train)

x_test_matrix <- as.matrix(x_test_scaled)
y_test_matrix <- as.matrix(y_test)

# Get number of features for the model's input layer
n_features <- ncol(x_train_matrix)
print(paste("Data converted to matrix. Ready for Keras with", n_features, "features."))


# Define the Keras (TensorFlow) Model

# FIXME: Not training correctly.  Check backwards step.

model <- keras_model_sequential() %>%
  # Input layer: Must match the number of features
  layer_dense(units = 64, activation = 'relu', input_shape = c(n_features)) %>%
  layer_dropout(rate = 0.4) %>%
  layer_dense(units = 32, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation = 'sigmoid')

print("Keras model defined.")
summary(model)

# Compile the Model (the learning process) here
model %>% compile(
  # 'binary_crossentropy' is the standard loss function for 0/1 classification
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam(learning_rate = 0.001),
  metrics = c('accuracy')
)

print("Keras model compiled.")


# Train the Model
print("Starting model training...")
history <- model %>% fit(
  x = x_train_matrix,
  y = y_train_matrix,
  epochs = 50,           # Number of times to see the full training data
  batch_size = 32,       # Number of samples per gradient update
  validation_split = 0.2, # Use 20% of training data for validation
  verbose = 2            # Show one line per epoch
)

print("Model training complete.")
print(plot(history))


# Predict and Evaluate the Model
# Predict probabilities on the test set
pred_probs <- model %>% predict(x_test_matrix)

# Convert probabilities to class labels (0 or 1) with 0.5 category cutoff
pred_classes_numeric <- ifelse(pred_probs > 0.5, 1, 0)

# Convert numeric predictions (0, 1) and test labels (0, 1) back to factors
pred_factor <- factor(
  pred_classes_numeric,
  levels = c(0, 1),
  labels = c(negative_class_name, positive_class_name)
)

y_test_factor <- factor(
  y_test, # y_test is already the numeric 0/1 vector
  levels = c(0, 1),
  labels = c(negative_class_name, positive_class_name)
)

# Use the same 'positive' class name
print(paste("Using", positive_class_name, "as the 'positive' class for evaluation."))

# Generate the confusion matrix
conf_matrix <- confusionMatrix(
  data = pred_factor,
  reference = y_test_factor,
  positive = positive_class_name
)

print("--- Model Evaluation on Test Set (TensorFlow/Keras) ---")
print(conf_matrix)