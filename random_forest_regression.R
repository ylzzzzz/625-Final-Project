library(randomForest)
library(caTools)
library(caret)

model_data = final_modeling_df_fully_cleaned

# Convert target variable to a factor
model_data$Alzheimers_Disease = as.factor(model_data$Alzheimers_Disease)
levels(model_data$Alzheimers_Disease) = make.names(levels(model_data$Alzheimers_Disease))
print(paste("Target variable 'Alzheimers_Disease' converted to factor with levels:", 
            paste(levels(model_data$Alzheimers_Disease), collapse = ", ")))


# Remove columns that are not predictors.
cols_to_remove = c("Sample_ID", "batch_id", "bio_group")
model_data = model_data[, !(names(model_data) %in% cols_to_remove)]

#print("Data dimensions after removing non-predictor columns:")
#print(dim(model_data))
#print("Column names of final modeling data:")
#print(names(model_data))


# Split Data into Training (70%) and Test Sets (30%)
split = sample.split(model_data$Alzheimers_Disease, SplitRatio = 0.70)

train_set = subset(model_data, split == TRUE)
test_set = subset(model_data, split == FALSE)

#print(paste("Training set size:", nrow(train_set)))
#print(paste("Test set size:", nrow(test_set)))

y_train = train_set$Alzheimers_Disease
x_train = train_set[, !(names(train_set) %in% c("Alzheimers_Disease"))]
y_test = test_set$Alzheimers_Disease
x_test = test_set[, !(names(test_set) %in% c("Alzheimers_Disease"))]

#print("Created x/y train/test splits to handle wide data.")
#print(paste("x_train dimensions:", paste(dim(x_train), collapse = " x ")))
#print(paste("y_train length:", length(y_train)))


# Train the Random Forest Model
rf_model = randomForest(
  x = x_train,
  y = y_train,
  ntree = 500,
  importance = TRUE,
  na.action = na.omit
)
print("Model training complete.")

# Print the model summary
print(rf_model)

# Plot variable importance (which predictors the model found most useful for predicting Alzheimer's)
varImpPlot(rf_model, main = "Variable Importance Plot")


# Predict and Evaluate the Model
predictions = predict(rf_model, newdata = x_test)

positive_class = levels(y_train)[2]
print(paste("Using", positive_class, "as the 'positive' class for evaluation."))

conf_matrix = confusionMatrix(
  data = predictions,
  reference = y_test,
  positive = positive_class
)

print("--- Model Evaluation on Test Set ---")
print(conf_matrix)

