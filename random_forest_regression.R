library(randomForest)
library(caTools)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

model_data = final_modeling_df_fully_cleaned

# Convert target variable to a factor
model_data$Alzheimers_Disease = as.factor(model_data$Alzheimers_Disease)
levels(model_data$Alzheimers_Disease) = make.names(levels(model_data$Alzheimers_Disease))
print(paste("Target variable 'Alzheimers_Disease' converted to factor with levels:", 
            paste(levels(model_data$Alzheimers_Disease), collapse = ", ")))


# Remove columns that are not predictors
cols_to_remove = c("Sample_ID", "batch_id", "bio_group")
model_data = model_data[, !(names(model_data) %in% cols_to_remove)]

# Split Data into Training (70%) and Test Sets (30%)
split = sample.split(model_data$Alzheimers_Disease, SplitRatio = 0.70)

train_set = subset(model_data, split == TRUE)
test_set = subset(model_data, split == FALSE)

y_train = train_set$Alzheimers_Disease
x_train = train_set[, !(names(train_set) %in% c("Alzheimers_Disease"))]
y_test = test_set$Alzheimers_Disease
x_test = test_set[, !(names(test_set) %in% c("Alzheimers_Disease"))]


# Train the Random Forest Model
time_taken = system.time({
  rf_model = randomForest(
    x = x_train,
    y = y_train,
    ntree = 1000,
    importance = TRUE,
    na.action = na.omit
  )
})

print("Model training complete.")
print(time_taken)

# Print the model summary
print(rf_model)

# Plot variable importance
varImpPlot(rf_model, main = "Variable Importance Plot")

png(filename = "RF_Base_VarImp.png", width = 800, height = 600)
varImpPlot(rf_model, main = "Variable Importance Plot")
dev.off()

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

# Get the final OOB Error rate from the model
oob_error = rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
oob_accuracy = 1 - oob_error

# Get the Test Accuracy from your Confusion Matrix
test_accuracy = conf_matrix$overall["Accuracy"]

# Compare them
diff = oob_accuracy - test_accuracy

print(paste("Training (OOB) Accuracy:", round(oob_accuracy, 4)))
print(paste("Test Set Accuracy:      ", round(test_accuracy, 4)))
print(paste("Difference:             ", round(diff, 4)))

if(diff > 0.05) {
  print("WARNING: Possible Overfitting detected (Gap > 5%)")
} else {
  print("Model generalizes well (Gap < 5%)")
}


################################################################################
# ROC Curve Plot


# Predict probabilities
pred_probs = predict(rf_model, newdata = x_test, type = "prob")

# Extract the probability of the "positive" class (e.g., "Alzheimers")
positive_class_prob = pred_probs[, "X1"]

# Create ROC object
roc_obj = roc(y_test, positive_class_prob)

# Plot using ggroc (part of pROC, works with ggplot layers)
p_roc = ggroc(roc_obj, color = "darkred", size = 1.2) +
  theme_minimal() +
  labs(
    title = paste("ROC Curve (AUC =", round(auc(roc_obj), 3), ")"),
    x = "Specificity (1 - False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey")

ggsave("RF_ROC_Curve.png", plot = p_roc, width = 6, height = 4)
################################################################################
# Feature importance plot
# Extract importance and convert to data frame
imp_df = as.data.frame(importance(rf_model))
imp_df$Variable = rownames(imp_df)

# Plot using ggplot 
# Reorder by MeanDecreaseAccuracy
p_imp = ggplot(imp_df, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for readable labels
  theme_minimal() +
  labs(
    title = "Variable Importance (Random Forest)",
    subtitle = "Predictors of Alzheimer's Disease",
    x = "Predictor",
    y = "Mean Decrease in Accuracy"
  )

ggsave("RF_Feature_Importance.png", plot = p_imp, width = 7, height = 5)

################################################################################
# Feature importance plot (Top 20 Only)

# Extract importance and convert to data frame
imp_df = as.data.frame(importance(rf_model))
imp_df$Variable = rownames(imp_df)

# Filter to just the Top 20 based on MeanDecreaseAccuracy
top_20_imp = imp_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(20)

# Plot using ggplot with the filtered data
p_imp = ggplot(top_20_imp, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for readable labels
  theme_minimal() +
  labs(
    title = "Top 20 Variable Importance (Random Forest)",
    subtitle = "Predictors of Alzheimer's Disease",
    x = "Predictor",
    y = "Mean Decrease in Accuracy"
  )

# Save the plot
ggsave("RF_Feature_Importance.png", plot = p_imp, width = 7, height = 5)
################################################################################
# Error plot
# Extract error rates
err_df = as.data.frame(rf_model$err.rate)
err_df$Trees = 1:nrow(err_df)

# Reshape for ggplot
err_long = pivot_longer(err_df, cols = -Trees, names_to = "ErrorType", values_to = "Error")

# Plot
p_err = ggplot(err_long, aes(x = Trees, y = Error, color = ErrorType)) +
  geom_line(size = 0.8) +
  theme_bw() +
  labs(
    title = "OOB and Class Error Rate vs. Number of Trees",
    y = "Error Rate",
    color = "Metric"
  ) +
  scale_color_manual(values = c("OOB" = "black", "Alzheimers_Disease.No" = "blue", "Alzheimers_Disease.Yes" = "red"))

ggsave("RF_Error_Rate.png", plot = p_err, width = 7, height = 4)
################################################################################
# Table statistics
# Number of Genes (Variables) Actually Selected
vars_used_count = varUsed(rf_model, count = TRUE)
num_genes_selected = sum(vars_used_count > 0)

print(paste("Total genes available:", length(vars_used_count)))
print(paste("Number of genes actually used in the forest:", num_genes_selected))


# Accuracy (Training vs Test)
# Training Accuracy
train_acc = 1 - rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]

# Test Accuracy
test_acc = conf_matrix$overall["Accuracy"]

print(paste("Training (OOB) Accuracy:", round(train_acc, 4)))
print(paste("Test Set Accuracy:      ", round(test_acc, 4)))


# RMSE (Root Mean Squared Error)
# Function to calculate RMSE for classification
calc_class_rmse = function(model, data, actual_outcomes, positive_label) {
  # Get probabilities for the positive class
  probs = predict(model, newdata = data, type = "prob")[, positive_label]
  
  # Convert actual class to 0 or 1 (1 = positive_label)
  # We use (actual == positive) to get TRUE/FALSE, then as.numeric to get 1/0
  actual_numeric = as.numeric(actual_outcomes == positive_label)
  
  # Calculate RMSE
  sqrt(mean((actual_numeric - probs)^2))
}

# This assumes the second level is the positive one
pos_class = levels(y_train)[2] 

# Calculate RMSE
train_rmse = calc_class_rmse(rf_model, x_train, y_train, pos_class)
test_rmse  = calc_class_rmse(rf_model, x_test, y_test, pos_class)

print(paste("Training (Probability) RMSE:", round(train_rmse, 4)))
print(paste("Test Set (Probability) RMSE:      ", round(test_rmse, 4)))