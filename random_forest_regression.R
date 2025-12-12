library(randomForest)
library(caTools)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(limma) # Required for removeBatchEffect

# Load Data
model_data = fread("./data/Alzheimers_Disease_cleandata_final.csv")

# Convert target variable to a factor immediately for the design matrix
model_data$Alzheimers_Disease = as.factor(model_data$Alzheimers_Disease)
levels(model_data$Alzheimers_Disease) = make.names(levels(model_data$Alzheimers_Disease))


# Perform Batch Correction (The TODO Step)
print("--- Starting Batch Correction ---")

# Identify metadata columns to exclude from the expression matrix
metadata_cols = c("Sample_ID", "batch_id", "bio_group", "Alzheimers_Disease")

# Extract the Gene Expression Matrix (Rows = Genes, Cols = Samples for limma)
# We subset the data to only gene columns, then transpose it
gene_expression = as.matrix(model_data[, !metadata_cols, with = FALSE])
gene_expression_t = t(gene_expression) 

# Create the Design Matrix (preserving the biological condition of interest)
design_mat = model.matrix(~ model_data$Alzheimers_Disease)

# Apply removeBatchEffect
# We protect the "Alzheimers_Disease" signal using the design matrix
batch_corrected_mat = removeBatchEffect(
  x = gene_expression_t,
  batch = model_data$batch_id,
  design = design_mat
)

# Transpose back to (Rows = Samples, Cols = Genes) for Machine Learning
batch_corrected_df = as.data.frame(t(batch_corrected_mat))

# Re-attach the target variable
batch_corrected_df$Alzheimers_Disease = model_data$Alzheimers_Disease

print("Batch correction complete. Replaced raw expression data with corrected residuals.")

# Overwrite model_data with the corrected version for the rest of the pipeline
model_data = batch_corrected_df


# Standard ML Setup (Splitting)
# Split Data into Training (70%) and Test Sets (30%)
set.seed(123) # Good practice to set seed for reproducibility
split = sample.split(model_data$Alzheimers_Disease, SplitRatio = 0.70)

train_set = subset(model_data, split == TRUE)
test_set = subset(model_data, split == FALSE)

y_train = train_set$Alzheimers_Disease
x_train = train_set[, !(names(train_set) %in% c("Alzheimers_Disease"))]
y_test = test_set$Alzheimers_Disease
x_test = test_set[, !(names(test_set) %in% c("Alzheimers_Disease"))]


# Train the Random Forest Model
print("--- Training Random Forest on Batch-Corrected Data ---")
time_taken = system.time({
  rf_model = randomForest(
    x = x_train,
    y = y_train,
    ntree = 200,
    importance = TRUE,
    na.action = na.omit
  )
})

print("Model training complete.")
print(time_taken)
print(rf_model)


# Variable Importance
# Save the base plot
png(filename = "RF_BatchCorrected_VarImp.png", width = 800, height = 600)
varImpPlot(rf_model, main = "Variable Importance (Batch Corrected)")
dev.off()

# Extract importance for ggplot
imp_df = as.data.frame(importance(rf_model))
imp_df$Variable = rownames(imp_df)

# Top 20 Plot
top_20_imp = imp_df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(20)

p_imp = ggplot(top_20_imp, aes(x = reorder(Variable, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 20 Variable Importance (Batch Corrected)",
    subtitle = "Predictors of Alzheimer's Disease",
    x = "Predictor",
    y = "Mean Decrease in Accuracy"
  )

ggsave("RF_BatchCorrected_Feature_Importance.png", plot = p_imp, width = 7, height = 5)


# Predictions & Evaluation
predictions = predict(rf_model, newdata = x_test)
positive_class = levels(y_train)[2] # Typically "X1" or "Yes"

conf_matrix = confusionMatrix(
  data = predictions,
  reference = y_test,
  positive = positive_class
)

print("--- Model Evaluation on Test Set ---")
print(conf_matrix)

# Overfitting Check
oob_error = rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
oob_accuracy = 1 - oob_error
test_accuracy = conf_matrix$overall["Accuracy"]
diff = oob_accuracy - test_accuracy

print(paste("Training (OOB) Accuracy:", round(oob_accuracy, 4)))
print(paste("Test Set Accuracy:      ", round(test_accuracy, 4)))
print(paste("Difference:             ", round(diff, 4)))

if(diff > 0.05) {
  print("WARNING: Possible Overfitting detected (Gap > 5%)")
} else {
  print("Model generalizes well (Gap < 5%)")
}


# ROC Curve
pred_probs = predict(rf_model, newdata = x_test, type = "prob")
# Dynamic selection of positive class column (usually the 2nd column)
positive_class_col = colnames(pred_probs)[2] 
positive_class_prob = pred_probs[, positive_class_col]

roc_obj = roc(y_test, positive_class_prob)

p_roc = ggroc(roc_obj, color = "darkgreen", size = 1.2) + # Changed color to distinguish from non-corrected
  theme_minimal() +
  labs(
    title = paste("ROC Curve (Batch Corrected) | AUC =", round(auc(roc_obj), 3)),
    x = "Specificity (1 - False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey")

ggsave("RF_BatchCorrected_ROC_Curve.png", plot = p_roc, width = 6, height = 4)


# Error Rate Plot
err_df = as.data.frame(rf_model$err.rate)
err_df$Trees = 1:nrow(err_df)
err_long = pivot_longer(err_df, cols = -Trees, names_to = "ErrorType", values_to = "Error")

p_err = ggplot(err_long, aes(x = Trees, y = Error, color = ErrorType)) +
  geom_line(size = 0.8) +
  theme_bw() +
  labs(
    title = "OOB and Class Error Rate vs. Trees (Batch Corrected)",
    y = "Error Rate",
    color = "Metric"
  )

ggsave("RF_BatchCorrected_Error_Rate.png", plot = p_err, width = 7, height = 4)


# RMSE Calculation
calc_class_rmse = function(model, data, actual_outcomes, positive_label) {
  probs = predict(model, newdata = data, type = "prob")[, positive_label]
  actual_numeric = as.numeric(actual_outcomes == positive_label)
  sqrt(mean((actual_numeric - probs)^2))
}

train_rmse = calc_class_rmse(rf_model, x_train, y_train, positive_class)
test_rmse  = calc_class_rmse(rf_model, x_test, y_test, positive_class)

print(paste("Training (Probability) RMSE:", round(train_rmse, 4)))
print(paste("Test Set (Probability) RMSE:       ", round(test_rmse, 4)))
