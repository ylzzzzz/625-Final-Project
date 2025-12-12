start_time <- Sys.time()
library(data.table)
mydata <- fread("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata_final.csv")
dim(mydata)
end_time <- Sys.time()
time_used <- end_time - start_time
print(time_used)

genedata_clean <- fread("~/Documents/DNA/2025Fall/Biostat625/Final/Batch_corrected_expression.csv")
clinical_data <- mydata[, 1:6]
# train data and test data
n <- nrow(mydata)
train_index <- sample(seq_len(n), size = 0.7 * n)

X_train <- genedata_clean[train_index,-1 ]
X_test  <- genedata_clean[-train_index, -1]
y_train <- unlist(clinical_data[train_index, 2])
y_test <- unlist(clinical_data[-train_index, 2])

# 检查划分结果
nrow(train_data)
nrow(test_data)

# 查看缺失值
missing_values <- is.na(test_data)
any(missing_values)
missing_values <- is.na(train_data)
any(missing_values)


#------------------lasso----------------------------

library(glmnet)
set.seed(1)
grid <- 10^seq(10, -2, length = 100)
# fit
start_time <- Sys.time()

lasso_mod <- glmnet(as.matrix( X_train), y_train, alpha = 1, lambda = grid,family = "binomial")

end_time <- Sys.time()
time_used <- end_time - start_time
print(time_used)

plot(lasso_mod)
# select parameter
cv_out_lasso <- cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
plot(cv_out_lasso)
bestlam_lasso <- cv_out_lasso$lambda.min
bestlam_lasso
# prediction
lasso_pred <- predict(lasso_mod, s = bestlam_lasso, newx = as.matrix(X_test),type = "response")
coef <- predict(lasso_mod , type = "coefficients", s = bestlam_lasso)
lasso_coef <- as.matrix(coef)
lasso_coef <- lasso_coef[-1, , drop = FALSE]
sum(lasso_coef!=0)
selected_features <- rownames(lasso_coef)[lasso_coef!=0]

library(pROC) 
roc_obj_lasso <- roc(response = as.factor(y_test),
                     predictor = as.numeric(lasso_pred))
roc_obj_lasso$auc
best_cutoff <- coords(roc_obj_lasso, "best", best.method = "youden")

pred_lasso <- ifelse(lasso_pred > best_cutoff$threshold, 1, 0)

# Calculate the confusion_mat in the training and test sets
library(caret)
confusion_mat <- confusionMatrix(as.factor(pred_lasso), as.factor(y_test), positive = "1")
print(confusion_mat)
caret::RMSE(lasso_pred, y_test)

# keep data
# objects_to_save <- pred_lasso
# file_path <- "~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_preddata_lasso.csv"
# save(objects_to_save, file=file_path)

LassoResult_After <- selected_features


