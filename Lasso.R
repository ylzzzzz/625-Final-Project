mydata <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata.csv")
gene_data <- mydata[, 7:ncol(mydata)]
clinical_data <- mydata[, 1:6]
# train data and test data
n <- nrow(mydata)
train_index <- sample(seq_len(n), size = 0.7 * n)

train_data <- mydata[train_index, ]
test_data  <- mydata[-train_index, ]
X_train <- train_data[, 7:ncol(mydata)]
X_test <- test_data[, 7:ncol(mydata)]
y_train <- train_data[, 2]
y_test <- test_data[, 2]

# 检查划分结果
nrow(train_data)
nrow(test_data)

# 查看缺失值
missing_values <- is.na(test_data)
any(missing_values)
missing_values <- is.na(train_data)
any(missing_values)


# deal with NA ------------------------------------------------------------
clean_NA <- function(X, max_na = 0.5, method = c("median", "mean")) {

  method <- match.arg(method)
  
  if (!is.matrix(X)) {
    X <- as.data.frame(X)
  }
  
  na_prop <- colMeans(is.na(X))
  too_many_na <- na_prop > max_na
  if (any(too_many_na)) {
    cat(sprintf("Removed %d variables with > %.0f%% missing values:\n", 
                sum(too_many_na), 100 * max_na))
    print(names(X)[too_many_na])
    X <- X[, !too_many_na, drop = FALSE]
  }
  

  cat("Imputing remaining NA values using", method, "...\n")
  for (j in seq_len(ncol(X))) {
    if (anyNA(X[[j]])) {
      if (method == "median") {
        fill_value <- median(X[[j]], na.rm = TRUE)
      } else if (method == "mean") {
        fill_value <- mean(X[[j]], na.rm = TRUE)
      }
      X[[j]][is.na(X[[j]])] <- fill_value
    }
  }
  
  X <- as.matrix(sapply(X, as.numeric))
  
  if (anyNA(X)) {
    stop("There are still NA values after imputation!")
  } else {
    cat("all missing values handled successfully.\n")
    cat(sprintf("Final dimensions: %d samples x %d features\n", nrow(X), ncol(X)))
  }
  
  return(X)
}

genedata_clean <- clean_NA(gene_data, max_na = 0.5, method = "median")
X_train_clean <- genedata_clean[train_index, ]
X_test_clean  <- genedata_clean[-train_index, ]
y_train <- train_data[, 2]
y_test <- test_data[, 2]

#------------------lasso----------------------------

library(glmnet)
set.seed(1)
grid <- 10^seq(10, -2, length = 100)
# fit
lasso_mod <- glmnet(X_train_clean, y_train, alpha = 1, lambda = grid,family = "binomial")
plot(lasso_mod)
# select parameter
cv_out_lasso <- cv.glmnet(X_train_clean, y_train, alpha = 1)
plot(cv_out_lasso)
bestlam_lasso <- cv_out_lasso$lambda.min
bestlam_lasso
# prediction
lasso_pred <- predict(lasso_mod, s = bestlam_lasso, newx = X_test_clean,type = "response")
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

#保存数据
objects_to_save <- pred_lasso
file_path <- "~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_preddata_lasso.csv"
save(objects_to_save, file=file_path)



