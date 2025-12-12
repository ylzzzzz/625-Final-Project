start_time <- Sys.time()

mydata <- read.csv("~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_cleandata_final.csv")
dim(mydata)
end_time <- Sys.time()
time_used <- end_time - start_time
print(time_used)

genedata_clean <- mydata[, 7:ncol(mydata)]
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


X_train_clean <- genedata_clean[train_index, ]
X_test_clean  <- genedata_clean[-train_index, ]
y_train <- train_data[, 2]
y_test <- test_data[, 2]

#------------------lasso----------------------------

library(glmnet)
set.seed(1)
grid <- 10^seq(10, -2, length = 100)
# fit
start_time <- Sys.time()

lasso_mod <- glmnet(X_train_clean, y_train, alpha = 1, lambda = grid,family = "binomial")

end_time <- Sys.time()
time_used <- end_time - start_time
print(time_used)

plot(lasso_mod)
# select parameter
cv_out_lasso <- cv.glmnet(as.matrix(X_train_clean), y_train, alpha = 1)
plot(cv_out_lasso)
bestlam_lasso <- cv_out_lasso$lambda.min
bestlam_lasso
# prediction
lasso_pred <- predict(lasso_mod, s = bestlam_lasso, newx = as.matrix(X_test_clean),type = "response")
coef <- predict(lasso_mod , type = "coefficients", s = bestlam_lasso)
lasso_coef <- as.matrix(coef)
lasso_coef <- lasso_coef[-1, , drop = FALSE]
sum(lasso_coef!=0)
selected_features <- rownames(lasso_coef)[lasso_coef!=0]


# SVM -------------------------------------------------------------------

library(caret)
set.seed(1)

X_train_afterlasso <- X_train_clean[,selected_features]
y_train_factor <- as.factor(y_train)

trainData1 <- data.frame(X_train_afterlasso,
                         y_train = as.factor(y_train))

tuneGrid <- expand.grid(C = c(0.1, 1, 10),
                        sigma = c(0.01, 0.1, 1))
control <- trainControl(method = "cv", number = 10)

start_time <- Sys.time()
svm_model <- train(
  y_train ~ ., 
  data = trainData1,
  method = "svmRadial",
  trControl = control,
  tuneGrid = tuneGrid
  #prob.model = TRUE
)# Use the tune function for parameter optimization

end_time <- Sys.time()

time_used <- end_time - start_time
print(time_used)

bestModel <- svm_model$bestTune

testData <- data.frame(X_test_clean[,selected_features],y_test = as.factor(y_test))
predictions <- predict(svm_model, X_test_clean[,selected_features]) # prediction
# predictions <- predict(svm_model, X_test_clean[,selected_features],type = "prob")

library(caret)
confusion_mat <- confusionMatrix(as.factor(predictions), as.factor(y_test), positive = "1")
print(confusion_mat)
#caret::RMSE(predictions$`1`, y_test) #RSME:0.3746263

# #保存数据
# pred_SVM <- as.data.frame(predictions)
# objects_to_save <- pred_SVM
# file_path <- "~/Documents/DNA/2025Fall/Biostat625/Final/Alzheimers_Disease_preddata_SVM.csv"
# save(objects_to_save, file=file_path)




