# Set seed for reproducibility
set.seed(123)

# Create a 70/30 split using base R
index <- sample(1:nrow(model_data), size = 0.7 * nrow(model_data))

train_data <- model_data[index, ]
test_data  <- model_data[-index, ]

print(paste("Training set dimensions:", paste(dim(train_data), collapse = "x")))
print(paste("Testing set dimensions:", paste(dim(test_data), collapse = "x")))# Set seed for reproducibility
set.seed(123)

# Create a 70/30 split using base R
index <- sample(1:nrow(model_data), size = 0.7 * nrow(model_data))

train_data <- model_data[index, ]
test_data  <- model_data[-index, ]

print(paste("Training set dimensions:", paste(dim(train_data), collapse = "x")))
print(paste("Testing set dimensions:", paste(dim(test_data), collapse = "x")))


# Fit the model
logit_model <- glm(Alzheimers_Disease ~ ., data = train_data, family = binomial)

# View summary (coeffecients, p-values, AIC)
summary(logit_model)


# Predict probabilities on the test set
predicted_probs <- predict(logit_model, newdata = test_data, type = "response")

# Convert probabilities to classes (using 0.5 threshold)
# Assuming the second level is the "Positive" case (e.g., has disease)
positive_class <- levels(model_data$Alzheimers_Disease)[2]
negative_class <- levels(model_data$Alzheimers_Disease)[1]

predicted_classes <- ifelse(predicted_probs > 0.5, positive_class, negative_class)

# Ensure they are factors with the same levels for the confusion matrix
predicted_classes <- factor(predicted_classes, levels = levels(test_data$Alzheimers_Disease))


library(caret) # Highly recommended for evaluation

# Create Confusion Matrix
conf_matrix <- confusionMatrix(predicted_classes, test_data$Alzheimers_Disease)
print(conf_matrix)