

########## Models : ##########################
library(caret)
library(e1071)       # For SVM
library(class)       # For knn (caret uses internally)
library(tidyverse)   # For data handling
library(rpart)
library(randomForest)
library(MASS)

# Let's start with models using the dataset after log transformation of features : 
data_model = data_log

## Train-Test split : 
set.seed(42)
train_index <- createDataPartition(data_model$koi_disposition, p = 0.8, list = FALSE)
train_data <- data_model[train_index, ]
test_data  <- data_model[-train_index, ]
##
## Linear Regression : 
lm_model <- lm(as.numeric(koi_disposition) ~ ., data = train_data)
lm_pred <- predict(lm_model, newdata = test_data)
lm_class_pred <- ifelse(lm_pred > 0.5, 1, 0)
table(Predicted = lm_class_pred, Actual = test_data$koi_disposition)

summary(lm_model)

## Logistic Regression : 
log_model <- glm(koi_disposition ~ ., data = train_data, family = binomial)
summary(log_model)
# Predict on test set (probabilities then classes)
log_probs <- predict(log_model, newdata = test_data, type = "response")
log_pred <- ifelse(log_probs > 0.5, 1, 0)
# Confusion matrix
table(Predicted = log_pred, Actual = test_data$koi_disposition)

## SVMs : 
svm_model <- svm(koi_disposition ~ ., data = train_data, probability = TRUE)
# Predict
svm_pred <- predict(svm_model, newdata = test_data)
# Confusion matrix
table(Predicted = svm_pred, Actual = test_data$koi_disposition)


## KNN : 
train_x <- scale(train_data[, !names(train_data) %in% "koi_disposition"])
test_x  <- scale(test_data[, !names(test_data) %in% "koi_disposition"], center = attr(train_x, "scaled:center"), scale = attr(train_x, "scaled:scale"))
train_y <- train_data$koi_disposition
test_y  <- test_data$koi_disposition
# Run kNN 
knn_pred <- knn(train = train_x, test = test_x, cl = train_y, k = 5)
# Confusion matrix
table(Predicted = knn_pred, Actual = test_y)


## Decision Tree : 
tree_model <- rpart(koi_disposition ~ ., data = train_data, method = "class")

# Predict on test set
tree_pred <- predict(tree_model, newdata = test_data, type = "class")

# Confusion matrix
table(Predicted = tree_pred, Actual = test_data$koi_disposition)


## Random-Forest : 
rf_model <- randomForest(koi_disposition ~ ., data = train_data, ntree = 500)

# Predict on test set
rf_pred <- predict(rf_model, newdata = test_data)

# Confusion matrix
table(Predicted = rf_pred, Actual = test_data$koi_disposition)

# variable importance
importance(rf_model)
varImpPlot(rf_model)


### Now with the dataset after PCA : 
data_model = as.data.frame(pca_capped$x)
data_model$koi_disposition = data_log$koi_disposition

## Train-Test split : 
set.seed(42)
train_index <- createDataPartition(data_model$koi_disposition, p = 0.8, list = FALSE)
train_data <- data_model[train_index, ]
test_data  <- data_model[-train_index, ]
##
## Linear Regression : 
lm_model <- lm(as.numeric(koi_disposition) ~ ., data = train_data)
lm_pred <- predict(lm_model, newdata = test_data)
lm_class_pred <- ifelse(lm_pred > 0.5, 1, 0)
table(Predicted = lm_class_pred, Actual = test_data$koi_disposition)

summary(lm_model)
## LDA : 
lda_model <- lda(koi_disposition ~ ., data = train_data)
# Predict on test data
lda_pred <- predict(lda_model, newdata = test_data)
# Predicted classes
lda_class <- lda_pred$class
# Confusion matrix
table(Predicted = lda_class, Actual = test_data$koi_disposition)

## QDA : 
qda_model <- qda(koi_disposition ~ ., data = train_data)

# Predict on test data
qda_pred <- predict(qda_model, newdata = test_data)

# Predicted classes
qda_class <- qda_pred$class

# Confusion matrix
table(Predicted = qda_class, Actual = test_data$koi_disposition)


## Logistic Regression : 
log_model <- glm(koi_disposition ~ ., data = train_data, family = binomial)
summary(log_model)
# Predict on test set (probabilities then classes)
log_probs <- predict(log_model, newdata = test_data, type = "response")
log_pred <- ifelse(log_probs > 0.5, 1, 0)
# Confusion matrix
table(Predicted = log_pred, Actual = test_data$koi_disposition)

## SVMs : 
svm_model <- svm(koi_disposition ~ ., data = train_data, probability = TRUE)
# Predict
svm_pred <- predict(svm_model, newdata = test_data)
# Confusion matrix
table(Predicted = svm_pred, Actual = test_data$koi_disposition)


## KNN : 
train_x <- scale(train_data[, !names(train_data) %in% "koi_disposition"])
test_x  <- scale(test_data[, !names(test_data) %in% "koi_disposition"], center = attr(train_x, "scaled:center"), scale = attr(train_x, "scaled:scale"))
train_y <- train_data$koi_disposition
test_y  <- test_data$koi_disposition
# Run kNN 
knn_pred <- knn(train = train_x, test = test_x, cl = train_y, k = 5)
# Confusion matrix
table(Predicted = knn_pred, Actual = test_y)


## Decision Tree : 
tree_model <- rpart(koi_disposition ~ ., data = train_data, method = "class")

# Predict on test set
tree_pred <- predict(tree_model, newdata = test_data, type = "class")

# Confusion matrix
table(Predicted = tree_pred, Actual = test_data$koi_disposition)


## Random-Forest : 
rf_model <- randomForest(koi_disposition ~ ., data = train_data, ntree = 500)

# Predict on test set
rf_pred <- predict(rf_model, newdata = test_data)

# Confusion matrix
table(Predicted = rf_pred, Actual = test_data$koi_disposition)

# variable importance
importance(rf_model)
varImpPlot(rf_model)
