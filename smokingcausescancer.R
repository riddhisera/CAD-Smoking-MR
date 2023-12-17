
selected_data <- subset(nongene_data, (current_tobacco_smoking == 1 | current_tobacco_smoking == 0) & 
                          (Smoking_status == 0 | Smoking_status == 2))
cad_0_subset <- subset(selected_data, CAD == 0)
cad_1_subset <- subset(selected_data, CAD == 1)
min_subset_size <- min(nrow(cad_0_subset), nrow(cad_1_subset))
cad_0_subset_sampled <- cad_0_subset[sample(1:nrow(cad_0_subset), min_subset_size), ]
cad_1_subset_sampled <- cad_1_subset[sample(1:nrow(cad_1_subset), min_subset_size), ]

# Combine the two subsets to get the final subset with an equal number of CAD values 0 and 1
equal_cad_subset <- rbind(cad_0_subset_sampled, cad_1_subset_sampled)

data_selected <- equal_cad_subset[c("smoking_current", "Diastolic_BP0", "Systolic_BP0", "PCE_Risk", "CAD")]
data_selected <- na.omit(data_selected)

library(caret)

# Normalizing or standardizing the numerical features
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
data_selected[c("Diastolic_BP0", "Systolic_BP0", "PCE_Risk")] <- as.data.frame(lapply(data_selected[c("Diastolic_BP0", "Systolic_BP0", "PCE_Risk")], normalize))

# Splitting the dataset into training and testing sets
set.seed(123) # for reproducible results
train_index <- sample(1:nrow(data_selected), 0.8 * nrow(data_selected))
train_data <- data_selected[train_index, ]
test_data <- data_selected[-train_index, ]

# Load necessary library
library(e1071)

# Train the SVM model
svm_model <- svm(CAD ~ ., data = train_data, type = 'C-classification', kernel = 'poly')

# Predict using the SVM model
svm_predictions <- predict(svm_model, test_data)



# Load necessary library
library(randomForest)

# Train the Random Forest model
rf_model <- randomForest(CAD ~ ., data = train_data, ntree = 1000)

# Predict using the Random Forest model
rf_predictions <- predict(rf_model, test_data)



# Calculate accuracy
svm_accuracy <- sum(test_data$CAD == svm_predictions) / nrow(test_data)
print(paste("Accuracy of SVM Model: ", svm_accuracy))


# Calculate accuracy
rf_accuracy <- sum(test_data$CAD == rf_predictions) / nrow(test_data)
print(paste("Accuracy of Random Forest Model: ", rf_accuracy))




normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

normalized_data <- as.data.frame(lapply(equal_cad_subset[factors], normalize))

standardize <- function(x) {
  return ((x - mean(x)) / sd(x))
}

standardized_data <- as.data.frame(lapply(equal_cad_subset, standardize))

## Converting numerical variable to factors
factors <- c("smoking_current", "Diastolic_BP0", "Systolic_BP0",
             "PCE_Risk")
equal_cad_subset[factors] <- lapply(equal_cad_subset[factors], factor)

library(glmnet)
library(e1071) 
library(caret)

## Selecting our input and output space
train_index <- createDataPartition(y = equal_cad_subset$CAD, p = 0.75, list = FALSE) 
X_train <- equal_cad_subset[train_index, -123]
y_train <- equal_cad_subset[train_index, 123]
X_test <- equal_cad_subset[-train_index, -123]
y_test <- equal_cad_subset[-train_index, 123]

# Linear
svm_linear <- svm(y_train ~ smoking_current + Diastolic_BP0 + Systolic_BP0 + PCE_Risk, data = X_train, kernel = "polynomial")
predictions_linear <- predict(svm_linear, X_test)
accuracy_linear <- sum(predictions_linear == y_test) / length(y_test)
pred_num_linear <- as.numeric(as.character(predictions_linear)) 
y_num_test <- as.numeric(as.character(y_test))
roc_linear <- roc(y_num_test, pred_num_linear)

print(paste("Accuracy (Linear Kernel):", accuracy_linear))

# Assuming your data is in a dataframe named df
# Logistic regression to predict CAD
model <- glm(CAD ~ Age + smoking_current + Diastolic_blood_pressure +
               Systolic_blood_pressure + Cholesterol + chol_ratio + Diabetes +
               Father_Heart_disease + Mother_Heart_disease + Siblings_Heart_disease +
               Beef_intake + Bread_intake + Cereal_intake +
               sleep_var, data = equal_cad_subset, family = binomial())

# Summary of the model
summary(model)



# Logistic Regression Model
model <- glm(CAD ~ smoking_current, data = equal_cad_subset, family = binomial())

# Summary of the model
summary(model)

# Predicting CAD
# This gives the probability of having CAD
predicted_probs <- predict(model, type = "response")

# Convert probabilities to a binary prediction based on a threshold (e.g., 0.5)
predicted_class <- ifelse(predicted_probs > 0.5, 1, 0)

# You can also evaluate the model's performance using appropriate metrics
# (like confusion matrix, accuracy, sensitivity, specificity, etc.)

# Assuming df$CAD is your actual outcomes and predicted_class are your model predictions

# Convert actual outcomes to the same format as your predictions if necessary
actual_class <- ifelse(equal_cad_subset$CAD == '1', 1, 0) # Adjust this based on how CAD is coded in your dataset

# Create a confusion matrix
confusion_matrix <- table(Predicted = predicted_class, Actual = actual_class)

# Calculate accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

# Print the accuracy
print(accuracy)

# Install and load necessary packages
if (!require("caret")) install.packages("caret")
if (!require("ggplot2")) install.packages("ggplot2")
library(caret)
library(ggplot2)

# Assuming predicted_class and actual_class are already defined as shown previously
confusion_matrix <- confusionMatrix(as.factor(predicted_class), as.factor(actual_class))

# Print confusion matrix
print(confusion_matrix)

# For a basic plot
fourfoldplot(confusion_matrix$table)


