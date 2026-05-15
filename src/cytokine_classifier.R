# Cytokine Stability Classification Pipeline
# Translational Immunology ML Project

library(tidyverse)
library(randomForest)
library(caret)

# Load cytokine classification tables
stable <- read_csv("results/stable_cytokines.csv") %>%
  mutate(class = "stable")

decreasing <- read_csv("results/decreasing_cytokines.csv") %>%
  mutate(class = "decreasing")

matrix_effect <- read_csv("results/matrix_effect_cytokines.csv") %>%
  mutate(class = "matrix_effect")

# Combine datasets
cytokine_data <- bind_rows(
  stable,
  decreasing,
  matrix_effect
)

# Keep numeric features only
numeric_data <- cytokine_data %>%
  select(where(is.numeric))

# Add classification labels
numeric_data$class <- cytokine_data$class

# Remove rows with missing values
numeric_data <- na.omit(numeric_data)

# Train/test split
set.seed(123)

train_index <- createDataPartition(
  numeric_data$class,
  p = 0.8,
  list = FALSE
)

train_data <- numeric_data[train_index, ]
test_data <- numeric_data[-train_index, ]

# Random Forest classifier
rf_model <- randomForest(
  class ~ .,
  data = train_data,
  importance = TRUE,
  ntree = 500
)

# Predictions
predictions <- predict(rf_model, test_data)

# Model evaluation
confusionMatrix(predictions, test_data$class)

# Variable importance
importance(rf_model)

# Save importance plot
png("figures/random_forest_feature_importance.png",
    width = 1200,
    height = 900)

varImpPlot(rf_model)

dev.off()
