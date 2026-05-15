# Matrix Prediction Random Forest Model
# Translational Immunology ML Project

library(tidyverse)
library(randomForest)
library(caret)
library(janitor)

# Load cleaned cytokine datasets
panel1 <- read_csv("data/data/Panel1_vert_clean.csv", show_col_types = FALSE)
panel2 <- read_csv("data/data/Panel2_vert_clean.csv", show_col_types = FALSE)
panel3 <- read_csv("data/data/Panel3_vert_clean.csv", show_col_types = FALSE)

# Combine all panels
all_data <- bind_rows(panel1, panel2, panel3)

# Check column names
print(colnames(all_data))

# Keep relevant columns
all_data <- all_data %>%
  select(
    Donor_ID,
    Matrix,
    Freeze_Thaw,
    Analyte,
    Cytokine_Concentration
  )

# Create log-transformed concentration values
all_data <- all_data %>%
  mutate(
    Log_Concentration = log10(Cytokine_Concentration + 1)
  )

# Convert long-format cytokine data into wide-format ML table
wide_data <- all_data %>%
  group_by(Donor_ID, Matrix, Freeze_Thaw, Analyte) %>%
  summarise(
    Log_Concentration = mean(Log_Concentration, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Analyte,
    values_from = Log_Concentration
  ) %>%
  janitor::clean_names()

# Prepare ML dataset
model_data <- wide_data %>%
  rename(Matrix = matrix) %>%
  select(-donor_id, -freeze_thaw) %>%
  mutate(Matrix = as.factor(Matrix))

# Remove missing values
model_data <- na.omit(model_data)

# Show matrix class balance
print("Matrix class counts:")
print(table(model_data$Matrix))

# Train/test split
set.seed(123)

train_index <- createDataPartition(
  model_data$Matrix,
  p = 0.8,
  list = FALSE
)

train_data <- model_data[train_index, ]
test_data  <- model_data[-train_index, ]

# Train random forest classifier
rf_model <- randomForest(
  Matrix ~ .,
  data = train_data,
  importance = TRUE,
  ntree = 500
)

# Print model summary
print(rf_model)

# Generate predictions
predictions <- predict(rf_model, test_data)

# Evaluate model
conf_matrix <- confusionMatrix(
  predictions,
  test_data$Matrix
)

print(conf_matrix)

# Feature importance
importance_df <- importance(rf_model)

print(importance_df)

# Save feature importance plot
png(
  "figures/matrix_prediction_feature_importance.png",
  width = 1200,
  height = 800
)

varImpPlot(
  rf_model,
  main = "Random Forest Feature Importance"
)

dev.off()

# Save prediction outputs
prediction_results <- data.frame(
  Actual = test_data$Matrix,
  Predicted = predictions
)

write_csv(
  prediction_results,
  "results/matrix_prediction_results.csv"
)

print("Model completed successfully.")
