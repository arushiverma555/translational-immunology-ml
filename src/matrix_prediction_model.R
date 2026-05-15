# Matrix Prediction Model
# Exploratory ML model to classify biological matrix type from cytokine profiles

library(tidyverse)
library(randomForest)
library(caret)

# Load processed cytokine data
panel1 <- read_csv("data/Panel1_vert_clean.csv", show_col_types = FALSE)
panel2 <- read_csv("data/Panel2_vert_clean.csv", show_col_types = FALSE)
panel3 <- read_csv("data/Panel3_vert_clean.csv", show_col_types = FALSE)

all_data <- bind_rows(panel1, panel2, panel3) %>%
  rename_with(~ str_replace_all(., "\\s+", "_")) %>%
  rename(
    Protein_Name = any_of(c("Protein_Name", "Analyte")),
    Cytokine_Concentration = any_of(c("Cytokine_Concentration", "Concentration", "Value")),
    Freeze_Thaw_Cycle = any_of(c("Freeze_Thaw_Cycle", "Freeze_Thaw"))
  ) %>%
  mutate(
    Matrix = as.factor(Matrix),
    Donor_ID = as.factor(Donor_ID),
    Freeze_Thaw_Cycle = as.factor(Freeze_Thaw_Cycle),
    Log_Concentration = log10(as.numeric(Cytokine_Concentration) + 1)
  )

# Convert from long format to wide format:
# each sample = one row, each cytokine = one feature
wide_data <- all_data %>%
  group_by(Donor_ID, Matrix, Freeze_Thaw_Cycle, Protein_Name) %>%
  summarise(Log_Concentration = mean(Log_Concentration, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Protein_Name,
    values_from = Log_Concentration
  )

# Keep matrix label + cytokine features
model_data <- wide_data %>%
  select(-Donor_ID, -Freeze_Thaw_Cycle)

# Remove cytokines with too much missing data
missing_rate <- colMeans(is.na(model_data %>% select(-Matrix)))
keep_features <- names(missing_rate[missing_rate < 0.5])

model_data <- model_data %>%
  select(Matrix, all_of(keep_features))

# Impute remaining missing values with feature medians
model_data <- model_data %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.), median(., na.rm = TRUE), .)
  ))

print("Matrix class counts:")
print(table(model_data$Matrix))

# Train/test split
set.seed(123)

train_index <- createDataPartition(
  model_data$Matrix,
  p = 0.75,
  list = FALSE
)

train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

# Random forest model
rf_model <- randomForest(
  Matrix ~ .,
  data = train_data,
  importance = TRUE,
  ntree = 500
)

# Evaluate model
predictions <- predict(rf_model, test_data)

print("Confusion matrix:")
print(confusionMatrix(predictions, test_data$Matrix))

print("Random forest model:")
print(rf_model)

# Save feature importance plot
png("figures/matrix_prediction_feature_importance.png",
    width = 1200,
    height = 900)

varImpPlot(rf_model, main = "Feature Importance for Matrix Prediction")

dev.off()

# Save predictions
prediction_results <- tibble(
  Actual_Matrix = test_data$Matrix,
  Predicted_Matrix = predictions
)

write_csv(prediction_results, "results/matrix_prediction_results.csv")

print("Saved outputs:")
print("figures/matrix_prediction_feature_importance.png")
print("results/matrix_prediction_results.csv")
