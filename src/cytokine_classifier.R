# Cytokine Stability Classification Pipeline
# Exploratory ML extension for cytokine freeze-thaw analysis

library(tidyverse)
library(randomForest)
library(janitor)

clean_results <- function(path, class_label) {
  read_csv(path, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(class = class_label) %>%
    mutate(across(
      everything(),
      ~ ifelse(. == "<0.001", "0.001", .)
    )) %>%
    mutate(across(
      matches("slope|percent|p_value|interaction"),
      ~ suppressWarnings(as.numeric(.))
    ))
}

stable <- clean_results("results/stable_cytokines.csv", "stable")
decreasing <- clean_results("results/decreasing_cytokines.csv", "decreasing")
matrix_effect <- clean_results("results/matrix_effect_cytokines.csv", "matrix_effect")

cytokine_data <- bind_rows(stable, decreasing, matrix_effect)

print("Class counts before modeling:")
print(table(cytokine_data$class))

model_data <- cytokine_data %>%
  select(where(is.numeric), class)

# Remove numeric columns that are completely empty
model_data <- model_data %>%
  select(where(~ !all(is.na(.))), class)

# Replace remaining missing numeric values with column medians
model_data <- model_data %>%
  mutate(across(
    where(is.numeric),
    ~ ifelse(is.na(.), median(., na.rm = TRUE), .)
  )) %>%
  mutate(class = as.factor(class))

print("Class counts after cleaning:")
print(table(model_data$class))

set.seed(123)

rf_model <- randomForest(
  class ~ .,
  data = model_data,
  importance = TRUE,
  ntree = 500
)

print(rf_model)
print(importance(rf_model))

png("figures/random_forest_feature_importance.png",
    width = 1200,
    height = 900)

varImpPlot(rf_model)

dev.off()
