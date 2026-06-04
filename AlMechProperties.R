##########################################################
# Data download and exploration
##########################################################

# Note: this process could take a couple of minutes

# Installation of required libraries
options(timeout = 120)
options(scipen = 999)
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(broom)) install.packages("broom", repos = "http://cran.us.r-project.org")
if(!require(caretEnsemble)) install.packages("caretEnsemble", repos = "http://cran.us.r-project.org")
library(tidyverse)
library(caret)
library(broom)
library(caretEnsemble)

# Aluminum alloy dataset:
# https://data.mendeley.com/public-api/zip/b6br4yk6r3/download/1
# Citation:
# Bhat, Ninad; Barnard, Amanda; Birbilis, Nick (2023), “Aluminium alloy dataset for supervised learning”, Mendeley Data, V1, doi: 10.17632/b6br4yk6r3.1
url = "https://data.mendeley.com/public-files/datasets/b6br4yk6r3/files/bc98e0d7-6c56-42dc-aa54-d689cbfe1230/file_downloaded"
# Read the file directly from the file connection
raw_alloy_data <- read.csv(url)

# Inspect the dataset
# Consists of an aluminum alloy index, a processing method, 25 elements, 3 mechanical properties, and a class based on the type of aluminum alloy
str(raw_alloy_data) 

# 9 classes corresponding to [Series]xxx alloy names generally indicative of families of alloy compositions
# 3 alloys are "outlier", but 8xxx series is defined as "Other" alloys that do not fit in other alloy series
# Considering that alloy series are already defined by composition and processing, class can be omitted for the model
# This also eliminates accounting for outliers
raw_alloy_data %>% group_by(class) %>% count()

# Select only the element columns for a grid of histograms 
elements <- c("Ag", "Al", "B", "Be", "Bi", "Cd", "Co", "Cr", "Cu", "Er", "Eu", "Fe", "Ga", "Li", "Mg", "Mn", "Ni", "Pb", "Sc", "Si", "Sn", "Ti", "V", "Zn", "Zr")
# Make a 5x5 tiled grid of histograms where each tile represents the distribution of an element's content in all alloys in the dataset
raw_alloy_data %>%
  select(all_of(elements)) %>% 
  pivot_longer(cols = everything(), names_to = "Element", values_to = "Percent") %>%
  ggplot(aes(x = Percent)) +
  geom_histogram(bins=8, fill="skyblue", color="white") +
  facet_wrap(~Element, nrow=5, ncol=5, scales="free") +
  theme(
    axis.text.x = element_text(
      angle = 90,          # Rotates text 90 degrees
      vjust = 0.5,         # Centers text vertically relative to tick
      hjust = 1,           # Aligns text to the right (useful for 45-90 degree angles)
      size = 8             # Shrinks font size (adjust as needed)
    )) +
  labs(title="Distributions of Compositions in Dataset by Element") +
  # Use log scales on y axis to help visualize distributions without being overly dominated by counts of alloying elements equal to 0 (most cases)
  scale_y_continuous( 
    breaks = c(1, 10, 100, 1000, 10000), 
    trans = scales::pseudo_log_trans(base = 10)
  )

# 10 processes, some of which are combinations
# Solutionized on its own has n=4 but is found in several other numerous combinations
raw_alloy_data %>% group_by(Processing) %>% count()




##########################################################
# Dataset cleaning
##########################################################

# Rename columns to names that are easier to work with 
# The "class" column name causes confusion for grouping and other exploration
# Renamed to "Series" as is common in the aluminum industry
alloy_data <- raw_alloy_data %>% rename(
  elongation = "Elongation....",
  tensile_strength = "Tensile.Strength..MPa.",
  yield_strength = "Yield.Strength..MPa.",
  series = "class")

# Simplify similar heat treatments and one-hot encode into 4 new columns
alloy_data <- alloy_data %>% mutate(
  solutionised = str_detect(Processing, "Solutionised"),
  naturally_aged = str_detect(Processing, "Natural"),
  artificially_aged = str_detect(Processing, "Artificial"),
  strain_hardened = str_detect(Processing, "Strain|Cold Work")
) 

# Change "outlier" values to "8"
alloy_data$series[alloy_data$series == "outlier"] <- "8"

# Defining balance as the difference between the sum of all elements as reported and a perfect 100%
alloy_data <- alloy_data %>% mutate( balance = 1 - Ag - Al - B - Be - Bi - Cd - Co - Cr - Cu - Er - Eu - Fe - Ga - Li - Mg - Mn - Ni - Pb - Sc - Si - Sn - Ti - V - Zn - Zr)

# Box-Cox Tranformation standard lambda meanings:
# -1 = inverse transform
# -0.5 = inverse square root transform
# 0 = log transform
# 0.5 = square root transform
# 1 = no transform
# 2 = square transform
lambdas <- c(-1, -0.5, 0, 0.5, 1, 2)

# Build a long data frame applying Box-Cox transformations with each lambda to the elongation column with null values removed
e_long <- data.frame(original = na.omit(alloy_data$elongation)) %>%
  cross_join(data.frame(Lambda = lambdas)) %>%
  mutate(
    Transformed = case_when(
      Lambda == 0 ~ log(original),
      TRUE        ~ (original^Lambda - 1) / Lambda
    )
  ) %>%
  group_by(Lambda) %>%
  mutate(
    P_Value = shapiro.test(Transformed)$p.value,
    # Create a clean label format for the facet headers that include the lambda used to construct each plot and the p-value of each against the null hypothesis that the distribution is normal
    Label = paste0(
      "Lambda = ", Lambda, "\n",
      "p = ", if_else(P_Value < 0.001, formatC(P_Value, format = "e", digits = 2), as.character(round(P_Value, 3)))
    )
  ) %>%
  ungroup()
# Plot grid of qq-plots for each Box-Cox transformation to select best approach for elongation
ggplot(e_long, aes(sample = Transformed)) +
  stat_qq() +
  stat_qq_line(color = "firebrick", linewidth = 0.8) +
  facet_wrap(~ Label, scales = "free", nrow = 2) +
  labs(title = "Box-Cox Lambda Grid Evaluation: Elongation",
       subtitle = "Higher p-values indicate closer fit to normal distribution",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 10)) # Make headers pop

# Build a long data frame applying Box-Cox transformations with each lambda to the tensile strength column with null values removed
ts_long <- data.frame(original = na.omit(alloy_data$tensile_strength)) %>%
  cross_join(data.frame(Lambda = lambdas)) %>%
  mutate(
    Transformed = case_when(
      Lambda == 0 ~ log(original),
      TRUE        ~ (original^Lambda - 1) / Lambda
    )
  ) %>%
  group_by(Lambda) %>%
  mutate(
    P_Value = shapiro.test(Transformed)$p.value,
    # Create a clean label format for the facet headers that include the lambda used to construct each plot and the p-value of each against the null hypothesis that the distribution is normal
    Label = paste0(
      "Lambda = ", Lambda, "\n",
      "p = ", if_else(P_Value < 0.001, formatC(P_Value, format = "e", digits = 2), as.character(round(P_Value, 3)))
    )
  ) %>%
  ungroup()
# Plot grid of qq-plots for each Box-Cox transformation to select best approach for tensile strength
ggplot(ts_long, aes(sample = Transformed)) +
  stat_qq() +
  stat_qq_line(color = "firebrick", linewidth = 0.8) +
  facet_wrap(~ Label, scales = "free", nrow = 2) +
  labs(title = "Box-Cox Lambda Grid Evaluation: Tensile Strength",
       subtitle = "Higher p-values indicate closer fit to normal distribution",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 10)) # Make headers pop

# Build a long data frame applying Box-Cox transformations with each lambda to the yield strength column with null values removed
ys_long <- data.frame(original = na.omit(alloy_data$yield_strength)) %>%
  cross_join(data.frame(Lambda = lambdas)) %>%
  mutate(
    Transformed = case_when(
      Lambda == 0 ~ log(original),
      TRUE        ~ (original^Lambda - 1) / Lambda
    )
  ) %>%
  group_by(Lambda) %>%
  mutate(
    P_Value = shapiro.test(Transformed)$p.value,
    # Create a clean label format for the facet headers that include the lambda used to construct each plot and the p-value of each against the null hypothesis that the distribution is normal
    Label = paste0(
      "Lambda = ", Lambda, "\n",
      "p = ", if_else(P_Value < 0.001, formatC(P_Value, format = "e", digits = 2), as.character(round(P_Value, 3)))
    )
  ) %>%
  ungroup()
# Plot grid of qq-plots for each Box-Cox transformation to select best approach for yield strength
ggplot(ys_long, aes(sample = Transformed)) +
  stat_qq() +
  stat_qq_line(color = "firebrick", linewidth = 0.8) +
  facet_wrap(~ Label, scales = "free", nrow = 2) +
  labs(title = "Box-Cox Lambda Grid Evaluation: Yield Strength",
       subtitle = "Higher p-values indicate closer fit to normal distribution",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 10)) # Make headers pop

# Square root transforms of mechanical properties
alloy_data <- alloy_data %>% mutate(sqrt_elongation = sqrt(elongation), sqrt_tensile_strength = sqrt(tensile_strength), sqrt_yield_strength = sqrt(yield_strength))




##########################################################
# Create train and test sets
##########################################################

# Data set that focuses only on predicting square roots of tensile strength from series, decomposed processing columns, and alloying elements
# Rows with null tensile_strength are omitted for tensile strength model development
ts_alloy_data <- alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_tensile_strength) %>% filter(!is.na(sqrt_tensile_strength))

# Data set that focuses only on predicting square roots of yield strength from series, decomposed processing columns, and alloying elements
# Rows with null yield_strength are omitted for yield strength model development
ys_alloy_data <- alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_yield_strength) %>% filter(!is.na(sqrt_yield_strength))

# Data set that focuses only on predicting square roots of elongation from series, decomposed processing columns, and alloying elements
# Rows with null elongation are omitted for elongation model development
e_alloy_data <- alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_elongation) %>% filter(!is.na(sqrt_elongation))


# Tensile strength model data split: test set will be 20% of ts_alloy_data subset with square root tensile strength data
set.seed(1, sample.kind="Rounding") # if using R 3.6 or later
ts_test_index <- createDataPartition(y = ts_alloy_data$sqrt_tensile_strength, times = 1, p = 0.2, list = FALSE)
ts_train <- ts_alloy_data[-ts_test_index,]
ts_test <- ts_alloy_data[ts_test_index,]
rm(ts_test_index)

# Yield strength model data split: test set will be 20% of ys_alloy_data subset with square root yield strength data
set.seed(1, sample.kind="Rounding") # if using R 3.6 or later
ys_test_index <- createDataPartition(y = ys_alloy_data$sqrt_yield_strength, times = 1, p = 0.2, list = FALSE)
ys_train <- ys_alloy_data[-ys_test_index,]
ys_test <- ys_alloy_data[ys_test_index,]
rm(ys_test_index)

# Elongation model data split: test set will be 20% of e_alloy_data subset with square root elongation data
set.seed(1, sample.kind="Rounding") # if using R 3.6 or later
e_test_index <- createDataPartition(y = e_alloy_data$sqrt_elongation, times = 1, p = 0.2, list = FALSE)
e_train <- e_alloy_data[-e_test_index,]
e_test <- e_alloy_data[e_test_index,]
rm(e_test_index)




##########################################################
# Tensile strength model
##########################################################

# Generalized linear model (GLM)
# Train with the ts_train data and predict against the ts_test data
ts_fit_glm <- train(sqrt_tensile_strength ~ ., 
                            method = "glm", 
                            data = ts_train)
ts_y_hat_glm <- predict(ts_fit_glm, newdata = ts_test)
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_glm_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_glm)
ggplot(data = ts_glm_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_glm)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: GLM Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Seed Model and Rsquared Results table with GLM results in a rounded format
format_df <- function(data) {
  data[] <- data[] %>% mutate(across(2:3, as.numeric))
  data[] <- lapply(data, function(x) {
    if(is.numeric(x)) {
      sprintf("%.3f", x)
    } else {
      x
    }
  })
  return(data)
}
ts_results_table <- data.frame(
  TensileStrengthModel = "GLM", 
  TrainRsquared = ts_fit_glm$results$Rsquared, # Rsq from train function
  TestRsquared = cor(ts_y_hat_glm, ts_test$sqrt_tensile_strength)^2) # Rsq from predicted vs actual on test set
format_df(ts_results_table) 



# K-nearest neighbors (KNN)
# Train with the ts_train data
ts_fit_knn <- train(sqrt_tensile_strength ~ ., 
                            method = "knn", 
                            data = ts_train,
                            tuneGrid = expand.grid(k = seq(3, 21, by = 1)))
# After tuning, best k is 3
best_k <- ts_fit_knn$bestTune$k
best_k
best_results <- ts_fit_knn$results[ts_fit_knn$results$k == best_k, ]
# Predict and evaluate against the ts_test data
ts_y_hat_knn <- predict(ts_fit_knn, newdata = ts_test)
rsq_ts_knn <- cor(ts_y_hat_knn, ts_test$sqrt_tensile_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_knn_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_knn)
ggplot(data = ts_knn_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_knn)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: KNN Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ts_results_table <- ts_results_table %>% rbind(c("KNN", best_results$Rsquared, rsq_ts_knn))
format_df(ts_results_table)



# Regression trees (rpart)
# Train with the ts_train data
ts_fit_rpart <- train(sqrt_tensile_strength ~ ., 
                              method = "rpart", 
                              data = ts_train,
                              tuneGrid = expand.grid(cp = seq(0, 0.05, len = 25)))
# After tuning, best cp is 0
best_cp <- ts_fit_rpart$bestTune$cp
best_cp
best_results <- ts_fit_rpart$results[ts_fit_rpart$results$cp == best_cp, ]
# Predict and evaluate against the ts_test data
ts_y_hat_rpart <- predict(ts_fit_rpart, newdata = ts_test)
rsq_ts_rpart <- cor(ts_y_hat_rpart, ts_test$sqrt_tensile_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_rpart_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_rpart)
ggplot(data = ts_rpart_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_rpart)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: Regression Trees Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ts_results_table <- ts_results_table %>% rbind(c("Regression Trees", best_results$Rsquared, rsq_ts_rpart))
format_df(ts_results_table)



# Random Forest (rf)
# Train with the ts_train data and predict against the ts_test data
ts_fit_rf <- train(sqrt_tensile_strength ~ ., 
                           method = "rf", 
                           data = ts_train,
                           tuneGrid = expand.grid(.mtry = c(2, 9, 16, 23, 31)))
# After tuning, best mtry is 16
best_mtry <- ts_fit_rf$bestTune$mtry
best_mtry
best_results <- ts_fit_rf$results[ts_fit_rf$results$mtry == best_mtry, ]
# Predict and evaluate against the ts_test data
ts_y_hat_rf <- predict(ts_fit_rf, newdata = ts_test)
rsq_ts_rf <- cor(ts_y_hat_rf, ts_test$sqrt_tensile_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_rf_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_rf)
ggplot(data = ts_rf_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_rf)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: Random Forest Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ts_results_table <- ts_results_table %>% rbind(c("Random Forest", best_results$Rsquared,  rsq_ts_rf))
format_df(ts_results_table)



# eXtreme Gradient Boosted Trees (xgbTree)
# Train with the ts_train data and predict against the ts_test data
ts_fit_xgbt <- train(sqrt_tensile_strength ~ ., 
                             method = "xgbTree", 
                             data = ts_train)
# After tuning, best nrounds is 150
best_nrounds <- ts_fit_xgbt$bestTune$nrounds 
best_nrounds
# After tuning, best max_depth is 3
best_max_depth <- ts_fit_xgbt$bestTune$max_depth 
best_max_depth
# After tuning, best eta is 0.3
best_eta <- ts_fit_xgbt$bestTune$eta 
best_eta
# After tuning, best gamma is 0
best_gamma <- ts_fit_xgbt$bestTune$gamma 
best_gamma
# After tuning, best colsample_bytree is 0.8
best_colsample_bytree  <- ts_fit_xgbt$bestTune$colsample_bytree 
best_colsample_bytree
# After tuning, best min_child_weight is 1
best_min_child_weight  <- ts_fit_xgbt$bestTune$min_child_weight 
best_min_child_weight
# After tuning, best subsample is 1
best_subsample <- ts_fit_xgbt$bestTune$subsample
best_subsample 
best_results <- ts_fit_xgbt$results[
  ts_fit_xgbt$results$nrounds == best_nrounds &
    ts_fit_xgbt$results$max_depth == best_max_depth &
    ts_fit_xgbt$results$eta == best_eta &
    ts_fit_xgbt$results$gamma == best_gamma &
    ts_fit_xgbt$results$colsample_bytree == best_colsample_bytree &
    ts_fit_xgbt$results$min_child_weight == best_min_child_weight &
    ts_fit_xgbt$results$subsample == best_subsample
  , ]
# Predict and evaluate against the ts_test data
ts_y_hat_xgbt <- predict(ts_fit_xgbt, newdata = ts_test)
rsq_ts_xgbt <- cor(ts_y_hat_xgbt, ts_test$sqrt_tensile_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_xgbt_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_xgbt)
ggplot(data = ts_xgbt_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_xgbt)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: xgbTrees Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ts_results_table <- ts_results_table %>% rbind(c("eXtreme Gradient Boosted Trees", best_results$Rsquared, rsq_ts_xgbt))
format_df(ts_results_table)



# Neural network (nnet)
# Train with the ts_train data and predict against the ts_test data
ts_fit_nnet <- train(sqrt_tensile_strength ~ ., 
                             method = "nnet", 
                             data = ts_train,
                             linout = TRUE,
                             preProcess = c("center", "scale"),
                             trace = FALSE)
# After tuning, best size is 3
best_size <- ts_fit_nnet$bestTune$size
best_size
# After tuning, best decay is 0
best_decay <- ts_fit_nnet$bestTune$decay
best_decay
best_results <- ts_fit_nnet$results[
  ts_fit_nnet$results$size == best_size &
    ts_fit_nnet$results$decay == best_decay
  , ]
# Predict and evaluate against the ts_test data
ts_y_hat_nnet <- predict(ts_fit_nnet, newdata = ts_test)
rsq_ts_nnet <- cor(ts_y_hat_nnet, ts_test$sqrt_tensile_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_nnet_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_nnet)
ggplot(data = ts_nnet_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_nnet)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: Neural Network Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ts_results_table <- ts_results_table %>% rbind(c("Neural Network", best_results$Rsquared, rsq_ts_nnet))
format_df(ts_results_table)



# Ensemble 1: All models tested
# Define training control for the base models
# Use 'method = "repeatedcv"' for robust out-of-fold predictions
ts_alg_list1 <- c("glm", "knn", "rpart", "rf", "xgbTree", "nnet")
my_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  savePredictions = "final"
)
# Train an ensemble with the ts_train data
ts_model_list1 <- caretEnsemble::caretList(
  sqrt_tensile_strength ~ .,
  data = ts_train,
  trControl = my_control,
  methodList = ts_alg_list1
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
ts_fit_ensemble1 <- caretEnsemble::caretStack(ts_model_list1, method = "glm")
# Summarize the meta-model
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
summary(ts_fit_ensemble1)



# Ensemble 2: Downselect to tree-based models with importance >0.10
# Train a tree-based ensemble with the ts_train data
ts_alg_list2 <- c("rf", "xgbTree")
ts_model_list2 <- caretEnsemble::caretList(
  sqrt_tensile_strength ~ .,
  data = ts_train,
  trControl = my_control,
  methodList = ts_alg_list2
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
ts_fit_ensemble2 <- caretEnsemble::caretStack(ts_model_list2, method = "glm")
summary(ts_fit_ensemble2)
# Predict and evaluate against the ts_test data
ts_y_hat_ensemble2 <- predict(ts_fit_ensemble2, newdata = ts_test)
ts_y_hat_ensemble2 <- as.vector(ts_y_hat_ensemble2$pred)
rsq_ts_ensemble2 <- cor(ts_y_hat_ensemble2, ts_test$sqrt_tensile_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ts_ensemble2_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_ensemble2)
ggplot(data = ts_ensemble2_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_ensemble2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Tensile Strength: Tree-Based Ensemble Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ts_results_table <- ts_results_table %>% rbind(c("Tree-Based Ensemble", ts_fit_ensemble2$ens_model$results$Rsquared, rsq_ts_ensemble2))
format_df(ts_results_table)



## [Optional] Ensemble 3: Downselect to tree-based models with importance >0.05 (including KNN)
## Train a tree-based ensemble with the ts_train data
#ts_alg_list3 <- c("knn", "rf", "xgbTree")
#ts_model_list3 <- caretEnsemble::caretList(
#  sqrt_tensile_strength ~ .,
#  data = ts_train,
#  trControl = my_control,
#  methodList = ts_alg_list3
#)
## Define the meta-model (e.g., a linear regression using "glm")
## The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
#ts_fit_ensemble3 <- caretEnsemble::caretStack(ts_model_list3, method = "glm")
#summary(ts_fit_ensemble3)
## Predict and evaluate against the ts_test data
#ts_y_hat_ensemble3 <- predict(ts_fit_ensemble3, newdata = ts_test)
#ts_y_hat_ensemble3 <- as.vector(ts_y_hat_ensemble3$pred)
#rsq_ts_ensemble3 <- cor(ts_y_hat_ensemble3, ts_test$sqrt_tensile_strength)^2
## Plot the predicted vs actual test values to visualize how close they are to being equal
#ts_ensemble3_plot <- data.frame(ts_test$sqrt_tensile_strength, ts_y_hat_ensemble3)
#ggplot(data = ts_ensemble3_plot, aes(x = ts_test$sqrt_tensile_strength, y = ts_y_hat_ensemble3)) +
#  geom_point() +
#  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
#  labs(
#    x = "Actual Values (Square Root)",
#    y = "Predicted Values (Square Root)",
#    title = "Tensile Strength: KNN + Tree-Based Ensemble Accuracy"
#  ) +
#  coord_equal() +
#  xlim(0, 40) + 
#  ylim(0, 40)
## Add new model Rsq to results table for comparison
#ts_results_table <- ts_results_table %>% rbind(c("KNN + Tree-Based Ensemble", ts_fit_ensemble3$ens_model$results$Rsquared, rsq_ts_ensemble3))
#format_df(ts_results_table)


# Best model for tensile strength: tree-based ensemble




##########################################################
# Yield strength model
##########################################################

# Generalized linear model (GLM)
# Train with the ys_train data and predict against the ys_test data
ys_fit_glm <- train(sqrt_yield_strength ~ ., 
                    method = "glm", 
                    data = ys_train)
ys_y_hat_glm <- predict(ys_fit_glm, newdata = ys_test)
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_glm_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_glm)
ggplot(data = ys_glm_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_glm)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: GLM Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Seed Model and Rsquared Results table with GLM results in a rounded format
ys_results_table <- data.frame(
  YieldStrengthModel = "GLM", 
  TrainRsquared = ys_fit_glm$results$Rsquared, # Rsq from train function
  TestRsquared = cor(ys_y_hat_glm, ys_test$sqrt_yield_strength)^2) # Rsq from predicted vs actual on test set
format_df(ys_results_table) 



# K-nearest neighbors (KNN)
# Train with the ys_train data
ys_fit_knn <- train(sqrt_yield_strength ~ ., 
                    method = "knn", 
                    data = ys_train,
                    tuneGrid = expand.grid(k = seq(3, 21, by = 1)))
# After tuning, best k is 3
best_k <- ys_fit_knn$bestTune$k
best_k
best_results <- ys_fit_knn$results[ys_fit_knn$results$k == best_k, ]
# Predict and evaluate against the ys_test data
ys_y_hat_knn <- predict(ys_fit_knn, newdata = ys_test)
rsq_ys_knn <- cor(ys_y_hat_knn, ys_test$sqrt_yield_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_knn_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_knn)
ggplot(data = ys_knn_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_knn)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: KNN Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ys_results_table <- ys_results_table %>% rbind(c("KNN", best_results$Rsquared, rsq_ys_knn))
format_df(ys_results_table)



# Regression trees (rpart)
# Train with the ys_train data
ys_fit_rpart <- train(sqrt_yield_strength ~ ., 
                      method = "rpart", 
                      data = ys_train,
                      tuneGrid = expand.grid(cp = seq(0, 0.05, len = 25)))
# After tuning, best cp is 0
best_cp <- ys_fit_rpart$bestTune$cp
best_cp
best_results <- ys_fit_rpart$results[ys_fit_rpart$results$cp == best_cp, ]
# Predict and evaluate against the ys_test data
ys_y_hat_rpart <- predict(ys_fit_rpart, newdata = ys_test)
rsq_ys_rpart <- cor(ys_y_hat_rpart, ys_test$sqrt_yield_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_rpart_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_rpart)
ggplot(data = ys_rpart_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_rpart)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: Regression Trees Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ys_results_table <- ys_results_table %>% rbind(c("Regression Trees", best_results$Rsquared, rsq_ys_rpart))
format_df(ys_results_table)



# Random Forest (rf)
# Train with the ys_train data and predict against the ys_test data
ys_fit_rf <- train(sqrt_yield_strength ~ ., 
                   method = "rf", 
                   data = ys_train,
                   tuneGrid = expand.grid(.mtry = c(2, 9, 16, 23, 31)))
# After tuning, best mtry is 16
best_mtry <- ys_fit_rf$bestTune$mtry
best_mtry
best_results <- ys_fit_rf$results[ys_fit_rf$results$mtry == best_mtry, ]
# Predict and evaluate against the ys_test data
ys_y_hat_rf <- predict(ys_fit_rf, newdata = ys_test)
rsq_ys_rf <- cor(ys_y_hat_rf, ys_test$sqrt_yield_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_rf_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_rf)
ggplot(data = ys_rf_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_rf)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: Random Forest Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ys_results_table <- ys_results_table %>% rbind(c("Random Forest", best_results$Rsquared, rsq_ys_rf))
format_df(ys_results_table)



# eXtreme Gradient Boosted Trees (xgbTree)
# Train with the ys_train data and predict against the ys_test data
ys_fit_xgbt <- train(sqrt_yield_strength ~ ., 
                     method = "xgbTree", 
                     data = ys_train)
# After tuning, best nrounds is 150
best_nrounds <- ys_fit_xgbt$bestTune$nrounds 
best_nrounds
# After tuning, best max_depth is 3
best_max_depth <- ys_fit_xgbt$bestTune$max_depth 
best_max_depth
# After tuning, best eta is 0.3
best_eta <- ys_fit_xgbt$bestTune$eta 
best_eta
# After tuning, best gamma is 0
best_gamma <- ys_fit_xgbt$bestTune$gamma 
best_gamma
# After tuning, best colsample_bytree is 0.8
best_colsample_bytree  <- ys_fit_xgbt$bestTune$colsample_bytree 
best_colsample_bytree
# After tuning, best min_child_weight is 1
best_min_child_weight  <- ys_fit_xgbt$bestTune$min_child_weight 
best_min_child_weight
# After tuning, best subsample is 1
best_subsample <- ys_fit_xgbt$bestTune$subsample
best_subsample 
best_results <- ys_fit_xgbt$results[
  ys_fit_xgbt$results$nrounds == best_nrounds &
    ys_fit_xgbt$results$max_depth == best_max_depth &
    ys_fit_xgbt$results$eta == best_eta &
    ys_fit_xgbt$results$gamma == best_gamma &
    ys_fit_xgbt$results$colsample_bytree == best_colsample_bytree &
    ys_fit_xgbt$results$min_child_weight == best_min_child_weight &
    ys_fit_xgbt$results$subsample == best_subsample
  , ]
# Predict and evaluate against the ys_test data
ys_y_hat_xgbt <- predict(ys_fit_xgbt, newdata = ys_test)
rsq_ys_xgbt <- cor(ys_y_hat_xgbt, ys_test$sqrt_yield_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_xgbt_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_xgbt)
ggplot(data = ys_xgbt_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_xgbt)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: xgbTrees Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ys_results_table <- ys_results_table %>% rbind(c("eXtreme Gradient Boosted Trees", best_results$Rsquared, rsq_ys_xgbt))
format_df(ys_results_table)



# Neural network (nnet)
# Train with the ys_train data and predict against the ys_test data
ys_fit_nnet <- train(sqrt_yield_strength ~ ., 
                     method = "nnet", 
                     data = ys_train,
                     linout = TRUE,
                     preProcess = c("center", "scale"),
                     trace = FALSE)
# After tuning, best size is 3
best_size <- ys_fit_nnet$bestTune$size
best_size
# After tuning, best decay is 0
best_decay <- ys_fit_nnet$bestTune$decay
best_decay
best_results <- ys_fit_nnet$results[
  ys_fit_nnet$results$size == best_size &
    ys_fit_nnet$results$decay == best_decay
  , ]
# Predict and evaluate against the ys_test data
ys_y_hat_nnet <- predict(ys_fit_nnet, newdata = ys_test)
rsq_ys_nnet <- cor(ys_y_hat_nnet, ys_test$sqrt_yield_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_nnet_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_nnet)
ggplot(data = ys_nnet_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_nnet)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: Neural Network Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ys_results_table <- ys_results_table %>% rbind(c("Neural Network", best_results$Rsquared, rsq_ys_nnet))
format_df(ys_results_table)



# Ensemble 1: All models tested
# Define training control for the base models
# Use 'method = "repeatedcv"' for robust out-of-fold predictions
ys_alg_list1 <- c("glm", "knn", "rpart", "rf", "xgbTree", "nnet")
my_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  savePredictions = "final"
)
# Train an ensemble with the ys_train data
ys_model_list1 <- caretEnsemble::caretList(
  sqrt_yield_strength ~ .,
  data = ys_train,
  trControl = my_control,
  methodList = ys_alg_list1
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
ys_fit_ensemble1 <- caretEnsemble::caretStack(ys_model_list1, method = "glm")
# Summarize the meta-model
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
summary(ys_fit_ensemble1)



# Ensemble 2: Downselect to tree-based models with importance >0.05
# Train a tree-based ensemble with the ys_train data
ys_alg_list2 <- c("rf", "xgbTree")
ys_model_list2 <- caretEnsemble::caretList(
  sqrt_yield_strength ~ .,
  data = ys_train,
  trControl = my_control,
  methodList = ys_alg_list2
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
ys_fit_ensemble2 <- caretEnsemble::caretStack(ys_model_list2, method = "glm")
summary(ys_fit_ensemble2)
# Predict and evaluate against the ys_test data
ys_y_hat_ensemble2 <- predict(ys_fit_ensemble2, newdata = ys_test)
ys_y_hat_ensemble2 <- as.vector(ys_y_hat_ensemble2$pred)
rsq_ys_ensemble2 <- cor(ys_y_hat_ensemble2, ys_test$sqrt_yield_strength)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
ys_ensemble2_plot <- data.frame(ys_test$sqrt_yield_strength, ys_y_hat_ensemble2)
ggplot(data = ys_ensemble2_plot, aes(x = ys_test$sqrt_yield_strength, y = ys_y_hat_ensemble2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Yield Strength: Tree-Based Ensemble Accuracy"
  ) +
  coord_equal() +
  xlim(0, 40) + 
  ylim(0, 40)
# Add new model Rsq to results table for comparison
ys_results_table <- ys_results_table %>% rbind(c("Tree-Based Ensemble", ys_fit_ensemble2$ens_model$results$Rsquared, rsq_ys_ensemble2))
format_df(ys_results_table)


# Best model for yield strength: tree-based ensemble




##########################################################
# Elongation model
##########################################################

# Generalized linear model (GLM)
# Train with the e_train data and predict against the e_test data
e_fit_glm <- train(sqrt_elongation ~ ., 
                    method = "glm", 
                    data = e_train)
e_y_hat_glm <- predict(e_fit_glm, newdata = e_test)
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_glm_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_glm)
ggplot(data = e_glm_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_glm)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: GLM Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)

# Seed Model and Rsquared/RMSE Results table with GLM results in a rounded format
e_results_table <- data.frame(
  ElongationModel = "GLM", 
  TrainRsquared = e_fit_glm$results$Rsquared, # Rsq from train function
  TestRsquared = cor(e_y_hat_glm, e_test$sqrt_elongation)^2) # Rsq from predicted vs actual on test set
format_df(e_results_table) 



# K-nearest neighbors (KNN)
# Train with the e_train data
e_fit_knn <- train(sqrt_elongation ~ ., 
                    method = "knn", 
                    data = e_train,
                    tuneGrid = expand.grid(k = seq(3, 21, by = 1)))
# After tuning, best k is 3
best_k <- e_fit_knn$bestTune$k
best_k
best_results <- e_fit_knn$results[e_fit_knn$results$k == best_k, ]
# Predict and evaluate against the e_test data
e_y_hat_knn <- predict(e_fit_knn, newdata = e_test)
rsq_e_knn <- cor(e_y_hat_knn, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_knn_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_knn)
ggplot(data = e_knn_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_knn)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: KNN Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("KNN", best_results$Rsquared, rsq_e_knn))
format_df(e_results_table)



# Regression trees (rpart)
# Train with the e_train data
e_fit_rpart <- train(sqrt_elongation ~ ., 
                      method = "rpart", 
                      data = e_train,
                      tuneGrid = expand.grid(cp = seq(0, 0.05, len = 25)))
# After tuning, best cp is 0
best_cp <- e_fit_rpart$bestTune$cp
best_cp
best_results <- e_fit_rpart$results[e_fit_rpart$results$cp == best_cp, ]
# Predict and evaluate against the e_test data
e_y_hat_rpart <- predict(e_fit_rpart, newdata = e_test)
rsq_e_rpart <- cor(e_y_hat_rpart, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_rpart_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_rpart)
ggplot(data = e_rpart_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_rpart)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: Regression Trees Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("Regression Trees", best_results$Rsquared, rsq_e_rpart))
format_df(e_results_table)



# Random Forest (rf)
# Train with the e_train data and predict against the e_test data
e_fit_rf <- train(sqrt_elongation ~ ., 
                   method = "rf", 
                   data = e_train,
                   tuneGrid = expand.grid(.mtry = c(2, 9, 16, 23, 31)))
# After tuning, best mtry is 16
best_mtry <- e_fit_rf$bestTune$mtry
best_mtry
best_results <- e_fit_rf$results[e_fit_rf$results$mtry == best_mtry, ]
# Predict and evaluate against the e_test data
e_y_hat_rf <- predict(e_fit_rf, newdata = e_test)
rsq_e_rf <- cor(e_y_hat_rf, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_rf_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_rf)
ggplot(data = e_rf_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_rf)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: Random Forest Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("Random Forest", best_results$Rsquared, rsq_e_rf))
format_df(e_results_table)



# eXtreme Gradient Boosted Trees (xgbTree)
# Train with the e_train data and predict against the e_test data
e_fit_xgbt <- train(sqrt_elongation ~ ., 
                     method = "xgbTree", 
                     data = e_train)
# After tuning, best nrounds is 150
best_nrounds <- e_fit_xgbt$bestTune$nrounds 
best_nrounds
# After tuning, best max_depth is 3
best_max_depth <- e_fit_xgbt$bestTune$max_depth 
best_max_depth
# After tuning, best eta is 0.3
best_eta <- e_fit_xgbt$bestTune$eta 
best_eta
# After tuning, best gamma is 0
best_gamma <- e_fit_xgbt$bestTune$gamma 
best_gamma
# After tuning, best colsample_bytree is 0.8
best_colsample_bytree  <- e_fit_xgbt$bestTune$colsample_bytree 
best_colsample_bytree
# After tuning, best min_child_weight is 1
best_min_child_weight  <- e_fit_xgbt$bestTune$min_child_weight 
best_min_child_weight
# After tuning, best subsample is 1
best_subsample <- e_fit_xgbt$bestTune$subsample
best_subsample 
best_results <- e_fit_xgbt$results[
  e_fit_xgbt$results$nrounds == best_nrounds &
    e_fit_xgbt$results$max_depth == best_max_depth &
    e_fit_xgbt$results$eta == best_eta &
    e_fit_xgbt$results$gamma == best_gamma &
    e_fit_xgbt$results$colsample_bytree == best_colsample_bytree &
    e_fit_xgbt$results$min_child_weight == best_min_child_weight &
    e_fit_xgbt$results$subsample == best_subsample
  , ]
# Predict and evaluate against the e_test data
e_y_hat_xgbt <- predict(e_fit_xgbt, newdata = e_test)
rsq_e_xgbt <- cor(e_y_hat_xgbt, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_xgbt_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_xgbt)
ggplot(data = e_xgbt_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_xgbt)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: xgbTrees Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("eXtreme Gradient Boosted Trees", best_results$Rsquared, rsq_e_xgbt))
format_df(e_results_table)



# Neural network (nnet)
# Train with the e_train data and predict against the e_test data
e_fit_nnet <- train(sqrt_elongation ~ ., 
                     method = "nnet", 
                     data = e_train,
                     linout = TRUE,
                     preProcess = c("center", "scale"),
                     trace = FALSE)
# After tuning, best size is 3
best_size <- e_fit_nnet$bestTune$size
best_size
# After tuning, best decay is 0
best_decay <- e_fit_nnet$bestTune$decay
best_decay
best_results <- e_fit_nnet$results[
  e_fit_nnet$results$size == best_size &
    e_fit_nnet$results$decay == best_decay
  , ]
# Predict and evaluate against the e_test data
e_y_hat_nnet <- predict(e_fit_nnet, newdata = e_test)
rsq_e_nnet <- cor(e_y_hat_nnet, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_nnet_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_nnet)
ggplot(data = e_nnet_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_nnet)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: Neural Network Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("Neural Network", best_results$Rsquared, rsq_e_nnet))
format_df(e_results_table)



# Ensemble 1: All models tested
# Define training control for the base models
# Use 'method = "repeatedcv"' for robust out-of-fold predictions
e_alg_list1 <- c("glm", "knn", "rpart", "rf", "xgbTree", "nnet")
my_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  savePredictions = "final"
)
# Train an ensemble with the e_train data
e_model_list1 <- caretEnsemble::caretList(
  sqrt_elongation ~ .,
  data = e_train,
  trControl = my_control,
  methodList = e_alg_list1
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
e_fit_ensemble1 <- caretEnsemble::caretStack(e_model_list1, method = "glm")
# Summarize the meta-model
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
summary(e_fit_ensemble1)



# Ensemble 2: Downselect to tree-based models with importance >0.20
# Train a tree-based ensemble with the e_train data
e_alg_list2 <- c("rf", "xgbTree")
e_model_list2 <- caretEnsemble::caretList(
  sqrt_elongation ~ .,
  data = e_train,
  trControl = my_control,
  methodList = e_alg_list2
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
e_fit_ensemble2 <- caretEnsemble::caretStack(e_model_list2, method = "glm")
summary(e_fit_ensemble2)
# Predict and evaluate against the e_test data
e_y_hat_ensemble2 <- predict(e_fit_ensemble2, newdata = e_test)
e_y_hat_ensemble2 <- as.vector(e_y_hat_ensemble2$pred)
rsq_e_ensemble2 <- cor(e_y_hat_ensemble2, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_ensemble2_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_ensemble2)
ggplot(data = e_ensemble2_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_ensemble2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: Tree-Based Ensemble Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("Tree-Based Ensemble", e_fit_ensemble2$ens_model$results$Rsquared, rsq_e_ensemble2))
format_df(e_results_table)



# Ensemble 3: Downselect to tree-based models with importance >0.10 (including KNN)
# Train a tree-based ensemble with the e_train data
e_alg_list3 <- c("knn", "rf", "xgbTree")
e_model_list3 <- caretEnsemble::caretList(
  sqrt_elongation ~ .,
  data = e_train,
  trControl = my_control,
  methodList = e_alg_list3
)
# Define the meta-model (e.g., a linear regression using "glm")
# The meta-model uses the predictions of the base models as input features and generates the relative importance of each model used in the ensemble model
e_fit_ensemble3 <- caretEnsemble::caretStack(e_model_list3, method = "glm")
summary(e_fit_ensemble3)
# Predict and evaluate against the e_test data
e_y_hat_ensemble3 <- predict(e_fit_ensemble3, newdata = e_test)
e_y_hat_ensemble3 <- as.vector(e_y_hat_ensemble3$pred)
rsq_e_ensemble3 <- cor(e_y_hat_ensemble3, e_test$sqrt_elongation)^2
# Plot the predicted vs actual test values to visualize how close they are to being equal
e_ensemble3_plot <- data.frame(e_test$sqrt_elongation, e_y_hat_ensemble3)
ggplot(data = e_ensemble3_plot, aes(x = e_test$sqrt_elongation, y = e_y_hat_ensemble3)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    x = "Actual Values (Square Root)",
    y = "Predicted Values (Square Root)",
    title = "Elongation: KNN + Tree-Based Ensemble Accuracy"
  ) +
  coord_equal() +
  xlim(0, 10) + 
  ylim(0, 10)
# Add new model Rsq to results table for comparison
e_results_table <- e_results_table %>% rbind(c("KNN + Tree-Based Ensemble", e_fit_ensemble3$ens_model$results$Rsquared, rsq_e_ensemble3))
format_df(e_results_table)




##########################################################
# Application: fill in blanks from original table
##########################################################

# Data set that focuses only on alloys that have at least one mechanical property null/missing
# Existing mechanical property data can be used on the incomplete table to assess Rsquared and root mean square error (RMSE)
empty_alloy_data <- alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_tensile_strength, sqrt_yield_strength, sqrt_elongation) %>% 
  filter(is.na(sqrt_tensile_strength) | is.na(sqrt_yield_strength) | is.na(sqrt_elongation))

# For predicting tensile strength, use ts_fit_ensemble2 (tree-based ensemble model) on existing data in the incomplete table
# For predicting yield strength, use ys_fit_ensemble2 (tree-based ensemble model)
# For predicting elongation, use e_fit_ensemble2 (tree-based ensemble model)
sqrt_tensile_strength_pred <- predict(ts_fit_ensemble2, newdata = empty_alloy_data)
sqrt_tensile_strength_pred <- as.vector(sqrt_tensile_strength_pred$pred)
sqrt_yield_strength_pred <- predict(ys_fit_ensemble2, newdata = empty_alloy_data)
sqrt_yield_strength_pred <- as.vector(sqrt_yield_strength_pred$pred)
sqrt_elongation_pred <- predict(e_fit_ensemble2, newdata = empty_alloy_data)
sqrt_elongation_pred <- as.vector(sqrt_elongation_pred$pred)
# Add square root predictions, their squares, and the squares of the existing properties in the incomplete table
empty_alloy_data <- empty_alloy_data %>% mutate(
  sqrt_tensile_strength_pred = sqrt_tensile_strength_pred, 
  sqrt_yield_strength_pred = sqrt_yield_strength_pred, 
  sqrt_elongation_pred = sqrt_elongation_pred, 
  tensile_strength_pred = sqrt_tensile_strength_pred^2, 
  yield_strength_pred = sqrt_yield_strength_pred^2, 
  elongation_pred = sqrt_elongation_pred^2, 
  tensile_strength = sqrt_tensile_strength^2, 
  yield_strength = sqrt_yield_strength^2, 
  elongation = sqrt_elongation^2) %>%
  select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_tensile_strength, sqrt_tensile_strength_pred, tensile_strength, tensile_strength_pred, sqrt_yield_strength, sqrt_yield_strength_pred, yield_strength, yield_strength_pred, sqrt_elongation, sqrt_elongation_pred, elongation, elongation_pred)

# Filter incomplete table into one table per property containing only existing data for evaluating Rsquared and RMSE
ts_remainder_data <- empty_alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_tensile_strength, sqrt_tensile_strength_pred, tensile_strength, tensile_strength_pred) %>% 
  filter(!is.na(sqrt_tensile_strength))
ys_remainder_data <- empty_alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_yield_strength, sqrt_yield_strength_pred, yield_strength, yield_strength_pred) %>% 
  filter(!is.na(sqrt_yield_strength))
e_remainder_data <- empty_alloy_data %>% select(series, solutionised, naturally_aged, artificially_aged, strain_hardened, all_of(elements), balance, sqrt_elongation, sqrt_elongation_pred, elongation, elongation_pred) %>% 
  filter(!is.na(sqrt_elongation))

# Rsquared/RMSE Results table with available data from the incomplete table in a rounded format
incomplete_holdout_results_table <- data.frame(
  Property = c("Tensile Strength", "Yield Strength", "Elongation"),
  # Rsquared formulas
  Rsquared = c(cor(ts_remainder_data$tensile_strength, ts_remainder_data$tensile_strength_pred)^2,
               cor(ys_remainder_data$yield_strength, ys_remainder_data$yield_strength_pred)^2,
               cor(e_remainder_data$elongation, e_remainder_data$elongation_pred)^2),
  # RMSE formulas
  RMSE = c(sqrt(mean((ts_remainder_data$tensile_strength - ts_remainder_data$tensile_strength_pred)^2)),
           sqrt(mean((ys_remainder_data$yield_strength - ys_remainder_data$yield_strength_pred)^2)),
           sqrt(mean((e_remainder_data$elongation - e_remainder_data$elongation_pred)^2)))
)
format_df(incomplete_holdout_results_table)




##########################################################
# Most important variables for models
##########################################################

# Weighted contribution of random forest to tensile strength ensemble model
ts_rf_wt <- as.numeric(ts_fit_ensemble2$ens_model$finalModel$coefficients[2])
# Weighted contribution of extreme gradient boosted trees to tensile strength ensemble model
ts_xgbt_wt <- as.numeric(ts_fit_ensemble2$ens_model$finalModel$coefficients[3])
# Ranked percent importance of the variables used in random forest component of tensile strength ensemble model for making predictions
ts_rf_importance_tbl <- varImp(ts_fit_ensemble2$models$rf)$importance %>%
  rownames_to_column(var = "Variable") %>%
  as.tibble() %>% 
  rename(rf_Overall = Overall) %>%
  arrange(Variable)
# Ranked percent importance of the variables used in extreme gradient boosted trees component of tensile strength ensemble model for making predictions
ts_xgbt_importance_tbl <- varImp(ts_fit_ensemble2$models$xgbTree)$importance %>%
  rownames_to_column(var = "Variable") %>%
  as.tibble() %>% 
  rename(xgbt_Overall = Overall) %>%
  arrange(Variable)
# Merge percent importance tables, average the contributions of each variable between the model components, and create an overall percent importance column normalized to the variable with the strongest importance
ts_importance_tbl <- merge(ts_rf_importance_tbl, ts_xgbt_importance_tbl) %>% 
  mutate(Overall = (rf_Overall * ts_rf_wt) + (xgbt_Overall * ts_xgbt_wt)) %>%
  arrange(desc(Overall))
ts_importance_tbl$Overall <- (ts_importance_tbl$Overall / ts_importance_tbl[1, 4]) * 100

# Weighted contribution of random forest to yield strength ensemble model
ys_rf_wt <- as.numeric(ys_fit_ensemble2$ens_model$finalModel$coefficients[2])
# Weighted contribution of extreme gradient boosted trees to yield strength ensemble model
ys_xgbt_wt <- as.numeric(ys_fit_ensemble2$ens_model$finalModel$coefficients[3])
# Ranked percent importance of the variables used in random forest component of yield strength ensemble model for making predictions
ys_rf_importance_tbl <- varImp(ys_fit_ensemble2$models$rf)$importance %>%
  rownames_to_column(var = "Variable") %>%
  as.tibble() %>% 
  rename(rf_Overall = Overall) %>%
  arrange(Variable)
# Ranked percent importance of the variables used in extreme gradient boosted trees component of yield strength ensemble model for making predictions
ys_xgbt_importance_tbl <- varImp(ys_fit_ensemble2$models$xgbTree)$importance %>%
  rownames_to_column(var = "Variable") %>%
  as.tibble() %>% 
  rename(xgbt_Overall = Overall) %>%
  arrange(Variable)
# Merge percent importance tables, average the contributions of each variable between the model components, and create an overall percent importance column normalized to the variable with the strongest importance
ys_importance_tbl <- merge(ys_rf_importance_tbl, ys_xgbt_importance_tbl) %>% 
  mutate(Overall = (rf_Overall * ys_rf_wt) + (xgbt_Overall * ys_xgbt_wt)) %>%
  arrange(desc(Overall))
ys_importance_tbl$Overall <- (ys_importance_tbl$Overall / ys_importance_tbl[1, 4]) * 100

# Weighted contribution of random forest to elongation ensemble model
e_rf_wt <- as.numeric(e_fit_ensemble2$ens_model$finalModel$coefficients[2])
# Weighted contribution of extreme gradient boosted trees to elongation ensemble model
e_xgbt_wt <- as.numeric(e_fit_ensemble2$ens_model$finalModel$coefficients[3])
# Ranked percent importance of the variables used in random forest component of elongation ensemble model for making predictions
e_rf_importance_tbl <- varImp(e_fit_ensemble2$models$rf)$importance %>%
  rownames_to_column(var = "Variable") %>%
  as.tibble() %>% 
  rename(rf_Overall = Overall) %>%
  arrange(Variable)
# Ranked percent importance of the variables used in extreme gradient boosted trees component of elongation ensemble model for making predictions
e_xgbt_importance_tbl <- varImp(e_fit_ensemble2$models$xgbTree)$importance %>%
  rownames_to_column(var = "Variable") %>%
  as.tibble() %>% 
  rename(xgbt_Overall = Overall) %>%
  arrange(Variable)
# Merge percent importance tables, average the contributions of each variable between the model components, and create an overall percent importance column normalized to the variable with the strongest importance
e_importance_tbl <- merge(e_rf_importance_tbl, e_xgbt_importance_tbl) %>% 
  mutate(Overall = (rf_Overall * e_rf_wt) + (xgbt_Overall * e_xgbt_wt)) %>%
  arrange(desc(Overall))
e_importance_tbl$Overall <- (e_importance_tbl$Overall / e_importance_tbl[1, 4]) * 100

# Table of top ten variables and respective strengths in tensile strength model
ts_top_10_table <- data.frame(
  Rank = seq(1, 10, 1),
  Variable = head(ts_importance_tbl$Variable, 10),
  Importance = head(ts_importance_tbl$Overall, 10)
)
ts_top_10_table

# Table of top ten variables and respective strengths in yield strength model
ys_top_10_table <- data.frame(
  Rank = seq(1, 10, 1),
  Variable = head(ys_importance_tbl$Variable, 10),
  Importance = head(ys_importance_tbl$Overall, 10)
)
ys_top_10_table

# Table of top ten variables and respective strengths in elongation model
e_top_10_table <- data.frame(
  Rank = seq(1, 10, 1),
  Variable = head(e_importance_tbl$Variable, 10),
  Importance = head(e_importance_tbl$Overall, 10)
)
e_top_10_table




##########################################################
# Application of models for composition adjustment
##########################################################

# Demonstration of applying the models to the first alloy in the alloy_data table
original_alloy <- alloy_data %>% slice(1)
new_alloy <- original_alloy %>% mutate( Al = Al - 0.01, Cu = Cu + 0.01 )
demo_table <- bind_rows(original_alloy, new_alloy)
demo_sqrt_tensile_strength_pred <- predict(ts_fit_ensemble2, newdata = demo_table)
demo_sqrt_tensile_strength_pred <- as.vector(demo_sqrt_tensile_strength_pred$pred)
demo_sqrt_yield_strength_pred <- predict(ys_fit_ensemble2, newdata = demo_table)
demo_sqrt_yield_strength_pred <- as.vector(demo_sqrt_yield_strength_pred$pred)
demo_sqrt_elongation_pred <- predict(e_fit_ensemble2, newdata = demo_table)
demo_sqrt_elongation_pred <- as.vector(demo_sqrt_elongation_pred$pred)
# Comparison of alloy mechanical properties before and after decreasing aluminum by 1% and increasing copper by 1%
demo_results <- data.frame(
  Alloy = c("Original", "New"),
  AlPercent = demo_table$Al,
  CuPercent = demo_table$Cu,
  TensileStrength = demo_sqrt_tensile_strength_pred^2,
  YieldStrength = demo_sqrt_yield_strength_pred^2,
  Elongation = demo_sqrt_elongation_pred^2
)
demo_results
