---
title: "Combined x radiomic model"
author: "Farid"
date: "`r Sys.Date()`"
output: html_document
---
## Loading packages and dataset
```{r}
library(tidyverse)
library(tidymodels)
library(ggpubr)
library(themis)
library(stacks)
library(finetune)
library(vip)
library(tidyposterior)
library(modeldata)
options(tidymodels.dark = TRUE)

library(readxl)
data_nonnormal <- read_excel("data_stro.xlsx", 
    col_types = c("text", "numeric", "text", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "text", "text", 
        "text", "text", "text"))

data_nonnormal <- mutate_if(data_nonnormal, is.character, factor)
```

## Comparing combined x radiomic model
```{r}
# selecting the key clinical predictors and the most important radiomic features.
data_comb <- data_nonnormal %>%
  select(estadiamento_cat, MTV, IDM, Esc_CAT)

# selecting the most important radiomic features.
data_rad <- data_nonnormal %>%
  select(MTV, IDM, Esc_CAT)

# Training and Testing split - combined model.
set.seed(12)
data_split_comb <- initial_split(data_comb, prop = 0.7, strata = Esc_CAT)
data_tr_comb <- training(data_split_comb)
data_te_comb <- testing(data_split_comb)

# Training and Testing split - radiomic model.
set.seed(12)
data_split_rad <- initial_split(data_rad, prop = 0.7, strata = Esc_CAT)
data_tr_rad <- training(data_split_rad)
data_te_rad <- testing(data_split_rad)

# Cross-validation folds - combined model.
set.seed(12)
data_folds_comb <- vfold_cv(data_tr_comb, 5, strata = Esc_CAT)

# Cross-validation folds - radiomic model.
set.seed(12)
data_folds_rad <- vfold_cv(data_tr_rad, 5, strata = Esc_CAT)

```

## Recipes
```{r}
# Recipes:
# First, the recipe() needs to be informed about the model being used and the training data.
# Then, variables that are too correlated with each other should be filtered out with step_corr().
# Any numeric variables with zero variance are then removed with step_zv().
# In the next step, numeric variables are normalized (centered and scaled) to account for differing scales, as the model being trained is sensitive to these differences (using step_normalize()).
# Lastly, the recipe() is prepped, meaning the steps are actually implemented with the training data, and the necessary parameters are estimated from the training data to apply this process to another dataset later.
rec_spec_comb <- recipe(Esc_CAT ~ ., data_tr_comb) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  step_dummy(all_factor(), -all_outcomes()) 

rec_spec_comb %>% prep() %>% juice() %>% glimpse()

rec_spec_rad <- recipe(Esc_CAT ~ ., data_tr_rad) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric())

rec_spec_rad %>% prep() %>% juice() %>% glimpse()

```

## Training the best model with the best hyperparameter values - combined.
```{r}
rf_model_comb <- rand_forest(
                  trees = 500,
                  min_n = 27
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res_comb <- fit_resamples(rf_model_comb,
                                   rec_spec_comb,
                                   resamples = data_folds_comb,
                                   control = control_resamples(save_pred = T))

```

## Training the best model with the best hyperparameter values - radiomic.
```{r}
rf_model_rad <- rand_forest(
                  trees = 500,
                  min_n = 2
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res_rad <- fit_resamples(rf_model_rad,
                                   rec_spec_rad,
                                   resamples = data_folds_rad,
                                   control = control_resamples(save_pred = T))

```

## Combined ROC_AUC
```{r}
rf_res_comb %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res_comb %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res_comb %>%
  collect_metrics()

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the training set - combined.
```{r}
library(probably)

thresholds <- seq(0.1, 0.9, by = 0.1)

threshold_data_tr_comb <- rf_res_comb %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_tr_comb

threshold_data_tr_comb <- threshold_data_tr_comb %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_tr_comb <- threshold_data_tr_comb %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_tr_comb

threshold_data_tr_comb %>%
  filter(.threshold == max_j_index_threshold_tr_comb)

```

## Radiomic ROC_AUC
```{r}
rf_res_rad %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res_rad %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res_rad %>%
  collect_metrics()

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the training set - Radiomic.
```{r}
threshold_data_tr_rad <- rf_res_rad %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_tr_rad

threshold_data_tr_rad <- threshold_data_tr_rad %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_tr_rad <- threshold_data_tr_rad %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_tr_rad

threshold_data_tr_rad %>%
  filter(.threshold == max_j_index_threshold_tr_rad)

```

## Mcnemar test: training set - comparing sensitivity
```{r}
pred_comb <- rf_res_comb %>% 
  collect_predictions() %>% 
  mutate(id = row_number()) %>%
  select(id, Esc_CAT, pred_comb = .pred_class)
  

pred_rad <- rf_res_rad %>% 
  collect_predictions() %>% 
  mutate(id = row_number()) %>%
  select(id, Esc_CAT, pred_rad = .pred_class)
  

combined_preds <- pred_comb %>%
  inner_join(pred_rad, by = c("id", "Esc_CAT"))

# filter comb = FN & rad = TP
fn_comb_not_rad <- combined_preds %>%
  dplyr::filter(Esc_CAT == 1, pred_comb == 0, pred_rad == 1)
# filter rad = FN & comb = TP
fn_rad_not_comb <- combined_preds %>%
  dplyr::filter(Esc_CAT == 1, pred_comb == 1, pred_rad == 0)

fn_comb_not_rad
fn_rad_not_comb

# Comparing sensitivity: combined x radiomic
b_tr <- 4 # number of positive instances in the training set misclassified as FN by the combined model but not by the radiomic model
c_tr <- 0 # number of positive instances in the training set misclassified as FN by the radiomic model but not by the combined model

# Contingency matrix
contingency_matrix <- matrix(c(0, b_tr, c_tr, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix, correct = TRUE) # continuity correction

```

## Mcnemar test: training set - comparing specificity
```{r}
# filter comb = FP and rad = TN
fp_comb_not_rad <- combined_preds %>%
  dplyr::filter(Esc_CAT == 0, pred_comb == 1, pred_rad == 0)
# filter rad = FP and comb = TN
fp_rad_not_comb <- combined_preds %>%
  dplyr::filter(Esc_CAT == 0, pred_comb == 0, pred_rad == 1)

fp_comb_not_rad
fp_rad_not_comb


# Comparing specificity: combined x radiomic
d_tr <- 0 # number of negative instances in the training set misclassified as FP by the combined model but not by the radiomic model
e_tr <- 3 # number of negative instances in the training set misclassified as FP by the radiomic model but not by the combined model

# Contingency matrix
contingency_matrix <- matrix(c(0, d_tr, e_tr, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix, correct = TRUE) # continuity correction

```

## Model comparison
Bayesian analysis was used here to answer the question: 'When looking at resampling results, are the differences between models real?' To address this, a model was created with the outcome being the resampling statistics (ROC_AUC).
```{r}
# comb x rad
model_comp <- inner_join(
  rf_res_comb %>% 
    collect_metrics(summarize = FALSE) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(id, comb = .estimate),
  rf_res_rad %>% 
    collect_metrics(summarize = FALSE) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(id, rad = .estimate)
)

set.seed(12)
class_post <- perf_mod(model_comp)
class_post

# Comparison visualization
autoplot(class_post) +
  theme_pubr() +
  theme(legend.position = "right") +
  scale_color_manual(values = c("comb" = "#FF8282", "rad" = "skyblue4"))

autoplot(contrast_models(class_post)) + 
  xlab("AUC difference (combined - radiomic)") +
  theme_pubr()
```

## The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.
```{r}
#rf.
final_wf_comb <- workflow() %>% 
  add_model(rf_model_comb) %>% 
  add_recipe(rec_spec_comb)

final_res_comb <- last_fit(final_wf_comb, data_split_comb) 

final_res_comb %>%
  collect_metrics

final_res_comb %>% 
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

set.seed(12)
final_res_comb_boot <- int_pctl(final_res_comb, alpha = 0.05)
final_res_comb_boot

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the testing set.
```{r}
threshold_data_te_comb <- final_res_comb %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_te_comb 

threshold_data_te_comb <- threshold_data_te_comb %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_te_comb <- threshold_data_te_comb %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_te_comb

threshold_data_te_comb %>%
  filter(.threshold == max_j_index_threshold_te_comb)

```

## The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.
```{r}
#rf.
final_wf_rad <- workflow() %>% 
  add_model(rf_model_rad) %>% 
  add_recipe(rec_spec_rad)

final_res_rad <- last_fit(final_wf_rad, data_split_comb) 

final_res_rad %>%
  collect_metrics

final_res_rad %>% 
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

set.seed(12)
final_res_rad_boot <- int_pctl(final_res_rad, alpha = 0.05)
final_res_rad_boot

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the testing set.
```{r}
threshold_data_te_rad <- final_res_rad %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_te_rad 

threshold_data_te_rad <- threshold_data_te_rad %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_te_rad <- threshold_data_te_rad %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_te_rad

threshold_data_te_rad %>%
  filter(.threshold == max_j_index_threshold_te_rad)

```

## DeLong´s test for AUC comparison: test set
```{r}
library(pROC)
# Extraindo previsões para o modelo radiômico
pred_rad <- final_res_rad %>% 
  collect_predictions() 

# Extraindo previsões para o modelo combinado
pred_comb <- final_res_comb %>% 
  collect_predictions() 

# Calculando as curvas ROC para ambos os modelos
roc_rad <- roc(pred_rad$Esc_CAT, pred_rad$.pred_1, levels = c("0", "1"))
roc_comb <- roc(pred_comb$Esc_CAT, pred_comb$.pred_1, levels = c("0", "1"),
                direction = "<")

# Realizando o teste DeLong para comparar as AUCs
set.seed(12)
test_result <- roc.test(roc_rad, roc_comb, method = "delong")

# Exibindo o resultado do teste
print(test_result)

```

## Mcnemar´s test: test set - comparing sensitivity
```{r}
# Combining models predictions
pred_comb_te <- final_res_comb %>% 
  collect_predictions() %>% 
  mutate(id = row_number()) %>%
  select(id, Esc_CAT, pred_comb_te = .pred_class)
  

pred_rad_te <- final_res_rad %>% 
  collect_predictions() %>% 
  mutate(id = row_number()) %>%
  select(id, Esc_CAT, pred_rad_te = .pred_class)
  

combined_preds_te <- pred_comb_te %>%
  inner_join(pred_rad_te, by = c("id", "Esc_CAT"))

# filter comb = FN & rad = TP
fn_comb_not_rad_te <- combined_preds_te %>%
  dplyr::filter(Esc_CAT == 1, pred_comb_te == 0, pred_rad_te == 1)
# filter rad = FN & comb = TP
fn_rad_not_comb_te <- combined_preds_te %>%
  dplyr::filter(Esc_CAT == 1, pred_comb_te == 1, pred_rad_te == 0)

fn_comb_not_rad_te
fn_rad_not_comb_te

# Comparing sensitivity: combined x radiomic
b <- 2 # number of positive instances in the test set misclassified as FN by the combined model but not by the radiomic model 
c <- 0 # number of positive instances in the test set misclassified as FN by the radiomic model but not by the combined model 

# Contingency matrix
contingency_matrix <- matrix(c(0, b, c, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix, correct = TRUE) # continuity correction

```

## Mcnemar´s test: test set - comparing specificity
```{r}

# filter comb = FP & rad = TN
fp_comb_not_rad_te <- combined_preds_te %>%
  dplyr::filter(Esc_CAT == 0, pred_comb_te == 1, pred_rad_te == 0)
# filter rad = FP & comb = TN
fp_rad_not_comb_te <- combined_preds_te %>%
  dplyr::filter(Esc_CAT == 0, pred_comb_te == 0, pred_rad_te == 1)

fp_comb_not_rad_te
fp_rad_not_comb_te

# Comparing specificity: combined x radiomic
d <- 0 # number of negative instances in the test set misclassified as FP by the combined model but not by the radiomic model
e <- 1 # number of negative instances in the test set misclassified as FP by the radiomic model but not by the combined model

# Contingency matrix
contingency_matrix <- matrix(c(0, d, e, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix, correct = TRUE) # continuity correction

```

