---
title: "Staging MTV and IDM - RF Combined model"
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
data_nonnormal <- read_excel("data_stro_NA.xlsx", 
    col_types = c("text", "text", "numeric", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", "text", 
        "text", "text", "text", "text", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "numeric", 
        "numeric", "numeric", "text", "text", 
        "text", "text", "text", "text", "text", 
        "text"))

data_nonnormal <- mutate_if(data_nonnormal, is.character, factor)

```

## Building the model for 5-PS prediction using Ann Arbor Staging, IDM and MTV (Combined model).
```{r}
# selecting the key clinical predictors and the most important radiomic features.
data_rf <- data_nonnormal %>%
  select(CL_Staging, Semantic_MTV, Texture_IDM, Outcome_5PS_CAT)

# Training and Testing split.
set.seed(12)
data_split <- initial_split(data_rf, prop = 0.7, strata = Outcome_5PS_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds <- vfold_cv(data_tr, 10, repeats = 5, strata = Outcome_5PS_CAT)

```

## Recipes
```{r}
# Recipes:
# First, the recipe() needs to be informed about the model being used and the training data.
# Then, variables that are too correlated with each other should be filtered out with step_corr().
# Any numeric variables with zero variance are then removed with step_zv().
# In the next step, numeric variables are normalized (centered and scaled) to account for differing scales, as the model being trained is sensitive to these differences (using step_normalize()).
# After that, the factor columns are converted into one or more numeric binary variables (0 and 1) for the levels present in the training data (using step_dummy()).
# Lastly, the recipe() is prepped, meaning the steps are actually implemented with the training data, and the necessary parameters are estimated from the training data to apply this process to another dataset later.
rec_spec <- recipe(Outcome_5PS_CAT ~ ., data_tr) %>%
  step_normalize(Semantic_MTV) %>%
  step_normalize(Texture_IDM) %>%
  step_dummy(all_factor(), -all_outcomes()) 

rec_spec %>% prep() %>% juice() %>% glimpse()
```

## Model specifications.
```{r}
# XGBoosting
xgb_spec <- boost_tree(
           trees = tune(),
           learn_rate = tune(),
           min_n = tune(),
           loss_reduction = tune(),
           tree_depth = tune(),
           sample_size = tune()
       ) %>%
           set_engine("xgboost", scale_pos_weight = 3.5) %>%
           set_mode("classification")

# Random Forest
rf_spec <- rand_forest(
                  trees = 500,
                  min_n = tune()
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")


# Neural networks
nn_spec <- mlp(
  hidden_units = tune(),
  epochs = tune(),
  activation = "softmax") %>%
  set_engine("nnet") %>%
  set_mode("classification")

nn_param <- nn_spec %>% 
   extract_parameter_set_dials() %>% 
   update(hidden_units = hidden_units(c(1, 2)),
          epochs = epochs(c(20,100)))

# Support vector machine
svm_spec <- svm_rbf(cost = tune(),
                    rbf_sigma = tune()) %>% 
   set_engine("kernlab") %>% 
   set_mode("classification")

svm_param <- svm_spec %>% 
   extract_parameter_set_dials() %>% 
   update(cost = cost(c(-10, 5), trans = transform_log2()),
          rbf_sigma = rbf_sigma(c(-2, 0), trans = transform_log10()))


# Logistic regression
lr_spec <- logistic_reg(
              ) %>%
                  set_engine("glm") %>%
                  set_mode("classification")

# Naive Bayes
library(discrim)
nb_spec <- naive_Bayes(
  smoothness = tune(),
  Laplace = tune()) %>%
    set_engine("klaR") %>%
    set_mode("classification")


```

## Model metrics and control.
```{r}
# Model metrics
model_metrics <- metric_set(roc_auc, bal_accuracy)

# Model control
model_control <-
   control_race(
      save_pred = TRUE,
      parallel_over = "everything",
      save_workflow = TRUE
   )
```

## Tuning hyperparameters.
```{r}
# Creating a set of workflows.
model_set <- workflow_set(
  preproc = list(regular = rec_spec),
  models = list(xgboost = xgb_spec, rf = rf_spec, nnet = nn_spec, svm = svm_spec, logreg = lr_spec, nbayes = nb_spec),
  cross = TRUE) %>%
  option_add(param_info = nn_param, id = "regular_nnet") %>%
  option_add(param_info = svm_param, id = "regular_svm")

# Tuning models using ANOVA-based racing optimization.
model_set <- model_set %>% 
  workflow_map("tune_race_anova",
               seed = 12,
               resamples = data_folds,
               grid = 25,
               metrics = model_metrics,
               control = model_control)

# Tuned model set results
model_set

```

## Evaluating model performance.
```{r}
autoplot(model_set)

autoplot(model_set, metric = "roc_auc", select_best = TRUE) + 
  ylab("AUROC") +
  theme_pubr() +
  theme(
    legend.position = "right"
  )
```

## Ranking model performance on resamples based on AUC.
```{r}
rank_results(model_set, rank_metric = "roc_auc", select_best = T) %>% 
  filter(.metric == "roc_auc")  %>%
  mutate(conf.low = mean - (std_err*1.96),
         conf.high = (std_err*1.96) + mean) %>%
  dplyr::select(wflow_id, mean, conf.low, conf.high)

# Ranking the best models for each classifier. It is important to check not only the model with the highest AUC but also the model with the highest AUC and the lowest standard error.

#xgb.
model_set %>% 
  extract_workflow_set_result("regular_xgboost") %>%
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_xgboost") %>% 
  show_best()

#rf.
model_set %>% 
  extract_workflow_set_result("regular_rf") %>% 
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_rf") %>% 
  show_best()

#nn.
model_set %>% 
  extract_workflow_set_result("regular_nnet") %>% 
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_nnet") %>% 
  show_best()

#svm.
model_set %>% 
  extract_workflow_set_result("regular_svm") %>% 
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_svm") %>% 
  show_best()

#lr.
model_set %>%
  extract_workflow_set_result("regular_logreg") %>%
  show_best()

#nb.
model_set %>%
  extract_workflow_set_result("regular_nbayes") %>%
  collect_metrics()

model_set %>%
  extract_workflow_set_result("regular_nbayes") %>%
  show_best()

```

## Model comparison
Bayesian analysis was used here to answer the question: 'When looking at resampling results, are the differences between models real?' To address this, a model was created with the outcome being the resampling statistics (ROC_AUC).
```{r}
set.seed(12)
class_post_wfset <- perf_mod(model_set)

autoplot(class_post_wfset) + theme_pubr() + theme(legend.position = "right")

contrast_rf_nb <- contrast_models(class_post_wfset,
                                  list_1 = "regular_rf",
                                  list_2 = "regular_xgboost",
                                  seed = 12)
summary(contrast_rf_nb, size = 0.02) %>%
  dplyr::select(-size)

contrast_wfset <- contrast_models(class_post_wfset)
# Comparison visualization
library(ggrepel)
autoplot(class_post_wfset, type = "ROPE", size = 0.02) +
  geom_text_repel(aes(label = workflow)) +
  theme_pubr() +
  theme(legend.position = "none" )
  
```

## Training the best models with the best hyperparameter values.
```{r}

wf_rf <- workflow() %>% 
  add_model(rf_spec) %>% 
  add_recipe(rec_spec)

# Expanding the range of hyperparameter values for the best model.
rf_grid <- grid_regular(
  min_n(range = c(10,40)),
  levels = 20
 )

set.seed(12)
rf_res <- tune_grid(
     object = wf_rf,
     resamples = data_folds,
     grid = rf_grid,
     control = control_grid(save_pred = TRUE)
 )

rf_res %>%
  collect_metrics()

show_best(rf_res, metric = "roc_auc")

show_best(rf_res, metric = "accuracy")


```

#### Random Forest
```{r}
# Training the best model with the best hyperparameter values.

rf_model <- rand_forest(
                  trees = 500,
                  min_n = 33
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res_2 <- fit_resamples(rf_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

### ROC_AUC Random Forest
```{r}
# ROC_AUC

rf_res_2 %>%
  collect_predictions() %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()


rf_res_2 %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

rf_res_2 %>%
  collect_metrics()

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the training set - Random forest.
```{r}
library(probably)

thresholds <- seq(0.1, 0.9, by = 0.001)

threshold_data_tr <- rf_res_2 %>% 
  collect_predictions() %>%
  threshold_perf(Outcome_5PS_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_tr %>%
  arrange(.threshold)

threshold_data_tr <- threshold_data_tr %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_tr <- threshold_data_tr %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_tr

threshold_data_tr %>%
  filter(.threshold == max_j_index_threshold_tr)

threshold_rf_cl_rad_tr <- threshold_data_tr %>%
  arrange(.threshold) %>%
  filter(near(.threshold, 0.297))

```

## Confusion matrix
```{r}
pred_rf_res_2 <- rf_res_2 %>% 
  collect_predictions() %>%
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.297, 1, 0))
pred_rf_res_2$pred_class_youden <- as.factor(pred_rf_res_2$pred_class_youden)

pred_rf_res_2$Outcome_5PS_CAT <- fct_recode(pred_rf_res_2$Outcome_5PS_CAT, "1-3" = "0", "4-5" = "1")
pred_rf_res_2$pred_class_youden <- fct_recode(pred_rf_res_2$pred_class_youden, "1-3" = "0", "4-5" = "1")

pred_rf_res_2 %>%
  conf_mat(Outcome_5PS_CAT, pred_class_youden) %>%
  autoplot("heatmap") +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 12, face = "bold"))
```

## RF variable importance.
```{r}
set.seed(12)
rf_model %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Outcome_5PS_CAT ~ .,
    data = juice(prep(rec_spec))
  ) %>%
  vip(geom = "col") + theme_pubr()

```

## The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.
```{r}
#rf.
final_wf_rf <- workflow() %>% 
  add_model(rf_model) %>% 
  add_recipe(rec_spec)

set.seed(12)
final_res_rf <- last_fit(final_wf_rf, data_split) 

final_res_rf %>%
  collect_metrics

final_res_rf %>% 
  collect_predictions() %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

# Bootstrapping for 95% confidence interval.
set.seed(12)
final_res_rf_boot <- int_pctl(final_res_rf, alpha = 0.05)
final_res_rf_boot

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the testing set.
```{r}
threshold_data_te <- final_res_rf %>% 
  collect_predictions() %>%
  threshold_perf(Outcome_5PS_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_te %>%
  arrange(.threshold)

threshold_data_te <- threshold_data_te %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_te <- threshold_data_te %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_te

threshold_cl_rad_te <- threshold_data_te %>%
  arrange(.threshold) %>%
  filter(near(.threshold, 0.240) | near(.threshold, 0.297))
```

## Confusion matrix
```{r}
pred_test_res_rf <- final_res_rf %>% 
  collect_predictions() %>%
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.297, 1, 0))
pred_test_res_rf$pred_class_youden <- as.factor(pred_test_res_rf$pred_class_youden)

pred_test_res_rf$Outcome_5PS_CAT <- fct_recode(pred_test_res_rf$Outcome_5PS_CAT, "1-3" = "0", "4-5" = "1")
pred_test_res_rf$pred_class_youden <- fct_recode(pred_test_res_rf$pred_class_youden, "1-3" = "0", "4-5" = "1")

pred_test_res_rf %>%
  conf_mat(Outcome_5PS_CAT, pred_class_youden) %>%
  autoplot("heatmap") +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 12, face = "bold"))

```
