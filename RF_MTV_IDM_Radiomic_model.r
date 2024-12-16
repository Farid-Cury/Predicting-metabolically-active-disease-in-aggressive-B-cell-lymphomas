---
title: "RF MTV & IDM - Radiomic model"
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
library(stacks)
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

# Changing the variables of type "character" to "factor".
data_nonnormal <- mutate_if(data_nonnormal, is.character, factor)
```

## Building the model for 5-PS prediction using IDM and MTV (Radiomic model).
```{r}
# selecting the most important radiomic features.
data_rf <- data_nonnormal %>%
  select(MTV, IDM, Esc_CAT)

# Training and Testing split.
set.seed(12)
data_split <- initial_split(data_rf, prop = 0.7, strata = Esc_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds <- vfold_cv(data_tr, 5, strata = Esc_CAT)

```
## Recipes
```{r}
# Recipes:
# First, the recipe() needs to be informed about the model being used and the training data.
# Then, variables that are too correlated with each other should be filtered out with step_corr().
# Any numeric variables with zero variance are then removed with step_zv().
# In the next step, numeric variables are normalized (centered and scaled) to account for differing scales, as the model being trained is sensitive to these differences (using step_normalize()).
# Lastly, the recipe() is prepped, meaning the steps are actually implemented with the training data, and the necessary parameters are estimated from the training data to apply this process to another dataset later.
rec_spec <- recipe(Esc_CAT ~ ., data_tr) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric())

rec_spec %>% prep() %>% juice() %>% glimpse()
```

## Model specifications.
```{r}
# XGBoosting
xgb_spec <- boost_tree(
           trees = 500,
           learn_rate = tune(),
           min_n = tune()
       ) %>%
           set_engine("xgboost", scale_pos_weight = tune()) %>%
           set_mode("classification")

# Random forest
rf_spec <- rand_forest(
                  trees = 500,
                  min_n = tune()
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

# Neural networks
nn_spec <- mlp(
  hidden_units = 2,
  epochs = 20,
  activation = "softmax") %>%
  set_engine("nnet") %>%
  set_mode("classification")

# Support vector machine
svm_spec <- svm_rbf(
  cost = 0.5
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

# Logistic regression
lr_spec <- logistic_reg(
              ) %>%
                  set_engine("glm") %>%
                  set_mode("classification")

# Naive Bayes
library(discrim)
nb_spec <- naive_Bayes(
  smoothness = tune(),
  Laplace = tune() %>%
    set_engine("klaR") %>%
    set_mode("classification")
)

```
## Model metrics and control.
```{r}
# Model metrics
model_metrics <- metric_set(roc_auc, bal_accuracy)

# Model control
model_control <- control_stack_grid()
```

## Tuning hyperparameters.
```{r}
# Creating a set of workflows.
model_set <- workflow_set(
  preproc = list(regular = rec_spec),
  models = list(xgboost = xgb_spec, randomForest = rf_spec, nnet = nn_spec, svm = svm_spec, logreg = lr_spec, nbayes = nb_spec),
  cross = TRUE
)

# Tuning models using ANOVA-based racing optimization.
set.seed(12)
model_set <- model_set %>% 
  workflow_map("tune_race_anova", resamples = data_folds, metrics = model_metrics)

# Tuned model set results
model_set

```
## Evaluating model performance.
```{r}
autoplot(model_set)
autoplot(model_set, select_best = TRUE) + theme_minimal()

```
## Evaluating model performance.
```{r}
autoplot(model_set)

autoplot(model_set, metric = "roc_auc", select_best = TRUE) + 
  ylab("AUC ROC") +
  theme_pubr() +
  theme(
    legend.position = "right"
  )

```
## Ranking model performance on resamples based on AUC.
```{r}
rank_results(model_set, rank_metric = "roc_auc") %>% 
  filter(.metric == "roc_auc")

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
  extract_workflow_set_result("regular_randomForest") %>% 
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_randomForest") %>% 
  show_best()

model_set %>% 
  extract_workflow_set_result("regular_randomForest") %>% 
  plot_race() + ylab("AUC") + theme_minimal()

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
Random forest and naive Bayes classifiers showed similar results.

## Training the best models with the best hyperparameter values.
```{r}

wf_rf <- workflow() %>% 
  add_model(rf_spec) %>% 
  add_recipe(rec_spec)

# Expanding the range of hyperparameter values for the best model.
rf_grid <- grid_latin_hypercube(
  min_n(range = c(20,40)),
  size = 20
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
#### Random forest
```{r}
# Training the best model with the best hyperparameter values.

rf_model <- rand_forest(
                  trees = 500,
                  min_n = 2
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res_2 <- fit_resamples(rf_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

#### Naive Bayes
```{r}
# Training the best model with the best hyperparameter values.

nb_model <- naive_Bayes(
  smoothness = 1.4458423,
  Laplace = 0.8565532 %>%
    set_engine("klaR") %>%
    set_mode("classification")
)

set.seed(12)
nb_res <- fit_resamples(nb_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```
## Random forest ROC_AUC
```{r}
# ROC_AUC 

rf_res_2 %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

rf_res_2 %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

rf_res_2 %>%
  collect_metrics()

```
## Finding the cut-point that gives the highest Youden’s index for the predictions of the training set - Random forest.
```{r}
library(probably)

thresholds <- seq(0.1, 0.9, by = 0.1)

threshold_data_tr <- rf_res_2 %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_tr

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

```

## Naive Bayes ROC_AUC
```{r}
nb_res %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

nb_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

nb_res %>%
  collect_metrics()

```
## Finding the cut-point that gives the highest Youden’s index for the predictions of the training set - Naive Bayes.
```{r}
threshold_data_tr_nb <- nb_res %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_tr_nb

threshold_data_tr_nb <- threshold_data_tr_nb %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_tr_nb <- threshold_data_tr_nb %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_tr_nb

threshold_data_tr_nb %>%
  filter(.threshold == max_j_index_threshold_tr_nb)

```

## Model comparison
Bayesian analysis was used here to answer the question: 'When looking at resampling results, are the differences between models real?' To address this, a model was created with the outcome being the resampling statistics (ROC_AUC).
```{r}
# rf x nb
model_comp <- inner_join(
  rf_res_2 %>% 
    collect_metrics(summarize = FALSE) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(id, rf = .estimate),
  nb_res %>% 
    collect_metrics(summarize = FALSE) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(id, nb = .estimate)
)

class_post <- perf_mod(model_comp)

# Comparison visualization
autoplot(class_post) + theme_minimal()
autoplot(contrast_models(class_post)) + theme_minimal()

```

## RF variable importance.
```{r}
rf_model %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Esc_CAT ~ .,
    data = juice(prep(rec_spec))
  ) %>%
  vip(geom = "point") + theme_minimal()

```

## The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.
```{r}
#rf.
final_wf_rf <- workflow() %>% 
  add_model(rf_model) %>% 
  add_recipe(rec_spec)

final_res_rf <- last_fit(final_wf_rf, data_split) 

final_res_rf %>%
  collect_metrics

final_res_rf %>% 
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() 

# Bootstrapping for 95% confidence interval
set.seed(12)
final_res_rf_boot <- int_pctl(final_res_rf, alpha = 0.05)
final_res_rf_boot

```

## Finding the cut-point that gives the highest Youden’s index for the predictions of the testing set.
```{r}
threshold_data_te <- final_res_rf %>% 
  collect_predictions() %>%
  threshold_perf(Esc_CAT, .pred_1, event_level = "second", thresholds)

threshold_data_te 

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

threshold_data_te %>%
  filter(.threshold == max_j_index_threshold_te)

```

## Building an API. 
```{r}
library(vetiver)

model_api <- extract_workflow(final_res_rf) %>%
  vetiver_model("ABCL_response_rf")

library(plumber)

# class probabilities (remove the "#" in the code to run the API locally).
pr() %>%
  vetiver_api(model_api, type = "prob") #  %>% pr_run() 

```

## Plotting the training and testing ROC curves
```{r}
# ROC - Training set
roc_rf_res_2 <- rf_res_2 %>%
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "Training")  

# ROC - Testing set
roc_final_res_rf <- final_res_rf %>%
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "Testing")  

# Combining ROC curves into one df
roc_data <- bind_rows(roc_rf_res_2, roc_final_res_rf)

# Plotting both ROCs in the same plot
ggplot(roc_data, aes(x = 1 - specificity,
                     y = sensitivity,
                     color = Model)) +
  geom_abline(lty = 1,
              color = "black",
              size = 0.5) +  
  geom_path(alpha = 0.8,
            size = 1.2) +  
  labs(x = "1 - Specificity",
    y = "Sensitivity") +
  scale_color_manual(values = c("Training" = "skyblue4",
                                "Testing" = "firebrick4"),
                     name = NULL) +  
  coord_equal() +
  theme_pubr()
```

