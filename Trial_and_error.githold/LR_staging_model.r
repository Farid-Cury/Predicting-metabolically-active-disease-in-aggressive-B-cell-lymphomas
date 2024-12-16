---
title: "Staging model"
author: "Farid"
date: "`r Sys.Date()`"
output: html_document
---

```{r}
# Trying a model with only Ann-Arbor staging as predictor variable.

library(tidymodels)
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

# Changing the variables of type "character" to "factor".
data_nonnormal <- mutate_if(data_nonnormal, is.character, factor)


# Selecting Ann-Arbor staging as predictor variable.
data_stag <- select(data_nonnormal, Esc_CAT, estadiamento_cat)


# Training and Testing split.
set.seed(12)
data_split <- initial_split(data_stag, prop = 0.7, strata = Esc_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds <- vfold_cv(data_tr, 5, strata = Esc_CAT)

```

```{r}
# Recipes:
# First, the recipe() needs to be informed about the model being used and the training data.
# After that, the factor columns are converted into one or more numeric binary variables (0 and 1) for the levels present in the training data (using step_dummy()).
# Lastly, the recipe() is prepped, meaning the steps are actually implemented with the training data, and the necessary parameters are estimated from the training data to apply this process to another dataset later.
rec_spec <- recipe(Esc_CAT ~ ., data_tr) %>%
  step_dummy(all_nominal_predictors()) 


rec_spec %>% prep() %>% juice() %>% glimpse()
```

```{r}
# Model specifications.

xgb_spec <- boost_tree(
           trees = 500,
           learn_rate = tune()
       ) %>%
           set_engine("xgboost", scale_pos_weight = tune()) %>%
           set_mode("classification")

rf_spec <- rand_forest(
                  trees = 500,
                  min_n = tune()
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")
                  
nn_spec <- mlp(
  hidden_units = 2,
  epochs = 20,
  activation = "softmax") %>%
  set_engine("nnet") %>%
  set_mode("classification")

lr_spec <- logistic_reg(
              ) %>%
                  set_engine("glm") %>%
                  set_mode("classification")

svm_spec <- svm_rbf(
  cost = 0.5
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

library(discrim)
nb_spec <- naive_Bayes(
  smoothness = tune(),
  Laplace = tune() %>%
    set_engine("klaR") %>%
    set_mode("classification")
)

```

```{r}
# Model metrics and control.
model_metrics <- metric_set(roc_auc, accuracy)

model_control <- control_stack_grid()
```

```{r}
# Tuning hyperparameters.

model_set <- workflow_set(
  preproc = list(regular = rec_spec),
  models = list(xgboost = xgb_spec, randomForest = rf_spec, nnet = nn_spec, svm = svm_spec, logreg = lr_spec, nbayes = nb_spec),
  cross = TRUE
)

set.seed(12)
model_set <- model_set %>% 
  workflow_map("tune_race_anova", resamples = data_folds, metrics = model_metrics)

model_set

```
```{r}
# Evaluating model performance.
autoplot(model_set)
autoplot(model_set, select_best = TRUE)

```
```{r}
# Ranking models.
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

#nnet.
model_set %>% 
  extract_workflow_set_result("regular_nnet") %>% 
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_nnet") %>% 
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
```{r}
# All the models had the same result.
# Choosing the simplest model.

lr_model <- lr_spec


#lr
set.seed(12)
lr_res <- fit_resamples(lr_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

```{r}
# ROC_AUC

#lr
lr_res %>% 
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot()

lr_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot()

lr_res %>% 
  collect_predictions() %>%
  conf_mat(Esc_CAT, .pred_class) %>% 
  autoplot(type = "heatmap")


```




```{r}
# The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.


#lr.
final_wf_lr <- workflow() %>% 
  add_model(lr_model) %>% 
  add_recipe(rec_spec)

final_res_lr <- last_fit(final_wf_lr, data_split) 

final_res_lr %>%
  collect_metrics

final_res_lr %>% 
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot()

final_res_lr %>% 
  collect_predictions() %>%
  conf_mat(Esc_CAT, .pred_class) %>% 
  autoplot(type = "heatmap")

set.seed(12)
final_res_lr_boot <- int_pctl(final_res_lr, alpha = 0.05)
final_res_lr_boot

```




