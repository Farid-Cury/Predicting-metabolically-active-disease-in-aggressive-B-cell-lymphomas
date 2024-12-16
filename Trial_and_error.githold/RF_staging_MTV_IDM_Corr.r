---
Title: "RF_staging_MTV_IDM_Corr"
Author: "Farid"
date: "10/08/2024"
output: html_document
---
```{r}
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

# selecting the most important radiomic feature in lasso.
data_lasso <- data_nonnormal %>%
  select(estadiamento_cat, IDM, MTV, Correlation, Esc_CAT)

# Changing the variables of type "character" to "factor".
data_lasso <- mutate_if(data_lasso, is.character, factor)

# Training and Testing split.
set.seed(12)
data_split <- initial_split(data_lasso, prop = 0.7, strata = Esc_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds <- vfold_cv(data_tr, 5, strata = Esc_CAT)

```

```{r}
#Recipes:
# Quantitative Variables:  Z-score normalization, remove all zero variance variables, remove correlated variables. 
# Downsample: Class imbalance correction through downsampling of cases from the majority class.
# Step_dummy: One-hot encoding of categorical variables.
rec_spec <- recipe(Esc_CAT ~ ., data_tr) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
  step_dummy(all_factor(), -all_outcomes())

rec_spec %>% prep() %>% juice() %>% glimpse()
```

```{r}
# Models specifications.
xgb_spec <- boost_tree(
           trees = 500,
           learn_rate = tune(),
           min_n = tune()
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

svm_spec <- svm_rbf(
  cost = 0.5
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

lr_spec <- logistic_reg(
              ) %>%
                  set_engine("glm") %>%
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

#lgb.
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
  plot_race()

#nnet.
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

```{r}
# Training the best model with the best hyperparameter values.

rf_model <- rand_forest(
                  trees = 500,
                  min_n = 29
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res_2 <- fit_resamples(rf_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

```{r}
# ROC_AUC and confusion matrix.

rf_res_2 %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot()

rf_res_2 %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot()

rf_res_2 %>%
  collect_metrics()


```

```{r}
# The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.

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
  autoplot()

set.seed(12)
final_res_rf_boot <- int_pctl(final_res_rf, alpha = 0.05)
final_res_rf_boot

```

