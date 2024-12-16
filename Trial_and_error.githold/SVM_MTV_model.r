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

# Building the model for 5-PS prediction using Ann Arbor Staging, IDM and MTV.

# selecting the key clinical predictors and the most important radiomic features in rf.
data_rf <- data_nonnormal %>%
  select(MTV, Esc_CAT)

# Changing the variables of type "character" to "factor".
data_rf <- mutate_if(data_rf, is.character, factor)

# Training and Testing split.
set.seed(12)
data_split <- initial_split(data_rf, prop = 0.7, strata = Esc_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds <- vfold_cv(data_tr, 5, strata = Esc_CAT)

```

```{r}
# Recipes:
# First, the recipe() needs to be informed about the model being used and the training data.
# Any numeric variables with zero variance are then removed with step_zv().
# In the next step, numeric variables are normalized (centered and scaled) to account for differing scales, as the model being trained is sensitive to these differences (using step_normalize()).
# Lastly, the recipe() is prepped, meaning the steps are actually implemented with the training data, and the necessary parameters are estimated from the training data to apply this process to another dataset later.
rec_spec <- recipe(Esc_CAT ~ ., data_tr) %>%
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric())
   


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
library(stacks)
model_metrics <- metric_set(roc_auc, bal_accuracy)

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
autoplot(model_set, select_best = TRUE) + theme_minimal()

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


```{r}
# Training the best model with the best hyperparameter values.

svm_model <- svm_rbf(
  cost = 0.5
  ) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

set.seed(12)
svm_res <- fit_resamples(svm_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

```{r}
# ROC_AUC and confusion matrix.

svm_res %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

svm_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

svm_res %>%
  collect_metrics()

```

```{r}
# Finding the cut-point that gives the highest Youden’s index for the predictions of the training set.
library(probably)

thresholds <- seq(0.1, 0.9, by = 0.1)

threshold_data_tr <- svm_res %>% 
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



```{r}
# The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.

#rf.
final_wf_svm <- workflow() %>% 
  add_model(svm_model) %>% 
  add_recipe(rec_spec)

final_res_svm <- last_fit(final_wf_svm, data_split) 

final_res_svm %>%
  collect_metrics

final_res_svm %>% 
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

set.seed(12)
final_res_svm_boot <- int_pctl(final_res_svm, alpha = 0.05)
final_res_svm_boot

```

```{r}
# Finding the cut-point that gives the highest Youden’s index for the predictions of the testing set.
threshold_data_te <- final_res_svm %>% 
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



