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
data_nonnormal <- read_excel("data_nonnormal.xlsx", 
     col_types = c("text", "text", "text", 
         "text", "text", "text", "text", "text", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "text", "text", "text", "numeric", 
         "text", "text", "text", "text"))

data_nonnormal <- mutate_if(data_nonnormal, is.character, factor)

# selecting only radiomic features as predictors
data_nonnormal_rad <- data_nonnormal %>%
  dplyr::select(SUVmax:Esc_CAT)

```

```{r}

# Training and testing split.
set.seed(12)
data_split <- initial_split(data_nonnormal_rad, prop = 0.7, strata = Esc_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds <- vfold_cv(data_tr, 5, strata = Esc_CAT)

```

```{r}
data_tr_rad <- data_tr %>%
  dplyr::select(SUVmax:SumVar)

library(corrplot)
corrplot(cor(data_tr_rad),
   method = "color", 
   addCoef.col="black", 
   order = "AOE", 
   number.cex=0.3,
   tl.cex=0.7,
   tl.col = "black")

library(tidyverse)
library(corrr)
correlate(data_tr_rad, method = "pearson")

cor_mat <- data_tr_rad %>%
  correlate() %>%    # Create correlation data frame (cor_df)
  rearrange()  # rearrange by correlations
fashion(cor_mat) 
network_plot(cor_mat, min_cor = 0.2)

# rplot(cor_mat, colors = c("#526c9d", "#556898", "#586592", "#5c618d", "#5f5d87", "#625982", "#65567d", "#695277", "#6c4e72", "#6f4a6c", "#724767", "#764361", "#793f5c", "#7c3c57", "#7f3851", "#83344c", "#863046", "#892d41", "#8c293c", "#902536", "#932231", "#961e2b", "#991a26", "#9d1620", "#a0131b", "#a30f16", "#a60b10", "#aa070b", "#ad0405", "#b00000")) # <- correlation matrix using corrr package.

# t-test betweent radiomic features using corrr package.
calc_ttest_p_value <- function(vec_a, vec_b){
  t.test(vec_a, vec_b)$p.value
}
colpair_map(data_tr_rad, calc_ttest_p_value)

```

```{r}
#Recipes:
# Quantitative Variables: The quantitative variables have already been standardized using the Z-score method (therefore, the "step_normalize" step is not necessary).
# Downsample: Class imbalance correction through downsampling of cases from the majority class.
# Step_dummy: One-hot encoding of categorical variables.

rec_spec <- recipe(Esc_CAT ~ ., data_tr) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
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
                  
lasso_spec <- logistic_reg(
  penalty = tune(),
  mixture = 1
  ) %>%
  set_engine("glmnet")

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
  models = list(xgboost = xgb_spec, randomForest = rf_spec, lasso = lasso_spec),
  cross = TRUE
)

set.seed(28)
model_set <- model_set %>% 
  workflow_map("tune_race_anova", resamples = data_folds, metrics = model_metrics)

model_set

```

```{r}
# Evaluating model performance.
autoplot(model_set) + theme_minimal()
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

#lasso.
model_set %>% 
  extract_workflow_set_result("regular_lasso") %>% 
  collect_metrics()

model_set %>%
  extract_workflow_set_result("regular_lasso") %>%
  show_best()

```

```{r}
# Training the best model with the best hyperparameter values.

rf_model <- rand_forest(
                  trees = 500,
                  min_n = 23
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(28)
rf_res <- fit_resamples(rf_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

```{r}
# Training the best model with the best hyperparameter values.

lasso_model <- logistic_reg(
  penalty = 0.034629959,
  mixture = 1
  ) %>%
  set_engine("glmnet")

set.seed(28)
lasso_res <- fit_resamples(lasso_model,
                                   rec_spec,
                                   resamples = data_folds,
                                   control = control_resamples(save_pred = T))

```

```{r}
# ROC_AUC and confusion matrix.

rf_res %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res %>%
  collect_metrics()

```

```{r}
# ROC_AUC and confusion matrix.

lasso_res %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

lasso_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

lasso_res %>%
  collect_metrics()

```

```{r}
# Model comparison

# rf x lasso
model_comp <- inner_join(
  rf_res %>% 
    collect_metrics(summarize = FALSE) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(id, rf = .estimate),
  lasso_res %>% 
    collect_metrics(summarize = FALSE) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::select(id, lasso = .estimate)
)

class_post <- perf_mod(model_comp)

autoplot(class_post) + theme_minimal()
autoplot(contrast_models(class_post)) + theme_minimal()

```

```{r}
rf_model %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Esc_CAT ~ .,
    data = juice(prep(rec_spec))
  ) %>%
  vip(geom = "point") + theme_minimal()

```

```{r}

final_wf_lasso <- workflow() %>% 
  add_model(lasso_model) %>% 
  add_recipe(rec_spec)

final_wf_lasso %>%
  fit(data_tr) %>%
  pull_workflow_fit() %>%
  tidy()

final_wf_lasso %>%
  fit(data_tr) %>%
  pull_workflow_fit() %>%
  vi()

library(tidyverse)
final_wf_lasso %>%
  fit(data_tr) %>%
  pull_workflow_fit() %>%
  vi() %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL)

```

```{r}
library(forcats)
data_tr$Esc_CAT <- fct_recode(data_tr$Esc_CAT, "responder" = "0", "nonresponder" = "1")

library(ggpubr)
ggdensity(data_tr, x = "IDM",
          add = "median", rug = TRUE, color = "Esc_CAT",
          fill = "Esc_CAT", palette = "Pastel1",
          xlab = "IDM", legend.title = "5-PS") + theme_minimal()

ggdensity(data_tr, x = "MTV",
          add = "median", rug = TRUE, color = "Esc_CAT",
          fill = "Esc_CAT", palette = "Pastel1",
          xlab = "MTV", legend.title = "5-PS") + theme_minimal()

ggdensity(data_tr, x = "Correlation",
          add = "median", rug = TRUE, color = "Esc_CAT",
          fill = "Esc_CAT", palette = "Pastel1",
          xlab = "Correlation", legend.title = "5-PS") + theme_minimal()

```
# Note that IDM and MTV values are very different between responders and non-responders (good separability between the groups), unlike the correlation values.

```{r}
data_ordinal <- read_excel("data_ordinal.xlsx", 
     col_types = c("text", "numeric", "text", 
         "text", "text", "text", "text", "text", 
         "text", "text", "text", "text", "text", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "numeric", "numeric", "numeric", 
         "text", "text"))

data_ordinal$Esc_CAT <- fct_recode(data_ordinal$Esc_CAT, "responder" = "0", "nonresponder" = "1")
data_ordinal$estadiamento_cat <- fct_recode(data_ordinal$estadiamento_cat, "limited" = "0", "advanced" = "1")
View(data_ordinal)

# Training and Testing split.
set.seed(12)
data_split <- initial_split(data_ordinal, prop = 0.7, strata = Esc_CAT)
data_tr_ord <- training(data_split)
data_te_ord <- testing(data_split)

library(ggridges)
library(viridis)
ggplot(data_tr_ord) + geom_density_ridges(
  aes(x = IDM,
      y = score_deauville,
      group = interaction(score_deauville,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.5) +
  xlab("Inverse Difference Momment") + 
  ylab("Deauville 5-point score") +
  labs(fill = "Ann Arbor Staging") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  theme_minimal()
  
ggplot(data_tr_ord) + geom_density_ridges(
  aes(x = MTV,
      y = score_deauville,
      group = interaction(score_deauville,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.5) +
  xlab("MTV") + 
  ylab("Deauville 5-point score") +
  labs(fill = "Ann Arbor Staging") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  theme_minimal() 
```
# Note that non-responders tend to have lower IDM values. The opposite occurs for MTV, where non-responders tend to have higher MTV values.

```{r}
# HK2 expression and most important radiomic features.
ggplot(data_tr_ord) + geom_col(
  aes(x = IDM,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("IDM") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "IDM levels according to treatment response, stratified by HK2 staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v50_expres_hkii_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggplot(data_tr_ord) + geom_col(
  aes(x = Correlation,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("Correlation") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "Correlation levels according to treatment response, stratified by HK2 staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v50_expres_hkii_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggplot(data_tr_ord) + geom_col(
  aes(x = MTV,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("MTV") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "MTV levels according to treatment response, stratified by HK2 staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v50_expres_hkii_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

```
```{r}
# LDHA expression and most important radiomic features.
ggplot(data_tr_ord) + geom_col(
  aes(x = IDM,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("IDM") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "IDM levels according to treatment response, stratified by LDHA staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v56_expres_ldh5_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggplot(data_tr_ord) + geom_col(
  aes(x = Correlation,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("Correlation") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "Correlation levels according to treatment response, stratified by LDHA staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v56_expres_ldh5_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggplot(data_tr_ord) + geom_col(
  aes(x = MTV,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("MTV") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "MTV levels according to treatment response, stratified by LDHA staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v56_expres_ldh5_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

```
```{r}
# MCT4 expression and most important radiomic features.
ggplot(data_tr_ord) + geom_col(
  aes(x = IDM,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("IDM") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "IDM levels according to treatment response, stratified by MCT4 staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v74_expres_mct4_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggplot(data_tr_ord) + geom_col(
  aes(x = Correlation,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("Correlation") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "Correlation levels according to treatment response, stratified by MCT4 staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v74_expres_mct4_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggplot(data_tr_ord) + geom_col(
  aes(x = MTV,
      y = Esc_CAT,
      group = interaction(Esc_CAT,estadiamento_cat),
      fill = estadiamento_cat), alpha = 0.8) +
  xlab("MTV") + 
  ylab("Categorized 5-PS") +
  labs(fill = "Ann Arbor Staging") +
  labs(title = "MTV levels according to treatment response, stratified by MCT4 staining score") +
  scale_fill_brewer(palette = "Pastel1", direction = -1) +
  facet_wrap(vars(v74_expres_mct4_tumo)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

```

```{r}
data_tr$Esc_CAT <- fct_recode(data_tr$Esc_CAT, "0" = "responder", "1" = "nonresponder")


# selecting the most important radiomic feature in rf.
data_rf <- data_nonnormal %>%
  select(estadiamento_cat, MTV, IDM, Esc_CAT)

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
#Recipes:
# Quantitative Variables:  Z-score normalization, remove all zero variance variables, remove correlated variables. 
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
model_metrics <- metric_set(roc_auc, accuracy)

model_control <- control_stack_grid()
```

```{r}
# Tuning hyperparameters.
model_set <- workflow_set(
  preproc = list(regular = rec_spec),
  models = list(xgboost = xgb_spec, randomForest = rf_spec, svm = svm_spec, logreg = lr_spec, nbayes = nb_spec),
  cross = TRUE
)

set.seed(800)
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
  plot_race() + ylab("AUC") + theme_minimal()

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

set.seed(75)
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
                  min_n = 34
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(543)
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
  autoplot() + theme_minimal()

rf_res_2 %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

rf_res_2 %>%
  collect_metrics()

```

```{r}
# Finding the cut-point that gives the highest Youden’s index for the predictions of the training set.
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
  autoplot() + theme_minimal()

set.seed(543)
final_res_rf_boot <- int_pctl(final_res_rf, alpha = 0.05)
final_res_rf_boot

```
# The model maintained its performance, even with random initialization parameters.

```{r}
# Finding the cut-point that gives the highest Youden’s index for the predictions of the testing set.
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
