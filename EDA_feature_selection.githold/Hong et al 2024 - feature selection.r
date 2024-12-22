---
title: "Hong 2024 - Feature selection"
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
library(ggrepel)
library(corrr)
library(gghighlight)
library(ggridges)
library(readxl)
options(tidymodels.dark = TRUE)

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

# Converting "character" variables into "factor".
data_nonnormal <- mutate_if(data_nonnormal, is.character, factor)
```

## Selecting only the training set for feature selection.
```{r}
# Training split.
set.seed(12)
data_split <- initial_split(data_nonnormal, prop = 0.7, strata = Esc_CAT)
data_tr_all <- training(data_split)

```

## Dichotomising expression scores acording to optimal cutpoints
```{r}
data_tr_all$v50_expres_hkii_tumo <- as.numeric(data_tr_all$v50_expres_hkii_tumo)
data_tr_all$v56_expres_ldh5_tumo <- as.numeric(data_tr_all$v56_expres_ldh5_tumo)
data_tr_all$v74_expres_mct4_tumo <- as.numeric(data_tr_all$v74_expres_mct4_tumo)


data_tr_cp <- data_tr_all %>%
  dplyr::mutate(HK2tu_opt_cp = case_when(
    v50_expres_hkii_tumo <= 4 ~ 0,
    v50_expres_hkii_tumo >= 5 ~ 1
  ),
  LDHAtu_opt_cp = case_when(
    v56_expres_ldh5_tumo <= 2 ~ 0,
    v56_expres_ldh5_tumo >= 3 ~ 1
  ),
  MCT4tu_opt_cp = case_when(
    v74_expres_mct4_tumo <= 2 ~ 0,
    v74_expres_mct4_tumo >= 3 ~ 1
  ))

data_tr_cp$HK2tu_opt_cp <- as.factor(data_tr_cp$HK2tu_opt_cp)
data_tr_cp$LDHAtu_opt_cp <- as.factor(data_tr_cp$LDHAtu_opt_cp)
data_tr_cp$MCT4tu_opt_cp <- as.factor(data_tr_cp$MCT4tu_opt_cp)
```

## Selecting candidates for feature selection.
```{r}
data_tr_fs <- data_tr_cp %>%
  select(v3_sexo:CAT_glut1_Tu_imput, CAT_cd147_Tu_imput, CAT_CAIX_Tumor, GLUT1_stroma:Esc_CAT, HK2tu_opt_cp:MCT4tu_opt_cp)
```

## Relocating columns to facilitate feature selection.
```{r}
data_tr_fs <- data_tr_fs %>%
  dplyr::relocate(c(HK2tu_opt_cp, LDHAtu_opt_cp, MCT4tu_opt_cp),
                  .after = v15_doenca_vol_bulky)
```

## Recipes - preprocessing steps
```{r}
# Recipes 1 - dummy and normalization preprocessing steps
rec_spec <- recipe(Esc_CAT ~ ., data = data_tr_fs) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>% # encoding categorical features
  step_normalize(SUVmax:MTV) %>% # normalizing semantic features
  step_normalize(Entropy:SumVar) # normalizing texture features

# Prepping the recipe.
rec_prep <- prep(rec_spec)

# Applying to dataset.
data_tr_fs_baked <- bake(rec_prep, new_data = NULL)

```

```{r}
# Transforming the dataset into long format
data_tr_fs_long <- data_tr_fs_baked %>%
  pivot_longer(
    cols = -Esc_CAT,            
    names_to = "variable",   
    values_to = "value"      
  )

```

## Step 1 - Univariate logistic regression
```{r}
# performing univariate binary logistic regression in all variables.
resultList <- data_tr_fs_long %>%
  group_by(variable) %>%
  group_map(
    function(.x, .y) {
      model <- glm(
        Esc_CAT ~ value,  
        data = .x, 
        family = "binomial"
      )
  
      tidy(model) %>%
        mutate(variable = .y)  
    }
  )

# biding results into one df.
step_1_df <- bind_rows(resultList)

# Removing rows with "X0".
step_1_df_cleaned_1 <- step_1_df %>%
  filter(!grepl("X0", variable$variable))

# Removing rows with "(Intercept)"
step_1_df_cleaned_2 <- step_1_df %>%
  filter(!grepl("X0", variable$variable) & !grepl("(Intercept)", term))

glimpse(step_1_df_cleaned_2)
```


## Step 2 - Screening of features with p.value <= 0.5
```{r}

step_2_df <- step_1_df_cleaned_2 %>%
  dplyr::mutate(relevance_screening = case_when(
    p.value > 0.5 ~ 0,
    p.value <= 0.5 ~ 1
  ))

step_2_df_filt <- step_2_df %>%
  filter(relevance_screening == 1)
```


## Step 3 - Correlation or VIF. In this case, Spearman correlation analysis was performed.
```{r}
# Filtering only numeric features (radiomic features) for correlation analysis.
step_2_df_filt_rad <- step_2_df_filt %>%
  filter(!grepl("X1", variable$variable))
```

```{r}
# Correlation analysis - Spearman.
data_tr_corr_semantic <- data_tr_fs %>%
  dplyr::select(MTV, TLG)

data_tr_corr_texture <- data_tr_fs %>%
  dplyr::select(ASM,
                ClusP,
                Clussh,
                Contrast,
                Correlation,
                DifEnt,
                DifVar,
                Entropy,
                Homogeneity,
                IDM,
                MC1,
                MC2,
                MaxP,
                SumEnt
                )

# Correlation matrix and network plot
cor_mat_semantic <- data_tr_corr_semantic %>%
  correlate(method = "spearman") %>%    # Create correlation data frame (cor_df)
  rearrange()  # rearrange by correlations
fashion(cor_mat_semantic)
network_plot(cor_mat_semantic, min_cor = 0.2)

# Correlation matrix and network plot
cor_mat_texture <- data_tr_corr_texture %>%
  correlate(method = "spearman") %>%    # Create correlation data frame (cor_df)
  rearrange()  # rearrange by correlations
fashion(cor_mat_texture) 
network_plot(cor_mat_texture, min_cor = 0.2)

any_over_90 <- function(x) any(x >= .8, na.rm = TRUE)
cor_mat_texture %>% 
  focus_if(any_over_90, mirror = TRUE) %>%
  shave() %>% 
  rplot(shape = 20, print_cor = T)

```
```{r}
# transforming the correlation df into a long format df.
data_tr_corr_texture_long <- data_tr_corr_texture %>%
  correlate(method = "spearman") %>%
  stretch() %>%
  mutate(Correlation = abs(r)) %>% # obtaining the absolute value of correlation.
  select(-r)

# Creating the columns p.value_x and p.value_y to compare p-values of highly correlated features (|r| ≥ 0.8) and retain the one with the smallest p-value.
data_tr_corr_texture_long_p <- data_tr_corr_texture_long %>% # **pesquisar uma melhor forma de fazer isso**
  dplyr::mutate(p.value_x = case_when(
    x == "ASM" ~ 0.31918769,
    x == "ClusP" ~ 0.04876238,
    x == "Clussh" ~ 0.47347491,
    x == "Contrast" ~ 0.36196240,
    x == "Correlation" ~ 0.29025037,
    x == "DifEnt" ~ 0.26134897,
    x == "DifVar" ~ 0.33231555,
    x == "Entropy" ~ 0.22299324,
    x == "Homogeneity" ~ 0.35326135,
    x == "IDM" ~ 0.06926849,
    x == "MC1" ~ 0.33017587,
    x == "MC2" ~ 0.45219128,
    x == "MaxP" ~ 0.17738158,
    x == "SumEnt" ~ 0.04702447
  ),
  p.value_y = case_when(
    y == "ASM" ~ 0.31918769,
    y == "ClusP" ~ 0.04876238,
    y == "Clussh" ~ 0.47347491,
    y == "Contrast" ~ 0.36196240,
    y == "Correlation" ~ 0.29025037,
    y == "DifEnt" ~ 0.26134897,
    y == "DifVar" ~ 0.33231555,
    y == "Entropy" ~ 0.22299324,
    y == "Homogeneity" ~ 0.35326135,
    y == "IDM" ~ 0.06926849,
    y == "MC1" ~ 0.33017587,
    y == "MC2" ~ 0.45219128,
    y == "MaxP" ~ 0.17738158,
    y == "SumEnt" ~ 0.04702447
  ))

# creating columns min_p.value and variable_to_keep.
data_tr_corr_texture_long_p_1 <- data_tr_corr_texture_long_p %>% 
  dplyr::mutate(min_p.value = pmin(p.value_x, p.value_y), # Smallest p.value between correlated features.
    variable_to_keep = ifelse(p.value_x <= p.value_y, x, y) # keeping the variable with the smallest p.value.
  )

# Filtering highly correlated pairs (|r| ≥ 0.8).
high_corr_pairs_1 <- data_tr_corr_texture_long_p_1 %>%
  filter(Correlation >= 0.8) %>% 
  arrange(desc(Correlation), min_p.value) 

high_corr_pairs_clean_1 <- high_corr_pairs_1 %>%
  filter(row_number() %% 2 != 0)

# Visualizing results.
glimpse(high_corr_pairs_clean_1)

```

```{r}
# filtering highly correlated features with the smallest p-value --> step 3.1
data_tr_corr_texture_long_p_2 <- data_tr_corr_texture_long_p_1 %>%
  filter(x %in% c("MC1", "DifEnt", "DifVar", "Homogeneity", "MaxP", "Entropy", "SumEnt") &
         y %in% c("MC1", "DifEnt", "DifVar", "Homogeneity", "MaxP", "Entropy", "SumEnt"))

# Filtering highly correlated pairs.
high_corr_pairs_2 <- data_tr_corr_texture_long_p_2 %>%
  filter(Correlation >= 0.8) %>% 
  arrange(desc(Correlation), min_p.value) 

high_corr_pairs_clean_2 <- high_corr_pairs_2 %>%
  filter(row_number() %% 2 != 0)

# Visualizing results.
glimpse(high_corr_pairs_clean_2)

```

```{r}
# filtering highly correlated features with the smallest p-value --> step 3.2
data_tr_corr_texture_long_p_3 <- data_tr_corr_texture_long_p_2 %>%
  filter(x %in% c("MC1", "DifEnt", "DifVar", "MaxP", "SumEnt") &
         y %in% c("MC1", "DifEnt", "DifVar", "MaxP", "SumEnt"))

high_corr_pairs_3 <- data_tr_corr_texture_long_p_3 %>%
  filter(Correlation >= 0.8) %>% 
  arrange(desc(Correlation), min_p.value) 

high_corr_pairs_clean_3 <- high_corr_pairs_3 %>%
  filter(row_number() %% 2 != 0)

# Visualizar os resultados
glimpse(high_corr_pairs_clean_3)

```

```{r}
# filtering highly correlated features with the smallest p-value --> step 3.3
data_tr_corr_texture_long_p_4 <- data_tr_corr_texture_long_p_3 %>%
  filter(x %in% c("MC1", "DifEnt", "SumEnt") &
         y %in% c("MC1", "DifEnt", "SumEnt"))

high_corr_pairs_4 <- data_tr_corr_texture_long_p_4 %>%
  filter(Correlation >= 0.8) %>% 
  arrange(desc(Correlation), min_p.value) 

high_corr_pairs_clean_4 <- high_corr_pairs_4 %>%
  filter(row_number() %% 2 != 0)

glimpse(high_corr_pairs_clean_4)
```

```{r}
print(high_corr_pairs_clean_1) # all highly correlated pairs must be removed, except DifEnt (highly correlated feature with the smallest p-value (step 3.3)).
```


```{r}
# selecting features from step 2 (p<0.5)
data_tr_p <- data_tr_fs %>% 
  select(ASM,
  CAIX_stroma,
  ClusP,
  Clussh,
  Contrast,
  Correlation,
  DifEnt,
  DifVar,
  Entropy,
  GLUT1_stroma,
  HK2_stroma,
  HK2tu_opt_cp,
  Homogeneity,
  IDM,
  MC1,
  MC2,
  MCT4_stroma,
  MTV,
  MaxP,
  SumEnt,
  TLG,
  estadiamento_cat,
  v10_ldh_elevada,
  v12_envo_sitios_extra,
  v15_doenca_vol_bulky,
  v3_sexo,
  Esc_CAT)

# removing highly correlated features identified at step 3.
data_tr_corr <- data_tr_p %>%
  select(-TLG, -Homogeneity, -Contrast, -DifVar, -ASM, -Entropy, -MaxP, -SumEnt)

# relocating cols to facilitate step 4
data_tr_corr <- data_tr_corr %>%
  dplyr::relocate(c(GLUT1_stroma, HK2_stroma, HK2tu_opt_cp),
                  .after = CAIX_stroma)

# Recipes 2: preprocessing steps for features selected at the end of step 3.
rec_corr <- recipe(Esc_CAT ~ ., data = data_tr_corr) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_normalize(MTV) %>% # normalizing semantic features
  step_normalize(ClusP:MC2) %>% # normalizing texture features
  step_corr(ClusP:MC2, threshold = 0.8, method = "spearman")
```

```{r}
# Prepping the recipe.
rec_prep_corr <- prep(rec_corr)

# Applying to the dataset.
data_tr_corr_baked <- bake(rec_prep_corr, new_data = NULL)

glimpse(data_tr_corr_baked)
```

## Step 4 - Selecting features with a p.value <= 0.2
```{r}
# Creating the column relevance_cor.
step_4_df_cor <- step_2_df_filt %>%
  dplyr::mutate(relevance_cor = case_when(
    p.value > 0.2 ~ 0,
    p.value <= 0.2 ~ 1
  ))

# Filtering rows with a p.value <= 0.2 (relevance_cor = 1).
step_4_df_cor_filt <- step_4_df_cor %>%
  filter(relevance_cor == 1)

# Visualizing variables with a p.value <= 0.2
print(step_4_df_cor_filt$variable)

```

## Step 5 - Final feature selection using 3 embedded feature selection algorithms (RF, LASSO, EN).
```{r}
# selecting features from the unprocessed initial df.
data_tr_final <- data_tr_fs %>%
  select(estadiamento_cat, CAIX_stroma, MCT4_stroma, HK2tu_opt_cp,
         ClusP, IDM, MTV, Esc_CAT) # selected features at the end of step 4 (not highly correlated & p.value <= 0.2)

glimpse(data_tr_final)

# Cross-validation folds.
set.seed(12)
data_folds_vi <- vfold_cv(data_tr_final, 5, strata = Esc_CAT)
```

```{r}
# Recipes 3: Preprocessing steps for features selected at the end of step 4
rec_final <- recipe(Esc_CAT ~ ., data = data_tr_final) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_normalize(MTV) %>% # normalizing semantic features
  step_normalize(ClusP:IDM) %>% # normalizing texture features
  step_corr(ClusP:IDM, threshold = 0.8, method = "spearman")

rec_final %>% prep() %>% juice() %>% glimpse()
```

## Model specifications for variable importance analysis.
```{r}
# Random Forest
rf_spec <- rand_forest(
                  trees = 500,
                  min_n = tune()
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

# LASSO              
lasso_spec <- logistic_reg(
  penalty = tune(),
  mixture = 1
  ) %>%
  set_engine("glmnet")

# Elastic net              
en_spec <- logistic_reg(
  penalty = tune(),
  mixture = tune()
  ) %>%
  set_engine("glmnet")

```

## Model metrics and control.
```{r}
# Model metrics.
model_metrics <- metric_set(roc_auc, bal_accuracy)

# Model control.
model_control <- control_stack_grid()
```

## Tuning hyperparameters.
```{r}
# Creating a set of workflows.
model_set <- workflow_set(
  preproc = list(regular = rec_final),
  models = list(randomForest = rf_spec, lasso = lasso_spec, en = en_spec),
  cross = TRUE
)

# Tuning models using ANOVA-based racing optimization.
set.seed(12)
model_set <- model_set %>% 
  workflow_map("tune_race_anova",
               resamples = data_folds_vi,
               metrics = model_metrics)

# Tuned model set results.
model_set

```

## Evaluating model performance.
```{r}
autoplot(model_set)
autoplot(model_set,
         select_best = TRUE) +
  theme_minimal()

```

## Ranking model performance on resamples based on AUC.
```{r}
rank_results(model_set, rank_metric = "roc_auc") %>% 
  filter(.metric == "roc_auc")

# Ranking the best models for each classifier. It is important to check not only the model with the highest AUC but also the model with the highest AUC and the lowest standard error.

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

#en.
model_set %>% 
  extract_workflow_set_result("regular_en") %>% 
  collect_metrics()

model_set %>%
  extract_workflow_set_result("regular_en") %>%
  show_best()

```

## Training the models with the best hyperparameter values.
#### Random Forest
```{r}
# Training the best model with the best hyperparameter values.

rf_model <- rand_forest(
                  trees = 500,
                  min_n = 39
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res <- fit_resamples(rf_model,
                                   rec_final,
                                   resamples = data_folds_vi,
                                   control = control_resamples(save_pred = T))

```

#### LASSO
```{r}
# Training the best model with the best hyperparameter values.

lasso_model <- logistic_reg(
  penalty = 7.851534e-04,
  mixture = 1
  ) %>%
  set_engine("glmnet")

set.seed(12)
lasso_res <- fit_resamples(lasso_model,
                                   rec_final,
                                   resamples = data_folds_vi,
                                   control = control_resamples(save_pred = T))

```


#### en
```{r}
# Training the best model with the best hyperparameter values.

en_model <- logistic_reg(
  penalty = 1.024186e-02,
  mixture = 0.1860880
  ) %>%
  set_engine("glmnet")

set.seed(12)
en_res <- fit_resamples(en_model,
                                   rec_final,
                                   resamples = data_folds_vi,
                                   control = control_resamples(save_pred = T))

```

## ROC_AUC Random Forest
```{r}
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

## ROC_AUC Lasso
```{r}
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

## ROC_AUC en
```{r}
en_res %>%
  collect_predictions() %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

en_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

en_res %>%
  collect_metrics()
```

## Random Forest variable importance.
```{r}
rf_model_imp <- rf_model %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Esc_CAT ~ ., data = juice(prep(rec_final)))

importance_data <- vi(rf_model_imp) %>%
  arrange(desc(Importance)) %>%
  mutate(highlight = if_else(row_number() <= 2, "Top 2", "Others"))  

ggplot(importance_data, aes(x = reorder(Variable, Importance), y = Importance, fill = highlight)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Top 2" = "#FF8282", "Others" = "#797A7B")) +
  theme_pubr() +
  labs(x = NULL,
       title = "Random Forest Variable Importance") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)
    )

```

## LASSO corfficients and variable importance.
```{r}
final_wf_lasso <- workflow() %>% 
  add_model(lasso_model) %>% 
  add_recipe(rec_final)

final_wf_lasso %>%
  fit(data_tr_final) %>%
  pull_workflow_fit() %>%
  tidy()

final_wf_lasso %>%
  fit(data_tr_final) %>%
  pull_workflow_fit() %>%
  vi()

final_wf_lasso %>%
  fit(data_tr_final) %>%
  pull_workflow_fit() %>%
  vi() %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL,
       title = "LASSO Variable Importance") +
  theme_pubr() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")
```

## En corfficients and variable importance.
```{r}
final_wf_en <- workflow() %>% 
  add_model(en_model) %>% 
  add_recipe(rec_final)

final_wf_en %>%
  fit(data_tr_final) %>%
  pull_workflow_fit() %>%
  tidy()

final_wf_en %>%
  fit(data_tr_final) %>%
  pull_workflow_fit() %>%
  vi()

final_wf_en %>%
  fit(data_tr_final) %>%
  pull_workflow_fit() %>%
  vi() %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL,
       title = "en Variable Importance") +
  theme_pubr() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")
```


