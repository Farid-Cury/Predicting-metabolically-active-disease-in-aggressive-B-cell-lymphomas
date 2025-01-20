---
title: "Hong 2024 - feature selection - clin-lab and rad"
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

## Selecting only the training set for exploratory data analysis
```{r}
# Training split.
set.seed(12)
data_split <- initial_split(data_nonnormal, prop = 0.7, strata = Outcome_5PS_CAT)
data_tr_all <- training(data_split)
```


## Selecting clinical-laboratory and radiomic features.
```{r}
data_tr_fs <- data_tr_all %>%
  select(CL_Sex:CL_Bulky_disease, Semantic_SUVmax:Outcome_5PS_CAT)
```


## Recipes - preprocessing steps
```{r}
# Recipes 1 - dummy and normalization preprocessing steps
rec_spec <- recipe(Outcome_5PS_CAT ~ ., data = data_tr_fs) %>%
  step_dummy(all_nominal_predictors()) %>% # encoding categorical features
  step_normalize(Semantic_SUVmax:Semantic_MTV) %>% # normalizing semantic features
  step_normalize(Texture_Entropy:Texture_SumVar) 

# Prepping the recipe.
rec_prep <- prep(rec_spec)

# Applying to dataset.
data_tr_fs_baked <- bake(rec_prep, new_data = NULL)

```

```{r}
# Transforming the dataset into long format
data_tr_fs_long <- data_tr_fs_baked %>%
  pivot_longer(
    cols = -Outcome_5PS_CAT,            
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
        Outcome_5PS_CAT ~ value,  
        data = .x, 
        family = "binomial"
      )
  
      tidy(model, conf.int = TRUE) %>%
        mutate(variable = .y)  
    }
  )

# biding results into one df.
step_1_df <- bind_rows(resultList)

# Removing rows with "(Intercept)"
step_1_df_cleaned <- step_1_df %>%
  filter(!grepl("(Intercept)", term))%>%
  mutate(OR = exp(estimate),
        conf.low.or = exp(conf.low),
        conf.high.or = exp(conf.high)
        ) %>%
  select(variable, OR, conf.low.or, conf.high.or, p.value)

View(step_1_df_cleaned)

```

## Step 2 - Screening of features with p.value <= 0.5
```{r}

step_2_df <- step_1_df_cleaned %>%
  dplyr::mutate(relevance_screening = case_when(
    p.value > 0.5 ~ 0,
    p.value <= 0.5 ~ 1
  ))

step_2_df_filt <- step_2_df %>%
  filter(relevance_screening == 1)

View(step_2_df_filt)
```


## Step 3 - Correlation or VIF. In this case, Spearman correlation analysis was performed.
```{r}
# Filtering only numeric features (radiomic features) for correlation analysis.
step_2_df_filt_rad <- step_2_df_filt %>%
  filter(!grepl("X1", variable$variable))
print(step_2_df_filt_rad$variable, n = 16)
```

### Correlation analysis & plots - semantic and texture-based features
```{r}
# Correlation analysis - Spearman.
data_tr_corr_semantic <- data_tr_fs %>%
  dplyr::select(Semantic_MTV,				
                Semantic_TLG)

data_tr_corr_texture <- data_tr_fs %>%
  dplyr::select(Texture_ASM,				
                Texture_ClusP,				
                Texture_Clussh,				
                Texture_Contrast,
                Texture_Correlation,
                Texture_DifEnt,				
                Texture_DifVar,
                Texture_Entropy,				
                Texture_Homogeneity,				
                Texture_IDM,				
                Texture_MC1,				
                Texture_MC2,				
                Texture_MaxP,				
                Texture_SumEnt
                )

# Correlation matrix and network plot
cor_mat_semantic <- data_tr_corr_semantic %>%
  correlate(method = "spearman")  
fashion(cor_mat_semantic)
network_plot(cor_mat_semantic, min_cor = 0.2)


any_over_90 <- function(x) any(x >= .8, na.rm = TRUE)
cor_mat_semantic %>% 
  focus_if(any_over_90, mirror = TRUE) %>%
  shave() %>% 
  rplot(shape = 20, print_cor = T) +
  theme(legend.position = "none")


# Correlation matrix and network plot
cor_mat_texture <- data_tr_corr_texture %>%
  correlate(method = "spearman") %>%    # Create correlation data frame (cor_df)
  rearrange()  # rearrange by correlations
fashion(cor_mat_texture) 
network_plot(cor_mat_texture, min_cor = 0.2)


# Correlation matrix for highly correlated features
cor_mat_texture %>% 
  focus_if(any_over_90, mirror = TRUE) %>%
  shave() %>% 
  rplot(shape = 20, print_cor = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

```

### Semantic -  Creating columns with p.values.
```{r}
data_tr_corr_semantic_long <- data_tr_corr_semantic %>%
  correlate(method = "spearman") %>%
  stretch() %>%
  mutate(Correlation = abs(r)) %>% # obtaining the absolute value of correlation.
  select(-r)

data_tr_corr_semantic_long_p <- data_tr_corr_semantic_long %>% # **pesquisar uma melhor forma de fazer isso**
  dplyr::mutate(p.value_x = case_when(
    x == "Semantic_MTV" ~ 0.01926608,
    x == "Semantic_TLG" ~ 0.07673238,
  ),
  p.value_y = case_when(
    y == "Semantic_MTV" ~ 0.01926608,
    y == "Semantic_TLG" ~ 0.07673238,
  )) %>% 
  dplyr::mutate(min_p.value = pmin(p.value_x, p.value_y), # Smallest p.value between correlated features.
    variable_to_keep = ifelse(p.value_x <= p.value_y, x, y) # keeping the variable with the smallest p.value.
  )
```
### Semantic -  Comparing p-values of highly correlated features (|r| ≥ 0.8) and retaining the one with the smallest p-value.
```{r}
filter_high_corr <- function(data_tr_corr_semantic_long_p) {
  # Lista inicial de variáveis a manter
  variables_to_keep <- unique(c(data_tr_corr_semantic_long_p$x, data_tr_corr_semantic_long_p$y))
  
  repeat {
    # Filtrar pares altamente correlacionados
    high_corr_pairs <- data_tr_corr_semantic_long_p %>%
      filter(Correlation >= 0.8) %>%
      arrange(desc(Correlation)) # Ordenar por |r|
    
    # Se não houver mais pares altamente correlacionados, parar
    if (nrow(high_corr_pairs) == 0) break

    # Selecionar o par com maior correlação
    pair <- high_corr_pairs[1, ]

    # Determinar variável a remover
    var_to_remove <- ifelse(pair$p.value_x > pair$p.value_y, pair$x, pair$y)

    # Atualizar a lista de variáveis a manter
    variables_to_keep <- setdiff(variables_to_keep, var_to_remove)

    # Remover todos os pares que incluem a variável eliminada
    data_tr_corr_semantic_long_p <- data_tr_corr_semantic_long_p %>%
      filter(x != var_to_remove & y != var_to_remove)
  }
  
  return(variables_to_keep)
}

# Executar o filtro
final_variables_semantic <- filter_high_corr(data_tr_corr_semantic_long_p) 
print(final_variables_semantic)

```

### Texture -  Creating columns with p.values.
```{r}
# Transforming the correlation df into a long format df.
data_tr_corr_texture_long <- data_tr_corr_texture %>%
  correlate(method = "spearman") %>%
  stretch() %>%
  mutate(Correlation = abs(r)) %>% # obtaining the absolute value of correlation.
  select(-r)

# Creating the columns p.value_x and p.value_y to compare p-values of highly correlated features (|r| ≥ 0.8) and retain the one with the smallest p-value.
data_tr_corr_texture_long_p <- data_tr_corr_texture_long %>% # **pesquisar uma melhor forma de fazer isso**
  dplyr::mutate(p.value_x = case_when(
    x == "Texture_ASM" ~ 0.31918769,
    x == "Texture_ClusP" ~ 0.04876238,
    x == "Texture_Clussh" ~ 0.47347491,
    x == "Texture_Contrast" ~ 0.36196240,
    x == "Texture_Correlation" ~ 0.29025037,
    x == "Texture_DifEnt" ~ 0.26134897,
    x == "Texture_DifVar" ~ 0.33231555,
    x == "Texture_Entropy" ~ 0.22299324,
    x == "Texture_Homogeneity" ~ 0.35326135,
    x == "Texture_IDM" ~ 0.06926849,
    x == "Texture_MC1" ~ 0.33017587,
    x == "Texture_MC2" ~ 0.45219128,
    x == "Texture_MaxP" ~ 0.17738158,
    x == "Texture_SumEnt" ~ 0.04702447
  ),
  p.value_y = case_when(
    y == "Texture_ASM" ~ 0.31918769,
    y == "Texture_ClusP" ~ 0.04876238,
    y == "Texture_Clussh" ~ 0.47347491,
    y == "Texture_Contrast" ~ 0.36196240,
    y == "Texture_Correlation" ~ 0.29025037,
    y == "Texture_DifEnt" ~ 0.26134897,
    y == "Texture_DifVar" ~ 0.33231555,
    y == "Texture_Entropy" ~ 0.22299324,
    y == "Texture_Homogeneity" ~ 0.35326135,
    y == "Texture_IDM" ~ 0.06926849,
    y == "Texture_MC1" ~ 0.33017587,
    y == "Texture_MC2" ~ 0.45219128,
    y == "Texture_MaxP" ~ 0.17738158,
    y == "Texture_SumEnt" ~ 0.04702447
  )) %>% 
  dplyr::mutate(min_p.value = pmin(p.value_x, p.value_y), # Smallest p.value between correlated features.
    variable_to_keep = ifelse(p.value_x <= p.value_y, x, y) # keeping the variable with the smallest p.value.
  )

```

```{r}
filter_high_corr <- function(data_tr_corr_texture_long_p) {
  # Lista inicial de variáveis a manter
  variables_to_keep <- unique(c(data_tr_corr_texture_long_p$x, data_tr_corr_texture_long_p$y))
  
  repeat {
    # Filtrar pares altamente correlacionados
    high_corr_pairs <- data_tr_corr_texture_long_p %>%
      filter(Correlation >= 0.8) %>%
      arrange(desc(Correlation)) # Ordenar por |r|
    
    # Se não houver mais pares altamente correlacionados, parar
    if (nrow(high_corr_pairs) == 0) break

    # Selecionar o par com maior correlação
    pair <- high_corr_pairs[1, ]

    # Determinar variável a remover
    var_to_remove <- ifelse(pair$p.value_x > pair$p.value_y, pair$x, pair$y)

    # Atualizar a lista de variáveis a manter
    variables_to_keep <- setdiff(variables_to_keep, var_to_remove)

    # Remover todos os pares que incluem a variável eliminada
    data_tr_corr_texture_long_p <- data_tr_corr_texture_long_p %>%
      filter(x != var_to_remove & y != var_to_remove)
  }
  
  return(variables_to_keep)
}

# Executar o filtro
final_variables <- filter_high_corr(data_tr_corr_texture_long_p) # final_variables contem tanto as variáveis que devem ser mantidas por não serem altamente correlacionadas (Texture_IDM, Texture_ClusP e Texture_Clussh), quanto as variáveis que estavam correlacionadas (abs(r) >= 0.8) com outra variável, porém, apresentaram um menor valor de p (Texture_DifVar e Texture_SumEnt --> apesar dessas duas variáveis estarem correlacionadas com outras variáveis que foram removidas, elas não estão correlacionadas entre si).
print(final_variables)
```

```{r}
# selecting features from step 2 (p<0.5) and removing highly correlated features identified at step 3.
data_tr_p <- data_tr_fs %>% 
  select(
  # texture-based features
  Texture_ClusP, Texture_Clussh, Texture_Correlation, 
  Texture_DifEnt, Texture_IDM, Texture_MC2, Texture_SumEnt,
  # semantic features
  Semantic_MTV,
  # clinical-laboratory features
  CL_Staging,
  CL_High_LDH,
  CL_Extranodal_site_involvement,
  CL_Bulky_disease,
  CL_Sex,
  Outcome_5PS_CAT)

# Recipes 2: preprocessing steps for features selected at the end of step 3.
rec_corr <- recipe(Outcome_5PS_CAT ~ ., data = data_tr_p) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(Semantic_MTV) %>% # normalizing semantic features
  step_normalize(Texture_ClusP:Texture_SumEnt) %>% # normalizing texture features
  step_corr(Texture_ClusP:Texture_SumEnt, threshold = 0.8, method = "spearman")
```

```{r}
# Prepping the recipe.
rec_prep_corr <- prep(rec_corr)

# Applying to the dataset.
rec_corr_baked <- bake(rec_prep_corr, new_data = NULL)

glimpse(rec_corr_baked)
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
data_tr_final <- data_tr_p %>%
  select(CL_Staging, 
         Semantic_MTV,
         Texture_ClusP,
         Texture_IDM,
         Texture_SumEnt,
         Outcome_5PS_CAT) # selected features at the end of step 4 (not highly correlated & p.value <= 0.2)

glimpse(data_tr_final)

# Cross-validation folds.
set.seed(12)
data_folds_vi <- vfold_cv(data_tr_final, 10, repeats = 5, strata = Outcome_5PS_CAT)
```

```{r}
# Recipes 3: Preprocessing steps for features selected at the end of step 4
rec_final <- recipe(Outcome_5PS_CAT ~ ., data = data_tr_final) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(Semantic_MTV) %>% # normalizing semantic features
  step_normalize(Texture_ClusP:Texture_SumEnt) %>% # normalizing texture features
  step_corr(Texture_ClusP:Texture_SumEnt, threshold = 0.8, method = "spearman")

rec_final %>% prep() %>% juice() %>% glimpse()
```

## Model specifications for variable importance analysis.
```{r}
# Random Forest
rf_spec <- rand_forest(
                  trees = 1000,
                  mtry = tune(),
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
model_control <- control_race(
      save_pred = TRUE,
      parallel_over = "everything",
      save_workflow = TRUE
   )
```

## Tuning hyperparameters.
```{r}
# Creating a set of workflows.
model_set <- workflow_set(
  preproc = list(regular = rec_final),
  models = list(rf = rf_spec, lasso = lasso_spec, en = en_spec),
  cross = TRUE
)

# Tuning models using ANOVA-based racing optimization.
model_set <- model_set %>% 
  workflow_map("tune_race_anova",
               seed = 12,
               resamples = data_folds_vi,
               grid = 25,
               metrics = model_metrics,
               control = model_control)

# Tuned model set results.
model_set

```

## Evaluating model performance.
```{r}
autoplot(model_set)
autoplot(model_set,
         select_best = TRUE) +
  geom_text_repel(aes(label = wflow_id), nudge_x = 0.4, nudge_y = 0) +
  theme_pubr() +
  border() +
  theme(legend.position = "none")

```

## Ranking model performance on resamples based on AUC.
```{r}
rank_results(model_set, select_best = TRUE, rank_metric = "roc_auc") %>%
  filter(.metric == "roc_auc") %>%
  mutate(conf.low = mean - (std_err*1.96),
         conf.high = (std_err*1.96) + mean) %>%
  dplyr::select(wflow_id, mean, conf.low, conf.high)

# Ranking the best models for each classifier. It is important to check not only the model with the highest AUC but also the model with the highest AUC and the lowest standard error.

#rf.
model_set %>% 
  extract_workflow_set_result("regular_rf") %>% 
  collect_metrics()

model_set %>% 
  extract_workflow_set_result("regular_rf") %>% 
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
                  trees = 1000,
                  mtry = 3,
                  min_n = 8
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
  penalty = 0.004997903,
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
  penalty = 0.02940719,
  mixture = 0.5495942
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
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

rf_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

rf_res %>%
  collect_metrics()

```

## ROC_AUC Lasso
```{r}
lasso_res %>%
  collect_predictions() %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

lasso_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

lasso_res %>%
  collect_metrics()
```

## ROC_AUC en
```{r}
en_res %>%
  collect_predictions() %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

en_res %>% 
  collect_predictions() %>% 
  group_by(id) %>% 
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  autoplot() + theme_minimal()

en_res %>%
  collect_metrics()
```

## Finding the cut-point that gives the highest Youden’s index for the predictions - Random forest.
```{r}
library(probably)

thresholds <- seq(0.1, 0.9, by = 0.001)

threshold_data_tr <- rf_res %>% 
  collect_predictions() %>%
  threshold_perf(Outcome_5PS_CAT, .pred_1, event_level = "second", thresholds)

threshold_save <- threshold_data_tr %>%
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

```

## Confusion matrix
```{r}

pred_rf_res <- rf_res %>% 
  collect_predictions() %>%
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.257, 1, 0))
pred_rf_res$pred_class_youden <- as.factor(pred_rf_res$pred_class_youden)

pred_rf_res$Outcome_5PS_CAT <- fct_recode(pred_rf_res$Outcome_5PS_CAT, "1-3" = "0", "4-5" = "1")
pred_rf_res$pred_class_youden <- fct_recode(pred_rf_res$pred_class_youden, "1-3" = "0", "4-5" = "1")

pred_rf_res %>%
  conf_mat(Outcome_5PS_CAT, pred_class_youden) %>%
  autoplot("heatmap") +
  theme(axis.ticks = element_blank(),
    axis.text = element_text(size = 12, face = "bold"))
```


## Random Forest variable importance.
```{r}
set.seed(12)
rf_model_imp <- rf_model %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Outcome_5PS_CAT ~ ., data = juice(prep(rec_final)))

importance_data <- vi(rf_model_imp) %>%
  arrange(desc(Importance)) %>%
  mutate(highlight = if_else(row_number() <= 3, "Top 3", "Others"))  

ggplot(importance_data, aes(x = reorder(Variable, Importance), y = Importance, fill = highlight)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Top 3" = "#FF8282", "Others" = "#797A7B")) +
  theme_pubr() +
  labs(x = NULL,
       y = "Importance",
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

set.seed(12)
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

set.seed(12)
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
       title = "Elastic Net Variable Importance") +
  theme_pubr() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")
```
















