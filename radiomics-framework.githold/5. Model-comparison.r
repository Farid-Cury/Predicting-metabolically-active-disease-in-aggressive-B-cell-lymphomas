---
title: "Model comparison - Comb x rad x clin"
author: "Farid"
date: "`r Sys.Date()`"
output: html_document
---
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
## clinical-radiomic model
roc_data_rf_cl_rad_tr <- read_xlsx("roc_data_rf_cl_rad_tr.xlsx") %>%
  mutate(across(c(.pred_class, Outcome_5PS_CAT), as.factor))
roc_data_rf_cl_rad_te <- read_xlsx("roc_data_rf_cl_rad_te.xlsx") %>%
  mutate(across(c(.pred_class, Outcome_5PS_CAT), as.factor))
## Radiomic model
roc_data_rf_rad_tr <- read_xlsx("roc_data_rf_rad_tr.xlsx") %>%
  mutate(across(c(.pred_class, Outcome_5PS_CAT), as.factor))
roc_data_rf_rad_te <- read_xlsx("roc_data_rf_rad_te.xlsx") %>%
  mutate(across(c(.pred_class, Outcome_5PS_CAT), as.factor))
## Clinical model
roc_data_lr_cl_tr <- read_xlsx("roc_data_lr_cl_tr.xlsx") %>%
  mutate(across(c(.pred_class, Outcome_5PS_CAT), as.factor))
roc_data_lr_cl_te <- read_xlsx("roc_data_lr_cl_te.xlsx") %>%
  mutate(across(c(.pred_class, Outcome_5PS_CAT), as.factor))

```

## Comparing AUROCs - training set
#### Combined x radiomic
```{r}
library(pROC)

# Combined x radiomic ROC curves
roc_comb_tr <- roc(roc_data_rf_cl_rad_tr$Outcome_5PS_CAT, roc_data_rf_cl_rad_tr$.pred_1, levels = c("0", "1"))
roc_rad_tr <- roc(roc_data_rf_rad_tr$Outcome_5PS_CAT, roc_data_rf_rad_tr$.pred_1, levels = c("0", "1"))
roc_cl_tr <- roc(roc_data_lr_cl_tr$Outcome_5PS_CAT, roc_data_lr_cl_tr$.pred_1, levels = c("0", "1"))

# DeLong´s test
set.seed(12)
test_result_tr_1 <- roc.test(roc_comb_tr, roc_rad_tr, method = "delong", conf.level = 0.95)

# Test results
print(test_result_tr_1)

```

#### Combined x clinical
```{r}
# Combined x clinical ROC curves
# DeLong´s test
set.seed(12)
test_result_tr_2 <- roc.test(roc_comb_tr, roc_cl_tr, method = "delong", conf.level = 0.95)

# Test results
print(test_result_tr_2)
```

#### Radiomic x clinical
```{r}
# Radiomic x clinical ROC curves
# DeLong´s test
set.seed(12)
test_result_tr_3 <- roc.test(roc_rad_tr, roc_cl_tr, method = "delong", conf.level = 0.95)

# Test results
print(test_result_tr_3)
```

## Plotting the ROC curves - training set
```{r}
# Combined model
prob_comb_tr <- roc_data_rf_cl_rad_tr %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "combined_tr")  
  
# Radiomic model
prob_rad_tr <- roc_data_rf_rad_tr %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "radiomic_tr")  

# Clinical model
prob_cl_tr <- roc_data_lr_cl_tr %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "clinical_tr")  
  

# Combining ROC curves into one df
combined_rocs_tr <- bind_rows(prob_comb_tr, prob_rad_tr, prob_cl_tr)

# Plotting both ROCs in the same plot
ggplot(combined_rocs_tr, aes(x = 1 - specificity,
                     y = sensitivity,
                     color = Model)) +
  geom_abline(lty = 1,
              color = "black",
              size = 0.5) +  
  geom_path(alpha = 0.7,
            size = 1.2) +  
  labs(x = "1 - Specificity",
    y = "Sensitivity") +
  scale_color_manual(values = c("combined_tr" = "skyblue4",
                                "radiomic_tr" = "firebrick4",
                                "clinical_tr" = "#8a5008"),
                     name = "Model") +  
  coord_equal() +
  theme_pubr() +
  theme(legend.position = "right")

```


## Combining models predictions - training set
```{r}
# Combined model
pred_comb_tr <- roc_data_rf_cl_rad_tr %>% 
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.297, 1, 0)) %>%
  mutate(pred_class_youden = as.factor(pred_class_youden)) %>%
  mutate(id = row_number()) %>%
  dplyr::select(id, Outcome_5PS_CAT, pred_comb_tr = pred_class_youden)
  
# Radiomic model
pred_rad_tr <- roc_data_rf_rad_tr %>% 
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.395, 1, 0)) %>%
  mutate(pred_class_youden = as.factor(pred_class_youden)) %>%
  mutate(id = row_number()) %>%
  dplyr::select(id, Outcome_5PS_CAT, pred_rad_tr = pred_class_youden)

# Clinical model
pred_cl_tr <- roc_data_lr_cl_tr %>% 
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.200, 1, 0)) %>%
  mutate(pred_class_youden = as.factor(pred_class_youden)) %>%
  mutate(id = row_number()) %>%
  dplyr::select(id, Outcome_5PS_CAT, pred_cl_tr = pred_class_youden)
  
# Combining predictions
combined_preds_tr <- pred_comb_tr %>%
  full_join(pred_rad_tr, by = c("id", "Outcome_5PS_CAT")) %>%
  full_join(pred_cl_tr, by = c("id", "Outcome_5PS_CAT"))

# Visualizar o resultado
combined_preds_tr

```

```{r}
# Sensibilidade e especificidade para o modelo combinado
metrics_comb_tr <- pred_comb_tr %>%
  conf_mat(Outcome_5PS_CAT, pred_comb_tr) %>%
  summary(event_level = "second") %>%
  filter(.metric %in% c("sens", "spec")) %>%
  mutate(Model = "Combined_tr")

# Sensibilidade e especificidade para o modelo radiômico
metrics_rad_tr <- pred_rad_tr %>%
  conf_mat(Outcome_5PS_CAT, pred_rad_tr) %>%
  summary(event_level = "second") %>%
  filter(.metric %in% c("sens", "spec")) %>%
  mutate(Model = "Radiomic_tr")

# Sensibilidade e especificidade para o modelo clínico
metrics_cl_tr <- pred_cl_tr %>%
  conf_mat(Outcome_5PS_CAT, pred_cl_tr) %>%
  summary(event_level = "second") %>%
  filter(.metric %in% c("sens", "spec")) %>%
  mutate(Model = "Clinical_tr")

# Combinar as métricas em um único data frame
metrics_all_tr <- bind_rows(metrics_comb_tr, metrics_rad_tr, metrics_cl_tr)

```

```{r}
library(ggrepel)
ggplot(metrics_all_tr, aes(x = Model, y = .estimate, color = .metric, group = .metric)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = round(.estimate, 3)),  # Adiciona os valores arredondados
            vjust = -1, nudge_y = 0.02, size = 3) +         # Ajusta a posição e o tamanho do texto
  labs(
    x = "Model",
    y = "Metric",
    color = "Metric"
  ) +
  scale_color_manual(values = c("sens" = "skyblue4", "spec" = "firebrick4")) +
  theme_pubr() +
  theme(
    legend.position = "right"
  )

```

#### Combined x radiomic model
```{r}
# filter comb = FN & rad = TP
fn_comb_not_rad_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_tr == 0, pred_rad_tr == 1)
# filter rad = FN & comb = TP
fn_rad_not_comb_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_tr == 1, pred_rad_tr == 0)

fn_comb_not_rad_tr
fn_rad_not_comb_tr
```

Comparing sensitivity: combined x radiomic
```{r}
b_tr_1 <- 6 # number of positive instances in the training set misclassified as FN by the combined model but not by the radiomic model
c_tr_1 <- 5 # number of positive instances in the training set misclassified as FN by the radiomic model but not by the combined model

# Contingency matrix
contingency_matrix_1 <- matrix(c(0, b_tr_1, c_tr_1, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_1, correct = TRUE) # continuity correction

```

```{r}
# filter comb = FP and rad = TN
fp_comb_not_rad_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_tr == 1, pred_rad_tr == 0)
# filter rad = FP and comb = TN
fp_rad_not_comb_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_tr == 0, pred_rad_tr == 1)

fp_comb_not_rad_tr
fp_rad_not_comb_tr
```

Comparing specificity: combined x radiomic
```{r}
d_tr_2 <- 7 # number of negative instances in the training set misclassified as FP by the combined model but not by the radiomic model
e_tr_2 <- 7 # number of negative instances in the training set misclassified as FP by the radiomic model but not by the combined model

# Contingency matrix
contingency_matrix_2 <- matrix(c(0, d_tr_2, e_tr_2, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_2, correct = TRUE) # continuity correction

```

#### Combined x clinical model
```{r}
# filter comb = FN & rad = TP
fn_comb_not_cl_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_tr == 0, pred_cl_tr == 1)
# filter rad = FN & comb = TP
fn_rad_not_comb_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_tr == 1, pred_cl_tr == 0)

fn_comb_not_cl_tr
fn_rad_not_comb_tr
```

Comparing sensitivity: combined x clinical
```{r}
b_tr_3 <- 16 # number of positive instances in the training set misclassified as FN by the combined model but not by the clinical model
c_tr_3 <- 0 # number of positive instances in the training set misclassified as FN by the clinical model but not by the combined model

# Contingency matrix
contingency_matrix_3 <- matrix(c(0, b_tr_3, c_tr_3, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_3, correct = TRUE) # continuity correction

```

```{r}

# filter comb = FP and cl = TN
fp_comb_not_cl_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_tr == 1, pred_cl_tr == 0)
# filter cl = FP and comb = TN
fp_cl_not_comb_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_tr == 0, pred_cl_tr == 1)

fp_comb_not_cl_tr
fp_cl_not_comb_tr
```

Comparing specificity: combined x clinical
```{r}
d_tr_4 <- 5 # number of negative instances in the training set misclassified as FP by the combined model but not by the clinical model
e_tr_4 <- 73 # number of negative instances in the training set misclassified as FP by the clinical model but not by the combined model

# Contingency matrix
contingency_matrix_4 <- matrix(c(0, d_tr_4, e_tr_4, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_4, correct = TRUE) # continuity correction

```

#### Radiomic x clinical model
```{r}
# filter rad = FN & cl = TP
fn_rad_not_cl_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_rad_tr == 0, pred_cl_tr == 1)
# filter cl = FN & rad = TP
fn_cl_not_rad_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_rad_tr == 1, pred_cl_tr == 0)

fn_rad_not_cl_tr
fn_cl_not_rad_tr
```

Comparing sensitivity: radiomic x clinical
```{r}
b_tr_5 <- 20 # number of positive instances in the training set misclassified as FN by the radiomic model but not by the clinical model
c_tr_5 <- 5 # number of positive instances in the training set misclassified as FN by the clinical model but not by the radiomic model

# Contingency matrix
contingency_matrix_5 <- matrix(c(0, b_tr_5, c_tr_5, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_5, correct = TRUE) # continuity correction

```

```{r}

# filter rad = FP and cl = TN
fp_rad_not_cl_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_rad_tr == 1, pred_cl_tr == 0)
# filter cl = FP and rad = TN
fp_cl_not_rad_tr <- combined_preds_tr %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_rad_tr == 0, pred_cl_tr == 1)

fp_rad_not_cl_tr
fp_cl_not_rad_tr
```

Comparing specificity: radiomic x clinical
```{r}
d_tr_6 <- 12 # number of negative instances in the training set misclassified as FP by the radiomic model but not by the clinical model
e_tr_6 <- 80 # number of negative instances in the training set misclassified as FP by the clinical model but not by the radiomic model

# Contingency matrix
contingency_matrix_6 <- matrix(c(0, d_tr_6, e_tr_6, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_6, correct = TRUE) # continuity correction

```

## Comparing AUROCs - test set
#### Combined x radiomic
```{r}
# Combined x radiomic ROC curves
roc_comb <- roc(roc_data_rf_cl_rad_te$Outcome_5PS_CAT, roc_data_rf_cl_rad_te$.pred_1, levels = c("0", "1"))
roc_rad <- roc(roc_data_rf_rad_te$Outcome_5PS_CAT, roc_data_rf_rad_te$.pred_1, levels = c("0", "1"))
roc_cl <- roc(roc_data_lr_cl_te$Outcome_5PS_CAT, roc_data_lr_cl_te$.pred_1, levels = c("0", "1"))

# DeLong´s test
set.seed(12)
test_result_1 <- roc.test(roc_comb, roc_rad, method = "delong", conf.level = 0.95)

# Test results
print(test_result_1)


```

#### Combined x clinical
```{r}
# Combined x clinical ROC curves
# DeLong´s test
set.seed(12)
test_result_2 <- roc.test(roc_comb, roc_cl, method = "delong", conf.level = 0.95)

# Test results
print(test_result_2)
```

#### Radiomic x clinical
```{r}
# Radiomic x clinical ROC curves
# DeLong´s test
set.seed(12)
test_result_3 <- roc.test(roc_rad, roc_cl, method = "delong", conf.level = 0.95)

# Test results
print(test_result_3)
```

## Plotting the ROC curves - test set
```{r}
# Combined model
prob_comb_te <- roc_data_rf_cl_rad_te %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "combined_te")  
  
# Radiomic model
prob_rad_te <- roc_data_rf_rad_te %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "radiomic_te")  

# Clinical model
prob_cl_te <- roc_data_lr_cl_te %>%
  roc_curve(Outcome_5PS_CAT, .pred_1, event_level = "second") %>%
  mutate(Model = "clinical_te")  
  

# Combining ROC curves into one df
combined_rocs_te <- bind_rows(prob_comb_te, prob_rad_te, prob_cl_te)

# Plotting both ROCs in the same plot
ggplot(combined_rocs_te, aes(x = 1 - specificity,
                     y = sensitivity,
                     color = Model)) +
  geom_abline(lty = 1,
              color = "black",
              size = 0.5) +  
  geom_path(alpha = 0.7,
            size = 1.2) +  
  labs(x = "1 - Specificity",
    y = "Sensitivity") +
  scale_color_manual(values = c("combined_te" = "skyblue4",
                                "radiomic_te" = "firebrick4",
                                "clinical_te" = "#8a5008"),
                     name = "Model") +  
  coord_equal() +
  theme_pubr() +
  theme(legend.position = "right")

```

## Combining predictions - test set
```{r}
# Combined model
pred_comb_te <- roc_data_rf_cl_rad_te %>% 
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.297, 1, 0)) %>%
  mutate(pred_class_youden = as.factor(pred_class_youden)) %>%
  mutate(id = row_number()) %>%
  dplyr::select(id, Outcome_5PS_CAT, pred_comb_te = pred_class_youden)
  
# Radiomic model
pred_rad_te <- roc_data_rf_rad_te %>% 
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.395, 1, 0)) %>%
  mutate(pred_class_youden = as.factor(pred_class_youden)) %>%
  mutate(id = row_number()) %>%
  dplyr::select(id, Outcome_5PS_CAT, pred_rad_te = pred_class_youden)

# Clinical model
pred_cl_te <- roc_data_lr_cl_te %>% 
  mutate(pred_class_youden = ifelse(.pred_1 >= 0.200, 1, 0)) %>%
  mutate(pred_class_youden = as.factor(pred_class_youden)) %>%
  mutate(id = row_number()) %>%
  dplyr::select(id, Outcome_5PS_CAT, pred_cl_te = pred_class_youden)
  
# Combining predictions
combined_preds_te <- pred_comb_te %>%
  full_join(pred_rad_te, by = c("id", "Outcome_5PS_CAT")) %>%
  full_join(pred_cl_te, by = c("id", "Outcome_5PS_CAT"))

# Visualizar o resultado
combined_preds_te

```

```{r}
# Sensibilidade e especificidade para o modelo combinado
metrics_comb_te <- pred_comb_te %>%
  conf_mat(Outcome_5PS_CAT, pred_comb_te) %>%
  summary(event_level = "second") %>%
  filter(.metric %in% c("sens", "spec")) %>%
  mutate(Model = "Combined_te")

# Sensibilidade e especificidade para o modelo radiômico
metrics_rad_te <- pred_rad_te %>%
  conf_mat(Outcome_5PS_CAT, pred_rad_te) %>%
  summary(event_level = "second") %>%
  filter(.metric %in% c("sens", "spec")) %>%
  mutate(Model = "Radiomic_te")

# Sensibilidade e especificidade para o modelo clínico
metrics_cl_te <- pred_cl_te %>%
  conf_mat(Outcome_5PS_CAT, pred_cl_te) %>%
  summary(event_level = "second") %>%
  filter(.metric %in% c("sens", "spec")) %>%
  mutate(Model = "Clinical_te")

# Combinar as métricas em um único data frame
metrics_all_te <- bind_rows(metrics_comb_te, metrics_rad_te, metrics_cl_te)

```

```{r}

ggplot(metrics_all_te, aes(x = Model, y = .estimate, color = .metric, group = .metric)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = round(.estimate, 3)),  # Adiciona os valores arredondados
            box.padding = 2,
            vjust = -1,
            nudge_y = 0.04,
            size = 3) + 
  labs(
    x = "Model",
    y = "Metric",
    color = "Metric"
  ) +
  scale_color_manual(values = c("sens" = "skyblue4", "spec" = "firebrick4")) +
  theme_pubr() +
  theme(
    legend.position = "right",
  )

```


#### Combined x radiomic model
```{r}
# filter comb = FN & rad = TP
fn_comb_not_rad_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_te == 0, pred_rad_te == 1)
# filter rad = FN & comb = TP
fn_rad_not_comb_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_te == 1, pred_rad_te == 0)

fn_comb_not_rad_te
fn_rad_not_comb_te
```

Comparing sensitivity: combined x radiomic
```{r}

b_te_1 <- 0 # number of positive instances in the test set misclassified as FN by the combined model but not by the radiomic model 
c_te_1 <- 2 # number of positive instances in the test set misclassified as FN by the radiomic model but not by the combined model 

# Contingency matrix
contingency_matrix_7 <- matrix(c(0, b_te_1, c_te_1, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_7, correct = TRUE) # continuity correction

```

```{r}
# filter comb = FP & rad = TN
fp_comb_not_rad_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_te == 1, pred_rad_te == 0)
# filter rad = FP & comb = TN
fp_rad_not_comb_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_te == 0, pred_rad_te == 1)

fp_comb_not_rad_te
fp_rad_not_comb_te

```

Comparing specificity: combined x radiomic
```{r}
d_te_2 <- 2 # number of negative instances in the test set misclassified as FP by the combined model but not by the radiomic model
e_te_2 <- 2 # number of negative instances in the test set misclassified as FP by the radiomic model but not by the combined model

# Contingency matrix
contingency_matrix <- matrix(c(0, d_te_2, e_te_2, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix, correct = TRUE) # continuity correction

```

#### Combined x clinical
```{r}
# filter comb = FN & rad = TP
fn_comb_not_cl_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_te == 0, pred_cl_te == 1)
# filter rad = FN & comb = TP
fn_rad_not_comb_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_comb_te == 1, pred_cl_te == 0)

fn_comb_not_cl_te
fn_rad_not_comb_te
```

Comparing sensitivity: combined x clinical
```{r}
b_te_3 <- 1 # number of positive instances in the training set misclassified as FN by the combined model but not by the clinical model
c_te_3 <- 0 # number of positive instances in the training set misclassified as FN by the clinical model but not by the combined model

# Contingency mateix
contingency_matrix_3 <- matrix(c(0, b_te_3, c_te_3, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_3, correct = TRUE) # continuity correction

```

```{r}
# filter comb = FP and cl = TN
fp_comb_not_cl_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_te == 1, pred_cl_te == 0)
# filter cl = FP and comb = TN
fp_cl_not_comb_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_comb_te == 0, pred_cl_te == 1)

fp_comb_not_cl_te
fp_cl_not_comb_te
```

Comparing specificity: combined x clinical
```{r}
d_te_4 <- 0 # number of negative instances in the training set misclassified as FP by the combined model but not by the clinical model
e_te_4 <- 8 # number of negative instances in the training set misclassified as FP by the clinical model but not by the combined model

# Contingency mateix
contingency_matrix_4 <- matrix(c(0, d_te_4, e_te_4, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_4, correct = TRUE) # continuity correction

```

#### Radiomic x clinical model
```{r}
# filter rad = FN & cl = TP
fn_rad_not_cl_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_rad_te == 0, pred_cl_te == 1)
# filter cl = FN & rad = TP
fn_cl_not_rad_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 1, pred_rad_te == 1, pred_cl_te == 0)

fn_rad_not_cl_te
fn_cl_not_rad_te
```

Comparing sensitivity: radiomic x clinical
```{r}
b_te_5 <- 3 # number of positive instances in the training set misclassified as FN by the radiomic model but not by the clinical model
c_te_5 <- 0 # number of positive instances in the training set misclassified as FN by the clinical model but not by the radiomic model

# Contingency mateix
contingency_matrix_5 <- matrix(c(0, b_te_5, c_te_5, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_5, correct = TRUE) # continuity correction
```

```{r}
# filter rad = FP and cl = TN
fp_rad_not_cl_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_rad_te == 1, pred_cl_te == 0)
# filter cl = FP and rad = TN
fp_cl_not_rad_te <- combined_preds_te %>%
  dplyr::filter(Outcome_5PS_CAT == 0, pred_rad_te == 0, pred_cl_te == 1)

fp_rad_not_cl_te
fp_cl_not_rad_te
```

Comparing specificity: radiomic x clinical
```{r}
d_te_6 <- 2 # number of negative instances in the training set misclassified as FP by the radiomic model but not by the clinical model
e_te_6 <- 10 # number of negative instances in the training set misclassified as FP by the clinical model but not by the radiomic model

# Contingency mateix
contingency_matrix_6 <- matrix(c(0, d_te_6, e_te_6, 0), nrow = 2)

# McNemar.test()
mcnemar.test(contingency_matrix_6, correct = TRUE) # continuity correction

```



