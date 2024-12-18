---
title: "EDA & feature selection"
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
library(corrplot)
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

## Selecting only the training set for exploratory data analysis.
```{r}
# Training split.
set.seed(12)
data_split <- initial_split(data_nonnormal, prop = 0.7, strata = Esc_CAT)
data_tr_all <- training(data_split)

```

## Normality tests
```{r}
# Normality test (shapiro-wilk).
data_tr_all %>%
  select(SUVmax:SumVar) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))

# qq plots.
data_tr_all %>%
  pivot_longer(SUVmax:SumVar, names_to = "variable", values_to = "value") %>%
  ggplot(aes(sample = value)) +
  stat_qq(geom = "point", alpha = 0.2) +
  stat_qq_line(color = "black") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "QQ Plots of Radiomic Features", 
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") + 
  theme_pubr() +
  theme(
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),  
    axis.title = element_text(size = 10),         
    axis.text = element_text(size = 8)             
  )
```

## Association tests - 5-PS
```{r}
### 1: Association tests; 2: plot

## 1. Association tests

# changing col names
data_tr_all <- data_tr_all %>%
  rename("Staging" = estadiamento_cat)

# Defining the significance level
significance_level <- 0.05
log_significance_level <- -log10(significance_level)  # Value for the cutoff line (for the plot)

# Function for association tests
p_values_rad <- data_tr_all %>%
  summarise(across(v3_sexo:SumVar, function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ Esc_CAT, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ Esc_CAT, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, Esc_CAT)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

## 2. Creating a plot to visualize the p-values from association tests
# Transforming the dataframe to long format for visualization
p_values_long_rad <- p_values_rad %>%
  pivot_longer(everything(),              # Converts all columns to key-value pairs, creating a long-format dataframe
               names_to = "Variable",     # Renames the key column as "Variable"
               values_to = "p_value") %>% # Renames the value column as "p_value"
  mutate(
    significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),  # Labels each p-value as significant or not
    label = ifelse(p_value < significance_level, Variable, "")  # Adds labels for significant variables
  )

# Plotting p-values from association tests
# manh_plot <-
ggplot(p_values_long_rad, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +  # Creates a bar for each variable with height based on -log10(p-value)
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  # Adds a horizontal cutoff line for significance level
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  # Customizes bar colors by significance
  xlab("Variable") +  # X-axis label
  ylab("-log10(p-value)") +  # Y-axis label
  theme_pubr() +  # Applies a clean, publication-friendly theme
  theme(
    axis.text.x = element_blank(),   # Removes text on the X-axis to reduce clutter
    legend.position = "right",       # Positions the legend on the right side
    axis.ticks.x = element_blank()   # Removes X-axis ticks
  ) +  
  labs(fill = "") +  # Clears the legend title
  geom_text_repel(             # Adds labels using ggrepel to avoid overlap
    aes(label = label),        # Adds variable names as labels for significant p-values
    vjust = 0.8,               # Adjusts vertical spacing from the bars
    box.padding = 0.1,         # Increases padding around each label box
    point.padding = 0.1,       # Adds spacing between the label and bar
    angle = 90,                # Rotates the labels for better visibility
    nudge_y = 0.1,             # Nudges labels slightly upward from the bar
    nudge_x = -0.2,            # Shifts labels slightly leftward
    segment.color = "transparent",  # Removes the line connecting labels to points
    size = 3                   # Sets a smaller font size for labels
  ) 
# ggsave("manhattan_plot.png", plot = manh_plot, width = 15, height = 13, units = "cm", dpi = 600)

# changing col names
data_tr_all <- data_tr_all %>%
  rename("estadiamento_cat" = Staging)

```

## Correlation analysis
```{r}
# Correlation analysis.
data_tr_corr <- data_tr_all %>%
  dplyr::select(SUVmax:SumVar)

# Correlation matrix
corrplot(cor(data_tr_corr),
   method = "color",
   type = "lower",
   addCoef.col="black", 
   order = "AOE", 
   number.cex=0.3,
   tl.cex=0.7,
   tl.col = "black"
   )

# Network plot
correlate(data_tr_corr, method = "pearson")
cor_mat <- data_tr_corr %>%
  correlate() %>%    # Create correlation data frame (cor_df)
  rearrange()  # rearrange by correlations
fashion(cor_mat) 
network_plot(cor_mat, min_cor = 0.2)

```

## Univariate binary logistic regression between clinical-laboratory variables and glycolysis-related proteins was performed to select key predictors of categorized 5-PS.
```{r}
# Univariate logistic regression - 5-PS.
set.seed(12)
glm_fit_sexo <- 
  logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v3_sexo, data = data_tr_all)
tidy(glm_fit_sexo)

set.seed(12)
glm_fit_idade <- 
  logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v6_idade, data = data_tr_all)
tidy(glm_fit_idade)

set.seed(12)
glm_fit_grupo_etario <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v7_grupo_etario, data = data_tr_all)
tidy(glm_fit_grupo_etario)

set.seed(12)
glm_fit_estadiamento <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ estadiamento_cat, data = data_tr_all)
tidy(glm_fit_estadiamento)

set.seed(12)
glm_fit_ldh_elevada <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v10_ldh_elevada, data = data_tr_all)
tidy(glm_fit_ldh_elevada)

set.seed(12)
glm_fit_envo_sitios_extra <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v12_envo_sitios_extra, data = data_tr_all)
tidy(glm_fit_envo_sitios_extra)

set.seed(12)
glm_fit_pres_sintom_b <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v14_pres_sintom_b, data = data_tr_all)
tidy(glm_fit_pres_sintom_b)

set.seed(12)
glm_fit_doenca_vol_bulky <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v15_doenca_vol_bulky, data = data_tr_all)
tidy(glm_fit_pres_sintom_b)

set.seed(12)
glm_fit_glut1 <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_glut1_Tu_imput, data = data_tr_all)
tidy(glm_fit_glut1)

set.seed(12)
glm_fit_hkii <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_hkii_tu_imput, data = data_tr_all)
tidy(glm_fit_hkii)

set.seed(12)
glm_fit_ldh5 <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_ldh5_tu_imput, data = data_tr_all)
tidy(glm_fit_ldh5)

set.seed(12)
glm_fit_cd147 <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_cd147_Tu_imput, data = data_tr_all)
tidy(glm_fit_cd147)

set.seed(12)
glm_fit_CAIX <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_CAIX_Tumor, data = data_tr_all)
tidy(glm_fit_CAIX)

set.seed(12)
glm_fit_MCT4 <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_MCT4_tumor, data = data_tr_all)
tidy(glm_fit_MCT4)

set.seed(12)
glm_fit_GLUT1_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ GLUT1_stroma, data = data_tr_all)
tidy(glm_fit_GLUT1_stro)

set.seed(12)
glm_fit_HK2_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ HK2_stroma, data = data_tr_all)
tidy(glm_fit_HK2_stro)

set.seed(12)
glm_fit_LDHA_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ LDHA_stroma, data = data_tr_all)
tidy(glm_fit_LDHA_stro)

set.seed(12)
glm_fit_CAIX_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAIX_stroma, data = data_tr_all)
tidy(glm_fit_CAIX_stro)

set.seed(12)
glm_fit_MCT4_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ MCT4_stroma, data = data_tr_all)
tidy(glm_fit_MCT4_stro)

set.seed(12)
glm_fit_CD147_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CD147_stroma, data = data_tr_all)
tidy(glm_fit_CD147_stro)

set.seed(12)
glm_fit_GLUT3_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ GLUT3_stroma, data = data_tr_all)
tidy(glm_fit_GLUT3_stro)

```
Only staging was significant in the univariate logistic regression. No p-values < 0.100 were observed.

## Visualization of univariate logistic regression results.
```{r}
# Creating a list of logistic regression models.
glm_models <- list(
  sex = glm_fit_sexo,
  age = glm_fit_idade,
  age_group = glm_fit_grupo_etario,
  staging = glm_fit_estadiamento,
  LDH_group = glm_fit_ldh_elevada,
  extranodal_inv = glm_fit_envo_sitios_extra,
  B_symptoms = glm_fit_pres_sintom_b,
  bulky_disease = glm_fit_doenca_vol_bulky,
  HK2tu = glm_fit_hkii,
  LDH5tu = glm_fit_ldh5,
  MCT4tu = glm_fit_MCT4,
  GLUT1stro = glm_fit_GLUT1_stro,
  HK2stro = glm_fit_HK2_stro,
  LDHAstro = glm_fit_LDHA_stro,
  CAIXstro = glm_fit_CAIX_stro,
  MCT4stro = glm_fit_MCT4_stro,
  CD147stro = glm_fit_CD147_stro,
  GLUT3stro = glm_fit_GLUT3_stro
)

# Combining and tidying model results.
glm_tidy_combined <- bind_rows(lapply(names(glm_models), function(name) {
  tidy(glm_models[[name]], conf.int = TRUE) %>%
    mutate(model = name)
}))

# Adding custom labels for the staging model.
glm_tidy_combined$custom_label <- ifelse(glm_tidy_combined$model == "staging" & glm_tidy_combined$term == "(Intercept)", 
                                         "intercept", 
                                         ifelse(glm_tidy_combined$model == "staging" & glm_tidy_combined$term != "(Intercept)", 
                                                "staging advanced", 
                                                ""))

# Setting up a ggplot visualization.
  ggplot(glm_tidy_combined, aes(x = model, y = estimate, color = term)) +
  geom_point(position = position_dodge(width = 0.5)) +  
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.5)) +  
  geom_hline(yintercept = 0, color = "black") +   
  labs(
    x = "Model",  
    y = "log-odds"
  ) +
  theme_pubr() +
  ylim(-5.5, 5.5) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  gghighlight(
    model == "staging",
    unhighlighted_params = list(color = "#797A7B"),
    use_direct_label = FALSE  
  ) +
  scale_color_manual(values = c("skyblue4", "#FF8282")) +
  geom_text(
    data = subset(glm_tidy_combined, model == "staging" & term == "(Intercept)"),
    aes(label = "intercept"),
    vjust = 1.2, hjust = 1.2,
    size = 3,
    color = "skyblue4", fontface = "italic"
  ) +
  geom_label(
    data = subset(glm_tidy_combined, model == "staging" & term != "(Intercept)"),
    aes(label = "advanced stage"),
    vjust = -0.7, hjust = 1,
    size = 3,
    color = "#FF8282", fontface = "bold",
    label.size = 0.3,
    label.r = unit(0.2, "lines"),
    fill = "white"
  )

```

## Radiomic feature selection - Building the models for variable importance analysis.
```{r}
# Selecting only radiomic features as predictors.
data_nonnormal_rad <- data_nonnormal %>%
  dplyr::select(SUVmax:Esc_CAT)
```

## Splitting the dataset into training and testing sets
```{r}
# Training and testing split.
set.seed(12)
data_split <- initial_split(data_nonnormal_rad, prop = 0.7, strata = Esc_CAT)
data_tr_rad <- training(data_split)
data_te_rad <- testing(data_split)

# Cross-validation folds.
set.seed(12)
data_folds_vi <- vfold_cv(data_tr_rad, 5, strata = Esc_CAT)

```

## Recipes
```{r}
# Recipes:
# First, the recipe() needs to be informed about the model being used and the training data.
# Then, variables that are too correlated with each other should be filtered out with step_corr().
# Any numeric variables with zero variance are then removed.
# As the final step, numeric variables are normalized (centered and scaled) to account for differing scales, as the model being trained is sensitive to these differences.
# Lastly, the recipe() is prepped, meaning the steps are actually implemented with the training data, and the necessary parameters are estimated from the training data to apply this process to another dataset later.

rec_spec <- recipe(Esc_CAT ~ ., data_tr_rad) %>%
  step_corr(all_numeric_predictors(), threshold = 0.8) %>%
  step_zv(all_numeric()) %>%
  step_normalize(all_numeric())

rec_spec %>% prep() %>% juice() %>% glimpse()

```

## Model specification for variable importance analysis.
```{r}
# XGBoosting
xgb_spec <- boost_tree(
           trees = 500,
           learn_rate = tune(),
           min_n = tune()
       ) %>%
           set_engine("xgboost", scale_pos_weight = tune()) %>%
           set_mode("classification")

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

```

## Model metrics and control.
```{r}
# Model metrics.
model_metrics <- metric_set(roc_auc, accuracy)

# Model control.
model_control <- control_stack_grid()
```

## Tuning hyperparameters.
```{r}
# Creating a set of workflows.
model_set <- workflow_set(
  preproc = list(regular = rec_spec),
  models = list(xgboost = xgb_spec, randomForest = rf_spec, lasso = lasso_spec),
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
Random forest and LASSO classifiers showed similar results.

## Training the best models with the best hyperparameter values.
#### Random Forest
```{r}
# Training the best model with the best hyperparameter values.

rf_model <- rand_forest(
                  trees = 500,
                  min_n = 23
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")

set.seed(12)
rf_res <- fit_resamples(rf_model,
                                   rec_spec,
                                   resamples = data_folds_vi,
                                   control = control_resamples(save_pred = T))

```

#### LASSO
```{r}
# Training the best model with the best hyperparameter values.

lasso_model <- logistic_reg(
  penalty = 0.034629959,
  mixture = 1
  ) %>%
  set_engine("glmnet")

set.seed(12)
lasso_res <- fit_resamples(lasso_model,
                                   rec_spec,
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

## ROC_AUC LASSO
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

## Model comparison
Bayesian analysis was used here to answer the question: 'When looking at resampling results, are the differences between models real?' To address this, a model was created with the outcome being the resampling statistics (ROC_AUC).
```{r}
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

# Visualizing model comparison.
autoplot(class_post) + theme_minimal()
autoplot(contrast_models(class_post)) + theme_minimal()

```

## Random Forest variable importance.
```{r}
rf_model_imp <- rf_model %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(Esc_CAT ~ ., data = juice(prep(rec_spec)))

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
  add_recipe(rec_spec)

final_wf_lasso %>%
  fit(data_tr_rad) %>%
  pull_workflow_fit() %>%
  tidy()

final_wf_lasso %>%
  fit(data_tr_rad) %>%
  pull_workflow_fit() %>%
  vi()

final_wf_lasso %>%
  fit(data_tr_rad) %>%
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

## Density plots for the most important radiomic features.
```{r}
data_tr_rad$Esc_CAT <- fct_recode(data_tr_rad$Esc_CAT, "presence of metabolically active disease" = "0", "absence of metabolically active disease" = "1")

# dI <-
ggdensity(data_tr_rad, x = "IDM",
          add = "median", rug = TRUE, color = "Esc_CAT",
          fill = "Esc_CAT", palette = "Pastel1",
          xlab = "IDM", legend.title = "5-PS") + theme_pubr()
# ggsave("density_IDM.png", plot = dI, width = 12, height = 10, units = "cm", dpi = 600)

# dM <-
ggdensity(data_tr_rad, x = "MTV",
          add = "median", rug = TRUE, color = "Esc_CAT",
          fill = "Esc_CAT", palette = "Pastel1",
          xlab = "MTV", legend.title = "5-PS") + theme_pubr()
# ggsave("density_MTV.png", plot = dM, width = 12, height = 10, units = "cm", dpi = 600)

ggdensity(data_tr_rad, x = "Correlation",
          add = "median", rug = TRUE, color = "Esc_CAT",
          fill = "Esc_CAT", palette = "Pastel1",
          xlab = "Correlation", legend.title = "5-PS") + theme_pubr()

```
IDM and MTV values show notable differences between patients with and without hypermetabolic lesions on EoT PET, demonstrating good separability between the groups. In contrast, the correlation values show no such distinction.

## Dot plot with group means and reference line.
```{r}
# Calculate the overall mean of MTV
mean_MTV <- mean(data_tr_all$MTV)

ggplot(data_tr_all) + 
  geom_point(
    aes(x = MTV,
        y = score_deauville,
        color = score_deauville),
    alpha = 0.8,
    show.legend = FALSE
  ) +
  # Vertical line at the mean of MTV
  geom_vline(
    xintercept = mean_MTV,
    color = "grey27", size = 0.9
  ) +
  # Mean points of MTV for each score_deauville group with lines connecting to the overall mean
  stat_summary(
    aes(x = MTV, y = score_deauville, color = score_deauville),
    fun = mean,
    geom = "point",
    shape = 16,
    size = 5,
    alpha = 1,
    show.legend = FALSE
  ) +
  stat_summary(
    aes(x = MTV, y = score_deauville, color = score_deauville),
    fun = mean,
    geom = "segment",
    xend = mean_MTV,
    size = 2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  xlab("MTV") + 
  ylab("5-PS") +
  scale_color_viridis_d() +
  theme_minimal()

# Calculate the overall mean of IDM
mean_IDM <- mean(data_tr_all$IDM)

ggplot(data_tr_all) + 
  geom_point(
    aes(x = IDM,
        y = score_deauville,
        color = score_deauville),
    alpha = 0.8,
    show.legend = FALSE
  ) +
  # Vertical line at the mean of IDM
  geom_vline(
    xintercept = mean_IDM,
    color = "grey27", size = 0.9
  ) +
  # Mean points of IDM for each score_deauville group with lines connecting to the overall mean
  stat_summary(
    aes(x = IDM, y = score_deauville, color = score_deauville),
    fun = mean,
    geom = "point",
    shape = 16,
    size = 5,
    alpha = 1,
    show.legend = FALSE
  ) +
  stat_summary(
    aes(x = IDM, y = score_deauville, color = score_deauville),
    fun = mean,
    geom = "segment",
    xend = mean_IDM,
    size = 2,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  xlab("IDM") + 
  ylab("5-PS") +
  scale_color_viridis_d() +
  theme_minimal()
```
Patients with metabolically active disease tend to exhibit lower IDM values. Conversely, these patients typically show higher MTV values.

## The last_fit argument fits the model on the training set and evaluates the model on the testing set. The testing set can only be used once to avoid data leakage.
```{r}
#rf.
final_wf_lasso <- workflow() %>% 
  add_model(lasso_model) %>% 
  add_recipe(rec_spec)

final_res_lasso <- last_fit(final_wf_lasso, data_split) 

final_res_lasso %>%
  collect_metrics

final_res_lasso %>% 
  collect_predictions() %>%
  roc_curve(Esc_CAT, .pred_1, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

set.seed(12)
final_res_lasso_boot <- int_pctl(final_res_lasso, alpha = 0.05)
final_res_lasso_boot

```

