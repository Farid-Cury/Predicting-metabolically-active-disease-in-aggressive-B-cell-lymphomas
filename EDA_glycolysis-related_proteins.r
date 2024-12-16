---
title: "EDA glycolysis-related proteins"
author: "Farid"
date: "`r Sys.Date()`"
output: html_document
---
## Loading packages and dataset
```{r}
library(tidyverse)
library(ggpubr)
library(tidymodels)

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
```

## Selecting only the training set for exploratory data analysis
```{r}
# Training split.
set.seed(12)
data_split <- initial_split(data_nonnormal, prop = 0.7, strata = Esc_CAT)
data_tr_all <- training(data_split)

```

## HK2 expression score optimal cutpoint 
```{r}
library(cutpointr)
data_tr_all$Esc_CAT <- fct_recode(data_tr_all$Esc_CAT, "1-3" = "0", "presence of 4-5" = "1")
data_tr_all$v50_expres_hkii_tumo <- as.numeric(data_tr_all$v50_expres_hkii_tumo)

cutpoint_HK2 <- cutpointr(data_tr_all, v50_expres_hkii_tumo, Esc_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_HK2)
plot(cutpoint_HK2)

```

## LDHA expression score optimal cutpoint 
```{r}
data_tr_all$v56_expres_ldh5_tumo <- as.numeric(data_tr_all$v56_expres_ldh5_tumo)

cutpoint_LDHA <- cutpointr(data_tr_all, v56_expres_ldh5_tumo, Esc_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_LDHA)
plot(cutpoint_LDHA) 

```

## MCT4 expression score optimal cutpoint 
```{r}
data_tr_all$v74_expres_mct4_tumo <- as.numeric(data_tr_all$v74_expres_mct4_tumo)

cutpoint_MCT4 <- cutpointr(data_tr_all, v74_expres_mct4_tumo, Esc_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_MCT4)
plot(cutpoint_MCT4) 

```

## dichotomising expression scores acording to optimal cutpoints
```{r}
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


```{r}
data_tr_cp$Esc_CAT <- fct_recode(data_tr_cp$Esc_CAT, "0" = "1-3", "1" = "4-5")
```

## Association tests & plots - 5-PS
```{r}
## 1. Association tests

# loading package
library(ggrepel)

# changing col names
data_tr_cp <- data_tr_cp %>%
  rename("Staging" = estadiamento_cat)

# Defining the significance level
significance_level <- 0.05
log_significance_level <- -log10(significance_level)  # Value for the cutoff line (for the plot)

# Function for association tests
p_values_rad <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp, MCT4tu_opt_cp), function(.x) {
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

```

## Univariate logistic regression was performed to select key clinical-laboratory features and glycolysis-related protein expression *at their optimal cutpoints* as predictors of categorized 5-PS.
```{r}
set.seed(12)
glm_fit_HK2_cp <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ HK2tu_opt_cp, data = data_tr_cp)
tidy(glm_fit_HK2_cp)

set.seed(12)
glm_fit_LDHA_cp <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ LDHAtu_opt_cp, data = data_tr_cp)
tidy(glm_fit_LDHA_cp)

set.seed(12)
glm_fit_MCT4_cp <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ MCT4tu_opt_cp, data = data_tr_cp)
tidy(glm_fit_MCT4_cp)

set.seed(12)
glm_fit_sexo <- 
  logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v3_sexo, data = data_tr_cp)
tidy(glm_fit_sexo)

set.seed(12)
glm_fit_idade <- 
  logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v6_idade, data = data_tr_cp)
tidy(glm_fit_idade)

set.seed(12)
glm_fit_grupo_etario <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v7_grupo_etario, data = data_tr_cp)
tidy(glm_fit_grupo_etario)

set.seed(12)
glm_fit_estadiamento <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ Staging, data = data_tr_cp)
tidy(glm_fit_estadiamento)

set.seed(12)
glm_fit_ldh_elevada <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v10_ldh_elevada, data = data_tr_cp)
tidy(glm_fit_ldh_elevada)

set.seed(12)
glm_fit_envo_sitios_extra <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v12_envo_sitios_extra, data = data_tr_cp)
tidy(glm_fit_envo_sitios_extra)

set.seed(12)
glm_fit_pres_sintom_b <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v14_pres_sintom_b, data = data_tr_cp)
tidy(glm_fit_pres_sintom_b)

set.seed(12)
glm_fit_doenca_vol_bulky <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ v15_doenca_vol_bulky, data = data_tr_cp)
tidy(glm_fit_pres_sintom_b)

set.seed(12)
glm_fit_glut1_CAT <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_glut1_Tu_imput, data = data_tr_cp)
tidy(glm_fit_glut1_CAT)

set.seed(12)
glm_fit_hkii_CAT <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_hkii_tu_imput, data = data_tr_cp)
tidy(glm_fit_hkii_CAT)

set.seed(12)
glm_fit_ldh5_CAT <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_ldh5_tu_imput, data = data_tr_cp)
tidy(glm_fit_ldh5_CAT)

set.seed(12)
glm_fit_cd147_CAT <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_cd147_Tu_imput, data = data_tr_cp)
tidy(glm_fit_cd147_CAT)

set.seed(12)
glm_fit_CAIX_CAT <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_CAIX_Tumor, data = data_tr_cp)
tidy(glm_fit_CAIX_CAT)

set.seed(12)
glm_fit_MCT4_CAT <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAT_MCT4_tumor, data = data_tr_cp)
tidy(glm_fit_MCT4_CAT)

set.seed(12)
glm_fit_GLUT1_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ GLUT1_stroma, data = data_tr_cp)
tidy(glm_fit_GLUT1_stro)

set.seed(12)
glm_fit_HK2_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ HK2_stroma, data = data_tr_cp)
tidy(glm_fit_HK2_stro)

set.seed(12)
glm_fit_LDHA_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ LDHA_stroma, data = data_tr_cp)
tidy(glm_fit_LDHA_stro)

set.seed(12)
glm_fit_CAIX_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CAIX_stroma, data = data_tr_cp)
tidy(glm_fit_CAIX_stro)

set.seed(12)
glm_fit_MCT4_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ MCT4_stroma, data = data_tr_cp)
tidy(glm_fit_MCT4_stro)

set.seed(12)
glm_fit_CD147_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ CD147_stroma, data = data_tr_cp)
tidy(glm_fit_CD147_stro)

set.seed(12)
glm_fit_GLUT3_stro <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ GLUT3_stroma, data = data_tr_cp)
tidy(glm_fit_GLUT3_stro)

```

## Visualization of univariate logistic regression results.
```{r}
# Creating a list of logistic regression models.
glm_models <- list(
  HK2tu_opt_cp = glm_fit_HK2_cp,
  LDHAtu_opt_cp = glm_fit_LDHA_cp,
  MCT4tu_opt_cp = glm_fit_MCT4_cp,
  sex = glm_fit_sexo,
  age = glm_fit_idade,
  age_group = glm_fit_grupo_etario,
  staging = glm_fit_estadiamento,
  LDH_group = glm_fit_ldh_elevada,
  extranodal_inv = glm_fit_envo_sitios_extra,
  B_symptoms = glm_fit_pres_sintom_b,
  bulky_disease = glm_fit_doenca_vol_bulky,
  HK2tu_CAT = glm_fit_hkii_CAT,
  LDH5tu_CAT = glm_fit_ldh5_CAT,
  MCT4tu_CAT = glm_fit_MCT4_CAT,
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
library(gghighlight)
glm_tidy_combined$custom_label <- ifelse(glm_tidy_combined$model == "staging" & glm_tidy_combined$term == "(Intercept)", 
                                         "intercept", 
                                         ifelse(glm_tidy_combined$model == "staging" & glm_tidy_combined$term != "(Intercept)", 
                                                "staging advanced", 
                                                ""))

# Setting up a ggplot visualization.
#lr_all_plot <-
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
    model == "staging" | model == "HK2tu_opt_cp",
    unhighlighted_params = list(color = "#797A7B"),
    use_direct_label = FALSE  # Desativa rótulos automáticos
  ) +
  scale_color_manual(values = c("skyblue4", "#FF8282", "#FF8282")) +
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
  ) +
  geom_text(
    data = subset(glm_tidy_combined, model == "HK2tu_opt_cp" & term == "(Intercept)"),
    aes(label = "intercept"),
    vjust = 1.5, hjust = 0,
    size = 3,
    color = "skyblue4", fontface = "italic"
  ) +
  geom_label(
    data = subset(glm_tidy_combined, model == "HK2tu_opt_cp" & term != "(Intercept)"),
    aes(label = "positive HK2tu"),
    vjust = -0.5, hjust = -0.1,
    size = 3,
    color = "#FF8282", fontface = "bold",
    label.size = 0.3,
    label.r = unit(0.2, "lines"),
    fill = "white"
  )
#ggsave("lr_all_plot.png", plot = lr_all_plot, width = 18, height = 10, units = "cm")

```

## Multivariate logistic regression - using HK2tu + staging to predict 5-PS
```{r}
set.seed(12)
glm_fit_HK2cp_staging <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Esc_CAT ~ HK2tu_opt_cp + Staging, data = data_tr_cp)
tidy(glm_fit_HK2cp_staging)

```

## Exploring associations between radiomic features and glycolysis-related protein expression - tumor compartment.
```{r}
# HK2
p_values_hk2tu <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, LDHAtu_opt_cp, MCT4tu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ HK2tu_opt_cp, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ HK2tu_opt_cp, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, HK2tu_opt_cp)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_hk2tu <- p_values_hk2tu %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_hk2tu, aes(x = Variable, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +  
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  
  xlab("Variable") +  
  ylab("-log10(p-value)") +  
  theme_pubr() +  
  theme(
    axis.text.x = element_blank(),   
    legend.position = "right",       
    axis.ticks.x = element_blank()   
  ) +  
  labs(fill = "") +  
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,               
    box.padding = 0.1,         
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = -0.2,            
    segment.color = "transparent",  
    size = 3                   
  ) 

# LDHA
p_values_ldhatu <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, MCT4tu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ LDHAtu_opt_cp, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ LDHAtu_opt_cp, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, LDHAtu_opt_cp)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_ldhatu <- p_values_ldhatu %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ldhatu, aes(x = Variable, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +  
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  
  xlab("Variable") +  
  ylab("-log10(p-value)") +  
  theme_pubr() +  
  theme(
    axis.text.x = element_blank(),   
    legend.position = "right",       
    axis.ticks.x = element_blank()   
  ) +  
  labs(fill = "") +  
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,               
    box.padding = 0.1,         
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = -0.2,
    segment.color = "transparent",  
    size = 3                   
  ) 

# MCT4
p_values_mct4tu <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ MCT4tu_opt_cp, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ MCT4tu_opt_cp, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, MCT4tu_opt_cp)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_mct4tu <- p_values_mct4tu %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_mct4tu, aes(x = Variable, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +  
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  
  xlab("Variable") +  
  ylab("-log10(p-value)") +  
  theme_pubr() +  
  theme(
    axis.text.x = element_blank(),   
    legend.position = "right",       
    axis.ticks.x = element_blank()   
  ) +  
  labs(fill = "") +  
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,               
    box.padding = 0.1,         
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = -0.2,            
    segment.color = "transparent",  
    size = 3                   
  ) 

```

## Exploring associations between radiomic features and glycolisis-related protein expression - stromal compartment.
```{r}
# GLUT1 stroma
p_values_glut1stro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ GLUT1_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ GLUT1_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, GLUT1_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_glut1stro <- p_values_glut1stro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_glut1stro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# HK2 stroma
p_values_hk2stro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ HK2_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ HK2_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, HK2_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_hk2stro <- p_values_hk2stro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_hk2stro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# LDHA stroma
p_values_ldhastro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ LDHA_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ LDHA_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, LDHA_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_ldhastro <- p_values_ldhastro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ldhastro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# CAIX stroma
p_values_ca9stro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ CAIX_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ CAIX_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, CAIX_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_ca9stro <- p_values_ca9stro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ca9stro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# MCT4 stroma
p_values_mct4stro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ MCT4_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ MCT4_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, MCT4_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_mct4stro <- p_values_mct4stro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_mct4stro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# CD147 stroma
p_values_cd147stro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ CD147_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ CD147_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, CD147_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_cd147stro <- p_values_cd147stro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_cd147stro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# GLUT3 stroma
p_values_glut3stro <- data_tr_cp %>%
  summarise(across(c(v3_sexo:SumVar, HK2tu_opt_cp, LDHAtu_opt_cp), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ GLUT3_stroma, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ GLUT3_stroma, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, GLUT3_stroma)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))

# Transforming the dataframe to long format for visualization
p_values_long_glut3stro <- p_values_glut3stro %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_glut3stro, aes(x = Variable, y = -log10(p_value), fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Hides the text on the X-axis
  labs(fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

```

## Extranodal site involvement predicting CAIX stromal expression
```{r}
set.seed(12)
glm_fit_CAIXstro_extranod <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(CAIX_stroma ~ v12_envo_sitios_extra, data = data_tr_cp)
tidy(glm_fit_CAIXstro_extranod)

glm_tidy_CAIXstro_extranod <- tidy(glm_fit_CAIXstro_extranod, conf.int = TRUE)

ggplot(glm_tidy_CAIXstro_extranod, aes(x = term, y = estimate)) +
  geom_point() +  
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +  
  geom_hline(yintercept = 0, color = "black") +   
  labs(
       x = "Variable",
       y = "log-odds") +
  theme_pubr()

```

## Co-expression analysis - score
In this analysis, scores were created to evaluate the co-expression of proteins in the tumor tissue, stroma, and their combination. Below is a detailed methodology:

1. Tumoral Co-expression Score
A 7-level score (0 to 6) was developed to assess protein expression in the tumor.
The score was calculated by summing the number of positively expressed proteins in the tumor:
- 0 proteins expressed: Score = 0
- 1 protein expressed: Score = 1
- 2 proteins expressed: Score = 2
- And so on, up to a maximum of 6 proteins.

The proteins considered were:
GLUT1, CD147, CAIX, HK2, LDHA, MCT4 (all evaluated in the tumor).

2. Stromal Co-expression Score
An 8-level score (0 to 7) was developed to assess protein expression in the stroma.
The calculation followed the same principle as the tumor expression:
- 0 proteins expressed: Score = 0
- 1 protein expressed: Score = 1
- 2 proteins expressed: Score = 2
- And so on, up to a maximum of 7 proteins.

The proteins included were:
GLUT1, GLUT3, CD147, CAIX, HK2, LDHA, MCT4 (all evaluated in the stroma).

3. Combined Tumor and Stroma Co-expression Score
A 14-level score (0 to 13) was developed to assess the combined co-expression between the tumor and stroma.
The sum considered protein expression in both compartments (tumor and stroma):
- 0 proteins expressed in the tumor and 0 in the stroma: Score = 0
- 1 protein expressed in the tumor and 0 in the stroma, or 0 in the tumor and 1 in the stroma: Score = 1
- 1 protein expressed in the tumor and 1 in the stroma: Score = 2
- And so on, up to a maximum of 13 proteins expressed in total.

The proteins included were the same as those in the individual analyses:
Tumor: GLUT1, CD147, CAIX, HK2, LDHA, MCT4.
Stroma: GLUT1, GLUT3, CD147, CAIX, HK2, LDHA, MCT4.
```{r}
data_tr_coexp <- data_tr_cp %>%
  dplyr::mutate(
    # Tumor coexpression score
    tu_coex_score = rowSums(across(c(CAT_glut1_Tu_imput, CAT_cd147_Tu_imput, CAT_CAIX_Tumor, HK2tu_opt_cp, LDHAtu_opt_cp, MCT4tu_opt_cp)) == 1),
    tu_coex_score = factor(tu_coex_score, levels = 0:6),
    # Stroma coexpression score
    stro_coex_score = rowSums(across(c(GLUT1_stroma, GLUT3_stroma, CD147_stroma, CAIX_stroma, HK2_stroma, LDHA_stroma, MCT4_stroma)) == 1),
    stro_coex_score = factor(stro_coex_score, levels = 0:7),
    # Tumor and stroma combined score
    tu_stro_coex_score = rowSums(across(c(
      CAT_glut1_Tu_imput, CAT_cd147_Tu_imput, CAT_CAIX_Tumor, HK2tu_opt_cp, LDHAtu_opt_cp, MCT4tu_opt_cp, 
      GLUT1_stroma, GLUT3_stroma, CD147_stroma, CAIX_stroma, HK2_stroma, LDHA_stroma, MCT4_stroma)) == 1),
    tu_stro_coex_score = factor(tu_stro_coex_score, levels = 0:13)
  )

str(data_tr_coexp)

```

## Tumor co-expression score visualization
```{r}
# Reclassify the levels of 'Esc_CAT' to make them more interpretable
data_tr_coexp$Esc_CAT <- fct_recode(data_tr_coexp$Esc_CAT, "1-3" = "0", "4-5" = "1")

# Count the frequency of combinations between 'tu_coex_score' and 'Esc_CAT'
data_tr_coexp_tu_freq <- data_tr_coexp %>%
  dplyr::count(tu_coex_score, Esc_CAT, name = "freq")

# Create a scatter plot using ggplot2
ggplot(data_tr_coexp_tu_freq, aes(x = tu_coex_score, y = Esc_CAT, color = Esc_CAT, size = freq)) +
  geom_point() +
  xlab("Tumor expression score") +
  ylab("5-PS") +
  scale_fill_manual(values = c("grey60", "grey40", "grey20")) +
  theme_pubr() +
  border() +
  labs(color = "5-PS") +
  scale_size_continuous(name = "Frequency", range = c(3, 10)) +
  theme(legend.position = "right")
  
```

## Stroma co-expression score visualization
```{r}
# Count the frequency of combinations between 'stro_coex_score' and 'Esc_CAT'
data_tr_coexp_stro_freq <- data_tr_coexp %>%
  dplyr::count(stro_coex_score, Esc_CAT, name = "freq")

# Create a scatter plot using ggplot2
ggplot(data_tr_coexp_stro_freq, aes(x = stro_coex_score, y = Esc_CAT, color = Esc_CAT, size = freq)) +
  geom_point() +
  xlab("Stroma expression score") +
  ylab("5-PS") +
  scale_fill_manual(values = c("grey60", "grey40", "grey20")) +
  theme_pubr() +
  border() +
  labs(color = "5-PS") +
  scale_size_continuous(name = "Frequency", range = c(3, 10)) +
  theme(legend.position = "right")
  
```

## Tumor & stroma co-expression score visualization
```{r}
# Count the frequency of combinations between 'tu_stro_coex_score' and 'Esc_CAT'
data_tr_coexp_tu_stro_freq <- data_tr_coexp %>%
  dplyr::count(tu_stro_coex_score, Esc_CAT, name = "freq")

# Create a scatter plot using ggplot2
ggplot(data_tr_coexp_tu_stro_freq, aes(x = tu_stro_coex_score, y = Esc_CAT, color = Esc_CAT, size = freq)) +
  geom_point() +
  xlab("Tumor & stroma expression score") +
  ylab("5-PS") +
  scale_fill_manual(values = c("grey60", "grey40", "grey20")) +
  theme_pubr() +
  border() +
  labs(color = "5-PS") +
  scale_size_continuous(name = "Frequency", range = c(3, 10)) +
  theme(legend.position = "right")
  
```

## Combined visualization
```{r}
data_tr_long_scores <- data_tr_coexp %>%
  dplyr::select(tu_coex_score, stro_coex_score,tu_stro_coex_score, Esc_CAT) %>%
  pivot_longer(cols = c(tu_coex_score:tu_stro_coex_score), names_to = "Variable", values_to = "Expression") %>%
  group_by(Variable, Expression)

data_tr_coexp_scores_freq <- data_tr_long_scores %>%
  dplyr::count(Expression, Esc_CAT, name = "freq")

ggplot(data_tr_coexp_scores_freq, aes(x = Expression, y = Esc_CAT, size = freq, color = Esc_CAT)) +
  geom_point() +
  xlab("Expression Score") +
  ylab("5-PS") +
  facet_wrap(vars(Variable), scales = "free_x") +
  scale_fill_manual(values = c("grey60", "grey40", "grey20")) +
  theme_pubr() +
  border() +
  labs(color = "5-PS") +
  scale_size_continuous(name = "Frequency", range = c(3, 10)) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8))
```

## Association tests between co-expression scores and 5-PS
```{r}
data_tr_coexp$tu_coex_score <- as.numeric(data_tr_coexp$tu_coex_score)
data_tr_coexp$stro_coex_score <- as.numeric(data_tr_coexp$stro_coex_score)
data_tr_coexp$tu_stro_coex_score <- as.numeric(data_tr_coexp$tu_stro_coex_score)

p_values_coexp <- data_tr_coexp %>%
  summarise(across(c(tu_coex_score:tu_stro_coex_score), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ Esc_CAT, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ Esc_CAT, alternative = "two.sided")$p.value
      }
    } 
    return(p_value)
  }))
p_values_coexp

```

## Correlation analysis between coexpression scores and radiomic features.
```{r}
data_tr_corr <- data_tr_coexp %>%
  dplyr::select(SUVmax:SumVar | tu_coex_score:tu_stro_coex_score)

data_tr_corr$tu_coex_score <- as.numeric(data_tr_corr$tu_coex_score)
data_tr_corr$stro_coex_score <- as.numeric(data_tr_corr$stro_coex_score)
data_tr_corr$tu_stro_coex_score <- as.numeric(data_tr_corr$tu_stro_coex_score)

library(corrr)
correlate(data_tr_corr, method = "pearson")

cor_mat <- data_tr_corr %>%
  correlate() %>%    # Create correlation data frame (cor_df)
  rearrange()  # rearrange by correlations
fashion(cor_mat) 

# p <-
network_plot(cor_mat, min_cor = 0.2)
# ggsave("network_plot.png", plot = p, width = 12, height = 12, units = "cm", dpi = 600) # <- best parameters for saving the network plot.

```
## Co-expression analysis - evaluating the significantly associated co-expressions
This co-expression analysis was based on the results of association tests evaluating the expression of specific proteins. Only proteins whose expression showed statistically significant associations were included in the scores.

1. Tumor Co-expression Associated Score (tu_coexp_assoc)
Tumor-related proteins significantly associated: HK2, LDHA, and MCT4.
If any of these proteins (HK2, LDHA, or MCT4) is not expressed (value = 0), the score is set to 0 (absence of associated co-expression).
If any of these proteins is expressed (value = 1), the score is set to 1 (presence of associated co-expression).

2. Stroma Co-expression Associated Score (stro_coexp_assoc)
CAIX and CD147 were significantly associated.
If either CAIX or CD147 is not expressed (value = 0), the score is set to 0 (absence of associated co-expression).
If either CAIX or CD147 is expressed (value = 1), the score is set to 1 (presence of associated co-expression).

3. Combined Tumor and Stroma Co-expression Associated Score (tu_stro_coexp_assoc)
Represents the co-expression of one significantly associated protein from the tumor compartment (LDHA) and one from the stroma (MCT4).
If either MCT4 (stroma) or LDHA (tumor) is not expressed (value = 0), the score is set to 0 (absence of associated co-expression).
If either MCT4 (stroma) or LDHA (tumor) is expressed (value = 1), the score is set to 1 (presence of associated co-expression).
```{r}
data_tr_coexp_assoc <- data_tr_coexp %>%
  dplyr::mutate(
    tu_coexp_assoc = case_when(
      HK2tu_opt_cp == 0 | LDHAtu_opt_cp == 0 | MCT4tu_opt_cp == 0 ~ 0,
      HK2tu_opt_cp == 1 | LDHAtu_opt_cp == 1 | MCT4tu_opt_cp == 1 ~ 1
    ),
    stro_coexp_assoc = case_when(
      CAIX_stroma == 0 | CD147_stroma == 0 ~ 0,
      CAIX_stroma == 1 | CD147_stroma == 1 ~ 1
    ),
    tu_stro_coexp_assoc = case_when(
      MCT4_stroma == 0 | LDHAtu_opt_cp == 0 ~ 0,
      MCT4_stroma == 1 | LDHAtu_opt_cp == 1 ~ 1
    )
  )

View(data_tr_coexp_assoc)
```

## Visualization of association tests between significantly associated tumor co-expressions and clinical-laboratory and radiomic features.
```{r}
data_tr_coexp_assoc <- data_tr_coexp_assoc %>%
  mutate(across(tu_coexp_assoc:tu_stro_coexp_assoc, as.factor))

p_values_tu_coexp_assoc <- data_tr_coexp_assoc %>%
  summarise(across(c(v3_sexo:v15_doenca_vol_bulky, SUVmax:Esc_CAT), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ tu_coexp_assoc, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ tu_coexp_assoc, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, tu_coexp_assoc)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))
p_values_tu_coexp_assoc

# Transforming the dataframe to long format for visualization
p_values_long_tu_coexp_assoc <- p_values_tu_coexp_assoc %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_tu_coexp_assoc, aes(x = Variable, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +  
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  
  xlab("Variable") +  
  ylab("-log10(p-value)") +  
  theme_pubr() +  
  theme(
    axis.text.x = element_blank(),   
    legend.position = "right",       
    axis.ticks.x = element_blank()   
  ) +  
  labs(fill = "") +  
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,               
    box.padding = 0.1,         
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = -0.2,            
    segment.color = "transparent",  
    size = 3                   
  ) 

```

## Visualization of association tests between significantly associated stroma co-expressions and clinical-laboratory and radiomic features.
```{r}
p_values_stro_coexp_assoc <- data_tr_coexp_assoc %>%
  summarise(across(c(v3_sexo:v15_doenca_vol_bulky, SUVmax:Esc_CAT), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ stro_coexp_assoc, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ stro_coexp_assoc, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, stro_coexp_assoc)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))
p_values_stro_coexp_assoc

# Transforming the dataframe to long format for visualization
p_values_long_stro_coexp_assoc <- p_values_stro_coexp_assoc %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_stro_coexp_assoc, aes(x = Variable, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +  
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  
  xlab("Variable") +  
  ylab("-log10(p-value)") +  
  theme_pubr() +  
  theme(
    axis.text.x = element_blank(),   
    legend.position = "right",       
    axis.ticks.x = element_blank()   
  ) +  
  labs(fill = "") +  
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,               
    box.padding = 0.1,         
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = -0.2,            
    segment.color = "transparent",  
    size = 3                   
  )

```

## Visualization of association tests between significantly associated tumor and stroma co-expressions and clinical-laboratory and radiomic features.
```{r}
p_values_tu_stro_coexp_assoc <- data_tr_coexp_assoc %>%
  summarise(across(c(v3_sexo:v15_doenca_vol_bulky, SUVmax:Esc_CAT), function(.x) {
    if (is.numeric(.x)) {  
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        p_value <- wilcox.test(.x ~ tu_stro_coexp_assoc, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        p_value <- t.test(.x ~ tu_stro_coexp_assoc, alternative = "two.sided")$p.value
      }
    } else {  
      contingency_table <- table(.x, tu_stro_coexp_assoc)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    return(p_value)
  }))
p_values_tu_stro_coexp_assoc

# Transforming the dataframe to long format for visualization
p_values_long_tu_stro_coexp_assoc <- p_values_tu_stro_coexp_assoc %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "p_value") %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_tu_stro_coexp_assoc, aes(x = Variable, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +  
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +  
  xlab("Variable") +  
  ylab("-log10(p-value)") +  
  theme_pubr() +  
  theme(
    axis.text.x = element_blank(),   
    legend.position = "right",       
    axis.ticks.x = element_blank()   
  ) +  
  labs(fill = "") +  
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,               
    box.padding = 0.1,         
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = -0.2,            
    segment.color = "transparent",  
    size = 3                   
  )
```

