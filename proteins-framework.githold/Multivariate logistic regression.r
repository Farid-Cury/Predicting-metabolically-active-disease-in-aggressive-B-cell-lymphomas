---
title: "Glycolysis-related proteins framework"
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

## Removing NAs.
```{r}
data_all <- recipe(Outcome_5PS_CAT ~., data = data_nonnormal) %>%
  step_naomit(all_nominal_predictors()) %>%
  prep() %>%
  bake(new_data = NULL) %>%
  glimpse()
```

### Correlation analysis & plots - semantic and texture-based features x tumor expression score
```{r}
# Correlation analysis - Spearman.
data_all <- data_all %>%
  dplyr::mutate(across(c(Tumor_expres_hk2,
                         Tumor_expres_ldha,
                         Tumor_expres_mct4,
                         Tumor_expres_glut1,
                         Tumor_expres_ca9,
                         Tumor_expres_cd147), 
                       ~ as.numeric(as.character(.))))

data_corr_all <- data_all %>%
  dplyr::select(Tumor_expres_hk2:Tumor_expres_cd147, Semantic_SUVmax:Texture_SumVar)

# Correlation matrix and network plot
cor_mat_tumor_score <- data_corr_all %>%
  correlate(method = "spearman") 

# Focusing on co-expression score associations with radiomic features
cor_focus_tumor_score <- cor_mat_tumor_score %>%
  focus(Tumor_expres_hk2:Tumor_expres_cd147)

# Transforming the dataset into long format
cor_long_tumor_score <- cor_focus_tumor_score %>%
  pivot_longer(cols = -term, names_to = "variable", values_to = "correlation") %>%
  mutate(term = reorder(term, correlation))  # Ordenando pelas correlações

# Bar plot
ggplot(cor_long_tumor_score, aes(x = term, y = correlation, fill = variable)) +
  geom_col(position = "dodge") +  
  coord_flip() + 
  labs(
    x = "Radiomic features",
    y = "r (Spearman)",
    fill = "Variable",
  ) +
  scale_fill_viridis_d() +
  theme_pubr() +
  theme(legend.position = "right")

```

## HK2 expression score optimal cutpoint 
```{r}
library(cutpointr)
data_all$Outcome_5PS_CAT <- fct_recode(data_all$Outcome_5PS_CAT, "1-3" = "0", "4-5" = "1")

cutpoint_HK2 <- cutpointr(data_all, Tumor_expres_hk2, Outcome_5PS_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_HK2)

dist_hk2 <- plot_x(cutpoint_HK2) +
  theme_pubr() +
  border() +
  theme(legend.position = "none",
        title = element_blank()) +
  coord_cartesian(xlim = c(-0.5, 6.5)) +
  coord_cartesian(ylim = c(0, 60))

roc_hk2 <- plot_roc(cutpoint_HK2) +
  geom_path(alpha = 0.5, size = 1.2, color = "skyblue4") +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() +
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
  
  
```

## LDHA expression score optimal cutpoint 
```{r}
cutpoint_LDHA <- cutpointr(data_all, Tumor_expres_ldha, Outcome_5PS_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_LDHA)

dist_ldha <- plot_x(cutpoint_LDHA) +
  theme_pubr() +
  border() +
  theme(legend.position = "none",
        title = element_blank()) +
  coord_cartesian(xlim = c(-0.5, 6.5)) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

roc_ldha <- plot_roc(cutpoint_LDHA) +
  geom_path(alpha = 0.5, size = 1.2, color = "firebrick4") +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() +
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
```

## MCT4 expression score optimal cutpoint 
```{r}
cutpoint_MCT4 <- cutpointr(data_all, Tumor_expres_mct4, Outcome_5PS_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_MCT4)

dist_mct4 <- plot_x(cutpoint_MCT4) +
  theme_pubr() +
  border() +
  theme(legend.position = "none",
        title = element_blank()) +
  coord_cartesian(xlim = c(-0.5, 6.5)) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

roc_mct4 <- plot_roc(cutpoint_MCT4) +
  geom_path(alpha = 0.5, size = 1.2, color = "darkgreen") +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() +
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
```

## GLUT1 expression score optimal cutpoint 
```{r}
cutpoint_GLUT1 <- cutpointr(data_all, Tumor_expres_glut1, Outcome_5PS_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_GLUT1)
dist_glut1 <- plot_x(cutpoint_GLUT1) +
  theme_pubr() +
  border() +
  theme(legend.position = "none",
        title = element_blank()) +
  coord_cartesian(xlim = c(-0.5, 6.5)) +
  coord_cartesian(ylim = c(0, 60))

roc_glut1 <- plot_roc(cutpoint_GLUT1) +
  geom_path(alpha = 0.5, size = 1.2, color = "grey10") +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() +
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
```

## CA9 expression score optimal cutpoint 
```{r}
cutpoint_CA9 <- cutpointr(data_all, Tumor_expres_ca9, Outcome_5PS_CAT, 
                method = maximize_metric, metric = youden)
summary(cutpoint_CA9)

dist_ca9 <- plot_x(cutpoint_CA9) +
  theme_pubr() +
  border() +
  theme(legend.position = "none",
        title = element_blank()) +
  coord_cartesian(xlim = c(-0.5, 6.5)) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

roc_ca9 <- plot_roc(cutpoint_CA9) +
  geom_path(alpha = 0.5, size = 1.2, color = "deeppink") +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() +
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
```

## CD147 expression score optimal cutpoint 
```{r}
cutpoint_CD147 <- cutpointr(data_all, Tumor_expres_cd147, Outcome_5PS_CAT, 
                method = maximize_metric, metric = sum_sens_spec)
summary(cutpoint_CD147)

dist_cd147 <- plot_x(cutpoint_CD147) +
  theme_pubr() +
  border() +
  theme(legend.position = "none",
        title = element_blank()) +
  scale_x_continuous(breaks = seq(0, 6)) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

roc_cd147 <- plot_roc(cutpoint_CD147) +
  geom_path(alpha = 0.5, size = 1.2, color = "#BA9A66") +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr() +
  theme(title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 
```
```{r}
library(patchwork)

cutpoint_rocs <- roc_hk2 + roc_ldha + roc_mct4 + roc_glut1 + roc_ca9 + roc_cd147
print(cutpoint_rocs)

cutpoint_dist <- dist_hk2 + dist_ldha + dist_mct4 + dist_glut1 + dist_ca9 + dist_cd147
print(cutpoint_dist)

```


## Dichotomising expression scores acording to optimal cutpoints
```{r}
data_cp <- data_all %>%
  dplyr::mutate(Tumor_HK2 = case_when(
    Tumor_expres_hk2 <= 4 ~ 0,
    Tumor_expres_hk2 >= 5 ~ 1
  ),
  Tumor_LDHA = case_when(
    Tumor_expres_ldha <= 3 ~ 0,
    Tumor_expres_ldha >= 4 ~ 1
  ),
  Tumor_MCT4 = case_when(
    Tumor_expres_mct4 <= 5 ~ 0,
    Tumor_expres_mct4 == 6 ~ 1
  ),
  Tumor_GLUT1 = case_when(
    Tumor_expres_glut1 <= 1 ~ 0,
    Tumor_expres_glut1 >= 2 ~ 1
  ),
  Tumor_CA9 = case_when(
    Tumor_expres_ca9 <= 1 ~ 0,
    Tumor_expres_ca9 >= 2 ~ 1
  ),
  Tumor_CD147 = case_when(
    Tumor_expres_cd147 <= 1 ~ 0,
    Tumor_expres_cd147 >= 2 ~ 1)
  )

data_cp$Tumor_HK2 <- as.factor(data_cp$Tumor_HK2)
data_cp$Tumor_LDHA <- as.factor(data_cp$Tumor_LDHA)
data_cp$Tumor_MCT4 <- as.factor(data_cp$Tumor_MCT4)
data_cp$Tumor_GLUT1 <- as.factor(data_cp$Tumor_GLUT1)
data_cp$Tumor_CA9 <- as.factor(data_cp$Tumor_CA9)
data_cp$Tumor_CD147 <- as.factor(data_cp$Tumor_CD147)
data_cp$Outcome_5PS_CAT <- fct_recode(data_cp$Outcome_5PS_CAT, "0" = "1-3", "1" = "4-5")
```

## Selecting features for association tests.
```{r}
data_assoc <- data_cp %>%
  select(CL_Sex:CL_Bulky_disease, Stroma_GLUT1:Texture_SumVar, Outcome_5PS_CAT:Tumor_CD147) %>%
  relocate(Tumor_HK2:Tumor_CD147, .after = CL_Bulky_disease)
```

## Exploring associations between clinical-laboratory-radiomic features and glycolysis-related protein expression - tumor compartment.
```{r}
## 1. Association tests

# loading package
library(ggrepel)

# Defining the significance level
significance_level <- 0.05
log_significance_level <- -log10(significance_level)  # Value for the cutoff line (for the plot)

# HK2
p_values_hk2tu <- data_assoc %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Tumor_LDHA:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Tumor_HK2, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Tumor_HK2, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Tumor_HK2)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_hk2tu <- p_values_hk2tu %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_HK2tu",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_hk2tu)

# Creating columns "significant" and "not significant"
p_values_long_hk2tu <- results_table_hk2tu %>%
  select(Variable_HK2tu, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_HK2tu, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_hk2tu, aes(x = Variable_HK2tu, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of HK2 tumor expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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
p_values_ldhatu <- data_assoc %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Tumor_MCT4:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Tumor_LDHA, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Tumor_LDHA, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Tumor_LDHA)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_ldhatu <- p_values_ldhatu %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_LDHAtu",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_ldhatu)

# Creating columns "significant" and "not significant"
p_values_long_ldhatu <- results_table_ldhatu %>%
  select(Variable_LDHAtu, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_LDHAtu, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ldhatu, aes(x = Variable_LDHAtu, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of LDHA tumor expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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
p_values_mct4tu <- data_assoc %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Tumor_GLUT1:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Tumor_MCT4, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Tumor_MCT4, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Tumor_MCT4)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_mct4tu <- p_values_mct4tu %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_MCT4tu",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_mct4tu)

# Creating columns "significant" and "not significant"
p_values_long_mct4tu <- results_table_mct4tu %>%
  select(Variable_MCT4tu, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_MCT4tu, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_mct4tu, aes(x = Variable_MCT4tu, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of MCT4 tumor expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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

## Exploring associations between clinical-laboratory-radiomic features and glycolisis-related protein expression - stromal compartment.
```{r}
# GLUT1 stroma
p_values_glut1stro <- data_assoc %>%
  summarise(across(c(CL_Sex:Tumor_CD147, Stroma_HK2:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_GLUT1, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_GLUT1, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_GLUT1)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_glut1stro <- p_values_glut1stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_glut1stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_glut1stro)

# Creating columns "significant" and "not significant"
p_values_long_glut1stro <- results_table_glut1stro %>%
  select(Variable_glut1stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_glut1stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_glut1stro, aes(x = Variable_glut1stro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10)) +  # Hides the text on the X-axis
  labs(title = "Association of GLUT1 stromal expression with clinical, laboratory, and radiomic features",
    fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# HK2 stroma
p_values_hk2stro <- data_assoc %>%
  summarise(across(c(CL_Sex:Stroma_GLUT1, Stroma_LDHA:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_HK2, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_HK2, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_HK2)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_hk2stro <- p_values_hk2stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_hk2stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_hk2stro)

# Creating columns "significant" and "not significant"
p_values_long_hk2stro <- results_table_hk2stro %>%
  select(Variable_hk2stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_hk2stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_hk2stro, aes(x = Variable_hk2stro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10)) +  # Hides the text on the X-axis
  labs(title = "Association of HK2 stromal expression with clinical, laboratory, and radiomic features",
       fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# LDHA stroma
p_values_ldhastro <- data_assoc %>%
  summarise(across(c(CL_Sex:Stroma_HK2, Stroma_CA9:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_LDHA, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_LDHA, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_LDHA)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_ldhastro <- p_values_ldhastro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_ldhastro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_ldhastro)

# Creating columns "significant" and "not significant"
p_values_long_ldhastro <- results_table_ldhastro %>%
  select(Variable_ldhastro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_ldhastro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ldhastro, aes(x = Variable_ldhastro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_pubr() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        legend.position = "right") +  
  labs(title = "Association of LDHA stromal expression with clinical, laboratory, and radiomic features",
       fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,
    box.padding = 0.1,
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = 0,            
    segment.color = "transparent",  
    size = 3                   
  ) 

# CAIX stroma
p_values_ca9stro <- data_assoc %>%
  summarise(across(c(CL_Sex:Stroma_LDHA, Stroma_MCT4:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_CA9, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_CA9, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_CA9)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_ca9stro <- p_values_ca9stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_ca9stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_ca9stro)

# Creating columns "significant" and "not significant"
p_values_long_ca9stro <- results_table_ca9stro %>%
  select(Variable_ca9stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_ca9stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ca9stro, aes(x = Variable_ca9stro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_pubr() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "right",
        axis.ticks.x = element_blank()) +  # Hides the text on the X-axis
  labs(title = "Association of CA9 stromal expression with clinical, laboratory, and radiomic features",
       fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.8,
    box.padding = 0.1,
    angle = 45,
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = 0,            
    segment.color = "transparent",  
    size = 2                   
  ) 

# MCT4 stroma
p_values_mct4stro <- data_assoc %>%
  summarise(across(c(CL_Sex:Stroma_CA9, Stroma_CD147:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_MCT4, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_MCT4, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_MCT4)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_mct4stro <- p_values_mct4stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_mct4stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_mct4stro)

# Creating columns "significant" and "not significant"
p_values_long_mct4stro <- results_table_mct4stro %>%
  select(Variable_mct4stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_mct4stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_mct4stro, aes(x = Variable_mct4stro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "grey40") +  # Cutoff line
  scale_fill_manual(values = c("Significant" = "lightblue3", "Not Significant" = "grey70")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10)) +  # Hides the text on the X-axis
  labs(title = "Association of MCT4 stromal expression with clinical, laboratory, and radiomic features",
       fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(
    aes(label = label),
    vjust = 1,             # Aligns vertically for better positioning on the column
    box.padding = 0.1,     # Adjusts spacing between labels and bars
    point.padding = 0.1,
    max.overlaps = Inf
  )

# CD147 stroma
p_values_cd147stro <- data_assoc %>%
  summarise(across(c(CL_Sex:Stroma_MCT4, Stroma_GLUT3:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_CD147, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_CD147, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_CD147)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_cd147stro <- p_values_cd147stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_cd147stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_cd147stro)

# Creating columns "significant" and "not significant"
p_values_long_cd147stro <- results_table_cd147stro %>%
  select(Variable_cd147stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_cd147stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_cd147stro, aes(x = Variable_cd147stro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_pubr() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        legend.position = "right"
        ) +  # Hides the text on the X-axis
  labs(title = "Association of CD147 stromal expression with clinical, laboratory, and radiomic features",
       fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(             
    aes(label = label),        
    vjust = 1,
    box.padding = 0.1,
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = 0.1,            
    segment.color = "transparent",  
    size = 3                   
  ) 

# GLUT3 stroma
p_values_glut3stro <- data_assoc %>%
  summarise(across(c(CL_Sex:Stroma_CD147, Semantic_SUVmax:Texture_SumVar), .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_GLUT3, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_GLUT3, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_GLUT3)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_glut3stro <- p_values_glut3stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_glut3stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_glut3stro)

# Creating columns "significant" and "not significant"
p_values_long_glut3stro <- results_table_glut3stro %>%
  select(Variable_glut3stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_glut3stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_glut3stro, aes(x = Variable_glut3stro, y = -log10(p_value),
                                fill = significant)) +
  geom_col() +
  geom_hline(yintercept = log_significance_level, color = "#797A7B") +  
  scale_fill_manual(values = c("Significant" = "#FF8282", "Not Significant" = "#797A7B")) +
  xlab("Variable") +
  ylab("-log10(p-value)") +
  theme_pubr() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 10),
        axis.ticks.x = element_blank(),
        legend.position = "right"
        ) +  # Hides the text on the X-axis
  labs(title = "Association of GLUT3 stromal expression with clinical, laboratory, and radiomic features",
       fill = "") +
  # Adds labels using ggrepel to avoid overlap
  geom_text_repel(             
    aes(label = label),        
    vjust = 0.3,
    angle = 45,
    box.padding = 0.1,
    point.padding = 0.1,       
    nudge_y = 0.1,             
    nudge_x = 0.1,            
    segment.color = "transparent",  
    size = 2                   
  ) 
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
data_coexp <- data_assoc %>%
  dplyr::mutate(
    # Tumor co-expression score
    Tumor_coexp_score = rowSums(across(c(Tumor_GLUT1, Tumor_CD147, Tumor_CA9, Tumor_HK2, Tumor_LDHA, Tumor_MCT4)) == 1),
    Tumor_coexp_score = factor(Tumor_coexp_score, levels = 0:6),
    # Stroma co-expression score
    Stroma_coexp_score = rowSums(across(c(Stroma_GLUT1, Stroma_GLUT3, Stroma_CD147, Stroma_CA9, Stroma_HK2, Stroma_LDHA, Stroma_MCT4)) == 1),
    Stroma_coexp_score = factor(Stroma_coexp_score, levels = 0:7),
    # Tumor and stroma combined score
    Tumor_stroma_coexp_score = rowSums(across(c(
      Tumor_GLUT1, Tumor_CD147, Tumor_CA9, Tumor_HK2, Tumor_LDHA, Tumor_MCT4, 
      Stroma_GLUT1, Stroma_GLUT3, Stroma_CD147, Stroma_CA9, Stroma_HK2, Stroma_LDHA, Stroma_MCT4)) == 1),
    Tumor_stroma_coexp_score = factor(Tumor_stroma_coexp_score, levels = 0:13)
  ) %>%
  mutate(across(c(Tumor_coexp_score,
    Stroma_coexp_score,
    Tumor_stroma_coexp_score),
  ~ as.numeric(as.character(.))))

glimpse(data_coexp)
```

## Co-expression Analysis: Evaluating Significant Associations
This analysis focused on the co-expression of specific proteins based on association tests. Only proteins with statistically significant associations were included in the scoring system.

**1. Tumor Co-expression Associated Score (tu_coexp_assoc)**
Significantly associated tumor proteins: HK2, LDHA, and MCT4.
Scoring criteria:
If any of these proteins (HK2, LDHA, or MCT4) is not expressed (value = 0), the score is set to 0 (absence of associated co-expression).
If all of these proteins are expressed (value = 1), the score is set to 1 (presence of associated co-expression).

**2. Stroma Co-expression Associated Score (stro_coexp_assoc)**
Significantly associated stroma proteins: CAIX and CD147.
Scoring criteria:
If either CAIX or CD147 is not expressed (value = 0), the score is set to 0 (absence of associated co-expression).
If both CAIX and CD147 are expressed (value = 1), the score is set to 1 (presence of associated co-expression).

**3. Combined Tumor and Stroma Co-expression Associated Score (tu_stro_coexp_assoc)**
Significantly associated proteins: LDHA (tumor) and MCT4 (stroma).
Scoring criteria:
If either MCT4 (stroma) or LDHA (tumor) is not expressed (value = 0), the score is set to 0 (absence of associated co-expression).
If both MCT4 (stroma) and LDHA (tumor) are expressed (value = 1), the score is set to 1 (presence of associated co-expression).
```{r}
data_coexp_assoc <- data_coexp %>%
  dplyr::mutate(
    Stroma_ca9_cd147_coexp_assoc = case_when(
      Stroma_CA9 == 0 | Stroma_CD147 == 0 ~ 0,
      Stroma_CA9 == 1 & Stroma_CD147 == 1 ~ 1
    ),
    Tumor_stroma_ldha_coexp_assoc = case_when(
      Tumor_LDHA == 0 | Stroma_LDHA == 0 ~ 0,
      Tumor_LDHA == 1 & Stroma_LDHA == 1 ~ 1
    )
  ) %>%
  dplyr::mutate(across(
    c(
      Stroma_ca9_cd147_coexp_assoc,
      Tumor_stroma_ldha_coexp_assoc,
    ),
    as.factor
  ))

glimpse(data_coexp_assoc)

```

```{r}
data_coexp_all <- data_coexp_assoc %>%
  dplyr::mutate(
    Tumor_stroma_coexp_hk2 = case_when(
      Tumor_HK2 == 0 | Stroma_HK2 == 0 ~ 0,
      Tumor_HK2 == 1 & Stroma_HK2 == 1 ~ 1
    ),
    Tumor_stroma_coexp_mct4 = case_when(
      Tumor_MCT4 == 0 | Stroma_MCT4 == 0 ~ 0,
      Tumor_MCT4 == 1 & Stroma_MCT4 == 1 ~ 1
    )
  ) %>%
  dplyr::mutate(across(
    c(
      Tumor_stroma_coexp_hk2,
      Tumor_stroma_coexp_mct4,
    ),
    as.factor
  ))

glimpse(data_coexp_all)
```

## Relocating columns to facilitate feature selection
```{r}
data_fs <- data_coexp_all %>%
  relocate(Tumor_coexp_score:Tumor_stroma_coexp_mct4, .after = Stroma_GLUT3)
glimpse(data_fs)
```

### Correlation analysis & plots - semantic and texture-based features x co-expression score
```{r}
# Correlation analysis - Spearman.
data_corr_coexp_score <- data_fs %>%
  dplyr::select(Tumor_coexp_score:Tumor_stroma_coexp_score, Semantic_SUVmax:Texture_SumVar)

# Correlation matrix and network plot
cor_mat_coexp_score <- data_corr_coexp_score %>%
  correlate(method = "spearman") 

# Focusing on co-expression score associations with radiomic features
cor_focus <- cor_mat_coexp_score %>%
  focus(Tumor_coexp_score, Stroma_coexp_score, Tumor_stroma_coexp_score)

# Transforming the dataset into long format
cor_long <- cor_focus %>%
  pivot_longer(cols = -term, names_to = "variable", values_to = "correlation") %>%
  mutate(term = reorder(term, correlation))  # Ordenando pelas correlações

# Bar plot
ggplot(cor_long, aes(x = term, y = correlation, fill = variable)) +
  geom_col(position = "dodge") +  
  coord_flip() + 
  labs(
    x = "Radiomic features",
    y = "r (Spearman)",
    fill = "Variable",
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_pubr() +
  theme(legend.position = "right")
```

## Association of protein co-expression with clinical, laboratory, and radiomic features
```{r}
# CA9 CD147 stromal co-expression
p_values_ca9_cd147_stro <- data_fs %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Semantic_SUVmax:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Stroma_ca9_cd147_coexp_assoc, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Stroma_ca9_cd147_coexp_assoc, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Stroma_ca9_cd147_coexp_assoc)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_ca9_cd147_stro <- p_values_ca9_cd147_stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_ca9_cd147_stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_ca9_cd147_stro)

# Creating columns "significant" and "not significant"
p_values_long_ca9_cd147_stro <- results_table_ca9_cd147_stro %>%
  select(Variable_ca9_cd147_stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_ca9_cd147_stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ca9_cd147_stro, aes(x = Variable_ca9_cd147_stro, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of CA9 and CD147 stromal co-expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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

# LDHA tumor and stroma co-expression
p_values_ldha_tu_stro <- data_fs %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Semantic_SUVmax:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Tumor_stroma_ldha_coexp_assoc, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Tumor_stroma_ldha_coexp_assoc, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Tumor_stroma_ldha_coexp_assoc)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_ldha_tu_stro <- p_values_ldha_tu_stro %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_ldha_tu_stro",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_ldha_tu_stro)

# Creating columns "significant" and "not significant"
p_values_long_ldha_tu_stro <- results_table_ldha_tu_stro %>%
  select(Variable_ldha_tu_stro, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_ldha_tu_stro, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_ldha_tu_stro, aes(x = Variable_ldha_tu_stro, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of LDHA tumor and stroma co-expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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

# HK2 tumor and stroma co-expression
p_values_Tumor_stroma_coexp_hk2 <- data_fs %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Semantic_SUVmax:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Tumor_stroma_coexp_hk2, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Tumor_stroma_coexp_hk2, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Tumor_stroma_coexp_hk2)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_Tumor_stroma_coexp_hk2 <- p_values_Tumor_stroma_coexp_hk2 %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_Tumor_stroma_coexp_hk2",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_Tumor_stroma_coexp_hk2)

# Creating columns "significant" and "not significant"
p_values_long_Tumor_stroma_coexp_hk2 <- results_table_Tumor_stroma_coexp_hk2 %>%
  select(Variable_Tumor_stroma_coexp_hk2, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_Tumor_stroma_coexp_hk2, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_Tumor_stroma_coexp_hk2, aes(x = Variable_Tumor_stroma_coexp_hk2, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of HK2 tumor and stroma co-expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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

# mct4 tumor and stroma co-expression
p_values_Tumor_stroma_coexp_mct4 <- data_fs %>%
  summarise(across(c(CL_Sex:CL_Bulky_disease, Semantic_SUVmax:Texture_SumVar), 
                   .fns = function(.x) {
    if (is.numeric(.x)) {
      # Normality test (Shapiro-Wilk)
      shapiro_p <- shapiro.test(.x)$p.value  
      if (shapiro_p < significance_level) {  
        test <- "Mann-Whitney"
        p_value <- wilcox.test(.x ~ Tumor_stroma_coexp_mct4, exact = FALSE, alternative = "two.sided")$p.value
      } else {  
        test <- "t-test"
        p_value <- t.test(.x ~ Tumor_stroma_coexp_mct4, alternative = "two.sided")$p.value
      }
    } else {
      # Categorical variables
      contingency_table <- table(.x, Tumor_stroma_coexp_mct4)
      expected_values <- chisq.test(contingency_table, simulate.p.value = FALSE)$expected
      if (any(expected_values < 5) | (sum(expected_values < 5) / length(expected_values) > 0.2)) {
        test <- "Fisher"
        p_value <- fisher.test(contingency_table, alternative = "two.sided")$p.value
      } else {
        test <- "Chi-squared"
        p_value <- chisq.test(contingency_table, simulate.p.value = FALSE)$p.value
      }
    }
    # Return the results as a dataframe
    return(data.frame(test = test, p_value = p_value))
  }, .names = "{.col}")) 

# Restructure the results into a table
results_table_Tumor_stroma_coexp_mct4 <- p_values_Tumor_stroma_coexp_mct4 %>%
  pivot_longer(cols = everything(),
               names_to = "Variable_Tumor_stroma_coexp_mct4",
               values_to = "Result") %>%
  unnest(Result)

# Display the results
print(results_table_Tumor_stroma_coexp_mct4)

# Creating columns "significant" and "not significant"
p_values_long_Tumor_stroma_coexp_mct4 <- results_table_Tumor_stroma_coexp_mct4 %>%
  select(Variable_Tumor_stroma_coexp_mct4, p_value) %>%
  mutate(significant = ifelse(p_value < significance_level, "Significant", "Not Significant"),
         label = ifelse(p_value < significance_level, Variable_Tumor_stroma_coexp_mct4, ""))

# Creating the plot with ggrepel to adjust labels
ggplot(p_values_long_Tumor_stroma_coexp_mct4, aes(x = Variable_Tumor_stroma_coexp_mct4, y = -log10(p_value),
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
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 10)
  ) +  
  labs(title = "Association of MCT4 tumor and stroma co-expression with clinical, laboratory, and radiomic features",
    fill = "") +  
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

## Recipes - preprocessing steps
```{r}
# Recipes 1 - dummy and normalization preprocessing steps
rec_spec <- recipe(Outcome_5PS_CAT ~ ., data = data_fs) %>%
  step_dummy(all_nominal_predictors()) %>% # encoding categorical features
  step_normalize(Semantic_SUVmax:Semantic_MTV) %>% # normalizing semantic features
  step_normalize(Texture_Entropy:Texture_SumVar) %>% # normalizing texture features
  step_normalize(Tumor_coexp_score:Tumor_stroma_coexp_score)

# Prepping the recipe.
rec_prep <- prep(rec_spec)

# Applying to dataset.
data_fs_baked <- bake(rec_prep, new_data = NULL)

```

```{r}
# Transforming the dataset into long format
data_fs_long <- data_fs_baked %>%
  pivot_longer(
    cols = -Outcome_5PS_CAT,            
    names_to = "variable",   
    values_to = "value"      
  )

```

## Step 1 - Univariate logistic regression
```{r}
# performing univariate binary logistic regression in all variables.
resultList <- data_fs_long %>%
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
  filter(!grepl("(Intercept)", term)) %>%
  mutate(OR = exp(estimate),
        conf.low.or = exp(conf.low),
        conf.high.or = exp(conf.high)
        ) %>%
  select(variable, OR, conf.low.or, conf.high.or, p.value)

print(step_1_df_cleaned$variable)
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

glimpse(step_2_df_filt)
```


## Step 3 - Correlation or VIF. In this case, Spearman correlation analysis was performed.
```{r}
# Filtering only numeric features (radiomic features) for correlation analysis.
step_2_df_filt_rad <- step_2_df_filt %>%
  filter(!grepl("X1", variable$variable))
glimpse(step_2_df_filt_rad$variable)
```

### Correlation analysis & plots - semantic and texture-based features
```{r}
# Correlation analysis - Spearman.
data_corr_semantic <- data_fs %>%
  dplyr::select(Semantic_SUVmean, Semantic_MTV, Semantic_TLG)

data_corr_texture <- data_fs %>%
  dplyr::select(Texture_ASM, Texture_ClusP, Texture_Clussh, Texture_Contrast, Texture_DifEnt,
                Texture_DifVar, Texture_Entropy, Texture_Homogeneity, Texture_IDM, Texture_MC1,
                Texture_MaxP, Texture_SumEnt
                )

# Correlation matrix and network plot
cor_mat_semantic <- data_corr_semantic %>%
  correlate(method = "spearman") 
fashion(cor_mat_semantic)
network_plot(cor_mat_semantic, min_cor = 0.2) +
  theme(legend.position = "none")

any_over_90 <- function(x) any(x >= .8)
cor_mat_semantic %>%
  focus_if(any_over_90, mirror = TRUE) %>%
  shave() %>% 
  rplot(shape = 20, print_cor = T) +
  theme(legend.position = "none")

# Correlation matrix and network plot
cor_mat_texture <- data_corr_texture %>%
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
data_corr_semantic_long <- data_corr_semantic %>%
  correlate(method = "spearman") %>%
  stretch() %>%
  mutate(Correlation = abs(r)) %>% # obtaining the absolute value of correlation.
  select(-r)

data_corr_semantic_long_p <- data_corr_semantic_long %>% # **pesquisar uma melhor forma de fazer isso**
  dplyr::mutate(p.value_x = case_when(
    x == "Semantic_MTV" ~ 0.004904301,
    x == "Semantic_TLG" ~ 0.026852108,
    x == "Semantic_SUVmean" ~ 0.375860802
  ),
  p.value_y = case_when(
    y == "Semantic_MTV" ~ 0.004904301,
    y == "Semantic_TLG" ~ 0.026852108,
    y == "Semantic_SUVmean" ~ 0.375860802
  )) %>% 
  dplyr::mutate(min_p.value = pmin(p.value_x, p.value_y), # Smallest p.value between correlated features.
    variable_to_keep = ifelse(p.value_x <= p.value_y, x, y) # keeping the variable with the smallest p.value.
  )
```
### Semantic -  Comparing p-values of highly correlated features (|r| ≥ 0.8) and retaining the one with the smallest p-value.
```{r}
filter_high_corr <- function(data_corr_semantic_long_p) {
  # Lista inicial de variáveis a manter
  variables_to_keep <- unique(c(data_corr_semantic_long_p$x, data_corr_semantic_long_p$y))
  
  repeat {
    # Filtrar pares altamente correlacionados
    high_corr_pairs <- data_corr_semantic_long_p %>%
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
    data_corr_semantic_long_p <- data_corr_semantic_long_p %>%
      filter(x != var_to_remove & y != var_to_remove)
  }
  
  return(variables_to_keep)
}

# Executar o filtro
final_variables_semantic <- filter_high_corr(data_corr_semantic_long_p) 
print(final_variables_semantic)

```

### Texture -  Creating columns with p.values.
```{r}
# Transforming the correlation df into a long format df.
data_corr_texture_long <- data_corr_texture %>%
  correlate(method = "spearman") %>%
  stretch() %>%
  mutate(Correlation = abs(r)) %>% # obtaining the absolute value of correlation.
  select(-r)

# Creating the columns p.value_x and p.value_y to compare p-values of highly correlated features (|r| ≥ 0.8) and retain the one with the smallest p-value.
data_corr_texture_long_p <- data_corr_texture_long %>% # **pesquisar uma melhor forma de fazer isso**
  dplyr::mutate(p.value_x = case_when(
    x == "Texture_ASM" ~ 0.203483648,
    x == "Texture_ClusP" ~ 0.016303613,
    x == "Texture_Clussh" ~ 0.462990824,
    x == "Texture_Contrast" ~ 0.140060859,
    x == "Texture_DifEnt" ~ 0.103565488,
    x == "Texture_DifVar" ~ 0.110214606,
    x == "Texture_Entropy" ~ 0.179820952,
    x == "Texture_Homogeneity" ~ 0.194154265,
    x == "Texture_IDM" ~ 0.026275430,
    x == "Texture_MC1" ~ 0.182457258,
    x == "Texture_MaxP" ~ 0.194537763,
    x == "Texture_SumEnt" ~ 0.025060147
  ),
  p.value_y = case_when(
    y == "Texture_ASM" ~ 0.203483648,
    y == "Texture_ClusP" ~ 0.016303613,
    y == "Texture_Clussh" ~ 0.462990824,
    y == "Texture_Contrast" ~ 0.140060859,
    y == "Texture_DifEnt" ~ 0.103565488,
    y == "Texture_DifVar" ~ 0.110214606,
    y == "Texture_Entropy" ~ 0.179820952,
    y == "Texture_Homogeneity" ~ 0.194154265,
    y == "Texture_IDM" ~ 0.026275430,
    y == "Texture_MC1" ~ 0.182457258,
    y == "Texture_MaxP" ~ 0.194537763,
    y == "Texture_SumEnt" ~ 0.025060147
  )) %>% 
  dplyr::mutate(min_p.value = pmin(p.value_x, p.value_y), # Smallest p.value between correlated features.
    variable_to_keep = ifelse(p.value_x <= p.value_y, x, y) # keeping the variable with the smallest p.value.
  )

```

```{r}
filter_high_corr <- function(data_corr_texture_long_p) {
  # Lista inicial de variáveis a manter
  variables_to_keep <- unique(c(data_corr_texture_long_p$x, data_corr_texture_long_p$y))
  
  repeat {
    # Filtrar pares altamente correlacionados
    high_corr_pairs <- data_corr_texture_long_p %>%
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
    data_corr_texture_long_p <- data_corr_texture_long_p %>%
      filter(x != var_to_remove & y != var_to_remove)
  }
  
  return(variables_to_keep)
}

# Executar o filtro
final_variables <- filter_high_corr(data_corr_texture_long_p) # final_variables contem tanto as variáveis que devem ser mantidas por não serem altamente correlacionadas (Texture_IDM, Texture_ClusP e Texture_Clussh), quanto as variáveis que estavam correlacionadas (abs(r) >= 0.8) com outra variável, porém, apresentaram um menor valor de p (Texture_DifVar e Texture_SumEnt --> apesar dessas duas variáveis estarem correlacionadas com outras variáveis que foram removidas, elas não estão correlacionadas entre si).
print(final_variables)
```

```{r}
# selecting features from step 2 (p<0.5) and removing highly correlated features identified at step 3.
data_p <- data_fs %>% 
  select(
  # texture-based features
  Texture_ClusP, Texture_Clussh, 
  Texture_DifEnt, Texture_IDM,
  Texture_SumEnt,
  # semantic features
  Semantic_MTV, Semantic_SUVmean,
  # stromal expression
  Stroma_CA9,
  Stroma_GLUT1,
  Stroma_GLUT3,
  Stroma_HK2,
  Stroma_MCT4,
  Stroma_LDHA,
  Stroma_ca9_cd147_coexp_assoc,
  # tumoral expression
  Tumor_CA9,
  Tumor_HK2,
  Tumor_LDHA,
  Tumor_stroma_ldha_coexp_assoc,
  Tumor_coexp_score,
  Tumor_stroma_coexp_hk2,
  # clinical-laboratory features
  CL_Staging,
  CL_High_LDH,
  CL_Extranodal_site_involvement,
  CL_Bulky_disease,
  CL_Sex,
  CL_Age,
  Outcome_5PS_CAT)

# Recipes 2: preprocessing steps for features selected at the end of step 3.
rec_corr <- recipe(Outcome_5PS_CAT ~ ., data = data_p) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(Semantic_MTV:Semantic_SUVmean) %>% # normalizing semantic features
  step_normalize(Texture_ClusP:Texture_SumEnt) %>% # normalizing texture features
  step_corr(Semantic_MTV:Semantic_SUVmean, threshold = 0.8, method = "spearman") %>%
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
View(step_4_df_cor_filt$variable)

```

## Step 5 - Multivariate logistic regression.
```{r}
# selecting features from the unprocessed initial df.
data_final <- data_p %>%
  select(
    # texture-based features
    Texture_ClusP, 
    Texture_DifEnt,
    Texture_IDM,
    Texture_SumEnt,
    # semantic features
    Semantic_MTV,
    # stromal expression
    Stroma_CA9,
    Stroma_HK2,
    # tumoral expression
    Tumor_HK2,
    Tumor_coexp_score,
    # clinical-laboratory features
    CL_Staging,
    CL_High_LDH,
    CL_Bulky_disease,
    CL_Sex,
    Outcome_5PS_CAT) # selected features at the end of step 4 (not highly correlated & p.value <= 0.2)

glimpse(data_final)
```

```{r}
# Recipes 3: Preprocessing steps for features selected at the end of step 4
rec_multi <- recipe(Outcome_5PS_CAT ~ ., data = data_final) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_normalize(Semantic_MTV) %>% # normalizing semantic features
  step_normalize(Texture_ClusP:Texture_SumEnt) %>% # normalizing texture features
  step_normalize(Tumor_coexp_score) %>%
  step_corr(Texture_ClusP:Texture_SumEnt, threshold = 0.8, method = "spearman")

# Prepping the recipe.
rec_multi_prep <- prep(rec_multi)

# Applying to the dataset.
rec_multi_baked <- bake(rec_multi_prep, new_data = NULL)

glimpse(rec_multi_baked)
```

```{r}
# Multivariate logistic regression
set.seed(12)
glm_fit_multi <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Outcome_5PS_CAT ~., data = rec_multi_baked)
glm_fit_multi

# Adding a column indicating if p.value < 0.05
glm_tidy_multi <- tidy(glm_fit_multi, conf.int = TRUE) %>%
  mutate(significant = p.value < 0.05) %>%
  mutate(OR = exp(estimate),
         conf.low.or = exp(conf.low),
         conf.high.or = exp(conf.high)
         )

print(glm_tidy_multi)

glm_save_multi <- glm_tidy_multi %>%
  select(term, OR, conf.low.or, conf.high.or, p.value)
print(glm_save_multi)

```


## Checking multicolinearity with VIF
```{r}
library(car)
vif_1 <- as.data.frame(vif(glm_fit_multi$fit))
print(vif_1)
```

## Removing the variable with VIF > 4 and the highest p-value among those with VIF > 4, then refitting the model.
```{r}
rec_multi_baked_vif_1 <- rec_multi_baked %>%
  select(-Texture_ClusP)

# Multivariate logistic regression
set.seed(12)
glm_fit_multi_vif_1 <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Outcome_5PS_CAT ~., data = rec_multi_baked_vif_1)
glm_fit_multi_vif_1

# Adding a column indicating if p.value < 0.05
glm_tidy_multi_vif_1 <- tidy(glm_fit_multi_vif_1, conf.int = TRUE) %>%
  mutate(significant = p.value < 0.05) %>%
  mutate(OR = exp(estimate),
         conf.low.or = exp(conf.low),
         conf.high.or = exp(conf.high)
         )

print(glm_tidy_multi_vif_1)

glm_save_multi_vif_1 <- glm_tidy_multi_vif_1 %>%
  select(term, OR, conf.low.or, conf.high.or, p.value)
print(glm_save_multi_vif_1)

```
## Checking multicolinearity with VIF
```{r}
vif_2 <- as.data.frame(vif(glm_fit_multi_vif_1$fit))
print(vif_2)
```

A variable still has a VIF > 4. At this stage, the variable was removed, and the model was refitted.
```{r}
rec_multi_baked_vif_2 <- rec_multi_baked_vif_1 %>%
  select(-Texture_DifEnt)

# Multivariate logistic regression
set.seed(12)
glm_fit_multi_vif_2 <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Outcome_5PS_CAT ~., data = rec_multi_baked_vif_2)
glm_fit_multi_vif_2

# Adding a column indicating if p.value < 0.05
glm_tidy_multi_vif_2 <- tidy(glm_fit_multi_vif_2, conf.int = TRUE) %>%
  mutate(significant = p.value < 0.05) %>%
  mutate(OR = exp(estimate),
         conf.low.or = exp(conf.low),
         conf.high.or = exp(conf.high)
         )

print(glm_tidy_multi_vif_2)

glm_save_multi_vif_2 <- glm_tidy_multi_vif_2 %>%
  select(term, OR, conf.low.or, conf.high.or, p.value)
print(glm_save_multi_vif_2)

```

## Checking for multicolinearity
```{r}
vif_3 <- as.data.frame(vif(glm_fit_multi_vif_2$fit))
print(vif_3)
```
All variables have a VIF < 4.

## Visualizing the results
```{r}
# Visualizing the results with gghighlight
ggplot(glm_tidy_multi_vif_2, aes(x = term, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  gghighlight(significant, label_key = term, use_direct_label = FALSE) +
  labs(
    x = "Variable",
    y = "log-odds",
    title = "Features with p < 0.05 in multivariate logistic regression"
  ) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

## Selecting only significant features
```{r}
rec_multi_baked_significant <- rec_multi_baked_vif_2 %>%
  select(CL_Staging_X1, Stroma_HK2_X1, Tumor_HK2_X1, Texture_IDM, Outcome_5PS_CAT)

rec_multi_baked_significant <- as.data.frame(rec_multi_baked_significant)

# Multivariate logistic regression
set.seed(12)
glm_fit_multi_significant <- logistic_reg(mode = "classification") %>%
  set_engine(engine = "glm") %>% 
  fit(Outcome_5PS_CAT ~., data = rec_multi_baked_significant)
glm_fit_multi_significant

# Adding a column indicating if p.value < 0.05
glm_tidy_multi_significant <- tidy(glm_fit_multi_significant, conf.int = TRUE) %>%
  mutate(significant = p.value < 0.05) %>%
  mutate(OR = exp(estimate),
         conf.low.or = exp(conf.low),
         conf.high.or = exp(conf.high)
         )

print(glm_tidy_multi_significant)

glm_save_multi_significant <- glm_tidy_multi_significant %>%
  select(term, OR, conf.low.or, conf.high.or, p.value)
print(glm_save_multi_significant)

```

## Checking for multicolinearity
```{r}
vif(glm_fit_multi_significant$fit)
```

## Residuals
```{r}
glm_fit_multi_significant_base <- glm(
  Outcome_5PS_CAT ~ .,
  data = rec_multi_baked_significant,
  family = binomial
)

summary(glm_fit_multi_significant_base$residuals)
#residuals_lr_vif <- 
residualPlots(glm_fit_multi_significant_base, tests = F)
#png("residuals_lr_vif.png", width = 16, height = 12, units = "cm", res = 600)
#residualPlots(glm_fit_multi_significant_base, tests = FALSE)
#dev.off()

```

## McFadden`s pseudo R^2
Calculating the R^2 and the p value for the relationship between the predictors and categorized 5PS
```{r}
ll.null <- glm_fit_multi_significant_base$null.deviance/-2
ll.proposed <- glm_fit_multi_significant_base$deviance/-2
## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
(ll.null - ll.proposed) / ll.null
 
## chi-square value = 2*(LL(Proposed) - LL(Null))
## p-value = 1 - pchisq(chi-square value, df = 2-1)
1 - pchisq(2*(ll.proposed - ll.null), df=1)
1 - pchisq((glm_fit_multi_significant_base$null.deviance - glm_fit_multi_significant_base$deviance), df=1)
```


## Prediction and model metrics
```{r}
set.seed(12)
pred_data_significant <- rec_multi_baked_significant %>%
  mutate(predicted_prob_significant = predict(glm_fit_multi_significant$fit,
                                  newdata = rec_multi_baked_significant, type = "response"))
```

## Choosing the threshold that maximizes sensitivity and specificity
```{r}
library(probably)
thresholds <- seq(0.1, 0.9, by = 0.001)

threshold_data_significant <- pred_data_significant %>%
  threshold_perf(Outcome_5PS_CAT, predicted_prob_significant, event_level = "second", thresholds)

threshold_save <- threshold_data_significant %>%
  arrange(.threshold)

threshold_data_significant <- threshold_data_significant %>%
  filter(.metric != "distance") %>%
  mutate(group = case_when(
    .metric == "sens" | .metric == "spec" ~ "1",
    TRUE ~ "2"
  ))

max_j_index_threshold_significant <- threshold_data_significant %>%
  filter(.metric == "j_index") %>%
  filter(.estimate == max(.estimate)) %>%
  pull(.threshold)

max_j_index_threshold_significant

threshold_data_significant %>%
  filter(.threshold == max_j_index_threshold_significant)

```

Best threshold = 0.1
```{r}
# If predicted probability >= 0.1, predicted class = 1. Else, predicted class = 0.
pred_data_significant <- pred_data_significant %>%
  mutate(predicted_class_significant = ifelse(predicted_prob_significant >= 0.1, 1, 0))
pred_data_significant$predicted_class_significant <- as.factor(pred_data_significant$predicted_class_significant)
```

## Confusion matrix and ROC.
```{r}
pred_data_significant$Outcome_5PS_CAT <- fct_recode(pred_data_significant$Outcome_5PS_CAT, "1-3" = "0", "4-5" = "1")
pred_data_significant$predicted_class_significant <- fct_recode(pred_data_significant$predicted_class_significant, "1-3" = "0", "4-5" = "1")
pred_data_significant %>%
  conf_mat(Outcome_5PS_CAT, predicted_class_significant) %>%
  autoplot("heatmap") +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 12, face = "bold"))

pred_data_significant %>%
  roc_curve(Outcome_5PS_CAT, predicted_prob_significant, event_level = "second") %>%
  autoplot() +
  geom_abline(lty = 1, color = "black", size = 0.5) +
  theme_pubr()

 pred_data_significant %>%
  roc_auc(Outcome_5PS_CAT, predicted_prob_significant, event_level = "second")
```
