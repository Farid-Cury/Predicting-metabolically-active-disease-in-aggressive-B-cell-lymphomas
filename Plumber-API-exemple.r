---
title: "PLumber API exemple"
author: "Farid"
date: "`r Sys.Date()`"
output: html_document
---
## Loading packages and dataset
```{r}
install.packages(c("tidymodels", "tidyverse", "readxl", "ranger", "plumber", "vetiver"))
library(tidyverse)
library(tidymodels)
library(readxl)
options(tidymodels.dark = TRUE)

data_clusters <- read_excel("data_clusters_suv4.xlsx")

data_clusters$Outcome_5PS_CAT <- fct_recode(data_clusters$Outcome_5PS_CAT, "RC" = "0", "DA"= "1")

data_ipi_rad <- data_clusters %>%
  select(Estadiamento, 
        `GLCM_NormalisedInverseDifferenceMoment(IBSI:1QCO)`,
         Outcome_5PS_CAT)


```

```{r}
set.seed(42)
data_split <- initial_split(data_ipi_rad, prop = 0.8, strata = Outcome_5PS_CAT)
data_tr <- training(data_split)
data_te <- testing(data_split)
```

```{r}
rec_rf_ipi_rad <- recipe(Outcome_5PS_CAT ~ ., data_tr)

rec_rf_ipi_rad %>% prep() %>% juice() %>% glimpse()
```

```{r}
rf_model_ipi_rad <- rand_forest(
  trees = 1000
              ) %>%
                  set_engine("ranger") %>%
                  set_mode("classification")
```

```{r}
rf_wf_ipi_rad <- workflow() %>% 
  add_model(rf_model_ipi_rad) %>% 
  add_recipe(rec_rf_ipi_rad)

set.seed(12)
rf_last_fit_ipi_rad <- last_fit(rf_wf_ipi_rad, data_split) 
```

```{r}
# Building an API. 
library(vetiver)

model_api <- extract_workflow(rf_last_fit_ipi_rad) %>%
  vetiver_model("RF_IPI_Rad")

library(plumber)

# class probabilities (remove the "#" in the code to run the API locally).
 pr() %>%
  vetiver_api(model_api, type = "prob") %>%
   pr_run()
```
