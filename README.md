# Predicting Metabolically Active Disease in Aggressive B-Cell Lymphomas
### 📖Overview
This repository contains the code and resources used in a pilot study to develop predictive models for metabolically active disease in aggressive B-cell lymphomas (ABCL). The study integrates clinical-laboratory data and radiomic features extracted from baseline 18F-FDG PET scans. The analysis pipeline involves exploratory data analysis (EDA), feature selection, and machine learning workflows to evaluate the predictive potential of these features.

**Radiomic Features**: Five semantic and sixteen texture-based features were extracted using the FIJI software.

**Objective**: To explore various combinations of features and model architectures to predict disease activity.

# Repository Structure

### 🖥️Exploratory Data Analysis (EDA):
- Feature Selection: Identification of relevant clinical and radiomic features.
- Glycolysis-Related Proteins: Correlation between textural features and glycolysis-related protein expression.

### 🧠Machine Learning (ML) Models:

- Combined Models: Integration of clinical and radiomic data.
- Radiomic-Only Models: Predictive models based exclusively on radiomic features.
- Comparison Studies: Performance evaluation of combined vs. radiomic-only models.

### ✅Trial-and-Error Experiments:

#### Oversampling Techniques:
- SMOTENC for combined models and SMOTE for radiomic models.

#### Feature-Specific Models:
- MTV (Metabolic Tumor Volume).
- IDM (Inverse Difference Moment).
- Ann Arbor Staging.

#### Combined Features:
Analyses using radiomic features selected by different feature selection methods.

# Prerequisites and Installation
- R (Version >= 4.00)
- Packages:
   - tidyverse
   - tidymodels
   - ggpubr
   - themis
   - stacks
   - finetune
   - vip
   - tidyposterior
   - modeldata
   - ggrepel
   - corrplot
   - corrr
   - gghighlight
   - ggridges
   - cutpointr
   - readxl

 #### Installation
 - Install the packages in R:

 install.packages(c(
  "tidyverse",
  "tidymodels",
  "ggpubr",
  "themis",
  "stacks",
  "finetune",
  "vip",
  "tidyposterior",
  "modeldata",
  "ggrepel",
  "corrplot",
  "corrr",
  "gghighlight",
  "ggridges",
  "readxl",
  "cutpointr"
))

- Install the dataset (data_stro.xlsx)

