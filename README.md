# Predicting Metabolically Active Disease in Aggressive B-Cell Lymphomas
### 📖Overview
This repository contains the code and resources used in a study to develop predictive models for metabolically active disease in aggressive B-cell lymphomas (ABCL). The study integrates clinical-laboratory data and radiomic features extracted from baseline 18F-FDG PET scans. The analysis pipeline involves exploratory data analysis (EDA), feature selection, and machine learning workflows to evaluate the predictive potential of these features. Prediction of 2-year overall survival was also assessed using similar modeling strategies.

**Radiomic Features**: 112 IBSI-compliant radiomic features were extracted using the LIFEx software.

**Objective**: To explore various combinations of features and model architectures to predict disease activity.

# Repository Structure

- [HTML code folder in Google Drive](https://drive.google.com/drive/u/0/folders/1p9nfra71X9MXYBI6lSm_Lpv_ZU1jSel9)

### 🖥️Exploratory Data Analysis (EDA):
- Feature Selection: Identification of relevant radiomic features.
- Glycolysis-Related Proteins: Correlation between textural features and glycolysis-related protein expression.

### 🧠Machine Learning (ML) Models:

- Combined Models: Integration of clinical and radiomic data.
- Radiomic-Only Models: Predictive models based exclusively on radiomic features.
- Clinical-Only Models: Models including key clinical predictors.
- Comparison Studies: Performance evaluation of combined vs. radiomic-only models.

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
   - discrim
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
  "cutpointr",
  "discrim"
))
