# 01_load_and_subset.R
# Load OsteoCAR human primary tumor data and subset to tumor cells

library(qs)
library(Seurat)
library(tidyverse)

# Load data
patient_prim <- qread("data/patient_prim.qs")

# Verify
cat("Total cells:", ncol(patient_prim), "\n")
cat("Annotation levels:", colnames(patient_prim@meta.data), "\n")
table(patient_prim$Ann_Level0)

# Subset to tumor cells only
tumor_prim <- subset(patient_prim, Ann_Level0 == "Tumor")
cat("Tumor cells:", ncol(tumor_prim), "\n")
table(tumor_prim$Ann_Level3)

# Free memory
rm(patient_prim)
gc()

# Save subsetted object
qsave(tumor_prim, "data/tumor_prim.qs")