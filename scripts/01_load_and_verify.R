library(qs)
library(Seurat)
library(tidyverse)

# Load data
patient_prim <- qread("data/patient_prim.qs")
patient_mets <- qread("data/patient_mets.qs")

# Verify
ncol(patient_prim)
ncol(patient_mets)
colnames(patient_prim@meta.data)
table(patient_prim$Ann_Level0)
table(patient_prim$Ann_Level3)