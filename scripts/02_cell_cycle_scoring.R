# 02_cell_cycle_scoring.R
# Assign cell cycle phase to each tumor cell

library(qs)
library(Seurat)
library(tidyverse)

# Load subsetted tumor cells
tumor_prim <- qread("data/tumor_prim.qs")

# Cell cycle scoring using Seurat built-in gene lists
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

tumor_prim <- CellCycleScoring(
  tumor_prim,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# Verify
table(tumor_prim$Phase)
table(tumor_prim$Ann_Level3, tumor_prim$Phase)

# Save
qsave(tumor_prim, "data/tumor_prim.qs")