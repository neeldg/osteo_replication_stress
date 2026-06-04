# 04_scoring.R
# Score each tumor cell for TRSS, proliferation modules, and repstress

library(qs)
library(Seurat)
library(tidyverse)

# Load data and gene signatures
tumor_prim <- qread("data/tumor_prim.qs")
gene_signatures <- readRDS("data/gene_signatures.rds")

# Extract signatures
trss_weights <- gene_signatures$trss_weights
e2f_genes <- gene_signatures$e2f_genes
g2m_genes <- gene_signatures$g2m_genes
osteo_prolif_genes <- gene_signatures$osteo_prolif_genes
repstress_genes <- gene_signatures$repstress_genes

# ── TRSS weighted score ──────────────────────────────────────────────────────
expr_matrix <- GetAssayData(tumor_prim, layer = "data")
trss_expr <- expr_matrix[names(trss_weights), ]
trss_score <- as.vector(t(trss_expr) %*% trss_weights)
tumor_prim$TRSS_score <- trss_score
cat("TRSS score range:", range(trss_score), "\n")

# ── Proliferation module scores (Cam's gene sets) ───────────────────────────
tumor_prim <- AddModuleScore(
  tumor_prim,
  features = list(e2f_genes),
  name = "E2F_score"
)

tumor_prim <- AddModuleScore(
  tumor_prim,
  features = list(g2m_genes),
  name = "G2M_score"
)

tumor_prim <- AddModuleScore(
  tumor_prim,
  features = list(osteo_prolif_genes),
  name = "Osteo_prolif"
)

# ── Repstress signature score ────────────────────────────────────────────────
tumor_prim <- AddModuleScore(
  tumor_prim,
  features = list(repstress_genes),
  name = "Repstress"
)

# Verify
head(tumor_prim@meta.data[, c("TRSS_score", "E2F_score1",
                              "G2M_score1", "Osteo_prolif1",
                              "Repstress1")])

# Save
qsave(tumor_prim, "data/tumor_prim.qs")