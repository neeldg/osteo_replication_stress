# 03_gene_signatures.R
# Define all gene signatures used in the analysis:
# TRSS weights, Cam's proliferation gene sets, repstress signature

library(tidyverse)
library(msigdbr)

# ── TRSS weights (Jungk & Kschischo 2025, Cell Stress) ──────────────────────
# 17/21 genes present — 4 histone genes absent due to poly-A capture limitation
trss_weights <- c(
  ARID2 = 3.252, TOP3A = 2.630, TAF12 = 2.358,
  AOC2 = 2.063, BTG2 = 1.906, TFDP2 = 1.848,
  PHF8 = 1.729, ATF2 = 1.591, GLI1 = -0.085,
  POLD4 = -0.144, TUBB = -0.390, BCAM = -0.793,
  UFL1 = -1.210, CMPK2 = -1.329, FBXL7 = -1.626,
  TP53 = -2.693, PPP2R5A = -3.138
)

# ── Proliferation gene sets ───────────────────────────────────────────

# E2F targets (Hallmark, MSigDB)
e2f_genes <- msigdbr(species = "Homo sapiens",
                     collection = "H") %>%
  filter(gs_name == "HALLMARK_E2F_TARGETS") %>%
  pull(gene_symbol) %>%
  unique()

# G2M checkpoint (Hallmark, MSigDB)
g2m_genes <- msigdbr(species = "Homo sapiens",
                     collection = "H") %>%
  filter(gs_name == "HALLMARK_G2M_CHECKPOINT") %>%
  pull(gene_symbol) %>%
  unique()

# Osteoblast proliferation (GO:0033690)
osteo_prolif_genes <- msigdbr(species = "Homo sapiens",
                              collection = "C5",
                              subcollection = "GO:BP") %>%
  filter(gs_name == "GOBP_POSITIVE_REGULATION_OF_OSTEOBLAST_PROLIFERATION") %>%
  pull(gene_symbol) %>%
  unique()

# ── Repstress signature (Takahashi et al. 2022, Cancer Res Commun) ──────────
repstress_genes <- c(
  "AURKB", "CCNA2", "GINS1", "LIG3", "MTF2",
  "ORC6", "PRPS1", "SRSF1", "SUV39H1", "TNPO2",
  "GADD45G", "POLA1", "POLD4", "POLE4", "RFC5",
  "RMI1", "RRM1"
)

# ── Summary ─────────────────────────────────────────────────────────────────
cat("TRSS genes:", length(trss_weights), "\n")
cat("E2F targets:", length(e2f_genes), "\n")
cat("G2M checkpoint:", length(g2m_genes), "\n")
cat("Osteoblast proliferation:", length(osteo_prolif_genes), "\n")
cat("Repstress genes:", length(repstress_genes), "\n")

# Save gene signatures as a list
gene_signatures <- list(
  trss_weights = trss_weights,
  e2f_genes = e2f_genes,
  g2m_genes = g2m_genes,
  osteo_prolif_genes = osteo_prolif_genes,
  repstress_genes = repstress_genes
)

saveRDS(gene_signatures, "data/gene_signatures.rds")