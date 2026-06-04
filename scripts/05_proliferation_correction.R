# 05_proliferation_correction.R
# Regress out proliferation signal from TRSS score
# Uses Cam's gene sets: E2F targets, G2M checkpoint, osteoblast proliferation

library(qs)
library(Seurat)
library(tidyverse)

# Load data
tumor_prim <- qread("data/tumor_prim.qs")

# ── v2 correction using proliferation gene sets ───────────────────────
# Regress TRSS score on E2F, G2M, and osteoblast proliferation scores
# Residual = proliferation-independent replication stress signal

trss_corrected_v2 <- residuals(
  lm(TRSS_score ~ E2F_score1 + G2M_score1 + Osteo_prolif1,
     data = tumor_prim@meta.data)
)

tumor_prim$TRSS_corrected_v2 <- trss_corrected_v2

# ── Summary table ────────────────────────────────────────────────────────────
score_summary_v2 <- tumor_prim@meta.data %>%
  group_by(Ann_Level3) %>%
  summarise(
    n_cells = n(),
    TRSS_corrected_mean = mean(TRSS_corrected_v2),
    TRSS_corrected_median = median(TRSS_corrected_v2),
    Repstress_mean = mean(Repstress1)
  ) %>%
  arrange(desc(TRSS_corrected_mean))

print(score_summary_v2)

# Save summary table
write.csv(score_summary_v2,
          "results/tables/TRSS_corrected_v2_summary.csv",
          row.names = FALSE)

# Save scored object
qsave(tumor_prim, "data/tumor_prim.qs")