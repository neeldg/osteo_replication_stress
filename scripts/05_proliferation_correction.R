# 05_proliferation_correction.R
# Regress out proliferation signal from TRSS score
# Uses Cam's gene sets: E2F targets and G2M checkpoint
# Note: osteoblast proliferation gene set excluded per Cam's recommendation
# — GO:0033690 captures normal osteoblast biology, not malignant tumor 
#   cell proliferation programs

library(qs)
library(Seurat)
library(tidyverse)

# Load data
tumor_prim <- qread("data/tumor_prim.qs")

# ── Proliferation correction using E2F + G2M only ───────────────────────────
trss_corrected_v3 <- residuals(
  lm(TRSS_score ~ E2F_score1 + G2M_score1,
     data = tumor_prim@meta.data)
)

tumor_prim$TRSS_corrected_v3 <- trss_corrected_v3

# Also correct repstress for comparison
repstress_corrected <- residuals(
  lm(Repstress1 ~ E2F_score1 + G2M_score1,
     data = tumor_prim@meta.data)
)

tumor_prim$Repstress_corrected <- repstress_corrected

# ── Summary table ────────────────────────────────────────────────────────────
score_summary_v3 <- tumor_prim@meta.data %>%
  group_by(Ann_Level3) %>%
  summarise(
    n_cells = n(),
    TRSS_corrected_mean = mean(TRSS_corrected_v3),
    TRSS_corrected_median = median(TRSS_corrected_v3),
    Repstress_raw_mean = mean(Repstress1),
    Repstress_corrected_mean = mean(Repstress_corrected)
  ) %>%
  arrange(desc(TRSS_corrected_mean))

print(score_summary_v3)

# ── Save ─────────────────────────────────────────────────────────────────────
dir.create("results/tables", showWarnings = FALSE)

write.csv(score_summary_v3,
          "results/tables/TRSS_corrected_v3_summary.csv",
          row.names = FALSE)

qsave(tumor_prim, "data/tumor_prim.qs")