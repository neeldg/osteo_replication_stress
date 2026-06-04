# 06_statistics.R
# Statistical testing of TRSS corrected scores across tumor subpopulations

library(qs)
library(Seurat)
library(tidyverse)

# Load data
tumor_prim <- qread("data/tumor_prim.qs")

# ── Full dataset statistics ──────────────────────────────────────────────────
cat("=== Full dataset (all cell cycle phases) ===\n\n")

kruskal_full <- kruskal.test(TRSS_corrected_v3 ~ Ann_Level3,
                             data = tumor_prim@meta.data)
print(kruskal_full)

pairwise_full <- pairwise.wilcox.test(
  tumor_prim$TRSS_corrected_v3,
  tumor_prim$Ann_Level3,
  p.adjust.method = "BH"
)
print(pairwise_full)

# ── G1-only statistics ───────────────────────────────────────────────────────
cat("\n=== G1 cells only (non-dividing) ===\n\n")

tumor_g1 <- subset(tumor_prim, Phase == "G1")
cat("G1 cells:", ncol(tumor_g1), "\n")
table(tumor_g1$Ann_Level3)

g1_summary <- tumor_g1@meta.data %>%
  group_by(Ann_Level3) %>%
  summarise(
    n_cells = n(),
    TRSS_corrected_mean = mean(TRSS_corrected_v3),
    TRSS_corrected_median = median(TRSS_corrected_v3)
  ) %>%
  arrange(desc(TRSS_corrected_mean))

print(g1_summary)

kruskal_g1 <- kruskal.test(TRSS_corrected_v3 ~ Ann_Level3,
                           data = tumor_g1@meta.data)
print(kruskal_g1)

pairwise_g1 <- pairwise.wilcox.test(
  tumor_g1$TRSS_corrected_v3,
  tumor_g1$Ann_Level3,
  p.adjust.method = "BH"
)
print(pairwise_g1)

# ── Save ─────────────────────────────────────────────────────────────────────
write.csv(g1_summary,
          "results/tables/TRSS_G1_summary.csv",
          row.names = FALSE)

qsave(tumor_g1, "data/tumor_g1.qs")