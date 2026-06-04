# 06_statistics.R
# Statistical testing of TRSS corrected scores across tumor subpopulations
# Kruskal-Wallis + pairwise Wilcoxon with Benjamini-Hochberg correction

library(qs)
library(Seurat)
library(tidyverse)

# Load data
tumor_prim <- qread("data/tumor_prim.qs")

# ── Full dataset statistics ──────────────────────────────────────────────────
cat("=== Full dataset (all cell cycle phases) ===\n\n")

# Overall test
kruskal_full <- kruskal.test(TRSS_corrected_v2 ~ Ann_Level3,
                             data = tumor_prim@meta.data)
print(kruskal_full)

# Pairwise comparisons
pairwise_full <- pairwise.wilcox.test(
  tumor_prim$TRSS_corrected_v2,
  tumor_prim$Ann_Level3,
  p.adjust.method = "BH"
)
print(pairwise_full)

# ── G1-only statistics ───────────────────────────────────────────────────────
cat("\n=== G1 cells only (non-dividing) ===\n\n")

# Subset to G1
tumor_g1 <- subset(tumor_prim, Phase == "G1")
cat("G1 cells:", ncol(tumor_g1), "\n")
table(tumor_g1$Ann_Level3)

# G1 summary
g1_summary <- tumor_g1@meta.data %>%
  group_by(Ann_Level3) %>%
  summarise(
    n_cells = n(),
    TRSS_corrected_mean = mean(TRSS_corrected_v2),
    TRSS_corrected_median = median(TRSS_corrected_v2)
  ) %>%
  arrange(desc(TRSS_corrected_mean))

print(g1_summary)

# Overall test
kruskal_g1 <- kruskal.test(TRSS_corrected_v2 ~ Ann_Level3,
                           data = tumor_g1@meta.data)
print(kruskal_g1)

# Pairwise comparisons
pairwise_g1 <- pairwise.wilcox.test(
  tumor_g1$TRSS_corrected_v2,
  tumor_g1$Ann_Level3,
  p.adjust.method = "BH"
)
print(pairwise_g1)

# ── Save results ─────────────────────────────────────────────────────────────
write.csv(g1_summary,
          "results/tables/TRSS_G1_summary.csv",
          row.names = FALSE)

# Save G1 object for visualization
qsave(tumor_g1, "data/tumor_g1.qs")