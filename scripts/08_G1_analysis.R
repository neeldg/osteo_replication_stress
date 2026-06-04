# 08_G1_analysis.R
# Stratified analysis in G1 cells only
# Validates TRSS finding without regression assumptions
# G1 cells are definitionally non-dividing — no proliferation correction needed

library(qs)
library(Seurat)
library(tidyverse)

# Load data
tumor_prim <- qread("data/tumor_prim.qs")

# ── Subset to G1 cells ───────────────────────────────────────────────────────
tumor_g1 <- subset(tumor_prim, Phase == "G1")
cat("Total G1 cells:", ncol(tumor_g1), "\n")
table(tumor_g1$Ann_Level3)

# ── Summary table ────────────────────────────────────────────────────────────
g1_summary <- tumor_g1@meta.data %>%
  group_by(Ann_Level3) %>%
  summarise(
    n_cells = n(),
    TRSS_corrected_mean = mean(TRSS_corrected_v2),
    TRSS_corrected_median = median(TRSS_corrected_v2)
  ) %>%
  arrange(desc(TRSS_corrected_mean))

print(g1_summary)

# ── Statistical tests ────────────────────────────────────────────────────────
cat("\n=== Kruskal-Wallis test (G1 cells) ===\n")
kruskal_g1 <- kruskal.test(TRSS_corrected_v2 ~ Ann_Level3,
                           data = tumor_g1@meta.data)
print(kruskal_g1)

cat("\n=== Pairwise Wilcoxon (G1 cells, BH correction) ===\n")
pairwise_g1 <- pairwise.wilcox.test(
  tumor_g1$TRSS_corrected_v2,
  tumor_g1$Ann_Level3,
  p.adjust.method = "BH"
)
print(pairwise_g1)

# ── Save ─────────────────────────────────────────────────────────────────────
write.csv(g1_summary,
          "results/tables/TRSS_G1_summary.csv",
          row.names = FALSE)

qsave(tumor_g1, "data/tumor_g1.qs")

cat("G1 analysis complete.\n")