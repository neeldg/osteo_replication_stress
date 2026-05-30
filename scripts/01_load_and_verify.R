library(qs)
library(Seurat)
library(tidyverse)

# Load data
patient_prim <- qread("data/patient_prim.qs")

# Subset to tumor cells only
tumor_prim <- subset(patient_prim, Ann_Level0 == "Tumor")
rm(patient_prim)
gc()

# Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

tumor_prim <- CellCycleScoring(
  tumor_prim,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)

# TRSS weighted score
trss_weights <- c(
  ARID2 = 3.252, TOP3A = 2.630, TAF12 = 2.358,
  AOC2 = 2.063, BTG2 = 1.906, TFDP2 = 1.848,
  PHF8 = 1.729, ATF2 = 1.591, GLI1 = -0.085,
  POLD4 = -0.144, TUBB = -0.390, BCAM = -0.793,
  UFL1 = -1.210, CMPK2 = -1.329, FBXL7 = -1.626,
  TP53 = -2.693, PPP2R5A = -3.138
)

expr_matrix <- GetAssayData(tumor_prim, layer = "data")
trss_expr <- expr_matrix[names(trss_weights), ]
trss_score <- as.vector(t(trss_expr) %*% trss_weights)
tumor_prim$TRSS_score <- trss_score

# Mechanistic RS gene modules
tier1_genes <- c("RPA1", "RPA2", "RPA3", "CHEK1", "SMARCAL1",
                 "ZRANB3", "FANCD2", "SETX", "DHX9", "RRM2",
                 "TIPIN", "TIMELESS", "CLSPN")

tier2_genes <- c("ATR", "ATRIP", "TOPBP1", "HLTF", "RADX",
                 "FANCI", "FANCA", "RAD51", "BRCA1", "BRCA2",
                 "AQR", "TOP1")

tumor_prim <- AddModuleScore(
  tumor_prim,
  features = list(tier1_genes),
  name = "Mechanistic_T1"
)

tumor_prim <- AddModuleScore(
  tumor_prim,
  features = list(tier2_genes),
  name = "Mechanistic_T2"
)

# Cell cycle correction via regression
tumor_prim$TRSS_corrected <- residuals(
  lm(TRSS_score ~ S.Score + G2M.Score,
     data = tumor_prim@meta.data))

tumor_prim$Mech_T1_corrected <- residuals(
  lm(Mechanistic_T11 ~ S.Score + G2M.Score,
     data = tumor_prim@meta.data))

tumor_prim$Mech_T2_corrected <- residuals(
  lm(Mechanistic_T21 ~ S.Score + G2M.Score,
     data = tumor_prim@meta.data))

# Summary table
score_summary_corrected <- tumor_prim@meta.data %>%
  group_by(Ann_Level3) %>%
  summarise(
    n_cells = n(),
    TRSS_corrected_mean = mean(TRSS_corrected),
    Mech_T1_corrected_mean = mean(Mech_T1_corrected),
    Mech_T2_corrected_mean = mean(Mech_T2_corrected)
  ) %>%
  arrange(desc(TRSS_corrected_mean))

print(score_summary_corrected)

# Statistical tests
kruskal.test(TRSS_corrected ~ Ann_Level3,
             data = tumor_prim@meta.data)

pairwise.wilcox.test(tumor_prim$TRSS_corrected,
                     tumor_prim$Ann_Level3,
                     p.adjust.method = "BH")

# Visualization
tumor_prim@meta.data$Ann_Level3 <- factor(
  tumor_prim@meta.data$Ann_Level3,
  levels = c("Basal_Progenitor", "COMA", "Proliferative",
             "Fibrogenic", "Interactive", "MP_Progenitor")
)

p_trss <- ggplot(tumor_prim@meta.data,
                 aes(x = Ann_Level3,
                     y = TRSS_corrected,
                     fill = Ann_Level3)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, fill = "white",
               outlier.size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Replication Stress (TRSS) by Tumor Subpopulation",
       subtitle = "Corrected for cell cycle effects",
       x = "",
       y = "Corrected TRSS Score")

ggsave("results/figures/TRSS_corrected_by_subpopulation.pdf",
       plot = p_trss, width = 8, height = 6)

ggsave("results/figures/TRSS_corrected_by_subpopulation.png",
       plot = p_trss, width = 8, height = 6, dpi = 300)

print(p_trss)