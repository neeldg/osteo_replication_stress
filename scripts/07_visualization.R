# 07_visualization.R
# Generate all figures for primary tumor analysis:
# 1. Violin plot — corrected TRSS v3 all cells
# 2. Violin plot — G1 cells only
# 3. UMAP — corrected TRSS score
# 4. UMAP — subpopulation labels
# 5. Heatmap — TRSS gene expression across subpopulations

library(qs)
library(Seurat)
library(tidyverse)
library(pheatmap)

# Load data
tumor_prim <- qread("data/tumor_prim.qs")
tumor_g1 <- qread("data/tumor_g1.qs")
gene_signatures <- readRDS("data/gene_signatures.rds")

# ── Factor ordering by TRSS score ────────────────────────────────────────────
subpop_order <- c("Basal_Progenitor", "COMA", "Fibrogenic",
                  "Interactive", "MP_Progenitor", "Proliferative")

tumor_prim@meta.data$Ann_Level3 <- factor(
  tumor_prim@meta.data$Ann_Level3, levels = subpop_order)

tumor_g1@meta.data$Ann_Level3 <- factor(
  tumor_g1@meta.data$Ann_Level3, levels = subpop_order)

# ── Figure 1: Violin plot — all cells ────────────────────────────────────────
p_violin_all <- ggplot(tumor_prim@meta.data,
                       aes(x = Ann_Level3,
                           y = TRSS_corrected_v3,
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
       subtitle = "Corrected for E2F and G2M proliferation",
       x = "",
       y = "Corrected TRSS Score")

ggsave("results/figures/TRSS_corrected_v3_by_subpopulation.pdf",
       plot = p_violin_all, width = 8, height = 6)
ggsave("results/figures/TRSS_corrected_v3_by_subpopulation.png",
       plot = p_violin_all, width = 8, height = 6, dpi = 300)
print(p_violin_all)

# ── Figure 2: Violin plot — G1 cells only ────────────────────────────────────
p_violin_g1 <- ggplot(tumor_g1@meta.data,
                      aes(x = Ann_Level3,
                          y = TRSS_corrected_v3,
                          fill = Ann_Level3)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, fill = "white",
               outlier.size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Replication Stress (TRSS) in G1 Cells Only",
       subtitle = "Non-dividing cells — no proliferation correction needed",
       x = "",
       y = "Corrected TRSS Score")

ggsave("results/figures/TRSS_G1_only.pdf",
       plot = p_violin_g1, width = 8, height = 6)
ggsave("results/figures/TRSS_G1_only.png",
       plot = p_violin_g1, width = 8, height = 6, dpi = 300)
print(p_violin_g1)

# ── Figure 3: UMAP — corrected TRSS score ────────────────────────────────────
umap_coords <- as.data.frame(Embeddings(tumor_prim, reduction = "umap"))
umap_coords$TRSS_corrected_v3 <- tumor_prim$TRSS_corrected_v3
umap_coords$Ann_Level3 <- tumor_prim$Ann_Level3

p_umap_trss <- ggplot(umap_coords,
                      aes(x = umap_1, y = umap_2,
                          color = TRSS_corrected_v3)) +
  geom_point(size = 0.1, alpha = 0.6) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "lightgrey",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-5, 5),
    oob = scales::squish
  ) +
  theme_classic() +
  labs(title = "Corrected TRSS Score",
       subtitle = "Proliferation-independent replication stress",
       color = "TRSS\nScore",
       x = "UMAP 1",
       y = "UMAP 2")

ggsave("results/figures/UMAP_TRSS_corrected_v3.pdf",
       plot = p_umap_trss, width = 8, height = 6)
ggsave("results/figures/UMAP_TRSS_corrected_v3.png",
       plot = p_umap_trss, width = 8, height = 6, dpi = 300)
print(p_umap_trss)

# ── Figure 4: UMAP — subpopulation labels ────────────────────────────────────
p_umap_clusters <- DimPlot(
  tumor_prim,
  reduction = "umap",
  group.by = "Ann_Level3",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.1
) +
  theme(legend.position = "right") +
  labs(title = "Tumor Subpopulations")

ggsave("results/figures/UMAP_subpopulations.pdf",
       plot = p_umap_clusters, width = 8, height = 6)
ggsave("results/figures/UMAP_subpopulations.png",
       plot = p_umap_clusters, width = 8, height = 6, dpi = 300)
print(p_umap_clusters)

# ── Figure 5: Heatmap — TRSS gene expression ─────────────────────────────────
trss_gene_names <- names(gene_signatures$trss_weights)
expr_matrix <- GetAssayData(tumor_prim, layer = "data")
trss_expr <- expr_matrix[trss_gene_names, ]
trss_expr_df <- as.data.frame(t(as.matrix(trss_expr)))
trss_expr_df$Ann_Level3 <- tumor_prim$Ann_Level3

gene_summary <- trss_expr_df %>%
  group_by(Ann_Level3) %>%
  summarise(across(all_of(trss_gene_names), mean))

heatmap_matrix <- gene_summary %>%
  column_to_rownames("Ann_Level3") %>%
  as.matrix() %>%
  t()

heatmap_matrix <- heatmap_matrix[, subpop_order]

gene_weights_df <- data.frame(
  Weight = ifelse(gene_signatures$trss_weights > 0,
                  "Positive", "Negative"),
  row.names = names(gene_signatures$trss_weights)
)

annotation_colors <- list(
  Weight = c(Positive = "#B2182B", Negative = "#2166AC")
)

pdf("results/figures/TRSS_gene_heatmap.pdf", width = 8, height = 7)
pheatmap(
  heatmap_matrix,
  annotation_row = gene_weights_df,
  annotation_colors = annotation_colors,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "TRSS Gene Expression Across Tumor Subpopulations",
  fontsize = 10
)
dev.off()

cat("All figures saved to results/figures/\n")