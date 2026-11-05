# =============================================================================
# 04_pca.R
# PCA on VSN-normalised intensities (proteins with no missing values),
# centred, not scaled.
#
# Output:
#   results/figures/pca_plot.png
#   results/tables/pca_variance_explained.tsv
# =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

prep   <- readRDS(file.path(dir_processed, "prepared.rds"))
design <- prep$design
vsn_normalised <- readRDS(file.path(dir_processed, "vsn_normalised.rds"))

# PCA cannot handle NA: complete cases only
pca_mat <- vsn_normalised[complete.cases(vsn_normalised), ]
cat("Proteins used for PCA:", nrow(pca_mat), "\n")

pca      <- prcomp(t(pca_mat), center = TRUE, scale. = FALSE)
var_expl <- round(100 * summary(pca)$importance[2, 1:4], 1)
print(var_expl)
write.table(data.frame(PC = names(var_expl), pct_variance = var_expl),
            file.path(dir_tables, "pca_variance_explained.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pca_df <- as.data.frame(pca$x[, 1:4]) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(Group = factor(group_labels[as.character(assign_group(Sample))],
                               levels = group_labels))

hull_df <- pca_df %>%
  dplyr::group_by(Group) %>%
  dplyr::slice(chull(PC1, PC2)) %>%
  dplyr::ungroup()

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Group, label = Sample)) +
  geom_polygon(data = hull_df, aes(fill = Group, group = Group),
               alpha = 0.15, colour = NA, show.legend = FALSE) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  labs(x = paste0("PC1 (", var_expl[1], "%)"),
       y = paste0("PC2 (", var_expl[2], "%)"),
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave(file.path(dir_figures, "pca_plot.png"), p_pca,
       width = 6, height = 5, dpi = 300)
