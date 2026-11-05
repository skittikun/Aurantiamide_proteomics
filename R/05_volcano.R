# =============================================================================
# 05_volcano.R
# Volcano plots: log2FC vs -log10 nominal p. Dashed lines at |log2FC| = fc_cut
# and p = p_cut; top n_label significant proteins labelled per panel.
#
# Output:
#   results/figures/volcano_all_comparisons.png      all six contrasts
#   results/figures/volcano_main_comparisons.png     Fig. 3B-D
# =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

res <- readRDS(file.path(dir_processed, "limma_results.rds"))

volcano_df <- res$all_results_uniprot %>%
  dplyr::filter(!is.na(P_value), !is.na(log2FC)) %>%
  dplyr::mutate(
    neg_log10_p = -log10(P_value),
    gene = ifelse(is.na(Gene.Names) | Gene.Names == "",
                  ifelse(is.na(Entry), tomatch, Entry),
                  sub(" .*", "", Gene.Names)),
    status = dplyr::case_when(
      P_value < p_cut & log2FC >=  fc_cut ~ "Up",
      P_value < p_cut & log2FC <= -fc_cut ~ "Down",
      TRUE                                ~ "Not significant"),
    status = factor(status, levels = c("Up", "Down", "Not significant"))
  )

stopifnot(all(volcano_df$Comparison %in% comparison_levels))
print(table(volcano_df$Comparison, volcano_df$status))

plot_volcano <- function(df) {
  df <- df %>%
    dplyr::mutate(Comparison = factor(Comparison,
                                      levels = intersect(comparison_levels,
                                                         unique(Comparison)),
                                      labels = comparison_labels[
                                        intersect(comparison_levels,
                                                  unique(Comparison))]))
  label_df <- df %>%
    dplyr::filter(status != "Not significant") %>%
    dplyr::group_by(Comparison) %>%
    dplyr::slice_min(P_value, n = n_label) %>%
    dplyr::ungroup()

  ggplot(df, aes(x = log2FC, y = neg_log10_p, colour = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(p_cut),      linetype = "dashed", colour = "grey40") +
    geom_text_repel(data = label_df, aes(label = gene),
                    size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = c("Up" = "#D7301F",
                                   "Down" = "#2C7BB6",
                                   "Not significant" = "grey70")) +
    facet_wrap(~ Comparison, nrow = 1) +
    labs(x = expression(log[2]~fold~change),
         y = expression(-log[10]~italic(p-value)),
         colour = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          strip.background = element_rect(fill = "grey95"))
}

ggsave(file.path(dir_figures, "volcano_all_comparisons.png"),
       plot_volcano(volcano_df), width = 12, height = 4.5, dpi = 300)

# %in% (not ==) so every protein in the three main contrasts is kept
ggsave(file.path(dir_figures, "volcano_main_comparisons.png"),
       plot_volcano(dplyr::filter(volcano_df, Comparison %in% main_comparisons)),
       width = 12, height = 4.5, dpi = 300)
