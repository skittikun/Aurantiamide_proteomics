# =============================================================================
# 00_config.R
# Shared paths, sample design, comparison labels and thresholds.
# Sourced at the top of every analysis script.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
})

# ---- Paths (relative to the repository root) --------------------------------
dir_raw       <- "/Users/ksongsom/Library/CloudStorage/OneDrive-Personal/SideHustles/Mai/Aurantiamide/manuscripts/github_update20261105/data/raw"
dir_processed <- "/Users/ksongsom/Library/CloudStorage/OneDrive-Personal/SideHustles/Mai/Aurantiamide/manuscripts/github_update20261105/data/processed"
dir_results   <- "results"
dir_figures   <- file.path(dir_results, "figures")
dir_tables    <- file.path(dir_results, "tables")

for (d in c(dir_processed, dir_results, dir_figures, dir_tables)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# MaxQuant protein-group intensity table (deposited on Zenodo)
file_protein_groups <- file.path(dir_raw,
  "proteinGroups_intensity_p-Value-posthoc-annotation_Q_cutoff.xlsx")

# UniProt ID mapping snapshot (retrieved 12 October 2026; deposited on Zenodo)
file_uniprot_map <- file.path(dir_raw, "idmapping_2026_10_12.tsv")

# ---- Groups ------------------------------------------------------------------
group_levels <- c("control", "AB", "Ar1AB", "Ar10AB")

group_labels <- c(control = "Control",
                  AB      = "AB42",
                  Ar1AB   = "AB42 + Auranofin 1",
                  Ar10AB  = "AB42 + Auranofin 10")

group_cols <- c("Control"             = "grey40",
                "AB42"                = "#D7301F",
                "AB42 + Auranofin 1"  = "#FDAE61",
                "AB42 + Auranofin 10" = "#2C7BB6")

# Assign group from the sample name (independent of column order)
assign_group <- function(sample) {
  g <- dplyr::case_when(
    grepl("^SF",   sample) ~ "control",
    grepl("^Ar10", sample) ~ "Ar10AB",
    grepl("^Ar1-", sample) ~ "Ar1AB",
    grepl("^AB",   sample) ~ "AB"
  )
  factor(g, levels = group_levels)
}

# ---- Contrasts ---------------------------------------------------------------
comparison_levels <- c("AB_vs_control",
                       "Ar1AB_vs_control",
                       "Ar10AB_vs_control",
                       "Ar1AB_vs_AB",
                       "Ar10AB_vs_AB",
                       "Ar10AB_vs_Ar1AB")

comparison_labels <- c(AB_vs_control     = "AB42 vs Control",
                       Ar1AB_vs_control  = "AB42 + Auranofin 1 vs Control",
                       Ar10AB_vs_control = "AB42 + Auranofin 10 vs Control",
                       Ar1AB_vs_AB       = "AB42 + Auranofin 1 vs AB42",
                       Ar10AB_vs_AB      = "AB42 + Auranofin 10 vs AB42",
                       Ar10AB_vs_Ar1AB   = "Auranofin 10 vs Auranofin 1")

# Comparisons shown in the main figure (Fig. 3B-D)
main_comparisons <- c("AB_vs_control", "Ar1AB_vs_AB", "Ar10AB_vs_AB")

# ---- Thresholds --------------------------------------------------------------
p_cut       <- 0.001  # nominal (unadjusted) limma p value
fc_cut      <- 1      # |log2 fold change| (1 = 2-fold)
min_valid   <- 3      # valid values required in at least one group
n_heatmap   <- 50     # proteins shown in the heatmap
n_label     <- 10     # proteins labelled per volcano panel
