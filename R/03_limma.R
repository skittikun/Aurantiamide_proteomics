# =============================================================================
# 03_limma.R
# Differential abundance with limma on VSN-normalised, non-imputed data.
# Missing values are kept; each protein is fitted on its available values.
# Six pairwise contrasts; BH-adjusted p values are reported per contrast,
# and proteins are called significant at nominal p < p_cut and
# |log2FC| >= fc_cut.
#
# UniProt annotation: upload results/tables/to_uniprot.txt to
# https://www.uniprot.org/id-mapping (UniProtKB AC/ID -> UniProtKB) and save
# the TSV as data/raw/idmapping_YYYY_MM_DD.tsv. The snapshot used for the
# paper (idmapping_2026_10_12.tsv) is on Zenodo.
#
# Output:
#   data/processed/limma_results.rds
#   results/tables/all_normalized_noimpute_limma_uniprot.xlsx
#   results/tables/significant_normalized_noimpute_limma_uniprot.xlsx
#   results/tables/significant_noimpute_limma_summary.xlsx
#   results/tables/quick_summary.xlsx
# =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages(library(limma))

prep   <- readRDS(file.path(dir_processed, "prepared.rds"))
design <- prep$design
vsn_normalised <- readRDS(file.path(dir_processed, "vsn_normalised.rds"))
stopifnot(identical(colnames(vsn_normalised), design$sample))

# ---- 1. Model and contrasts --------------------------------------------------
groups <- design$group
design_matrix <- model.matrix(~ 0 + groups)
colnames(design_matrix) <- levels(groups)

contrasts <- makeContrasts(
  AB_vs_control     = AB     - control,
  Ar1AB_vs_control  = Ar1AB  - control,
  Ar10AB_vs_control = Ar10AB - control,
  Ar1AB_vs_AB       = Ar1AB  - AB,
  Ar10AB_vs_AB      = Ar10AB - AB,
  Ar10AB_vs_Ar1AB   = Ar10AB - Ar1AB,
  levels = design_matrix
)

fit        <- lmFit(vsn_normalised, design_matrix)
fit_ebayes <- eBayes(contrasts.fit(fit, contrasts))

# ---- 2. Long-format results --------------------------------------------------
to_long <- function(m, value) {
  as.data.frame(m) %>%
    tibble::rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Comparison", values_to = value)
}

all_results <- to_long(fit_ebayes$p.value, "P_value") %>%
  dplyr::group_by(Comparison) %>%
  dplyr::mutate(padj = p.adjust(P_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(to_long(fit_ebayes$coefficients, "coefficient"),
                   by = c("Protein", "Comparison")) %>%
  dplyr::mutate(
    log2FC      = round(coefficient, 2),
    fold_change = round(ifelse(coefficient > 0, 1, -1) * 2^abs(coefficient), 2),
    direction   = ifelse(coefficient > 0, "up", "down"),
    tomatch     = gsub(";.*", "", Protein)   # leading accession
  )

cat("Proteins passing BH FDR < 0.05 per contrast:\n")
print(all_results %>% dplyr::group_by(Comparison) %>%
        dplyr::summarise(n_fdr05 = sum(padj < 0.05, na.rm = TRUE)))

# ---- 3. UniProt annotation ---------------------------------------------------
writeLines(unique(all_results$tomatch),
           file.path(dir_tables, "to_uniprot.txt"))

uniprot_out <- read.csv(file_uniprot_map, sep = "\t")

all_results_uniprot <- merge(all_results, uniprot_out,
                             by.x = "tomatch", by.y = "From",
                             all.x = TRUE, all.y = FALSE)

# Nominal p < p_cut (no fold-change filter) for the full significant table
sig_results_uniprot <- all_results_uniprot %>%
  dplyr::filter(P_value < p_cut)

write.xlsx(all_results_uniprot,
           file.path(dir_tables, "all_normalized_noimpute_limma_uniprot.xlsx"))
write.xlsx(sig_results_uniprot,
           file.path(dir_tables, "significant_normalized_noimpute_limma_uniprot.xlsx"))

sig_results_summary <- sig_results_uniprot %>%
  dplyr::group_by(Protein) %>%
  dplyr::summarise(
    significant_comparison = paste(Comparison,    collapse = ", "),
    log2FC                 = paste(log2FC,        collapse = ", "),
    fold_change            = paste(fold_change,   collapse = ", "),
    direction              = paste(direction,     collapse = ", "),
    Protein.names          = paste(Protein.names, collapse = ", ")
  )
write.xlsx(sig_results_summary,
           file.path(dir_tables, "significant_noimpute_limma_summary.xlsx"))

# ---- 4. Counts per contrast (p < p_cut and |log2FC| >= fc_cut) ----------------
sig_summary <- all_results %>%
  dplyr::filter(!is.na(P_value), !is.na(log2FC)) %>%
  dplyr::mutate(status = dplyr::case_when(
    P_value < p_cut & log2FC >=  fc_cut ~ "Up",
    P_value < p_cut & log2FC <= -fc_cut ~ "Down",
    TRUE ~ "NS")) %>%
  dplyr::group_by(Comparison) %>%
  dplyr::summarise(
    Tested    = dplyr::n(),
    Up        = sum(status == "Up"),
    Down      = sum(status == "Down"),
    Total_sig = Up + Down,
    Pct_sig   = round(100 * Total_sig / Tested, 1),
    .groups   = "drop") %>%
  dplyr::mutate(Comparison = factor(Comparison, levels = comparison_levels)) %>%
  dplyr::arrange(Comparison)

print(sig_summary)
write.xlsx(sig_summary, file.path(dir_tables, "quick_summary.xlsx"))

saveRDS(list(fit = fit_ebayes,
             all_results = all_results,
             all_results_uniprot = all_results_uniprot),
        file.path(dir_processed, "limma_results.rds"))
