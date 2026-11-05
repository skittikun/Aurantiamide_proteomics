# =============================================================================
# 06_heatmap.R
# Heatmap of the top n_heatmap differentially abundant proteins
# (p < p_cut, |log2FC| >= fc_cut in at least one contrast; lowest p kept per
# protein). Un-scaled VSN intensities; missing values shown in white.
# Rows split by GO-keyword functional category and clustered within each
# category (Pearson distance, complete linkage). Columns ordered by group.
#
# Output:
#   results/figures/heatmap_function.png
#   results/tables/heatmap_top_proteins.xlsx
# =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

prep   <- readRDS(file.path(dir_processed, "prepared.rds"))
design <- prep$design
vsn_normalised <- readRDS(file.path(dir_processed, "vsn_normalised.rds"))
res <- readRDS(file.path(dir_processed, "limma_results.rds"))

# ---- 1. Select proteins ------------------------------------------------------
top_prot <- res$all_results_uniprot %>%
  dplyr::filter(P_value < p_cut, abs(log2FC) >= fc_cut) %>%
  dplyr::group_by(Protein) %>%
  dplyr::slice_min(P_value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::slice_min(P_value, n = n_heatmap, with_ties = FALSE) %>%
  dplyr::mutate(gene = dplyr::case_when(
    !is.na(Gene.Names) & Gene.Names != "" ~ word(Gene.Names, 1),
    !is.na(Entry)                         ~ Entry,
    TRUE                                  ~ tomatch))

# ---- 2. Functional category (first matching keyword group) -------------------
assign_function <- function(go) {
  go <- tolower(go)
  dplyr::case_when(
    is.na(go) | go == "" ~ "Other",
    str_detect(go, "mitochond|respirat|oxidative phosph|atp synth")              ~ "Mitochondria / energy",
    str_detect(go, "proteasom|ubiquitin|protein folding|chaperone|autophag")     ~ "Proteostasis",
    str_detect(go, "oxidative stress|reactive oxygen|antioxidant|glutathione")   ~ "Oxidative stress",
    str_detect(go, "inflammat|immune|cytokine|interferon")                       ~ "Immune / inflammation",
    str_detect(go, "translation|ribosom|rna processing|splicing")                ~ "RNA / translation",
    str_detect(go, "cytoskelet|actin|microtubule|cell adhesion")                 ~ "Cytoskeleton / adhesion",
    str_detect(go, "transport|vesicle|endocyt|exocyt|synap")                     ~ "Transport / vesicle",
    str_detect(go, "apopto|cell death")                                          ~ "Apoptosis",
    str_detect(go, "signal transduction|kinase|phosphorylation")                 ~ "Signalling",
    TRUE ~ "Other"
  )
}

top_prot <- top_prot %>%
  dplyr::mutate(Function = assign_function(Gene.Ontology..biological.process.))

write.xlsx(top_prot, file.path(dir_tables, "heatmap_top_proteins.xlsx"))

# ---- 3. Matrix and annotations -----------------------------------------------
hm_mat <- vsn_normalised[top_prot$Protein, , drop = FALSE]
rownames(hm_mat) <- top_prot$gene

sample_group <- factor(group_labels[as.character(assign_group(colnames(hm_mat)))],
                       levels = group_labels)
ord          <- order(sample_group)
hm_mat       <- hm_mat[, ord]
sample_group <- sample_group[ord]

fun_levels <- c(setdiff(sort(unique(top_prot$Function)), "Other"),
                intersect("Other", top_prot$Function))
fun_cols   <- setNames(brewer.pal(max(3, length(fun_levels)), "Set3")[seq_along(fun_levels)],
                       fun_levels)
row_fun    <- factor(top_prot$Function, levels = fun_levels)

int_col  <- colorRamp2(quantile(hm_mat, c(0.02, 0.98), na.rm = TRUE),
                       c("#F7FBFF", "#08306B"))
rng      <- quantile(hm_mat, c(0.02, 0.50, 0.98), na.rm = TRUE)
row_mean <- rowMeans(hm_mat, na.rm = TRUE)

col_ha <- HeatmapAnnotation(Group = sample_group,
                            col = list(Group = group_cols),
                            show_annotation_name = FALSE)

row_ha <- rowAnnotation(
  Function = row_fun,
  `Mean log2\nintensity` = anno_simple(row_mean, col = int_col),
  col = list(Function = fun_cols),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7),
  annotation_name_rot = 90
)

# ---- 4. Heatmap --------------------------------------------------------------
ht <- Heatmap(hm_mat,
              name = "log2 intensity\n(VSN)",
              col  = int_col,
              na_col = "white",
              top_annotation = col_ha,
              right_annotation = row_ha,
              column_split = sample_group,
              column_title_gp = gpar(fontsize = 8),
              column_title_rot = 0,
              row_split = row_fun,
              cluster_columns = FALSE,
              cluster_row_slices = FALSE,
              clustering_distance_rows = "pearson",
              clustering_method_rows = "complete",
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 9),
              column_names_rot = 90,
              row_title_rot = 0,
              row_title_gp = gpar(fontsize = 8),
              heatmap_legend_param = list(
                title_gp  = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 8),
                at = round(seq(rng[1], rng[3], length.out = 5), 1)),
              border = TRUE)

png(file.path(dir_figures, "heatmap_function.png"),
    width = 9, height = 10, units = "in", res = 300)
draw(ht, merge_legend = TRUE)
dev.off()
