# =============================================================================
# 07_qc_sample_correlation.R
# Sample-level QC on filtered log2 intensities (before normalisation):
# valid values, median, total intensity, pairwise correlation and a
# correlation-distance dendrogram.
#
# Output:
#   results/tables/qc_sample_summary.tsv
#   results/tables/qc_sample_correlation.tsv
#   results/figures/qc_sample_dendrogram.png
# =============================================================================

source("R/00_config.R")

prep <- readRDS(file.path(dir_processed, "prepared.rds"))
intensity_filtered <- prep$intensity_filtered

lg <- log2(intensity_filtered)

qc <- data.frame(sample    = colnames(lg),
                 n_valid   = colSums(!is.na(lg)),
                 median    = round(apply(lg, 2, median, na.rm = TRUE), 2),
                 total_raw = signif(colSums(intensity_filtered, na.rm = TRUE), 3))
print(qc)
write.table(qc, file.path(dir_tables, "qc_sample_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cor_mat <- cor(lg, use = "pairwise.complete.obs")
write.table(data.frame(sample = rownames(cor_mat), round(cor_mat, 2), check.names = FALSE),
            file.path(dir_tables, "qc_sample_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(dir_figures, "qc_sample_dendrogram.png"),
    width = 7, height = 5, units = "in", res = 300)
plot(hclust(as.dist(1 - cor_mat)), main = "Sample correlation (1 - r)",
     xlab = "", sub = "")
dev.off()
