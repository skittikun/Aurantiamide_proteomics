# =============================================================================
# 02_normalise.R
# Compare normalisation methods with NormalyzerDE, then apply VSN to the
# filtered (non-imputed) intensities.
#
# Output:
#   results/normalyzer/DN_Loess_VSN_Comparison/   NormalyzerDE report
#   data/processed/vsn_normalised.rds             VSN matrix (glog2 scale)
#   results/tables/vsn_normalised_intensity.tsv   same, as text
# =============================================================================

source("R/00_config.R")

suppressPackageStartupMessages({
  library(NormalyzerDE)
  library(limma)
  library(vsn)
})

prep <- readRDS(file.path(dir_processed, "prepared.rds"))
intensity_filtered <- prep$intensity_filtered

# ---- 1. Benchmark normalisation methods --------------------------------------
# NormalyzerDE log2-transforms the input; zeros are treated as missing.
dir_norm <- file.path(dir_results, "normalyzer")
dir.create(dir_norm, showWarnings = FALSE)

NormalyzerDE::normalyzer(
  jobName        = "DN_Loess_VSN_Comparison",
  designPath     = file.path(dir_processed, "design.tsv"),
  dataPath       = file.path(dir_processed, "intensity_raw.tsv"),
  outputDir      = dir_norm,
  noLogTransform = FALSE,
  zeroToNA       = TRUE,
  sampleColName  = "sample",
  groupColName   = "group"
)

# ---- 2. Candidate methods on the filtered matrix ------------------------------
# Cyclic loess kept for comparison only; VSN was selected.
loess_normalised <- limma::normalizeCyclicLoess(log2(intensity_filtered))

vsn_normalised <- NormalyzerDE::performVSNNormalization(intensity_filtered)
vsn_normalised <- as.matrix(as.data.frame(vsn_normalised))
rownames(vsn_normalised) <- rownames(intensity_filtered)
colnames(vsn_normalised) <- colnames(intensity_filtered)

# VSN output is already on a generalised log2 scale: do not log again.
saveRDS(vsn_normalised, file.path(dir_processed, "vsn_normalised.rds"))

write.table(data.frame(Protein = rownames(vsn_normalised), vsn_normalised,
                       check.names = FALSE),
            file.path(dir_tables, "vsn_normalised_intensity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
