# Run the full analysis from the repository root:
#   Rscript run_all.R

scripts <- c("R/01_prepare_data.R",
             "R/02_normalise.R",
             "R/03_limma.R",
             "R/04_pca.R",
             "R/05_volcano.R",
             "R/06_heatmap.R",
             "R/07_qc_sample_correlation.R")

for (s in scripts) {
  message("\n==== ", s, " ====")
  source(s, local = new.env())
}

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
