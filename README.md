# Aurantiamide proteomics analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21932776.svg)](https://doi.org/10.5281/zenodo.21932776)

Analysis code for Structure-based virtual screening identifies aurantiamide acetate as an IDO1 inhibitor with neuroprotective effects against Aβ42 toxicity. Archived at v1.0.0.

## Design
16 samples: AB, Ar1-AB, Ar10-AB, SF (control), n = 4 each.
17,803 protein groups (MaxQuant).

## Quick start
1. Download `data_archive.zip` from Zenodo (DOI above), unpack into `data/raw/`.
2. `renv::restore()`
3. `Rscript each R`

## Scripts
| Script | Purpose |
|---|---|
| `00_config.R` | Paths, seeds, thresholds, design helper. No absolute paths. |
| `01_impute_missforest.R` | Zeros -> NA, missForest imputation (seed 111). |
| `02_normalise_and_limma.R` | NormalyzerDE benchmark, cyclic loess vs VSN, VSN chosen; limma + 6 contrasts. |
| `03_pca.R` | PCA on raw (log2), imputed (log2), VSN (already glog2). |
| `04_enrichment_volcano.R` | Volcanoes, UpSet, topGO, clusterProfiler, KEGG, Enrichr. |


## Reproducibility notes
- missForest is stochastic; seed fixed in `00_config.R`. The exact imputed
  matrix used for the paper is deposited on Zenodo.
- VSN output is on a generalised log2 scale, so it is NOT logged again
  before PCA or limma. logFC values are already log2 differences.
- Enrichr is a live web service. `USE_LIVE_ENRICHR <- FALSE`
  replays the deposited `ENRICHR.rds` snapshot (queried 2025-02-08).

## Data availability
Raw MS: PRIDE PXD[TO BE ADDED]. Processed data + intermediates: Zenodo DOI above.
