# Aurantiamide proteomics analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21932775.svg)](https://doi.org/10.5281/zenodo.21932775)


Analysis code for *Structure-based virtual screening identifies aurantiamide acetate as an IDO1 inhibitor with neuroprotective effects against Aβ42 toxicity*. Archived at v[X.Y.Z].

## Design
12 samples: SF (control), AB, Ar1-AB, Ar10-AB; n = 3 each.
Input: MaxQuant protein-group intensity table.

## Quick start
1. Download `data_archive.zip` from Zenodo (DOI above) and unpack it into `data/raw/`.
   It contains `proteinGroups_intensity_p-Value-posthoc-annotation.xlsx` and
   `idmapping_2026_10_12.tsv`.
2. `renv::restore()`
3. From the repository root: `Rscript run_all.R`

## Scripts
| Script | Purpose |
|---|---|
| `R/00_config.R` | Paths, group labels and colours, contrasts, thresholds. Sourced by every script. |
| `R/01_prepare_data.R` | Read MaxQuant table, build design, zeros -> NA (no imputation), keep proteins with ≥ 3 valid values in at least one group. |
| `R/02_normalise.R` | NormalyzerDE benchmark; cyclic loess vs VSN; VSN chosen and applied to the filtered matrix. |
| `R/03_limma.R` | limma on VSN data with missing values kept; 6 contrasts; BH-adjusted p; UniProt annotation; result tables and counts per contrast. |
| `R/04_pca.R` | PCA on VSN data (proteins with no missing values; centred, not scaled). |
| `R/05_volcano.R` | Volcano plots for all six contrasts and for the three main contrasts (Fig. 3B–D). |
| `R/06_heatmap.R` | Top 50 proteins; un-scaled VSN intensities; GO-keyword functional categories (Fig. 3E). |
| `R/07_qc_sample_correlation.R` | Sample QC: valid values, medians, correlation, dendrogram. |
| `run_all.R` | Runs 01–07 in order and saves `sessionInfo()`. |

## Outputs
- `data/processed/` – intermediate `.rds` files and the design table.
- `results/tables/` – VSN matrix, limma results (`all_` and `significant_normalized_noimpute_limma_uniprot.xlsx`), `quick_summary.xlsx`, PCA variance, heatmap protein list, QC tables.
- `results/figures/` – PCA, volcano, heatmap, QC dendrogram.
- `results/normalyzer/` – NormalyzerDE normalisation report.

## Thresholds
Set in `R/00_config.R`:
- Differentially abundant: nominal limma p < 0.001 and |log2FC| ≥ 1.
- BH-adjusted p values are calculated per contrast and included in the result tables.

## Reproducibility notes
- No imputation. Zero intensities are treated as missing; limma fits each protein on its available values.
- VSN output is on a generalised log2 scale, so it is NOT logged again
  before PCA or limma. logFC values are already log2 differences.
- Groups are assigned from sample names, not column order, and checked with `stopifnot()`.
- UniProt annotation uses the deposited ID-mapping snapshot (`idmapping_2026_10_12.tsv`).
  To refresh it, upload `results/tables/to_uniprot.txt` to
  [UniProt ID mapping](https://www.uniprot.org/id-mapping) and update `file_uniprot_map` in `R/00_config.R`.
- Protein–protein interaction and GO enrichment were run on the STRING web server (v12.0) using the selected protein list.

## Data availability
Raw MS: PRIDE PXD082592. Processed data + intermediates: Zenodo DOI above.
