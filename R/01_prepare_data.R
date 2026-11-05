# =============================================================================
# 01_prepare_data.R
# Read the MaxQuant protein-group table, select samples, build the design,
# treat zero intensities as missing (no imputation) and apply the
# valid-value filter.
#
# Output:
#   data/processed/design.tsv              sample design (sample, group, replicate)
#   data/processed/intensity_raw.tsv       unfiltered intensities (for NormalyzerDE)
#   data/processed/prepared.rds            list(annot, design, intensity_all,
#                                               intensity_filtered)
# =============================================================================

source("/Users/ksongsom/Library/CloudStorage/OneDrive-Personal/SideHustles/Mai/Aurantiamide/manuscripts/github_update20261105/zenodo/github/R/00_config.R")

# ---- 1. Read and select samples ---------------------------------------------
origin_dt <- read.xlsx(file_protein_groups)

# Replicates N2-N4 are analysed and renumbered N1-N3
origin_dt <- origin_dt %>%
  dplyr::select(-matches("N1")) %>%
  dplyr::rename_with(~ .x %>%
                       str_replace("N2", "N1") %>%
                       str_replace("N3", "N2") %>%
                       str_replace("N4", "N3"))

intensity_cols <- grep("Intensity$", colnames(origin_dt), value = TRUE)
stopifnot(length(intensity_cols) == 12)

# ---- 2. Design (group assigned from sample name, not column order) ----------
design <- data.frame(
  sample    = sub("[-_]Intensity$", "", intensity_cols),
  stringsAsFactors = FALSE
)
design$group     <- assign_group(design$sample)
design$replicate <- as.integer(sub(".*N(\\d)$", "\\1", design$sample))

stopifnot(!any(is.na(design$group)),
          all(table(design$group) == 3))

# ---- 3. Intensity matrix -----------------------------------------------------
annot <- origin_dt %>% dplyr::select(Protein.IDs:Gene.names)

intensity_all <- as.matrix(origin_dt[, intensity_cols])
colnames(intensity_all) <- design$sample
rownames(intensity_all) <- origin_dt$Protein.IDs
stopifnot(identical(colnames(intensity_all), design$sample))

cat("Protein groups in input table:", nrow(intensity_all), "\n")

# Files for the NormalyzerDE benchmark (02_normalise.R)
write.table(design, file.path(dir_processed, "design.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(intensity_all, file.path(dir_processed, "intensity_raw.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Zeros -> NA (no imputation) ------------------------------------------
intensity_na <- intensity_all
intensity_na[intensity_na == 0] <- NA

# ---- 5. Valid-value filter ---------------------------------------------------
# Keep proteins with >= min_valid non-missing values in at least one group
valid_per_group <- sapply(levels(design$group), function(g)
  rowSums(!is.na(intensity_na[, design$group == g, drop = FALSE])))
keep <- rowSums(valid_per_group >= min_valid) >= 1

intensity_filtered <- intensity_na[keep, ]
cat("Protein groups after valid-value filter:", nrow(intensity_filtered), "\n")

saveRDS(list(annot              = annot,
             design             = design,
             intensity_all      = intensity_na,
             intensity_filtered = intensity_filtered),
        file.path(dir_processed, "prepared.rds"))
