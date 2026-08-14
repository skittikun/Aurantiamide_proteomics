library(GO.db)
library(reshape2)
library(dplyr)
library(ggrepel)
library('heatmaply')
#library(xlsx)
library(msqrob2)
library(limma)
library(org.Hs.eg.db) 
library(enrichplot)
#library(EnsDb.Hsapiens.v75)

library(clusterProfiler)
library(UpSetR)
library(topGO)
library(org.Hs.eg.db)
library(limma)
library(gridExtra)
library(emmeans)
library(ggplot2)
library(dplyr)
library(multcompView)
library(NormalyzerDE)
library(scales)
library(factoextra)
library(tidyr)
library(ggrepel)
library(ggplot2)
library(tidyr)
library(SWATH2stats)
library(dplyr)
library(iq)
library(openxlsx)
library(limma)
library(missForest)
library(limma)      # For cyclic Loess
library(vsn)        # For VSN normalization
library(preprocessCore) # For quantile normalization
library(NormalyzerDE)
library(enrichR)

listEnrichrSites()
setEnrichrSite("Enrichr") 
websiteLive <- TRUE
dbs <- listEnrichrDbs()

dtIN <- read.xlsx("proteinGroups_intensity_p-Value-posthoc-annotation.xlsx")
dt_info <- dtIN[, !grepl("Intensity", colnames(dtIN))]

dt_update <- readRDS("randomForest_imputed_intensity_dt.rds")
dt_intensity.imputed <- readRDS("randomForest_imputed_intensity.rds")
uniprot_out <- read.csv("Uniprot_retrieve_idmapping_2025_02_08.tsv", sep = "\t")
vsn_normalized_data <- readRDS("randomForest_imputed_intensity_normalised.rds")

#pca on raw reads
    str(dt_update)
    str(dt_intensity.imputed)
    str(vsn_normalized_data)
    
    ## =============================================================
    ## PCA for three proteomics matrices: raw, imputed, VSN
    ## Samples: AB, Ar1-AB, Ar10-AB, SF  x  N1-N4  (16 samples)
    ## =============================================================
    
    library(ggplot2)
    library(patchwork)   # optional, for combining plots
    
    ## ---------- 1. Extract the three matrices --------------------
    
    intensity_cols <- grep("[-_]Intensity$", colnames(dt_update), value = TRUE)
    stopifnot(length(intensity_cols) == 16)
    
    mat_raw     <- as.matrix(dt_update[, intensity_cols])
    rownames(mat_raw) <- dt_update$Protein.IDs
    
    mat_imputed <- as.matrix(dt_intensity.imputed)
    rownames(mat_imputed) <- dt_update$Protein.IDs
    
    mat_vsn     <- as.matrix(vsn_normalized_data)
    
    ## ---------- 2. Sample metadata from column names -------------
    
    make_meta <- function(cn) {
      clean <- sub("[-_]Intensity$", "", cn)
      group <- sub("-N[0-9]+$", "", clean)
      rep   <- regmatches(clean, regexpr("N[0-9]+$", clean))
      data.frame(
        sample    = cn,
        label     = clean,
        group     = factor(group, levels = c("AB", "Ar1-AB", "Ar10-AB", "SF")),
        replicate = rep,
        stringsAsFactors = FALSE
      )
    }
    
    meta <- make_meta(colnames(mat_vsn))
    print(meta)
    
    ## ---------- 3. PCA function ----------------------------------
    
    run_pca <- function(mat, meta, do_log = TRUE, label = "dataset",
                        min_var = 1e-8) {
      
      m <- mat[, meta$sample, drop = FALSE]   # enforce identical column order
      
      if (do_log) {
        m[m <= 0] <- NA                       # zeros are missing, not real signal
        m <- log2(m)
      }
      
      keep_complete <- stats::complete.cases(m)
      m <- m[keep_complete, , drop = FALSE]
      
      vars <- apply(m, 1, stats::var)
      m <- m[is.finite(vars) & vars > min_var, , drop = FALSE]
      
      message(sprintf("%s: %d of %d proteins retained (%.1f%%)",
                      label, nrow(m), nrow(mat), 100 * nrow(m) / nrow(mat)))
      
      pca <- stats::prcomp(t(m), center = TRUE, scale. = TRUE)
      
      pv <- 100 * pca$sdev^2 / sum(pca$sdev^2)
      
      scores <- data.frame(pca$x[, 1:min(4, ncol(pca$x))], meta,
                           dataset = label, row.names = NULL)
      
      scree <- data.frame(
        PC      = factor(seq_along(pv), levels = seq_along(pv)),
        percent = pv,
        cumul   = cumsum(pv),
        dataset = label
      )
      
      list(pca = pca, scores = scores, scree = scree,
           pv = pv, n_proteins = nrow(m), label = label)
    }
    
    ## ---------- 4. Run on all three ------------------------------
    
    datasets <- list(
      list(mat = mat_raw,     do_log = TRUE,  label = "Raw (log2)"),
      list(mat = mat_imputed, do_log = TRUE,  label = "Imputed (log2)"),
      list(mat = mat_vsn,     do_log = FALSE, label = "VSN normalised")
    )
    
    res <- lapply(datasets, function(d)
      run_pca(d$mat, meta, do_log = d$do_log, label = d$label))
    names(res) <- vapply(res, `[[`, character(1), "label")
    
    lapply(res, function(r) summary(r$pca)$importance[, 1:4])
    
    ## ---------- 5. Score plots -----------------------------------
    
    pal <- c("AB" = "#4C72B0", "Ar1-AB" = "#DD8452",
             "Ar10-AB" = "#55A868", "SF" = "#C44E52")
    
    plot_scores <- function(r, pcx = 1, pcy = 2, show_labels = TRUE) {
      df <- r$scores
      xv <- paste0("PC", pcx); yv <- paste0("PC", pcy)
      
      p <- ggplot(df, aes(.data[[xv]], .data[[yv]],
                          colour = group, shape = group)) +
        geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
        geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
        geom_point(size = 3.5, alpha = 0.9) +
        scale_colour_manual(values = pal) +
        labs(
          title    = r$label,
          subtitle = sprintf("%s proteins", format(r$n_proteins, big.mark = ",")),
          x = sprintf("PC%d (%.1f%%)", pcx, r$pv[pcx]),
          y = sprintf("PC%d (%.1f%%)", pcy, r$pv[pcy]),
          colour = "Group", shape = "Group"
        ) +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank())
      
      if (show_labels) {
        p <- p + geom_text(aes(label = replicate), size = 2.6,
                           vjust = -1.1, show.legend = FALSE)
      }
      p
    }
    
    score_plots <- lapply(res, plot_scores)
    
    # Side-by-side comparison of the three pipelines
    combined <- (score_plots[[1]] | score_plots[[2]] | score_plots[[3]]) +
      patchwork::plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = "PC1 vs PC2 across processing stages",
        theme = theme(plot.title = element_text(face = "bold"))
      )
    print(combined)
    
    ## ---------- 6. Scree plots -----------------------------------
    
    scree_all <- do.call(rbind, lapply(res, `[[`, "scree"))
    scree_all$dataset <- factor(scree_all$dataset, levels = names(res))
    
    p_scree <- ggplot(subset(scree_all, as.integer(PC) <= 10),
                      aes(PC, percent, fill = dataset)) +
      geom_col(position = position_dodge(0.8), width = 0.75) +
      labs(x = "Principal component", y = "Variance explained (%)",
           fill = NULL, title = "Scree comparison (first 10 PCs)") +
      theme_bw(base_size = 11) +
      theme(panel.grid.major.x = element_blank())
    print(p_scree)
    
    ## ---------- 7. PC3/PC4 to check for hidden structure ---------
    
    p34 <- lapply(res, plot_scores, pcx = 3, pcy = 4)
    print((p34[[1]] | p34[[2]] | p34[[3]]) +
            patchwork::plot_layout(guides = "collect"))
    
    ## ---------- 8. Outlier flagging (2 / 3 SD on PC1-PC2) --------
    
    flag_outliers <- function(r, nsd = 3) {
      s <- r$scores
      d <- sqrt(((s$PC1 - mean(s$PC1)) / stats::sd(s$PC1))^2 +
                  ((s$PC2 - mean(s$PC2)) / stats::sd(s$PC2))^2)
      data.frame(dataset = r$label, sample = s$label, group = s$group,
                 dist_sd = round(d, 2), outlier = d > nsd)
    }
    outliers <- do.call(rbind, lapply(res, flag_outliers))
    print(outliers[order(-outliers$dist_sd), ])
    
    ## ---------- 9. Optional: top loadings on PC1 (VSN) -----------
    
    top_loadings <- function(r, pc = 1, n = 20) {
      ld <- r$pca$rotation[, pc]
      ids <- names(sort(abs(ld), decreasing = TRUE))[1:n]
      idx <- match(ids, dt_update$Protein.IDs)
      data.frame(Protein.IDs = ids,
                 Gene.names  = dt_update$Gene.names[idx],
                 loading     = round(ld[ids], 4),
                 row.names   = NULL)
    }
    print(top_loadings(res[["VSN normalised"]], pc = 1, n = 20))
    
    ## ---------- 10. Save outputs ---------------------------------
    
    dir.create("output", showWarnings = FALSE)
    ggsave("pca_pc1_pc2_three_datasets.png", combined,
           width = 15, height = 5, dpi = 300)
    ggsave("pca_scree_comparison.png", p_scree,
           width = 9, height = 4.5, dpi = 300)
    write.csv(do.call(rbind, lapply(res, `[[`, "scores")),
              "pca_scores_all_datasets.csv", row.names = FALSE)
    write.csv(scree_all, "pca_variance_explained.csv", row.names = FALSE)
    
    
    
    