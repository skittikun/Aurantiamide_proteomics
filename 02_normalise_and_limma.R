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


# Example dataset with missing values
data(iris)
set.seed(111)
iris.na <- iris
for (i in 1:4) iris.na[sample(150, sample(20, 1)), i] <- NA

# Impute missing values using missForest
iris.imputed <- missForest(iris.na)$ximp

# View the imputed dataset
head(iris.imputed)

dtIN <- read.xlsx("proteinGroups_intensity_p-Value-posthoc-annotation.xlsx")

dt_info <- dtIN[, !grepl("Intensity", colnames(dtIN))]
dt_intensity.imputed <- readRDS("randomForest_imputed_intensity.rds")

head(dt_intensity.imputed)
head(dt_intensity[1:100,])

dt_update <- cbind(dt_info, dt_intensity.imputed)
head(dt_update)
dim(dt_update)

#saveRDS(dt_update, "randomForest_imputed_intensity_dt.rds")

dt_update <- readRDS("randomForest_imputed_intensity_dt.rds")

#to retrieve from Uniprot
  #write.table(unique(dt_update$Protein.IDs), "to_uniprot.txt", quote = F, row.names = F, col.names = F)
  uniprot_out <- read.csv("Uniprot_retrieve_idmapping_2025_02_08.tsv", sep = "\t")
#validate normalisation method
  library(limma)      # For cyclic Loess
  library(vsn)        # For VSN normalization
  library(preprocessCore) # For quantile normalization
  library(NormalyzerDE)   # For normalization evaluation
  
  intensity_data <- dt_intensity.imputed
  head(intensity_data)
  
  design <- data.frame(
    sample = gsub("-Intensity", "", colnames(intensity_data)),
    group = rep(c("AB", "Ar1AB", "Ar10AB", "control"), each = 4),
    replicates = rep(1:4, 4)
  )
  # 3. Apply different normalization methods

  # 4. Evaluate with NormalyzerDE
  setwd("")
  job_dir <- "normalization_evaluation"
  dir.create(job_dir)
  
  intensity_matrix <- as.matrix(intensity_data)
  colnames(intensity_matrix) <- gsub("\\-Intensity", "", colnames(intensity_matrix))

  write.table(design, file = "design_matrix.tsv", sep = "\t", row.names = FALSE, quote = F)
  write.table(intensity_matrix, file = "intensity_data.tsv", sep = "\t", quote = F, row.names = F)
  
  design$sample == colnames(intensity_matrix)
  
  # 2. Define paths
  design_path <- "design_matrix.tsv"
  data_path <- "intensity_data.tsv"
  job_dir <- "DN_analysis"
  
  # 3. Run NormalyzerDE pipeline
  norm_obj <- NormalyzerDE::normalyzer(
    jobName = "DN_Loess_VSN_Comparison",
    designPath = design_path,  # Path to design file
    dataPath = data_path,       # Path to data file
    outputDir = job_dir,
    noLogTransform = FALSE,     # Let NormalyzerDE handle log2 (matches paper's MaxQuant workflow[1])
    sampleColName = "sample",
    groupColName = "group"
  )
  
  
#decided to use VSN normalisation
# Perform normalization
  loess_normalized_data <- normalizeCyclicLoess(as.matrix((dt_intensity.imputed)))  
  head(loess_normalized_data[1:10,1:10])
  vsn_normalized_data <- performVSNNormalization(as.matrix(dt_intensity.imputed))
  head(vsn_normalized_data)
  vsn_normalized_data <- as.matrix(as.data.frame(vsn_normalized_data))
  row.names(vsn_normalized_data) <- dt_info$Protein.IDs
  
#perform statistical test
  groups <- factor(design$group, levels = c("AB", "Ar1AB", "Ar10AB", "control"))
  design_matrix <- model.matrix(~0 + groups)
  colnames(design_matrix) <- levels(groups)
  
  # Fit linear model
  fit <- lmFit(vsn_normalized_data, design_matrix)
  contrasts <- makeContrasts(
    # Treatment vs Control
    Ar1AB_vs_control = Ar1AB - control,
    Ar10AB_vs_control = Ar10AB - control,
    AB_vs_control = AB - control,
    
    # Treatment vs Treatment
    Ar10AB_vs_Ar1AB = Ar10AB - Ar1AB,
    Ar1AB_vs_AB = Ar1AB - AB,
    Ar10AB_vs_AB = Ar10AB - AB,
    
    levels = design_matrix
  )

  fit_contrasts <- contrasts.fit(fit, contrasts)
  fit_ebayes <- eBayes(fit_contrasts)  
  
  pval_dt <- as.data.frame(fit_ebayes$p.value)
  pval_dt_long <- pval_dt %>% 
    tibble::rownames_to_column("Protein") %>% 
    pivot_longer(
      cols = -Protein,
      names_to = "Comparison",
      values_to = "P_value"
    )
  
  coef_dt <- as.data.frame(fit_ebayes$coefficients)
  coef_dt_long <- coef_dt %>% 
    tibble::rownames_to_column("Protein") %>% 
    pivot_longer(
      cols = -Protein,
      names_to = "Comparison",
      values_to = "coefficient"
    )
  
  pval_dt_long$padj <- p.adjust(pval_dt_long$P_value, method = "fdr")
  
  sum(pval_dt_long$padj < 0.05) #non-adjust work
  
  # Filter for significant hits (FDR < 0.05)
  sig_results <- pval_dt_long[pval_dt_long$P_value < 0.05,]
  dim(sig_results)
  
  sig_results <- sig_results %>% left_join(coef_dt_long, by = c("Protein", "Comparison")) %>%
    dplyr::mutate(
      log2FC = round(coefficient, 2),  # Keep original signed values
      fold_change = round(ifelse(coefficient > 0, 1, -1)*2^abs(coefficient), 2),  # Magnitude only
      direction = ifelse(coefficient > 0, "up", "down"),  # Add direction
      tomatch = gsub("\\;.*", "", Protein)
    )
  sig_results_uniprot <- merge(sig_results, uniprot_out, by.x = "tomatch", by.y = "From", all.x = T, all.y = F)
  write.xlsx(sig_results_uniprot, "significant_imputed_limma_uniprot.xlsx")
  
  #with data
  sig_results_summary <- sig_results_uniprot %>%
    dplyr::group_by(Protein) %>%
    dplyr::summarise(
      significant_comparison = paste(Comparison, collapse = ", "),
      log2FC = paste(log2FC, collapse = ", "),  # Keep original signed values
      fold_change = paste(fold_change, collapse = ", "),  # Magnitude only
      direction = paste(direction, collapse = ", "),  # Add direction
      Protein.names = paste(Protein.names, collapse = ", ")
    )
  
  write.xlsx(sig_results_summary, "significant_imputed_limma_summary.xlsx")
  
  #
  head(sig_results_uniprot)
  comparison_list <- setNames(sig_results_uniprot$direction, 
                              sig_results_uniprot$Comparison)
  
  OUTLIST <- list()
  summary_dt <- sig_results_uniprot %>% dplyr::group_by(Comparison, direction) %>% dplyr::summarise(n = n())
  OUTLIST[["summary"]] <- summary_dt
  for(compareIN in unique(sig_results_uniprot$Comparison)){
    #compareIN = "Ar10AB_vs_Ar1AB"
    for(UPDO in unique(sig_results_uniprot$direction)){
      #UPDO = "up"
      OUTLIST[[paste(compareIN, UPDO, sep = "_")]] <- sig_results_uniprot %>% dplyr::filter(Comparison %in% compareIN, direction %in% UPDO)
    }
  }
  write.xlsx(OUTLIST, "significant_imputed_limma_uniprot_subset.xlsx")
  
 #upset plot of protein
  OUTLIST
  str(OUTLIST)
    str(sig_results_uniprot)
    sig_results_uniprot$group <- paste(sig_results_uniprot$Comparison, sig_results_uniprot$direction, sep = "_")
    
    # Generate binary matrix
    binary_group_matrix <- as.data.frame.matrix(
      table(
        unique(sig_results_uniprot[c("Protein", "group")])
      )
    )  
    png("upsetplots/group.png",
        width = 10, height = 6, units = "in", res = 400)
    upset(
      binary_group_matrix,
      nsets = ncol(binary_matrix),  # Show all sets
      order.by = "freq",           # Sort by intersection size
      text.scale = 1.2,
      mainbar.y.label = "Protein Intersections"
    )
    dev.off()
    
    #comparison
    binary_comparison_matrix <- as.data.frame.matrix(
      table(
        unique(sig_results_uniprot[c("Protein", "Comparison")])
      )
    )  
    png("upsetplots/comparison.png",
        width = 10, height = 6, units = "in", res = 400)
    upset(
      binary_comparison_matrix,
      nsets = ncol(binary_matrix),  # Show all sets
      order.by = "freq",           # Sort by intersection size
      text.scale = 1.2,
      mainbar.y.label = "Protein Intersections"
    )
    dev.off()
    
    #get table of intersection
    combinations <- unlist(lapply(1:ncol(binary_comparison_matrix), 
                                  function(x) combn(colnames(binary_comparison_matrix), x, simplify = FALSE)))
    
    # Extract proteins for each intersection
    intersection_list <- lapply(combinations, function(cols) {
      rownames(binary_comparison_matrix)[rowSums(binary_comparison_matrix[, cols, drop = FALSE]) == length(cols)]
    })
    
    # Format as dataframe
    intersection_df <- data.frame(
      Intersection = sapply(combinations, paste, collapse = "&"),
      Proteins = sapply(intersection_list, paste, collapse = ";")
    )
    head(intersection_df)
    
##go enrichment
  GO_list <- list()
  for(compareIN in unique(sig_results_uniprot$Comparison)){
    #compareIN = "Ar10AB_vs_Ar1AB"
    for(UPDO in unique(sig_results_uniprot$direction)){
      print(paste(compareIN, UPDO, sep = "_"))
        toenrich <- OUTLIST[[paste(compareIN, UPDO, sep = "_")]]
        str(  toenrich)
        
        uniprot_ids <- uniprot_out$From #toenrich$Protein
        
        # Map UniProt to Entrez IDs
        entrez_ids <- mapIds(org.Hs.eg.db,
                             keys = toenrich$Protein,
                             column = "ENTREZID",
                             keytype = "UNIPROT",
                             multiVals = "first")
        
        # 2. Create gene scores using ACTUAL ENTREZ IDS
        gene_scores <- setNames(rep(1, length(entrez_ids)), entrez_ids)
        matched_genes <- intersect(names(gene_scores), entrez_ids)
        gene_scores[matched_genes] <- toenrich$P_value
        
        # 3. Define significance properly
        topDiffGenes <- function(scores) {
          return(scores < 0.05) # Using raw p-value threshold
        }
        
        # 4. Create GO data with verified identifiers
        go_data <- new("topGOdata",
                       ontology = "BP",
                       allGenes = gene_scores,
                       geneSel = topDiffGenes,
                       nodeSize = 10,
                       annot = annFUN.org,
                       mapping = "org.Hs.eg.db",
                       ID = "entrez")
        
        result_fisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
        
        # KS test with elim algorithm to reduce redundancy
        result_ks <- runTest(go_data, algorithm = "elim", statistic = "ks")
      
        top_terms <- GenTable(go_data,
                              classicFisher = result_fisher,
                              elimKS = result_ks,
                              orderBy = "elimKS",
                              ranksOf = "classicFisher",
                              topNodes = 20,
                              numChar = 120)
        
        # Filter significant terms (e.g., elimKS < 0.05)
        significant_terms <- subset(top_terms, elimKS < 0.05)  
        
        GO_list[[paste(compareIN, UPDO, sep = "_")]] <- significant_terms
      }
  }
  GO_list
  write.xlsx(GO_list, "GO_bioprocess.xlsx")