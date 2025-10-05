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

# Perform normalization
loess_normalized_data <- normalizeCyclicLoess(as.matrix((dt_intensity.imputed)))  
head(loess_normalized_data[1:10,1:10])
vsn_normalized_data <- performVSNNormalization(as.matrix(dt_intensity.imputed))
head(vsn_normalized_data)
vsn_normalized_data <- as.matrix(as.data.frame(vsn_normalized_data))
row.names(vsn_normalized_data) <- dt_info$Protein.IDs

#saveRDS(vsn_normalized_data, "randomForest_imputed_intensity_normalised.rds")

################Stat-test#########################
    intensity_data <- dt_intensity.imputed
    head(intensity_data)
    
    design <- data.frame(
        sample = gsub("-Intensity", "", colnames(intensity_data)),
        group = rep(c("AB", "Ar1AB", "Ar10AB", "control"), each = 4),
        replicates = rep(1:4, 4)
    )
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
    
    
    
    design_matrix
    
    str(fit_contrasts$contrasts)
    
   
    #check foldchanges VSN-normalized data is already log2-scaled
        top_genes <- topTable(fit_ebayes, coef = "Ar1AB_vs_control", number = Inf)
        head(top_genes)
        coef_dt_long %>% dplyr::filter(Protein == "Q9Y2Y6")
        logFC <- top_genes$logFC
        vsn_normalized_data[rownames(vsn_normalized_data) %in% "Q9Y2Y6",grepl("Ar1-|SF", colnames(vsn_normalized_data))]
        trt <- (vsn_normalized_data[rownames(vsn_normalized_data) %in% "Q9Y2Y6",grepl("Ar1-", colnames(vsn_normalized_data))])
        cont <- (vsn_normalized_data[rownames(vsn_normalized_data) %in% "Q9Y2Y6",grepl("SF", colnames(vsn_normalized_data))])
        mean(log2(trt) - log2(cont))
        mean((trt) - (cont))
    
    sum(pval_dt_long$padj < 0.05) #non-adjust work
    
    # Filter for significant hits (FDR < 0.05)
    sig_results <- pval_dt_long[pval_dt_long$P_value < 0.05,]
    dim(sig_results)

################Volcano###########################
    allcombinations <- attributes(contrasts)$dimnames$Contrasts

    sig_tb <- list()
    uplist <- list()
    downlist <- list()
    updowntab <- c()
    vol_list <- list()
    for(i in 1:length(allcombinations)){
        print(i)
        logFC_thr <- 2
        SoftThreIN <- 0.05#/nrow(vol_in)
        
        pre_vol <- topTable(fit_ebayes, coef = allcombinations[i], number = Inf)
        pre_vol$uniprot <- gsub(";.*", "", rownames(pre_vol))
        vol_in <- merge(pre_vol, uniprot_out, by.x = "uniprot", by.y = "From", all.x = T, all.y = F)
        head(pre_vol)
        
        # add a column of NAs
        vol_in$diffexpressed <- "NO"
        # if log2logFC > 0.6 and pvalue < 0.05, set as "UP" 
        vol_in$diffexpressed[vol_in$logFC < -logFC_thr & vol_in$P.Value < SoftThreIN] <- "DOWN"
        # if log2logFC < -0.6 and pvalue < 0.05, set as "DOWN"
        vol_in$diffexpressed[vol_in$logFC > logFC_thr & vol_in$P.Value < SoftThreIN] <- "UP"
        #label gene of interest
        #add up down label
        #only few if more than 100 and show label only top 1%
        
        
        vol_in$datlabel <- ifelse(vol_in$diffexpressed %in% c("UP", "DOWN", "GOI") &
                                      (vol_in$logFC > -logFC_thr | vol_in$logFC < logFC_thr) &
                                      vol_in$P.Value < SoftThreIN, gsub(";.*", "", vol_in$Gene.Names), NA)  
        #}
        
        #vol_in$diffexpressed[(vol_in$Gene.Names %in% compileGOI$x)] <- "GOI"
        table(vol_in$diffexpressed)
        
        mycolors <- c("#77DD79","#A7D2E9",  "grey")
        names(mycolors) <- c("DOWN", "UP", "NO")
        names(mycolors) <- factor(names(mycolors), levels= c("DOWN", "UP", "NO"))
        
        p <- ggplot(data=vol_in, aes(x=logFC, y=-log10(P.Value), label = datlabel, col=diffexpressed)) + 
            geom_point()+
            theme_classic(base_size = 16)+
            geom_text_repel(
                size=2.0, color = "black",
                fontface = 'bold', color = 'white',
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            ) + geom_vline(xintercept=(c(-logFC_thr, logFC_thr)), col="gold") +
            geom_hline(yintercept=-log10(SoftThreIN), col="gold") + 
            theme(plot.title = element_text(hjust = 0.5)) + 
            ggtitle(allcombinations[i]) + scale_colour_manual(values = mycolors, 
                                                     labels = c("Down-regulated","Not DEGs","Up-regulated",   "Genes of Interest"),
                                                     name = "") +
            xlab("Log2 fold change") + ylab("-Log10 of odds")
        #p
        #ggsave(paste0("aurantiamide_plots/volcano/", allcombinations[i], ".png"),
        #       width = 12, height = 8, dpi = 300)
        
        ####summary table
        vol_sig <- union(
            vol_in[vol_in$logFC < -logFC_thr & vol_in$P.Value < SoftThreIN,]$Gene.Names,
            vol_in[vol_in$logFC > logFC_thr & vol_in$P.Value < SoftThreIN,]$Gene.Names
        )
        sig_tb[[i]] <- vol_sig
        names(sig_tb)[i] <- allcombinations[i]
        
        upeach <- gsub(";.*", "", rownames(vol_in[vol_in$logFC > logFC_thr & vol_in$P.Value < SoftThreIN,]))
        upeach <- upeach[!grepl("^NA.", upeach)]
        downeach <- gsub(";.*", "", rownames(vol_in[vol_in$logFC < -logFC_thr & vol_in$P.Value < SoftThreIN,]))
        downeach <- downeach[!grepl("^NA.", downeach)]
        
        
        uptab <- vol_in[vol_in$logFC > logFC_thr & vol_in$P.Value < SoftThreIN,] #%>% left_join(allprottab, by = "Gene.Names")
        downtab <- vol_in[vol_in$logFC < -logFC_thr & vol_in$P.Value < SoftThreIN,] #%>% left_join(allprottab, by = "Gene.Names")
        
        uptab <- cbind.data.frame(allcombinations[i], expression = paste0("up-regulated", " (n = ", nrow(uptab %>% dplyr::filter(!is.na(P.Value))), ")"), 
                                  uptab[c("uniprot", "Protein.names", "logFC", "P.Value")]) %>% arrange(P.Value, logFC) %>% 
            dplyr::filter(!is.na(P.Value))
        downtab <- cbind.data.frame(allcombinations[i], expression = paste0("down-regulated",  " (n = ", nrow(downtab %>% dplyr::filter(!is.na(P.Value))), ")")
                                    , downtab[c("uniprot", "Protein.names", "logFC", "P.Value")]) %>% arrange(P.Value, logFC) %>% 
            dplyr::filter(!is.na(P.Value))
        
        updowntab <- rbind.data.frame(updowntab, rbind.data.frame(uptab, downtab))
        
        upeach <- gsub(";.*", "", uptab$Protein.IDs)
        downeach <- gsub(";.*", "", downtab$Protein.IDs)
        
        uplist[[i]] <- upeach
        names(uplist)[i] <- paste0("up_", allcombinations[i])
        downlist[[i]] <- downeach
        names(downlist)[i] <- paste0("down_", allcombinations[i])
        
        vol_list[[i]] <- vol_in
        names(vol_list)[i] <- paste0("", allcombinations[i])
        
    }

saveRDS(vol_list, "aurantiamide_plots/volcano/vol_list.rds")


################KEGG###########################
    #convert Uniprot to Enzres ID
    hs <- org.Hs.eg.db
    toent <- function(IN){
        #IN <- vol_list$Ar1AB_vs_control$uniprot
        OUT <- AnnotationDbi::select(hs, 
                                     keys = IN,
                                     columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "UNIPROT"),
                                     keytype = "UNIPROT")
        return(OUT)  
    }
    
    allprot <- vol_list #list(uplist, downlist)
    
    ENTID <- list()
    for(k in 1:length(allprot)){
        #k = 1
        logFC_thr <- 2
        SoftThreIN <- 0.05#/nrow(vol_in)
        LISTIN <- allprot[[k]]
        LISTIN <- LISTIN[LISTIN$logFC < -logFC_thr & LISTIN$P.Value < SoftThreIN,]
        SET <- names(allprot)[k]
        ent_num <- toent(LISTIN$uniprot) #toent(LISTIN[[i]])$ENTREZID
        ENTID[[k]] <- ent_num
        names(ENTID)[k] <- SET
    }
    
    #ENTID <- vol_list
    GOtab <- list()
    EDO <- list()
    KEGG <- list()
    ENRICHR <- list()
    for(l in 1:length(ENTID)){
        pth <- 0.05
        qth <- 0.05
        #l <- 5
        #l <- 1
        print(l)
        names(ENTID)
        gene <- ENTID[[l]]
        head(gene)
        gene <- (gene$ENTREZID)
        #geneList <- na.omit(unique(unlist(ENTID[[l]])))#all genes
        geneList <- unique(unlist(ENTID))#all genes
        
        gene.df <- bitr(gene, fromType = "ENTREZID",
                        toType = c("ENSEMBL", "SYMBOL"),
                        OrgDb = org.Hs.eg.db)
        
        ggo <- groupGO(gene     = gene,
                       OrgDb    = org.Hs.eg.db,
                       ont      = "BP",
                       level    = 4,
                       readable = TRUE)
        
        ##GO over-representation
        pth <- 1
        qth <- 1
        ego <- enrichGO(gene          = gene.df$ENTREZID,#gene,
                        keyType = "ENTREZID",
                        #universe      = (geneList),
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "fdr",
                        pvalueCutoff  = pth,
                        qvalueCutoff  = qth,
                        readable      = TRUE)
        
        ego_df <- ego@result %>% 
            dplyr::filter(p.adjust < 0.05, qvalue < 1)
        
        edo <- pairwise_termsim(ego)
        
        
        NAMEIN <- names(ENTID)[l]
        NAMEIN <- gsub("\\|", "VS", NAMEIN)
        #ego
        
        print(ifelse(sum(ego@result$p.adjust < pth) == 0, paste0("FAIL_", names(ENTID)[l]), NA))
        
        #loop to increase p threshold
        #while(sum(ego@result$p.adjust < pth) == 0){
        #pth <- pth + 0.05
        #qth <- qth + 0.05
        #ego <- enrichGO(gene          = gene.df$ENTREZID,#gene,
        #                keyType = "ENTREZID",
        #                universe      = (geneList),
        #                OrgDb         = org.Hs.eg.db,
        #                ont           = "BP",
        #                pAdjustMethod = "fdr",
        #                pvalueCutoff  = pth,
        #                qvalueCutoff  = qth,
        #                readable      = TRUE)
        if(sum(ego@result$p.adjust < pth) == 0){
            
            kk <- enrichKEGG(gene         = gene,
                             organism     = 'hsa',
                             pvalueCutoff = pth)
            
        }else{
            NAMEIN <- paste(NAMEIN, "_pval=", pth, sep = "")
        }
        
        #barplot(ego, showCategory=8)
        #dotplot(ego)
        
        #if goplot error
        tryCatch(
            {
                #ego_filtered <- ego %>%
                #    dplyr::filter(ego$p.adjust < 0.05) %>%
                #    dplyr::slice_min(p.adjust, n = 20)
                
                #GOPLOT <- goplot(ego_filtered)
                
                GOPLOT <- clusterProfiler::goplot(ego) #go plot
                gg_go <- GOPLOT +  ggtitle(NAMEIN) + theme(plot.title = element_text(hjust = 0.5))
                ggplot2::ggsave(gg_go, filename = paste0("/cluster_profiler/GOplot/", NAMEIN, ".png"),
                                width = 12, height = 12, dpi = 300)
                
                
                CNET <- cnetplot(ego)
                gg_cnet <- CNET +  ggtitle(NAMEIN) + theme(plot.title = element_text(hjust = 0.5))
                ggplot2::ggsave(gg_cnet, filename = paste0("/cluster_profiler/cnet/", NAMEIN, ".png"),
                                width = 12, height = 12, dpi = 300)
                #enrichplot::emapplot(ego, showCategory=30) 
                
                
                ENRMAP <- emapplot(edo) #enrichmap
                ENRMAP +  ggtitle(NAMEIN) + theme(plot.title = element_text(hjust = 0.5))
                ggplot2::ggsave(paste0("/cluster_profiler/enrichment_map/", NAMEIN, ".png"),
                                width = 12, height = 12, dpi = 300) 
                
                dotplot(kk) +  ggtitle(NAMEIN) + theme(plot.title = element_text(hjust = 0.5)) + 
                    scale_y_discrete(name = "Pathway")
                ggplot2::ggsave(paste0("/cluster_profiler/keggpathway/", NAMEIN, ".png"),
                                width = 12, height = 12, dpi = 300) 
                
            }, error=function(e){}
        )
        
        #save data frame
        topcutoff <- 20
        if(pth == 0.05){
            eachGOtab <- ego@result[ego@result$p.adjust < pth,]
        }else{
            eachGOtab <- ego@result[1:topcutoff,]
        }
        GOtab[[l]] <- eachGOtab
        names(GOtab)[l] <- NAMEIN
        
        if(pth == 0.05){
            eachEDOtab <- edo@result[edo@result$p.adjust < pth,]
        }else{
            eachEDOtab <- edo@result[1:topcutoff,]
        }
        EDO[[l]] <- eachEDOtab 
        names(EDO)[l] <- NAMEIN
        
        
        KEGG[[l]] <- kk@result
        names(KEGG)[l] <- NAMEIN
        
        
        #Enrichr
        NAMEIN <- names(ENTID)[l]
        dbs <- listEnrichrDbs()
        if (is.null(dbs)) websiteLive <- FALSE
        if (websiteLive) head(dbs)
        dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015", "Human_Phenotype_Ontology", "KEGG_2021_Human")
        if (websiteLive) {
            enriched <- enrichr(gene.df$SYMBOL, dbs)
        }
        
        if (websiteLive) plotEnrich(enriched[[4]], showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title = paste0(NAMEIN, "_Human_Phenotype_Ontology"))
        ggplot2::ggsave(paste0("/cluster_profiler/enrichr/", NAMEIN, "_Human_Phenotype_Ontology", ".png"),
                        width = 12, height = 6, dpi = 300)
        
        if (websiteLive) plotEnrich((enriched[[3]]), showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title = paste0(NAMEIN, "_Biological_Process"))
        ggplot2::ggsave(paste0("/cluster_profiler/enrichr/", NAMEIN, "_Biological_Process", ".png"),
                        width = 12, height = 6, dpi = 300)
        
        
        
        if(!all(enriched[[5]]$Adjusted.P.value < 0.05)){
            enriched[[5]]$Adjusted.P.value = 0.03  
            if (websiteLive) plotEnrich(enriched[[5]], , showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title = paste0(NAMEIN, "_KEGG_all"))
            ggplot2::ggsave(paste0("/cluster_profiler/enrichr/", NAMEIN, "_KEGG_all", ".png"),
                            width = 12, height = 6, dpi = 300)
        }else{
            if (websiteLive) plotEnrich(enriched[[5]], , showTerms = 20, numChar = 100, y = "Count", orderBy = "P.value", title = paste0(NAMEIN, "_KEGG_Biological_Process"))
            ggplot2::ggsave(paste0("/cluster_profiler/enrichr/", NAMEIN, "_KEGG_Biological_Process", ".png"),
                            width = 12, height = 6, dpi = 300)
            
        }
        
        ENRICHR[[l]] <- enriched
        names(ENRICHR)[l] <- NAMEIN
        
    }
    
    #GOtab <- readRDS("/cluster_profiler/GOENRICH.rds")
    #saveRDS(GOtab, "/cluster_profiler/GOENRICH2.rds")
    #saveRDS(EDO, "/cluster_profiler/GOpairwise.rds")
    
    saveRDS(KEGG, "/cluster_profiler/KEGG.rds")
    saveRDS(ENRICHR, "/cluster_profiler/ENRICHR.rds")
    #openxlsx::write.xlsx(GOtab, file = "/cluster_profiler/GOENRICHpro.xlsx")
    write.xlsx(GOtab, file = "/cluster_profiler/GOENRICHpro.xlsx")
    write.xlsx(EDO, file = "/cluster_profiler/GOpairwise.xlsx")
    
