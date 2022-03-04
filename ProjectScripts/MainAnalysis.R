#BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")
packageF("pheatmap")
packageF("DESeq2")


AllData <- readRDS("Data/AllEstimatesAndMeta.Rds")
#IHC_Microglia <- read.table("Data/Microglia_IHC.txt", header = T, sep = "\t")

AllData %<>% mutate(across(contains("MGP"), rescale, c(0,1)))

#AllData <- merge(AllData, IHC_Microglia, by = "BioBankID", all.x = T, sort = F )

AllData_long <- pivot_longer(AllData, cols = matches("^T._"), names_sep = "_",
                             names_to = c("Tissue", "Matherial", "Method", "CellType"), values_to = "Estimate")



CellTypes <- unique(AllData_long$CellType)

RNA_MethodComparison <-   sapply(CellTypes, function(CellType){
  T1 = cor.test(formula(paste0("~T1_RNA_MGP_", CellType, "+T1_RNA_Bisque_", CellType)), data = AllData)
  T2 = cor.test(formula(paste0("~T2_RNA_MGP_", CellType, "+T2_RNA_Bisque_", CellType)), data = AllData)
  data.frame(CellType = CellType, T1 = T1$est, T2 = T2$est)
}, simplify = F) %>% rbindlist() %>% data.frame()




  
RNA_TissueComparison <- sapply(unique(AllData$Resequenced), function(TS){
  sapply(CellTypes, function(CellType){
    MGP = cor.test(formula(paste0("~T1_RNA_MGP_", CellType, "+T2_RNA_MGP_", CellType)), data = AllData %>% filter(Resequenced == TS))
    Bisque = cor.test(formula(paste0("~T1_RNA_Bisque_", CellType, "+T2_RNA_Bisque_", CellType)), data = AllData %>% filter(Resequenced == TS))
    data.frame(CellType = CellType, CorMGP = MGP$est, CorBisque = Bisque$est)
  }, simplify = F) %>% rbindlist() %>% data.frame() %>% mutate(SameTissue = TS)
}, simplify = F) %>% rbindlist() %>% data.frame()

RNA_TissueComparisonPlot <- sapply(c("MGP", "Bisque"), function(Method){
  PlotList <- sapply(CellTypes, function(CellType){
    ggplot(AllData, aes_string(paste0("T1_RNA_", Method, "_", CellType),
                               paste0("T2_RNA_", Method, "_", CellType),
                               color = "Resequenced")) +
      theme_bw() +
      labs(x = "Sample1", y = "Sample2", title = CellType) +
      #geom_abline(slope = 1, intercept = 0, lty = "dashed", color = "red") +
      geom_point() +
      scale_color_manual(values = c("grey50", "orange")) +
      stat_cor(method = "pearson", show.legend = F)
  }, simplify = F)
}, simplify = F)

ggarrange(ggarrange(plotlist = RNA_TissueComparisonPlot$MGP, nrow = 1, common.legend = T),
          ggarrange(plotlist = RNA_TissueComparisonPlot$Bisque, nrow = 1, common.legend = T),
          nrow = 2, labels = c("MGP", "Bisque"))

RNA_BS_Comparison <- sapply(c("Bisque", "MGP"), function(Method){
  sapply(CellTypes, function(CellType){
    Same_tissue = cor.test(formula(paste0("~T1_RNA_", Method, "_", CellType, "+T1_BS_MSP_", CellType)), data = AllData)
    Different_tissue = cor.test(formula(paste0("~T2_RNA_", Method, "_", CellType, "+ T1_BS_MSP_", CellType)), data = AllData %>% filter(Resequenced == "No"))
    data.frame(CellType = CellType, CorSameTissue = Same_tissue$est, CorDifferentTissue = Different_tissue$est)
  }, simplify = F) %>% rbindlist() %>% data.frame()
}, simplify = F)


RNA_H3K27_Comparison <- sapply(c("Bisque", "MGP"), function(Method){
  sapply(CellTypes, function(CellType){
    T1 = cor.test(formula(paste0("~T1_RNA_", Method, "_", CellType, " + T3_H3K27_MSP_", CellType)), data = AllData %>% filter(Resequenced == "No"))
    T2 = cor.test(formula(paste0("~T2_RNA_", Method, "_", CellType, " + T3_H3K27_MSP_", CellType)), data = AllData %>% filter(Resequenced == "No"))
    data.frame(CellType = CellType, CorT1 = T1$est, CorT2 = T2$est)
  }, simplify = F) %>% rbindlist() %>% data.frame()
}, simplify = F)

BS_H3K27_Comparison <- sapply(CellTypes, function(CellType){
  Cor = cor.test(formula(paste0("~T1_BS_MSP_", CellType, " + T3_H3K27_MSP_", CellType)), data = AllData %>% filter(Resequenced == "No"))
  data.frame(CellType = CellType, Cor = Cor$est)
}, simplify = F) %>% rbindlist() %>% data.frame()



AnnotationsDF <- data.frame(Name =  unique(sapply(grep("^T._", names(AllData), value = T), function(x){
  gsub(paste0(paste0("_", CellTypes), collapse = "|"), "", x)
  })))
AnnotationsDF$Tissue <- sapply(AnnotationsDF$Name, function(x){
  strsplit(x, "_")[[1]][1]
})

AnnotationsDF$Method <- sapply(AnnotationsDF$Name, function(x){
  strsplit(x, "_")[[1]][3]
})

AnnotationsDF$OmicsType <- sapply(AnnotationsDF$Name, function(x){
  strsplit(x, "_")[[1]][2]
})

AnnotationsDF$Fraction <- sapply(AnnotationsDF$OmicsType, function(x){
  if(x == "RNA"){
    "WholeTissue"
  } else {
    "Nuclear"
  }
})

rownames(AnnotationsDF) <- AnnotationsDF$Name
AnnotationsDF %<>% select(-Name) 

AnnoCol = list(Tissue = c("T1" = MoviePalettes$MoonRiseKingdomColors[2],
                          "T2" = MoviePalettes$MoonRiseKingdomColors[4],
                          "T3" = MoviePalettes$MoonRiseKingdomColors[10]),
               Method = c("CETS" = MoviePalettes$MoonRiseKingdomColors[5],
                          "MSP" = MoviePalettes$MoonRiseKingdomColors[1],
                          "MGP" = MoviePalettes$MoonRiseKingdomColors[4],
                          "Bisque" = MoviePalettes$MoonRiseKingdomColors[3]),
               Fraction = c("Nuclear" = MoviePalettes$BugsLife[2],
                            "WholeTissue" = MoviePalettes$BugsLife[4]),
               OmicsType = c("BS" = MoviePalettes$LittleShopOfHorrors[1],
                             "RNA" = MoviePalettes$LittleShopOfHorrors[5],
                             "H3K27" = MoviePalettes$LittleShopOfHorrors[9]))

#Look at the similarity
similarityPlotList <- sapply(CellTypes, function(CellType){
  Cor <- cor(AllData %>% filter(Resequenced == "No") %>%  select(matches(CellType)))
  rownames(Cor) <- sapply(rownames(Cor), function(x){
    gsub(paste0("_", CellType), "", x)
  })
  colnames(Cor) <- sapply(rownames(Cor), function(x){
    gsub(paste0("_", CellType), "", x)
  })
  diag(Cor) <- NA
  pheatmap(Cor, border_color = NA, na_col = "white", cluster_rows = T, cluster_cols = T,fontsize_number = 12, angle_col = 90,
           annotation_col = AnnotationsDF, annotation_colors = AnnoCol, cutree_rows = 2, cutree_cols = 2,
           scale = "none", method = "ward.D2", display_numbers = T, breaks=seq(-1, 1, length.out=100), main = CellType)
})




pheatmap(Cor, border_color = NA, na_col = "white", cluster_rows = T, cluster_cols = T,fontsize_number = 12,
         scale = "none", method = "ward.D2", display_numbers = T, breaks=seq(-1, 1, length.out=100), main = CellType)