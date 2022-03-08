#BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")
packageF("pheatmap")
packageF("DESeq2")

ResultsPath = "ResultsEstimateComparison"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

AllData <- readRDS("Data/AllEstimatesAndMeta.Rds")
#IHC_Microglia <- read.table("Data/Microglia_IHC.txt", header = T, sep = "\t")

AllData %<>% mutate(across(contains("MGP"), rescale, c(0,1)))

#AllData <- merge(AllData, IHC_Microglia, by = "BioBankID", all.x = T, sort = F )

AllData_long <- pivot_longer(AllData, cols = matches("^T._"), names_sep = "_",
                             names_to = c("Tissue", "Matherial", "Method", "CellType"), values_to = "Estimate")



CellTypes <- unique(AllData_long$CellType)

RNA_MethodComparison <-   sapply(CellTypes, function(CellType){
  T1 = cor.test(formula(paste0("~T1_RNA_MGP_", CellType, "+T1_RNA_Bisque_", CellType)), data = AllData)
  T2 = cor.test(formula(paste0("~T2_RNA_MGP_", CellType, "+T2_RNA_Bisque_", CellType)), data = AllData %>% filter(Resequenced == "No"))
  data.frame(CellType = CellType, TS1 = T1$est, TS2 = T2$est)
}, simplify = F) %>% rbindlist() %>% data.frame() %>% pivot_longer(-CellType,
                                                                   names_to = "TissueSample",
                                                                   values_to = "Cor") %>% data.frame()
cor.test(AllData$T1_BS_MSP_Neurons, AllData$T1_BS_CETS_Neurons)



P1 <- ggplot(RNA_MethodComparison, aes(TissueSample, CellType, fill = Cor)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "") +
  geom_tile() +
  geom_text(aes(label = round(Cor, 2)), color = "white") +
  scale_fill_gradient(low = MoviePalettes$GreenPallete[6], high = MoviePalettes$GreenPallete[3]) +
  facet_wrap(~TissueSample, scales = "free_x")


P2 <- ggplot(RNA_MethodComparison, aes(Cor)) +
  theme_minimal() +
  labs(x = "MGP-Bisque correlation") +
  geom_density()

ggarrange(P1, P2, heights = c(1.3, 1), nrow = 2)
ggsave(paste0(ResultsPath, "MethodcomparisonRNA.pdf"), device = "pdf",
       width = 10, height = 6, dpi = 300, useDingbats = F)


RNA_TissueComparisonPlot <- sapply(c("MGP", "Bisque"), function(Method){
  PlotList <- sapply(CellTypes, function(CellType){
    ggplot(AllData, aes_string(paste0("T1_RNA_", Method, "_", CellType),
                               paste0("T2_RNA_", Method, "_", CellType),
                               color = "Resequenced")) +
      theme_bw() +
      labs(x = "Sample1", y = "Sample2", title = CellType) +
      geom_point() +
      scale_color_manual(values = MoviePalettes$BugsLife[c(2,4)]) +
      stat_cor(method = "pearson", show.legend = F)
  }, simplify = F)
}, simplify = F)

ggarrange(ggarrange(plotlist = RNA_TissueComparisonPlot$MGP, nrow = 1, common.legend = T),
          ggarrange(plotlist = RNA_TissueComparisonPlot$Bisque, nrow = 1, common.legend = T),
          nrow = 2, labels = c("MGP", "Bisque"))



ggarrange(P3, P4, nrow = 2)

ggsave(paste0(ResultsPath, "TScomparisonRNA.pdf"), device = "pdf",
       width = 14, height = 6, dpi = 300, useDingbats = F)


RNA_TissueComparisonCorr <- sapply(c("MGP", "Bisque"), function(Method){
  PlotList <- sapply(CellTypes, function(CellType){
    sapply(unique(AllData$Resequenced), function(SampleType){
      Cor <- cor.test(formula(paste0("~T1_RNA_", Method, "_", CellType,
                                     " + T2_RNA_", Method, "_", CellType)),
                      data = AllData %>% filter(Resequenced == SampleType))
      TS <- if(SampleType == "Yes"){
        "Same sample"
      } else {
        "Different sample"
      }
      data.frame(Method = Method, SampleType = TS, CellType = CellType, Corr = Cor$estimate)
    }, simplify = F) %>% rbindlist()
  }, simplify = F) %>% rbindlist()
}, simplify = F) %>% rbindlist() %>% data.frame()

RNA_TissueComparisonCorr$Method <- factor(RNA_TissueComparisonCorr$Method,
                                          levels = unique(RNA_TissueComparisonCorr$Method))

RNA_TissueComparisonCorr$SampleType <- factor(RNA_TissueComparisonCorr$SampleType,
                                              levels = c("Same sample", "Different sample"))

ggplot(RNA_TissueComparisonCorr, aes(SampleType, CellType, fill = Corr)) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "") +
  geom_tile() +
  geom_text(aes(label = round(Corr, 2)), color = "white") +
  scale_fill_gradient(low = MoviePalettes$GreenPallete[6], high = MoviePalettes$GreenPallete[3]) +
  facet_wrap(~Method)



ggsave(paste0(ResultsPath, "TScomparisonRNA_CorTable.pdf"), device = "pdf",
       width = 10, height = 4, dpi = 300, useDingbats = F)

RNA_BS_Comparison <- sapply(c("Bisque", "MGP"), function(Method){
  sapply(CellTypes, function(CellType){
    Same_tissue = cor.test(formula(paste0("~T1_RNA_", Method, "_", CellType, "+T1_BS_MSP_", CellType)), data = AllData)
    Different_tissue = cor.test(formula(paste0("~T2_RNA_", Method, "_", CellType, "+ T1_BS_MSP_", CellType)), data = AllData %>% filter(Resequenced == "No"))
    data.frame(CellType = CellType, "Same sample" = Same_tissue$est, "Different sample" = Different_tissue$est)
  }, simplify = F) %>% rbindlist() %>%
    pivot_longer(-CellType,
                 names_to = "SampleType",
                 values_to = "Corr") %>% mutate(Method = Method) 
}, simplify = F) %>% rbindlist() %>%  data.frame()

RNA_BS_Comparison$SampleType <- factor(RNA_BS_Comparison$SampleType,
                                       levels = unique(RNA_BS_Comparison$SampleType))

RNA_BS_Comparison$Method <- factor(RNA_BS_Comparison$Method, levels = c("MGP", "Bisque"))
  
  
P3 <- ggplot(RNA_BS_Comparison, aes(SampleType, CellType, fill = Corr)) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "") +
  geom_tile() +
  geom_text(aes(label = round(Corr, 2)), color = "white") +
  scale_fill_gradient2(low = "white", mid = MoviePalettes$GreenPallete[6], high = MoviePalettes$GreenPallete[3], midpoint = 0) +
  facet_wrap(~Method)

P4 <- ggplot(RNA_BS_Comparison, aes(Corr, fill = SampleType)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(), legend.position = "bottom") +
  labs(x = "RNA-WGBS correlation") +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = MoviePalettes$BugsLife[c(4, 2)], name = "") +
  facet_wrap(~Method) 

ggarrange(P3, P4, nrow = 2)

ggsave(paste0(ResultsPath, "RNA_BScomparison.pdf"), device = "pdf",
       width = 10, height = 6, dpi = 300, useDingbats = F)

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

ggarrange(similarityPlotList)


pheatmap(Cor, border_color = NA, na_col = "white", cluster_rows = T, cluster_cols = T,fontsize_number = 12,
         scale = "none", method = "ward.D2", display_numbers = T, breaks=seq(-1, 1, length.out=100), main = CellType)