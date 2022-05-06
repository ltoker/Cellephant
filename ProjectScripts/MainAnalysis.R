#BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")
packageF("pheatmap")
packageF("DESeq2")
packageF("gridExtra")

ResultsPath = "ResultsEstimateComparison"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

# AllData <- readRDS("Data/AllEstimatesAndMeta.Rds")
# AllData %<>% mutate(across(contains("MGP"), rescale, c(0,1)))

BS_H3K27acEstimates <- readRDS("Data/BS_H3K27acEstimates.Rds")

RNAestimates_MGP_Bisque <- read.table("EstimatesResults/RNAestimatesMGP_Bisque.tsv",
                                      header = T, sep = "\t")


AllData <- merge(RNAestimates_MGP_Bisque, BS_H3K27acEstimates,
                 by = "Biobank_ID")


names(AllData) <- sapply(names(AllData), function(x){
  gsub("s$", "", x)
})

SuppTable2 <- AllData %>% select(Biobank_ID, Cohort, Group, Age, Sex, PMI, RIN, Resequenced)
write.table(SuppTable2, paste0(ResultsPath, "SuppTable2.tsv"), sep = "\t", row.names = F, col.names = T)


AllData_long <- pivot_longer(AllData, cols = matches("^T._"), names_sep = "_",
                             names_to = c("Tissue", "Matherial", "Method", "CellType"), values_to = "Estimate")



CellTypes <- unique(AllData_long$CellType)

BrainDecode <- sapply(list.files("Data", pattern = "CIBERSORT|dtangle"), function(File){
  temp <- read.table(paste0("Data/", File), sep = "\t", header = T)
  Method = strsplit(File, "_")[[1]][1]
  if(grepl("RNAseq1", File)){
    IDcol = "RNAseq_id_ParkOme1"
    Pref  = paste0("T1_RNA_", Method, "_")
  } else if(grepl("RNAseq2", File)){
    IDcol = "RNAseq_id_ParkOme2"
    Pref  = paste0("T2_RNA_", Method, "_")
  }
  rownames(temp) <- temp[,1]
  temp <- temp[,-1]
  names(temp) <- paste0(Pref, names(temp))
  
  temp$Biobank_ID <- AllData$Biobank_ID[match(rownames(temp), AllData[[IDcol]])]
  temp
}, simplify = F)

AllRNAestimates <- RNAestimates_MGP_Bisque

for(Val in names(BrainDecode)){
  AllRNAestimates <- merge(AllRNAestimates, BrainDecode[[Val]], by = "Biobank_ID")
}

AllRNAestimateslong <- pivot_longer(AllRNAestimates, cols = matches("^T._"), names_sep = "_",
                                    names_to = c("Tissue", "Matherial", "Method", "CellType"), values_to = "Estimate")
AllRNAestimateslong$Method <- factor(AllRNAestimateslong$Method, levels = c("Bisque", "dtangle", "CIBERSORT", "MGP"))
AllRNAestimateslong$CellType <- factor(AllRNAestimateslong$CellType, levels = unique(AllRNAestimateslong$CellType))

AllRNAestimateslong$Dataset <- sapply(AllRNAestimateslong$Tissue, function(x){
  if(x == "T1"){
    "RNAseq1"
  } else {
    "RNAseq2"
  }
})

######## Comparisson to IHC data and all method correlations #############

#Get IHC quantification of cortical samples from Patrick et el. 2020
IHC_Patrick <- sapply(c("astro", "endo", "microglia", "oligo", "neuro"), function(CellType){
  temp <- read.table(paste0("https://github.com/ellispatrick/CortexCellDeconv/raw/master/CellTypeDeconvAnalysis/Data/IHC.",
                                    CellType, ".txt"), sep = "\t", header = T) %>% t %>% data.frame()
  names(temp)[1] <- "Estimate"
  CellType <- if(CellType == "astro"){
    "Astrocyte"
  } else if (CellType == "endo"){
    "Endothelial"
  } else if (CellType == "microglia"){
    "Microglia"
  }
  else if (CellType == "oligo"){
    "Oligodendrocyte"
  }
  else if (CellType == "neuro"){
    "Neuron"
  }
  temp %>% mutate(SubjectID = rownames(.), CellType = CellType, Dataset = "Patrick_et.al", Method = "IHC" )
}, simplify = F) %>% rbindlist() %>% data.frame()

#Plotting comparison between IHC-based estimates from Patrick et al 2020 and across-cell based
#estimation approaches
IHC_RNA <- rbind(AllRNAestimateslong %>% filter(Method != "MGP") %>% mutate(RNAdataset = Dataset) %>%
                   select(Estimate, CellType, Dataset, RNAdataset, Method),
                 IHC_Patrick %>% select(-SubjectID) %>% mutate(RNAdataset = "RNAseq1"),
                 IHC_Patrick %>% select(-SubjectID) %>% mutate(RNAdataset = "RNAseq2"),
                 data.frame(Estimate = 0.5,  CellType = CellTypes, Dataset = "temp",  RNAdataset = "RNAseq1", Method = "bla"),
                 data.frame(Estimate = 0.5,  CellType = CellTypes, Dataset = "temp", RNAdataset = "RNAseq2",  Method = "bla")) %>% droplevels()

IHC_RNA$Method <- factor(IHC_RNA$Method, levels = levels(IHC_RNA$Method)[c(1:3, 5, 4)])

ggplot(IHC_RNA, aes(CellType, Estimate, fill = Method, color = Method)) +
  theme_bw() +
  labs(x = "", y = "Estimate (proportion of cells)") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#831F02", MoviePalettes$SpiritedAway[c(4,9)], "green", "black")) +
  scale_color_manual(values = c("#831F02", MoviePalettes$SpiritedAway[c(4,9)], "green", "black")) +
  geom_point(position=position_jitterdodge(jitter.width = 0.3), show.legend = F, size = 0.3) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), lty = "dashed") + 
  facet_wrap(~RNAdataset, nrow = 2) 
ggsave(paste0(ResultsPath, "RNA_IHC_estimates.pdf"), device = "pdf", width = 10, height = 5, useDingbats = F)

ggplot(IHC_RNA %>% filter(Dataset != "temp"), aes(Method , Estimate, fill = Method, color = Method)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =  1)) +
  labs(x = "", y = "Estimate (proportion of cells)") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#831F02", MoviePalettes$SpiritedAway[c(4,9)], "black")) +
  scale_color_manual(values = c("#831F02", MoviePalettes$SpiritedAway[c(4,9)], "black")) +
  geom_point(position=position_jitterdodge(jitter.width = 1), size = 0.3, show.legend = F) +
  facet_grid(CellType~Dataset, scales = "free_x", space = "free")
ggsave(paste0(ResultsPath, "RNA_IHC_estimates2.pdf"), device = "pdf", width = 7, height = 7, useDingbats = F)

AllRNAcomparison <- sapply(CellTypes, function(CellType){
  sapply(c("T1", "T2"), function(TS){
     Data <- AllRNAestimates %>% select(matches(CellType)) %>%
       select(matches(TS))
     cor(Data)
  }, simplify = F)
 
}, simplify = F)

AllRNAcomparisonDF <- sapply(names(AllRNAcomparison), function(CellType){
  temp <- AllRNAcomparison[[CellType]]
  sapply(names(temp), function(TS){
    temp2 <- temp[[TS]] %>% data.frame()
    diag(temp2) <- NA
    names(temp2) <- sapply(names(temp2), function(x){
      strsplit(x, "_")[[1]][[3]]
    })
    temp2 %<>%  mutate(CellType = CellType, Tissue = TS, Dataset = TS)
    temp2$Dataset <- sapply(temp2$Dataset, function(x){
      gsub("T", "RNAseq", x)
    })
    temp2
  }, simplify = F) %>% rbindlist() %>% data.frame()
}, simplify = F) %>% rbindlist() %>% data.frame() %>%
  pivot_longer(cols = c(1:4), names_to = "Method", values_to = "Correlation") %>% data.frame()

AllRNAcomparisonDF$CellType <- factor(AllRNAcomparisonDF$CellType, levels = unique(AllRNAcomparisonDF$CellType) )
AllRNAcomparisonDF$Method <- factor(AllRNAcomparisonDF$Method, levels = c("Bisque", "dtangle",
                                                                          "CIBERSORT", "MGP"))

ggplot(AllRNAcomparisonDF, aes(CellType, Correlation, color = Method, fill = Method)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  labs(x = "") +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.08), show.legend = F) +
  scale_color_manual(values = c("#831F02", MoviePalettes$SpiritedAway[c(4,9)], "#339656")) +
  scale_fill_manual(values = c("#831F02", MoviePalettes$SpiritedAway[c(4,9)], "#339656")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), lty = "dashed")
ggsave(paste0(ResultsPath, "AllRNA_estimatesCor.pdf"), device = "pdf", width = 10, height = 3, useDingbats = F)

  
##################################################################################################

RNA_MethodComparison <-   sapply(CellTypes, function(CellType){
  T1 = cor.test(formula(paste0("~T1_RNA_MGP_", CellType, "+T1_RNA_Bisque_", CellType)), data = AllData)
  T2 = cor.test(formula(paste0("~T2_RNA_MGP_", CellType, "+T2_RNA_Bisque_", CellType)), data = AllData %>% filter(Resequenced == "No"))
  data.frame(CellType = CellType, TS1 = T1$est, TS2 = T2$est)
}, simplify = F) %>% rbindlist() %>% data.frame() %>% pivot_longer(-CellType,
                                                                   names_to = "TissueSample",
                                                                   values_to = "Cor") %>% data.frame()
cor.test(AllData$T1_BS_MSP_Neuron, AllData$T1_BS_CETS_Neuron)



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
               Method = c("CETS" = MoviePalettes$SpiritedAway[2],
                          "MSP" = MoviePalettes$SpiritedAway[4],
                          "MGP" = MoviePalettes$SpiritedAway[6],
                          "Bisque" = MoviePalettes$SpiritedAway[7]),
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
  p <-  pheatmap(Cor, border_color = NA, na_col = "white", cluster_rows = T, cluster_cols = T,fontsize_number = 12, angle_col = 90,
                 annotation_row  = AnnotationsDF %>% select(Method, OmicsType),
                 annotation_col  = AnnotationsDF %>% select(-Method, -OmicsType),
                 annotation_colors = AnnoCol, cutree_rows = 2, cutree_cols = 2,
                 scale = "none", method = "ward.D2", display_numbers = T,
                 breaks=seq(-1, 1, length.out=100), main = CellType, filename = paste0(ResultsPath, CellType,
                                                                                       "AllEstimateCor.pdf"),
                 width = 8, height = 6)
  p$gtable
 
},simplify = F)


