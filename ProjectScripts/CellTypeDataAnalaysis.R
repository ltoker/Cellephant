if(!"devtools" %in% rownames(installed.packages())){
  BiocManager::install("devtools")
}

library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

source("ProjectScripts/ProjectFunctions.R")

plotMA = DESeq2::plotMA
packageF("org.Hs.eg.db")
packageF("tabulizer")
packageF("limma")
packageF("edgeR")
packageF("Rsubread")
packageF("cluster")
ResultsPath = "ResultsMarziReanalysis/"
packageF("ggpubr")
packageF("GGally")

ResultsPath = "ResultsEstimateComparison"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

ResultsPathMarzi = "ResultsMarziReanalysis/"


if("CellTypeSpecificCounts.Rda" %in% list.files(path = "Data")){
  load("Data/CellTypeSpecificCounts.Rda")
} else {
  #Import the peakset from Marzi et al. and create a SAF file
  temp <- read.table("Data/GSE102538_H3K27ac_EntorhinalCortex.bed.gz", header = F, sep = " ") %>%
    select(V4, V1, V2, V3)
  
  names(temp) <- c("GeneID" , "Chr",  "Start", "End")
  temp$Strand <- "*"
  
  #Quantified the H3K27ac reads from the NeuN+/- cellls based on Girdahal et al 2018 inside Marzi et al peaks
  temp2 <- featureCounts(files = list.files(path = "Data/CellTypesH3K27ac/BAMfiles/", full.names = T),nthreads = 5,
                         annot.inbuilt = "hg19",annot.ext = temp, isPairedEnd = T, countMultiMappingReads = F, countChimericFragments=F)
  
  SampleNames <- list.files(path = "Data/CellTypesH3K27ac/BAMfiles/", full.names = F)
  SampleNames <- sapply(SampleNames, function(x){
    x <- paste0(strsplit(x, "_")[[1]][c(7, 11)], collapse = "_")
    x <- gsub(".bam", "", x)
    if(grepl("\\+", x)){
      gsub("\\+", "pos", x)
    } else {
      gsub("-", "neg", x)
    }
  })
  
  CountsCells <- temp2$counts
  colnames(CountsCells) <- SampleNames
  
  #Annotate the peaks
  annoFileCollapsed <- GetGenomeAnno(genome = "hg19")
  PeakLocation <- temp %>% as(., "GRanges")
  seqlevelsStyle(PeakLocation) <- "UCSC"
  PeakLocation <- keepSeqlevels(PeakLocation, paste0("chr",c(1:22, "X", "Y")), pruning.mode = "coarse")
  
  PeakAnnoFile <- mergeByOverlaps(annoFileCollapsed, PeakLocation, maxgap = 0, type = "any", select = "all") %>% data.frame()
  names(PeakAnnoFile)[1:4] <- c("Region.CHR", "Region.START", "Region.END", "Region.width")
  PeakAnnoFile %<>% select(-matches("strand|_id|^id|^PeakName|annoFileCollapsed|^GeneID"))
  
  names(PeakAnnoFile)[grepl("PeakLocation", names(PeakAnnoFile))] <- c("Peak.CHR", "Peak.START", "Peak.END", "Peak.width", "PeakName")
  PeakAnnoFile %<>% mutate(Peak.Location = paste0(Peak.CHR, ":", Peak.START, "-", Peak.END))
  
  
  save(CountsCells, PeakAnnoFile, file = "Data/CellTypeSpecificCounts.Rda")
}
  
CellTypeMeta <- read.table("Data/CellTypesH3K27ac/MSSM_U01MH103392_EpiMap_Metadata_ChIPseq_August2016Release.csv",
                           header = T, sep = ",", comment.char = "!") %>% filter(HistoneMark == "H3K27ac", BrainRegion == "DLPFC")
CellTypeMeta2 <- read.table("Data/CellTypesH3K27ac/MSSM_U01MH103392_EpiMap_Metadata_clinical.csv",
                            header = T, sep = ",") 

CellTypeMeta <- merge(CellTypeMeta, CellTypeMeta2, by = "Individual_ID")

CellTypeMeta$Sample_ID2 <- sapply(as.character(CellTypeMeta$File_Name), function(x){
  x <- paste0(strsplit(x, "_")[[1]][c(7, 11)], collapse = "_")
  if(grepl("\\+", x)){
    gsub("\\+", "pos", x)
  } else {
    gsub("-", "neg", x)
  }
})

CellTypeMeta$CellType2 <- sapply(as.character(CellTypeMeta$CellType), function(x){
  if(grepl("\\+", x)){
    "Neuron"
  } else {
    "Glia"
  }
})

CellTypeMeta %<>% mutate(Hemisphere = tolower(Hemisphere))

#Run edgeR for the difference between the cell types, similar to analysis by Marzi et al.
#Keep peaks with cpm > 1 in at least two samples
keep <- rowSums(cpm(CountsCells) > 1) >= 2
CountsCells <- CountsCells[keep,]

#Arrance the order of samples in the metadata file
CellTypeMeta <- CellTypeMeta[match(colnames(CountsCells), CellTypeMeta$Sample_ID2),]

CellTypecountList<-DGEList(counts=CountsCells, group = CellTypeMeta$CellType2)
CellTypecountList <- calcNormFactors(CellTypecountList)


CellTypeModelMatrix <- model.matrix(as.formula("~Sex + Hemisphere + AgeDeath + pH + CellType2"), data = CellTypeMeta)

CellTypecountList <- estimateDisp(CellTypecountList, CellTypeModelMatrix)

fitTMM <- glmQLFit(CellTypecountList, CellTypeModelMatrix)
qlf_CellType <- glmQLFTest(fitTMM, coef = "CellType2Neuron")
CellType_results <- topTags(qlf_CellType, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.))

CellTypeResult_Anno <- CellType_results %>%
  AnnotDESeqResult(CountAnnoFile = PeakAnnoFile, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)


##### load AD analysis objects ############################
ResultsMarziCETs <- readRDS(paste0(ResultsPathMarzi, "MarziedgeR_CETs.Rds"))
ResultsMarziMSP <- readRDS(paste0(ResultsPathMarzi, "MarziedgeR_MSP.Rds"))
ResultsMarziMSPneuronal <- readRDS(paste0(ResultsPathMarzi, "MarziedgeR_MSPneuronAsFactor.Rds"))
ResultsMarziNoCellCorrection <- readRDS(paste0(ResultsPathMarzi, "MArziedgeR_NoCellCorretction.Rds"))
ResultsMarziCETsRandom <- readRDS(paste0(ResultsPathMarzi, "MarziedgeR_CETsRandom.Rds"))

###########################################################

SignifCETS <- ResultsMarziCETs %>% filter(FDR < 0.05, !duplicated(PeakName)) %>% select(PeakName, logCPM, logFC, PValue, FDR)
SignifMSPall <- ResultsMarziMSP %>% filter(FDR < 0.05, !duplicated(PeakName)) %>% select(PeakName, logCPM, logFC, PValue, FDR)
SignifMSPneuronF <- ResultsMarziMSPneuronal %>% filter(FDR < 0.05, !duplicated(PeakName)) %>% select(PeakName, logCPM, logFC, PValue, FDR)
SignifNoCorrect <- ResultsMarziNoCellCorrection %>% filter(FDR < 0.05, !duplicated(PeakName)) %>% select(PeakName, logCPM, logFC, PValue, FDR)
SignifCETsRandom <- ResultsMarziCETsRandom %>% filter(FDR < 0.05, !duplicated(PeakName)) %>% select(PeakName, logCPM, logFC, PValue, FDR)

SignifMarzi <- list(SignifNoCorrect, SignifCETsRandom, SignifCETS, SignifMSPneuronF, SignifMSPall)
names(SignifMarzi) <- c("NoCellCorrection", "CETS_Shuffled",  "CETs", "MSPneuronF", "MSPall")

CellType_resultsSignifMarzi <- sapply(names(SignifMarzi), function(Mod){
  temp <- CellType_results %>% filter(PeakName %in% SignifMarzi[[Mod]]$PeakName)
  temp$MethodSignif = Mod
  temp$DirectictionChange <- sapply(temp$PeakName, function(x){
    if(SignifMarzi[[Mod]] %>% filter(PeakName == x) %>% .$logFC < 0){
      "Hypoacetylated"
    } else if(SignifMarzi[[Mod]] %>% filter(PeakName == x) %>% .$logFC > 0){
      "Hyperacetylated"
    }
  })
  temp$logFC_AD <- SignifMarzi[[Mod]]$logFC[match(temp$PeakName, SignifMarzi[[Mod]]$PeakName)]
  temp
}, simplify = F) %>% rbindlist()


CellType_resultsSignifMarzi$DirectictionChange <- factor(CellType_resultsSignifMarzi$DirectictionChange, levels = c("Hypoacetylated", "Hyperacetylated"))
CellType_resultsSignifMarzi$MethodSignif <- factor(CellType_resultsSignifMarzi$MethodSignif, levels = c("NoCellCorrection", "CETS_Shuffled",  "CETs", "MSPneuronF", "MSPall"))



AllResults <- rbind(ResultsMarziNoCellCorrection %>% select(PeakName, logFC, FDR) %>% filter(!duplicated(PeakName)) %>% mutate(Method = "NoCellCorrection"),
                    ResultsMarziCETsRandom %>% select(PeakName, logFC, FDR) %>% filter(!duplicated(PeakName)) %>% mutate(Method = "CETS_Shuffled"),
                    ResultsMarziCETs %>% select(PeakName, logFC, FDR) %>% filter(!duplicated(PeakName)) %>% mutate(Method = "CETs"),
                    ResultsMarziMSPneuronal %>% select(PeakName, logFC, FDR) %>% filter(!duplicated(PeakName)) %>% mutate(Method = "MSPneuronF"),
                    ResultsMarziMSP %>% select(PeakName, logFC, FDR) %>% filter(!duplicated(PeakName)) %>% mutate(Method = "MSPall"))

AllResults <- merge(AllResults, CellType_results %>% select(PeakName, logFC, FDR), by = "PeakName", suffixes = c("_AD", "_Neurons"))
AllResults$Method <- factor(AllResults$Method, levels = c("NoCellCorrection", "CETS_Shuffled",  "CETs", "MSPneuronF", "MSPall"))

CorDF <- data.frame(Method = c("NoCellCorrection", "CETS_Shuffled",  "CETs", "MSPneuronF", "MSPall"),
                    x1 = -6, y1 = 2,
                    x2 = -6, y2 = 1.8)

CorDF$CorAll <- sapply(CorDF$Method, function(method){
  temp <- cor.test(~logFC_Neurons + logFC_AD, data = AllResults %>% filter(Method == method))
  round(temp$estimate, digits = 2)
})

CorDF$CorSignif <- sapply(CorDF$Method, function(method){
  temp <- cor.test(~logFC_Neurons + logFC_AD, data = AllResults %>% filter(Method == method, FDR_AD < 0.05))
  round(temp$estimate, digits = 2)
})

CorDF %<>% mutate(Text = paste0("'r'[italic('All peaks')] ", "*", "' = '", CorAll))
CorDF %<>% mutate(Text2 = paste0("'r'[italic('Significant peaks')] ", "*", "' = '", CorSignif))
CorDF$Method <- factor(CorDF$Method, levels = unique(CorDF$Method))

Plot1 <- ggplot(CellType_resultsSignifMarzi, aes(DirectictionChange, logFC)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),legend.position = c(0.85,0.4),
         legend.background = element_blank()) +
  labs(x = "", y = "logFC (Neurons vs. Glia)") +
  geom_violin(aes(fill = DirectictionChange)) +
  geom_boxplot(width = 0.2, outlier.shape = NA)+
  scale_fill_manual(values =  c("darkseagreen", "darkorchid4"), name = "") +
  geom_hline(yintercept = 0, color = "red") +
  facet_wrap(~MethodSignif, nrow = 2)


AllResults$Direction <- AllResults$logFC_AD
AllResults$Direction[AllResults$FDR_AD > 0.05] <- "NS"
AllResults$Direction[AllResults$FDR_AD < 0.05 & AllResults$logFC_AD < 0] <- "Down"
AllResults$Direction[AllResults$FDR_AD < 0.05 & AllResults$logFC_AD > 0] <- "Up"

Plot2 <- ggally_density(AllResults, aes_string("logFC_Neurons", "logFC_AD", fill = "..level..", color = "Direction", group = "Direction")) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "logFC (Neurons vs. Glia)", y = "logFC_AD") +
  scale_color_manual(values = c("orange", "white", "orange")) +
  geom_hline(yintercept = 0, color = "red") +
  geom_vline(xintercept = 0, color = "red") +
  geom_text(data = CorDF, aes(x = x1, y = y1, label = Text), parse = T, hjust = 0) +
  geom_text(data = CorDF, aes(x = x2, y = y2, label = Text2), parse = T, hjust = 0) +
  facet_wrap(~Method, nrow = 1)


Plot3 <- ggarrange(Plot1, Plot2, nrow = 1, widths = c(1,2.2), heights = c(1,1))
ggsave(paste0("MethodComparisonAll_new.pdf"), Plot3, device = "pdf", width = 14, height = 3, dpi = 300, path = ResultsPath, useDingbats = F)



####### Look at the cell-type specific enhancer regions based on Nott et al.
#Loading cell type specific promoter and enhancer data generated by Nott at al.
if(!"NOTT_2019.interactome.rda" %in% list.files("Data")){
  download.file("https://github.com/RajLabMSSM/echolocatoR/raw/master/data/NOTT_2019.interactome.rda",
                "Data/NOTT_2019.interactome.rda")
}
load("Data/NOTT_2019.interactome.rda")
names(NOTT_2019.interactome) <- sapply(names(NOTT_2019.interactome), function(x){
  gsub(" ", "_", x)
})


PeaksGRanges <- merge(PeakAnnoFile %>% filter(!duplicated(PeakName)) %>%
                        select(Peak.CHR, Peak.START, Peak.END, PeakName), AllResults, by = "PeakName") %>%
  select(c("Peak.CHR", "Peak.START", "Peak.END", names(AllResults))) %>%
  as("GRanges")

Overlaps <- sapply(names(NOTT_2019.interactome[-c(1:4)]), function(CellType){
  Type = NOTT_2019.interactome[[CellType]]
  temp <- as(Type, "GRanges")
  Overlap <- findOverlaps(temp,PeaksGRanges)
  PeaksGRanges[subjectHits(Overlap)] %>% data.frame() %>% mutate(CellType = CellType)
}, simplify = F) %>% rbindlist() %>% data.frame()

#Since the same peak can overlap with the several enhancer/promoter regions, filter to include
#each peak only once for each region/method/cell type
#Note - the same region can be defined as enahancer/promoter for multiple cell types
Overlaps %<>% mutate(FilterCol = paste(PeakName, Method, CellType, sep = ".")) %>%
  filter(!duplicated(FilterCol))


Overlaps$RegionType <- sapply(Overlaps$CellType, function(x){
  strsplit(x, "_")[[1]][2]
})

Overlaps$CellType <- sapply(Overlaps$CellType, function(x){
  strsplit(x, "_")[[1]][1]
})

UniqeRegions <- Overlaps %>% filter(Method == "CETs") %>% droplevels() %>%
  .$PeakName %>% table %>% data.frame()

names(UniqeRegions)[1] <- "PeakName"
Overlaps <- merge(Overlaps, UniqeRegions, by = "PeakName")
Overlaps$CellType <- factor(Overlaps$CellType, levels = c("Astrocyte", "Microglia", "Oligo", "Neuronal"))
Overlaps %<>% mutate(logFC_Glia = -1*logFC_Neurons) #For easier comparison with AD

Overlaps$DirectionAD <- sapply(Overlaps$logFC_AD, function(x){
  if(x < 0){
    "Hypoacetylated"
  } else {
    "Hyperacetylated"
  }
})

Overlaps$DirectionNeurons <- sapply(Overlaps$logFC_Neurons, function(x){
  if(x < 0){
    "Hypoacetylated"
  } else {
    "Hyperacetylated"
  }
})

Overlaps$DirectionGlia <- sapply(Overlaps$logFC_Glia, function(x){
  if(x < 0){
    "Hypoacetylated"
  } else {
    "Hyperacetylated"
  }
})

Prop <- Overlaps  %>% filter(Freq == 1, FDR_AD < 0.05) %>%
  group_by(Method, CellType, DirectionAD) %>% summarise(n = n()) %>%
  data.frame() %>% mutate(HypoDAR = "", HyperDAR = "", TotalDAR = NA)

for(x in unique(Prop$Method)){
  for(Direct in c("Hypoacetylated", "Hyperacetylated")){
    Col = gsub("acetylated", "DAR", Direct)
    Prop[Prop$Method==x & Prop$DirectionAD == Direct, Col] <- Prop %>% filter(Method == x, DirectionAD == Direct ) %>% .$n %>% sum 
  }
}

for(x in unique(Prop$Method)){
    Prop[Prop$Method==x, "TotalDAR"] <- Prop %>% filter(Method == x) %>% .$n %>% sum 
}  

Prop %<>% mutate(Proportion = n/TotalDAR,
                 Proportion = signif(Proportion, digits = 2))

Prop$ProportionText <- sapply(Prop$Proportion, function(x){
  x <- as.character(x)
  gsub("0*?0$", "", x)
})

Prop %<>% mutate(DirectionNumber = paste0(DirectionAD, "\n(", HypoDAR, HyperDAR, ")"))


PlotNottDAR <- ggplot(Prop, aes(DirectionNumber, Proportion, fill = CellType)) +
  theme_classic() +
  labs(x = "", y = "Proportion") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = MoviePalettes$SpiritedAway[c(5,7,9, 1)]) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  facet_wrap(~Method, scales = "free_x", nrow = 1)


PlotNottAll <- ggplot(Overlaps %>% filter(Freq == 1), aes(CellType, logFC_AD, fill = CellType)) +
  theme_classic() +
  labs(x = "") +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white") +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  scale_fill_manual(values = MoviePalettes$SpiritedAway[c(5,7,9, 1)]) +
  facet_wrap(~Method, nrow = 1)


ggarrange(PlotNottDAR, PlotNottAll, nrow = 2, common.legend = T)
ggsave(paste0(ResultsPath, "NottValidation.pdf"), device = "pdf", width = 12, height = 6, dpi = 300)


PropGlia_vs_Neuron <- Overlaps  %>% filter(Freq == 1, FDR_Neurons < 0.05, Method == "CETs") %>%
  group_by(Method, CellType, DirectionGlia) %>% summarise(n = n()) %>%
  data.frame() %>% mutate(HypoDAR = "", HyperDAR = "", TotalDAR = NA)

for(Direct in c("Hypoacetylated", "Hyperacetylated")){
  Col = gsub("acetylated", "DAR", Direct)
  PropGlia_vs_Neuron[PropGlia_vs_Neuron$DirectionGlia == Direct, Col] <- PropGlia_vs_Neuron %>%
    filter(DirectionGlia == Direct ) %>% .$n %>% sum 
}

PropGlia_vs_Neuron$TotalDAR <- PropGlia_vs_Neuron$n %>% sum 

PropGlia_vs_Neuron %<>% mutate(Proportion = n/TotalDAR,
                               Proportion = signif(Proportion, digits = 2))


PropGlia_vs_Neuron %<>% mutate(DirectionNumber = paste0(DirectionGlia, "\n(", HypoDAR, HyperDAR, ")"))


PlotNottGlia_vs_NeuronDAR <- ggplot(PropGlia_vs_Neuron, aes(DirectionNumber, Proportion, fill = CellType)) +
  theme_classic() +
  labs(x = "", y = "Proportion") +
  geom_bar(stat = "identity", show.legend = F) +
  scale_fill_manual(values = MoviePalettes$SpiritedAway[c(5,7,9, 1)]) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5))

PlotNottGlia_vs_NeuronAll <-  ggplot(Overlaps %>% filter(Method == "CETs", Freq == 1), aes(CellType,logFC_Glia, fill = CellType)) +
  theme_classic() +
  labs(x = "", y = "logFC glia vs. neurons") +
  geom_violin() +
  geom_boxplot(width = 0.2, fill = "white") +
  geom_hline(yintercept = 0, color = "red", lty = "dashed") +
  scale_fill_manual(values = MoviePalettes$SpiritedAway[c(5,7,9, 1)])




ggarrange(ggarrange(PlotNottGlia_vs_NeuronDAR, PlotNottGlia_vs_NeuronAll, ncol = 2, widths = c(1,1.5)),
          ggarrange(PlotNottDAR, PlotNottAll, nrow = 2, common.legend = T), nrow = 2, heights = c(1, 2))

ggsave(paste0(ResultsPath, "NottValidationCombined.pdf"), device = "pdf", width = 12, height = 8, dpi = 300)

write.table(Prop %>% select(-ProportionText, -DirectionNumber), paste0(ResultsPath, "ADoverlapNott.txt"),
            sep = "\t", row.names = F, col.names = T)

save.image(paste0(ResultsPath, "CellTypeComparison.Rdata"))
