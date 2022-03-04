if(!"devtools" %in% rownames(installed.packages())){
  BiocManager::install("devtools")
}

library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

source("ProjectScripts/ProjectFunctions.R")

plotMA = DESeq2::plotMA
packageF("org.Hs.eg.db")
if(!"tabulizer" %in% rownames(installed.packages())){
  install_github("ropensci/tabulizer")
}

library("tabulizer")
packageF("limma")
packageF("edgeR")
packageF("RColorBrewer")
packageF("ggpubr")
packageF("venn")
packageF("cluster")
packageF("patchwork")
packageF("tidyverse")
packageF("GGally")


ResultsPath = "ResultsMarziReanalysis"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = F)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")


CellTypePeakCountLoc = "Data/marzi_counts.tsv"
CountMatrixLoc = "Data/MarziPaperCounts.tsv"
Cohort = "Marzi"

MetadataSup <- extract_tables("Data/41593_2018_253_MOESM1_ESM.pdf", pages = c(27:30)) %>% lapply(data.frame) %>% rbindlist() #Reading the metadata from the supplementary table 
names(MetadataSup) <- c("SampleID", "Group", "BraakStage", "Age", "Sex", "NeuralProportion", "PMI", "Experiment")
MetadataSup <- MetadataSup[-c(1:3),]
MetadataSup$Age <- as.numeric(as.character(MetadataSup$Age))
MetadataSup$PMI <- round(as.numeric(as.character(MetadataSup$PMI))/60, digits = 1)
MetadataSup$NeuralProportionNumeric <- sapply(MetadataSup$NeuralProportion, function(x){
  if(grepl("%", x)){
    x = gsub("%", "", x) %>% as.character() %>% as.numeric()
    x/100
  } else {
    NA
  }
})

MetadataSup %<>% mutate(MergeColumn = paste(Group, Age, Sex, NeuralProportionNumeric, sep = "_"))

#Getting additional metadata
softDown("GSE102538", file = "Data/GSE102538.soft")
Metadata <- ReadSoft("Data/GSE102538.soft") %>% select(-matches("antib|Sample_source|Sample_platform|orga")) %>% data.frame()
names(Metadata) <- c("GSM", "SampleName", "Group", "Age", "Sex", "CETS")
Metadata$Age <- as.numeric(as.character(Metadata$Age))
Metadata$CETS <- as.numeric(as.character(Metadata$CETS))

Metadata %<>% arrange(CETS)
Metadata$Eno2 <- c(0.65, 0.99, 0.56, 0.48, 0.56, 0.55, 0.69, 0.26, 1.27, 0.88,
                      0.42, 0.51, 2.07, 0.84, 0.52, 1.30, 0.90, 0.80, 0.19, 1.16,
                      0.34, 0.69, 0.34, 0.77, 1.43, 0.72, 0.33, 0.36, 0.29, 0.59,
                      0.42, 0.83, 0.49, 1.53, 0.98, 1.31, 0.85, 1.61, 0.65, 1.11,
                      0.72, 1.26, 0.80, 0.48,  2.28, 1.17, NA)
Metadata$deltaCTEno2 <- sapply(Metadata$Eno2, function(x){
  if(!is.na(x)){
    -log2(x) 
  } else {
    NA
  }
}) 

Metadata$Group <- sapply(as.character(Metadata$Group), function(x) gsub("C", "Control", x))
Metadata$Group <- factor(Metadata$Group, levels = c("Control", "AD"))

Metadata$NeuralProportion <- round(Metadata$CETS, digits = 2)
Metadata %<>% mutate(MergeColumn = paste(Group, Age, Sex, NeuralProportion, sep = "_"))

Metadata <- merge(Metadata %>% select(-NeuralProportion), MetadataSup %>% select(SampleID, BraakStage, PMI, MergeColumn), by = "MergeColumn", all.x = T, sort = F )


#Following the steps in Marzi et al.
Metadata$CETS[is.na(Metadata$CETS)] <- median(Metadata$CETS, na.rm = T)
Metadata$CETSif <-  cut(Metadata$CETS, breaks=5, ordered_result = T)
Metadata$Agef <-  cut(Metadata$Age, breaks=5, ordered_result = T)
Metadata$BraakStage <- as.numeric(as.character(Metadata$BraakStage))
Metadata %<>% arrange(GSM)

# In the merging of the supplement metadada and GEO metada, could not differentiate between Control10 and Control 13.
# Since they have similar braak stages, but different PMI times (31/55) and PMI is not included in the analysis, randomly assign the PMI
Metadata %<>% filter(!duplicated(GSM))

rownames(Metadata) <- Metadata$GSM %>% as.character()
Metadata %<>% select(-MergeColumn)


##### Add relative cell proportion for based on differential NeuN positive and negative cell H3K27ac peaks ##########  
Metadata <- GetCellularProportions(Metadata = Metadata)
names(Metadata)[grepl("Microglia", names(Metadata))] <- sapply(names(Metadata)[grepl("Microglia", names(Metadata))], function(x){
  gsub("ivation", "", x)
})

Metadata$NeuNall_MSPif <-  cut(Metadata$NeuNall_MSP, breaks=5, ordered_result = T)
Metadata$Oligo_MSPif <-  cut(Metadata$Oligo_MSP, breaks=5, ordered_result = T)
Metadata$Microglia_MSPif <-  cut(Metadata$Microglia_MSP, breaks=5, ordered_result = T)

EstimateMelt <- Metadata %>% mutate(CETS = rescale(CETS, c(0,1)),
                                    deltaCTEno2 = rescale(deltaCTEno2, c(0,1))) %>%
  select(Group, CETS, deltaCTEno2, NeuNall_MSP, Agef, Sex) %>% gather(key = "Method", value = "Estimate", -Group, -Agef, -Sex)

EstimateStat <- sapply(unique(EstimateMelt$Method), function(method){
  data <- EstimateMelt %>% filter(Method == method)
  stat <- lm(Estimate~Group + Sex + Agef, data = data) %>% summary %>% .$coef %>% .[2,4]
  signif <- if(stat > 0.05){
    ""
  } else if(stat > 0.01){
    "*"
  } else if(stat > 0.001){
    "**"
  }
  data.frame(Method = method, Text = paste0("p = ", signif(stat, digits = 2), signif), x = 1.5, y = 1.2)
}, simplify = F) %>% rbindlist()



#Look at the differences in the different estimated of cell types
CellData <- Metadata %>% select(-matches("__|Total")) %>% gather(matches("MSP|CETS$"), key = "CellType", value = "MSP") 
CellData$CellType <- sapply(CellData$CellType, function(x) gsub("_MSP", "", x)) %>%
  factor(levels = c("Astrocyte", "Endothelial", "Microglia", "Microglia_act",
                    "Microglia_deact","Oligo", "OligoPrecursors",
                    "GabaVIPReln", "Pyramidal", "NeuNall","CETS"))

#Get confidence intervals for differences for all cell type MSPs
#First round
CellTypeStats <- sapply(levels(CellData$CellType), function(cellType){
  Data = CellData %>% filter(CellType == cellType)
  lm(MSP~Group + Sex + Agef, data = Data) 
}, simplify = F)

#Adjusting for differences in neurons, since they change across the groups
CellData2 <- Metadata %>% select(-matches("__|Total")) %>% gather(matches("MSP|CETS$"), -NeuNall_MSP, key = "CellType", value = "MSP") 
CellData2$CellType <- sapply(CellData2$CellType, function(x) gsub("_MSP", "", x)) %>%
  factor(levels = c("Astrocyte", "Endothelial", "Microglia", "Microglia_act",
                    "Microglia_deact","Oligo", "OligoPrecursors",
                    "GabaVIPReln", "Pyramidal","CETS"))

CellTypeStats2 <- sapply(levels(CellData2$CellType), function(cellType){
  Data = CellData2 %>% filter(CellType == cellType)
  lm(MSP~Group + Sex + Agef + NeuNall_MSP, data = Data) 
}, simplify = F)

#Adjusting for both neurons and miroglia (microglia is differential after adjustment for neurons)
CellData3 <- Metadata %>% select(-matches("__|Total")) %>% gather(matches("MSP|CETS$"), -NeuNall_MSP, -Microglia_MSP, key = "CellType", value = "MSP") 
CellData3$CellType <- sapply(CellData3$CellType, function(x) gsub("_MSP", "", x)) %>%
  factor(levels = c("Astrocyte", "Endothelial", "Microglia_act",
                    "Microglia_deact","Oligo", "OligoPrecursors",
                    "GabaVIPReln", "Pyramidal","CETS"))

CellTypeStats3 <- sapply(levels(CellData3$CellType), function(cellType){
  Data = CellData3 %>% filter(CellType == cellType)
  lm(MSP~Group + Sex + Agef + NeuNall_MSP + Microglia_MSP, data = Data) 
}, simplify = F)


GetConfInt <- function(StatData, Name){
  CellTypeStatSummary <- lapply(StatData, function(x){
    temp <- summary(x) %>% .$coef %>% .[2, c(1,4)]
    temp2 <- confint(x) %>% .[2,]
    out <- c(temp, temp2) %>% t
    colnames(out) <- c("Coeficient", "pValue", "Low", "High")
    out
  }) %>% do.call(rbind, .) %>% data.frame()
  CellTypeStatSummary$Signif <- sapply(CellTypeStatSummary$pValue, function(x){
    if(x < 0.05){
      "Yes"
    } else {
      "No"
    }
  })
  rownames(CellTypeStatSummary) <- names(StatData)
  CellTypeStatSummary %<>% mutate(CellType = factor(rownames(.), levels = rownames(.)),
                                  Model = Name)
  return(CellTypeStatSummary)
}

CellTypeStatSummaryNoMSPcorrection <- GetConfInt(CellTypeStats, Name = "Agef and Sex adjusted")
CellTypeStatSummaryNeuNallcorrection <- GetConfInt(CellTypeStats2, Name = "Agef, Sex and NeunAll adjusted")
CellTypeStatSummaryNeuNallMicrogliacorrection <- GetConfInt(CellTypeStats3, Name = "Agef, Sex, NeunAll and Microglia adjusted")

CellTypeStatSummaryCombined <- rbind(CellTypeStatSummaryNoMSPcorrection,
                                     CellTypeStatSummaryNeuNallcorrection,
                                     CellTypeStatSummaryNeuNallMicrogliacorrection)


MethodsBoxplot_plot <- ggplot(EstimateMelt, aes(Group, Estimate)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.box.spacing = unit(-0.4, "cm"),
        strip.text = element_text(face = "bold", size = 12)) +
  labs(x = "") +
  geom_boxplot(outlier.shape = NA, aes(fill = Group)) +
  geom_jitter(height = 0, width = 0.2) +
  scale_fill_manual(values = c("dodgerblue4" , "chocolate1")) +
  geom_text(data = EstimateStat, aes(x = x, y = y, label = Text))+
  facet_wrap(~Method, nrow = 1)

MSPsCI_plot <- ggplot(CellTypeStatSummaryCombined, aes(CellType, Coeficient, color = Signif )) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(face = "bold")) +
  labs(x = "", title = ) +
  geom_hline(yintercept = 0, color = "red", linetype = 2, size = 0.9) +
  geom_errorbar(aes(ymin = Low, ymax = High, size = Signif), show.legend = F) +
  geom_point(size = 2, show.legend = F) +
  scale_color_manual(values = c("black", "palevioletred3")) +
  scale_size_manual(values = c(1, 1.2))+
  facet_wrap(~Model, scales = "free", nrow = 3)


Plot <- ggarrange(MethodsBoxplot_plot, MSPsCI_plot, widths = c(1.1,1), labels = c("a", "b"))
ggsave(paste0("CellTypeDifferences", Cohort, ".pdf"), plot = Plot, scale = 1.3, device = "pdf", width = 10, height = 6, dpi = 300, useDingbats = F, path = ResultsPath)


############################# HTseq counts ######################################################################

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")

HTseqCounts <- read.table(CountMatrixLoc, header = T, sep = "\t")

names(HTseqCounts) <- sapply(names(HTseqCounts), function(x){
  x = gsub(".*bamfiles.0?", "", x)
  strsplit(x,"_")[[1]][1]
})


HTseqCounts %<>% mutate(Peak.Location = paste0("chr", CHR, ":", START, "-", END))
names(HTseqCounts)[2:4] <- c("CHR", "START", "END")

#Arrange the samples to match Metadata order
HTseqCounts %<>% select(c("Geneid", "CHR", "START", "END", "Strand",  "Length", as.character(Metadata$GSM)))

#Keep peaks with cpm > 1 in at least two samples
keep <- rowSums(cpm(HTseqCounts %>% select(matches("GSM"))) > 1) >= 2
HTseqCounts <- HTseqCounts[keep,]

HTseqCounts$MaxCount = apply(HTseqCounts %>% select(matches("GSM")), 1, max)
HTseqCounts %<>% mutate(NormalizedMaxCount = MaxCount/Length)  


AllCalledData <- GetCountMatrixHTseq(countsDF = HTseqCounts, meta = Metadata, MetaCol = c("GSM", "SampleName", "Group", "Age", "Agef", "Sex", "BraakStage", "PMI", "CETS", "CETSif", grep("Eno2|MSP", names(Metadata), value = "T")),
                                     MetaSamleCol = "GSM", countSampleRegEx = "GSM", OtherNormRegEx = "^C1orf43_|^CHMP2A_|^EMC7_|^GPI_|^PSMB2_|^PSMB4_|^RAB7A_|^REEP5_|^SNRPD3_|^VCP_|^VPS29")

ggplot(AllCalledData$SampleInfo, aes(Group, RiP_NormMeanRatioOrg)) +
  geom_boxplot() +
  geom_point()

############ Look at the variance explained by different factors  #########################
countMatrixFullAllCalled <- GetCollapsedMatrix(countsMatrixAnnot = AllCalledData$countsMatrixAnnot %>% filter(!duplicated(.$PeakName)), collapseBy = "PeakName",MetaSamleCol = "GSM",
                                               CorMethod = "pearson", countSampleRegEx = "GSM",MetaSamleIDCol = "GSM", groupCol = "Group",
                                               FilterBy = "", meta = AllCalledData$SampleInfo, title = paste0("Sample correlation, ", Cohort))

#Get the pvalues for associasion of each covariate with the first 5 PCs
PCAsamples <- prcomp(t(countMatrixFullAllCalled$CPMdata), scale. = T)
countMatrixFullAllCalled$Metadata %<>% mutate(PC1 = PCAsamples$x[,1],
                                              PC2 = PCAsamples$x[,2],
                                              PC3 = PCAsamples$x[,3],
                                              PC4 = PCAsamples$x[,4],
                                              PC5 = PCAsamples$x[,5]) 
VarExplained <- PCAsamples %>% summary() %>% .$importance %>%
  .[2, 1:sum(grepl("^PC", names(countMatrixFullAllCalled$Metadata)))]*100 


CovarPvalues <- sapply(grep("^PC", names(countMatrixFullAllCalled$Metadata), value = T), function(PC){
  temp <- lm(as.formula(paste0(PC, "~ Group + Age + Sex + CETS + NeuNall_MSP + Astrocyte_MSP + Microglia_MSP + Oligo_MSP + Eno2")),
             data = countMatrixFullAllCalled$Metadata) %>% summary
  temp$coefficients[-1,4]
}, simplify = F) %>% do.call(cbind, .) %>% data.frame()

names(CovarPvalues) <- paste0(names(CovarPvalues), "(", round(VarExplained, digits = 1), "%)")
CovarPvalues %<>% mutate(Variable = factor(rownames(CovarPvalues), levels = rownames(CovarPvalues)))
CovarPvaluesMelt <- gather(CovarPvalues, key = "PC", value = "pValue", -Variable)
CovarPvaluesMelt %<>% mutate(pPvalue = -log10(pValue)) 
CovarPvaluesMelt$Significant <- sapply(CovarPvaluesMelt$pValue, function(x){
  if(x < 0.05){
    "Yes"
  } else {
    "No"
  }
}) %>% factor(levels = c("No", "Yes"))

Plot  <- ggplot(CovarPvaluesMelt, aes(PC, Variable)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x.top = element_text(angle = 90, vjust=0.5, hjust=0), axis.text = element_text(size = 14)) +
  labs(x = "", y = "", title = "Association of variables with main PCs") +
  geom_tile(aes(fill = pPvalue), colour = "white") +
  scale_fill_gradient(high = "steelblue", low = "gray94", name = "-log10(p)") +
  geom_text(aes(label = signif(pValue, 2), colour = Significant), show.legend = F, size = 5) +
  scale_color_manual(values = c("black", "red"))
ggsave(paste0("AssociationWithPCs", Cohort, ".pdf"), plot = Plot, device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F, path = ResultsPath)



#Detect outliers
MedianCor <- apply(countMatrixFullAllCalled$SampleCor, 1, function(x) median(x, na.rm = TRUE))
Outlier <- MedianCor[MedianCor < (median(MedianCor) - 1.5*iqr(MedianCor))]


# ADgene <- read.table("Data/AD_implicated_genes.txt", header = T, sep = "\t")
# PDgene <- read.table("../ChIPseqPD_reproduce/GeneralResults/PDgeneStat.tsv", header = T, sep = "\t")


######### Run analyses (repeating the steps in Marzi et al.) ##################
ddsAD <- DESeqDataSetFromMatrix(countData = countMatrixFullAllCalled$countMatrix, colDat = countMatrixFullAllCalled$Metadata, design = ~Group)
ddsAD <- estimateSizeFactors(ddsAD)

ADcounts <- counts(ddsAD)

ADcounts<-as.data.frame(ADcounts)


ADcountList<-DGEList(counts=ADcounts, group=Metadata$Group)
ADcountList <- calcNormFactors(ADcountList)


MarziModelMatrix <- model.matrix(as.formula("~Agef + CETSif + Group"), data = Metadata)

ADcountList <- estimateDisp(ADcountList, MarziModelMatrix)


fitTMM <- glmQLFit(ADcountList, MarziModelMatrix)
qlf_group <- glmQLFTest(fitTMM, coef = "GroupAD")
group_results <- topTags(qlf_group, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.))

groupResult_Anno <- group_results %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)



#Heatmap of significant peaks (Attempt to reproduce Fig. 4 from Marzi et al.)
PlotPeakHeatmap <- function(data, meta, title){
  Neuron_CETsColor <- rev(brewer.pal(n = 5, "Greys"))
  names(Neuron_CETsColor) <- levels(meta$CETSif)
  
  AgefColor <- brewer.pal(n = 5, "Greens")
  names(AgefColor) <- levels(meta$Agef)
  
  annoCol = data.frame(Age = meta$Agef,
                       Sex = meta$Sex,
                       Neuron_CETs = meta$CETSif,
                       NeuNall_MSP = meta$NeuNall_MSP,
                       Oligo_MSP = meta$Oligo_MSP,
                       Microglia_MSP = meta$Microglia_MSP,
                       Braak = meta$BraakStage,
                       Group = meta$Group,
                       row.names = meta$GSM)
  annoColors = list(Group = c(Control = "dodgerblue4" , AD = "chocolate1"),
                    Sex = c(F = "indianred4", M = "cornflowerblue"),
                    Age = AgefColor,
                    Neuron_CETs = Neuron_CETsColor,
                    NeuNall_MSP = c("chartreuse4","gray97","maroon"),
                    Oligo_MSP = c("chartreuse4","gray97","maroon"),
                    Microglia_MSP = c("chartreuse4","gray97","maroon"),
                    Braak = c("azure", "darkorchid4"))
  Plot <- pheatmap(data, angle_col = 90, border_color = NA,
                   #color = colorRampPalette(c("darkblue", "gold2"))(999),
                   scale = "none",
                   show_rownames = F,
                   show_colnames = T,
                   annotation_col = annoCol,
                   annotation_colors = annoColors,
                   main = title,
                   filename = paste0(ResultsPath, title, ".pdf"), useDingbats = F, width = 10, height = 8)
  
}
SignifCETsPeaksHyper <- group_results %>% filter(FDR < 0.05, logFC > 0) %>% select(matches("PeakName|logFC|FDR"))
SignifCETsPeaksHypo <- group_results %>% filter(FDR < 0.05, logFC < 0) %>% select(matches("PeakName|logFC|FDR"))

CETSTMM <- cpm(ADcountList, log = T)

SignifCETSTMMhyper <- CETSTMM[rownames(CETSTMM) %in% SignifCETsPeaksHyper$PeakName,]
SignifCETSTMMhypo <- CETSTMM[rownames(CETSTMM) %in% SignifCETsPeaksHypo$PeakName,]

PlotPeakHeatmap(SignifCETSTMMhyper, Metadata, title = "CETs_Hyper")
PlotPeakHeatmap(SignifCETSTMMhypo, Metadata, title = "CETs_Hypo")



####### Repeat after adjusting for MSPs ##############
#First - similar to CETS
MarziModelMatrixMSP_f <- model.matrix(as.formula("~Agef + NeuNall_MSPif + Group"), data = countMatrixFullAllCalled$Metadata)

ADcountList2_f <- estimateDisp(ADcountList, MarziModelMatrixMSP_f)

fitTMM_MSP_f <- glmQLFit(ADcountList2_f, MarziModelMatrixMSP_f)
qlf_groupMSP_f <- glmQLFTest(fitTMM_MSP_f, coef = "GroupAD")
group_resultsMSP_f <- topTags(qlf_groupMSP_f, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.)) 


groupResult_MSPAnno_f <- group_resultsMSP_f %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)

#Now all relevant MSPs as continuous
MarziModelMatrixMSP <- model.matrix(as.formula("~Agef + NeuNall_MSP + Microglia_MSP + Oligo_MSP + Group"), data = countMatrixFullAllCalled$Metadata)

ADcountList2 <- estimateDisp(ADcountList, MarziModelMatrixMSP)

fitTMM_MSP <- glmQLFit(ADcountList2, MarziModelMatrixMSP)
qlf_groupMSP <- glmQLFTest(fitTMM_MSP, coef = "GroupAD")
group_resultsMSP <- topTags(qlf_groupMSP, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.)) 


groupResult_MSPAnno <- group_resultsMSP %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)


#Heatmap of significant peaks
SignifMSPsPeaksHyper <- group_resultsMSP %>% filter(FDR < 0.05, logFC > 0) %>% select(matches("PeakName|logFC|FDR"))
SignifMSPsPeaksHypo <- group_resultsMSP %>% filter(FDR < 0.05, logFC < 0) %>% select(matches("PeakName|logFC|FDR"))

MSPsTMM <- cpm(ADcountList2, log = T)

SignifMSPsTMMhyper <- MSPsTMM[rownames(MSPsTMM) %in% SignifMSPsPeaksHyper$PeakName,]
SignifMSPsTMMhypo <- MSPsTMM[rownames(MSPsTMM) %in% SignifMSPsPeaksHypo$PeakName,]

PlotPeakHeatmap(SignifMSPsTMMhyper, Metadata, title = "MSPs_Hyper")
PlotPeakHeatmap(SignifMSPsTMMhypo, Metadata, title = "MSPs_Hypo")

#########################################

#Merge results from both analyses
CompareResultsDF <- merge(group_results, group_resultsMSP, by = "PeakName", suffixes = c("_CETs", "_MSP"), sort = FALSE)
CompareResultsDF$MethodSignif <- apply(CompareResultsDF %>% select(FDR_CETs, FDR_MSP), 1, function(x){
  if(x[1] < 0.05 & x[2] < 0.05){
    "Both"
  } else if (x[1] > 0.05 & x[2] > 0.05) {
    "NS"
  } else if(x[1] < 0.05 & x[2] > 0.05 ){
    "CETs"
  } else if(x[1] > 0.05 & x[2] < 0.05){
    "MSP"
  }
})

CompareResultsDF$MethodSignif <- factor(CompareResultsDF$MethodSignif, levels = c("NS", "CETs", "MSP", "Both"))


#Combijne all four together


temp <- merge(group_resultsMSP_f %>% select(PeakName, logFC, PValue, FDR),
              group_resultsMSP %>% select(PeakName, logFC, PValue, FDR),
              by = "PeakName", suffixes = c("_MSPneuronF", "_MSPall"))


AllThreeCombinedSup <- merge(group_results %>% select(PeakName, logFC, PValue, FDR), temp, by = "PeakName")
names(AllThreeCombinedSup)[names(AllThreeCombinedSup) %in% c("logFC", "PValue", "FDR")] <- paste0(names(AllThreeCombinedSup)[names(AllThreeCombinedSup) %in% c("logFC", "PValue", "FDR")], "_CETs")

AllThreeCombined <- pivot_longer(AllThreeCombinedSup, col = -matches("PeakName|_CETs"),
                      names_to = c(".value", "Method"), names_pattern = "(.*)_(.*)") %>% data.frame()

AllThreeCombined$MethodSignif <- apply(AllThreeCombined %>% select(FDR_CETs, FDR, Method), 1, function(x){
  if(x[1] < 0.05 & x[2] < 0.05){
    "Both"
  } else if (x[1] > 0.05 & x[2] > 0.05) {
    "NS"
  } else if(x[1] < 0.05 & x[2] > 0.05 ){
    "CETs"
  } else if(x[1] > 0.05 & x[2] < 0.05){
    "Alternative"
  }
})

AllThreeCombined$Method <- factor(AllThreeCombined$Method, levels = c("CETs", "MSPneuronF", "MSPall"))
AllThreeCombined$MethodSignif <- factor(AllThreeCombined$MethodSignif, levels = c("NS", "CETs", "Alternative", "Both"))


CorDF <- data.frame(Method = c("MSPneuronF", "MSPall"),
                    x1 = 0.2, y1 = 2.5,
                    x2 = 0.2, y2 = 2.2)

CorDF$CorAll <- sapply(CorDF$Method, function(method){
  temp <- cor.test(~logFC_CETs + logFC, data = AllThreeCombined %>% filter(Method == as.character(method)))
  round(temp$estimate, digits = 2)
})

CorDF$CorSignif <- sapply(CorDF$Method, function(method){
  temp <- cor.test(~logFC_CETs + logFC, data = AllThreeCombined %>% filter(Method == as.character(method), FDR_CETs < 0.05))
  round(temp$estimate, digits = 2)
})

CorDF %<>% mutate(Text = paste0("'r'[italic('All peaks')] ", "*", "' = '", "*", CorAll))
CorDF %<>% mutate(Text2 = paste0("'r'[italic('CETs significant peaks')] ", "*", "' = '", "*", CorSignif))

ggplot(AllThreeCombined, aes(logFC_CETs, logFC)) +
  theme_classic() +
  theme(legend.background = element_blank()) +
  labs(x = "logFC_AD (CETs)", y = "logFC_AD (MSP)") +
  geom_bin2d(bins = 100) +
  #geom_point(aes(color = MethodSignif)) +
  geom_point(data = AllThreeCombined %>% filter(FDR_CETs < 0.05), color = "orange", alpha = 0.2) +
  # geom_point(data = AllFourCombined %>% filter(FDR_MSPneuronalAsFactor < 0.05, FDR_CETS > 0.05), color = "darkolivegreen4") +
  # geom_point(data = AllFourCombined %>% filter(FDR_MSPneuronalAsFactor < 0.05, FDR_CETS < 0.05), color = "purple") +
  #scale_fill_manual(values = c("grey80", "orange",  "darkolivegreen4", "purple")) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_text(data = CorDF, aes(x = x1, y = y1, label = Text), parse = T, hjust = 0) +
  geom_text(data = CorDF, aes(x = x2, y = y2, label = Text2), parse = T, hjust = 0) +
  facet_wrap(~Method, ncol = 2)

ggsave(paste0("MethodComparison2.pdf"), device = "pdf", width = 12, height = 4, dpi = 300, path = ResultsPath, useDingbats = F)

VennData <- AllThreeCombined %>% filter(MethodSignif != "NS")  %>%  select(FDR_CETs, FDR, Method, PeakName)
VennData$Signif_CETs <- sapply(VennData$FDR_CETs, function(x){
  if(x < 0.05){
    1
  } else {
    0
  }
})

VennData$Signif_Alternative <- sapply(VennData$FDR, function(x){
  if(x < 0.05){
    1
  } else {
    0
  }
})

VennData <- pivot_wider(VennData %>% select(-FDR), names_from = Method,
                        values_from =  Signif_Alternative,
                        names_prefix = "Signif_", values_fill = list(Signif_Alternative = 0)) %>% data.frame()


pdf(paste0(ResultsPath,"VennDiagram.pdf"),width = 6, height = 4.5, useDingbats = F)
venn(VennData  %>% select(matches("Signif")), box = F, opacity = 0.4, ilcs = 1.5, ilabels = T, zcolor = "style")
closeDev()



#Find Overlaps
SignifCETS <- group_results %>% filter(FDR < 0.05)
SignifMSP <- group_resultsMSP %>% filter(FDR < 0.05)



UniquePeakMSP <- SignifMSP %>% filter(!PeakName %in% SignifCETS$PeakName)
UniquePeakMSP <- groupResult_MSPAnno %>% filter(PeakName %in% UniquePeakMSP$PeakName) %>% select(PeakName, logFC, logCPM, PValue, FDR, symbol, Peak_Gene, GeneAnnoType, Peak.width, Peak.Location)

UniquePeakCETS <- SignifCETS %>% filter(!PeakName %in% SignifMSP$PeakName)


AllThreeCombinedSup <- merge(HTseqCounts %>% select(Geneid, CHR, START, END), AllThreeCombinedSup,
                             by.x = "Geneid", by.y = "PeakName") %>% data.frame %>% arrange(FDR_CETs)


########### Check the impact of cell correction ###############
#Run without cell correction
MarziModelMatrixNoCell <- model.matrix(as.formula("~Agef + Group"), data = countMatrixFullAllCalled$Metadata)

ADcountList2_No <- estimateDisp(ADcountList, MarziModelMatrixNoCell)

fitTMM_No <- glmQLFit(ADcountList2_No, MarziModelMatrixNoCell)
qlf_group_No <- glmQLFTest(fitTMM_No, coef = "GroupAD")
group_results_No <- topTags(qlf_group_No, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.)) 


groupResult_Anno_No <- group_results_No %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)


#Run with randomly assigning CETs inside a group 

MetaDummy <- countMatrixFullAllCalled$Metadata
MetaDummy$CETSif[MetaDummy$Group == "AD"] <- sample(MetaDummy$CETSif[MetaDummy$Group == "AD"], 24, replace = F)
MetaDummy$CETSif[MetaDummy$Group == "Control"] <- sample(MetaDummy$CETSif[MetaDummy$Group == "Control"], 23, replace = F)

MarziModelMatrixRandomCell <- model.matrix(as.formula("~Agef + CETSif + Group"), data = MetaDummy)

ADcountList2_Random <- estimateDisp(ADcountList, MarziModelMatrixRandomCell)

fitTMM_Random <- glmQLFit(ADcountList2_Random, MarziModelMatrixRandomCell)
qlf_group_Random <- glmQLFTest(fitTMM_Random, coef = "GroupAD")
group_results_Random <- topTags(qlf_group_Random, n = Inf) %>% data.frame() %>% mutate(PeakName = rownames(.)) 

groupResult_Anno_Random <- group_results_Random %>%
  AnnotDESeqResult(CountAnnoFile = AllCalledData$countsMatrixAnnot, by.x = "PeakName", by.y = "PeakName") %>% arrange(FDR)


temp <- merge(group_results_No %>% select(PeakName, logFC, PValue, FDR),
              group_results_Random %>% select(PeakName, logFC, PValue, FDR),
              by = "PeakName", suffixes = c("_NoCellAdjustment", "_Shuffled"))


CombinedAll <- merge(AllThreeCombinedSup, temp, by.x = "Geneid", by.y = "PeakName")
write.table(CombinedAll, paste0(ResultsPath, "Supplementary tableS1.tsv"), row.names = F, col.names = T, sep = "\t")


CombinedAll$CETs_Signif <- sapply(CombinedAll$FDR_CETs, function(x){
  if(x < 0.05){
    "CETs_Signif"
  } else {
    "CETS_NS"
  }
})

CombinedAll %<>% mutate(Direction = CETs_Signif)
CombinedAll$Direction[CombinedAll$logFC_CETs < 0 & CombinedAll$FDR_CETs < 0.05] <- "Down"
CombinedAll$Direction[CombinedAll$logFC_CETs > 0 & CombinedAll$FDR_CETs < 0.05] <- "Up"


CombinedAll_toPlot <- CombinedAll %>% select(logFC_CETs, logFC_NoCellAdjustment, logFC_Shuffled,
                                             logFC_MSPneuronF, logFC_MSPall, CETs_Signif, Direction) %>% mutate(CETs_Signif = factor(CETs_Signif))

CombinedAll_Sub <- rbind(CombinedAll_toPlot %>% filter(CETs_Signif == "CETs_Signif"), CombinedAll_toPlot %>% filter(CETs_Signif == "CETS_NS") %>%
                           .[sample(1:nrow(.), 10000, replace = F),]) 



Matrix <- ggpairs(CombinedAll_toPlot, columns = grep("logFC", names(CombinedAll_toPlot)),
                  lower = list(continuous = "density", mapping = aes_string(fill = "..level..", group = "Direction")),
                  upper = list(continuos = ggally_cor,  mapping = aes(color = CETs_Signif)),
                  diag = list(continuous = "blankDiag")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom")

for(j in 1:4){
  for(i in (j+1):5){
    Matrix[i,j] <- Matrix[i,j] +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#726F6F") +
      scale_color_manual(values = c("#2B5B83", "#FFA500"))
  }
}

for(i in 1:4){
  for(j in (i+1):5){
    Matrix[i,j] <- Matrix[i,j] +
      scale_color_manual(values = c("#2B5B83", "#FFA500"))
  }
}


Matrix2 <- ggpairs(CombinedAll_toPlot %>% filter(CETs_Signif == "CETs_Signif"), columns = grep("logFC", names(CombinedAll_toPlot)),
                   lower = list(continuous = "density", mapping = aes_string(fill = "..level..", group = "Direction")),
                   upper = list(continuos = ggally_cor,  mapping = aes_string(color = "CETs_Signif")),
                   diag = list(continuous = "blankDiag")) +
  theme_bw() +
  theme(panel.grid = element_blank()) 

for(j in 1:4){
  for(i in (j+1):5){
    Matrix2[i,j] <- Matrix2[i,j] +
      scale_fill_gradient(low = "orange", high = "white")
  }
}


ggsave(paste0("MethodComparisonAll_EachOther.pdf"), Matrix, device = "pdf", width = 12, height = 6, dpi = 300, path = ResultsPath, useDingbats = F)
ggsave(paste0("MethodComparisonAll_EachOther_SignifOnly.pdf"), Matrix2, device = "pdf", width = 12, height = 6, dpi = 300, path = ResultsPath, useDingbats = F)


CorAll_ns <- CombinedAll %>%  filter(FDR_CETs >= 0.05) %>%  select(matches("logFC")) %>% cor
CorAll_CETsSignif  <- CombinedAll %>% filter(FDR_CETs < 0.05) %>% select(matches("logFC")) %>% cor

rm(list = ls(pat = "^fit|ADcountList|AnnoFile|qlf", temp))
save.image(paste0(ResultsPath, Cohort, ".Rdata"))
saveRDS(groupResult_Anno, file = paste0(ResultsPath, Cohort, "edgeR_CETs.Rds"))
saveRDS(groupResult_MSPAnno, file = paste0(ResultsPath, Cohort, "edgeR_MSP.Rds"))
saveRDS(groupResult_MSPAnno_f, file = paste0(ResultsPath, Cohort, "edgeR_MSPneuronAsFactor.Rds"))
saveRDS(groupResult_Anno_No, file = paste0(ResultsPath, Cohort, "edgeR_NoCellCorretction.Rds"))
saveRDS(groupResult_Anno_Random, file = paste0(ResultsPath, Cohort, "edgeR_CETsRandom.Rds"))
saveRDS(CombinedAll, file = paste0(ResultsPath, Cohort, "AllfiveModels.Rds"))

closeDev()
