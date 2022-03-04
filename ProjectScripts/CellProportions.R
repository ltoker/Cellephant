source("/data/Rprojects/GeneralScripts/generalFunc.R")
source("ProjectScripts/ProjectFunctions.R")
plotMA = DESeq2::plotMA

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

annoFileCollapsed <- GetGenomeAnno(genome = "hg19")
BroadPeak <-  read.table("CellTypeH3K27ac/data/H3K27ac_consensus_signal_celltype_homogenate.saf", header = T, sep = "\t") %>% as(., "GRanges")

BroadPeak <-  mergeByOverlaps(annoFileCollapsed, BroadPeak, maxgap = 0, type = "any", select = "all") %>% data.frame()

HTseqCountsCellTypes <- read.table("CellTypeH3K27ac/data/H3K27ac_featureCounts_Celltype_homogenate.csv", header = T, sep = ",")
MetaCellTypes <- read.table("CellTypeH3K27ac/meta/MSSM_U01MH103392_EpiMap_Metadata_ChIPseq_August2016Release.csv", header = T, sep = ",")
MetaCellTypes2 <- read.table("CellTypeH3K27ac/meta/MSSM_U01MH103392_EpiMap_Metadata_clinical.csv", header = T, sep = ",")

MetaCellTypes <- merge(MetaCellTypes, MetaCellTypes2, by = "Individual_ID") %>% select(-Assay, -File_Name, -RunType, -Grant, -StudyName, -Organism) %<>% filter(HistoneMark == "H3K27ac", BrainRegion == "DLPFC")

#Filter the relevant samples from the count matrix

HTseqCountsCellTypes %<>% select(c("X", as.character(MetaCellTypes$Sample_ID)))
countMatrix <- as.matrix(HTseqCountsCellTypes[,-1])
rownames(countMatrix) <- as.character(HTseqCountsCellTypes$X)

MetaCellTypes$RiP <- apply(HTseqCountsCellTypes[,-1], 2, sum)
MetaCellTypes$CellType2 <- sapply(MetaCellTypes$CellType, function(x){
  if(grepl("-", x)){
    "Glia"
  } else {
    "Neuron"
  }
})

Model = as.formula("~CellType2 + ChromatinAmount + LibraryBatch + Sex + Hemisphere + AgeDeath + pH")

DESeqDS <- DESeqDataSetFromMatrix(countData = countMatrix, colData = MetaCellTypes, design = Model)
DESeqOut <- DESeq(DESeqDS)

DESeqOutCellTypeResults <- results(DESeqOut, name = "CellType2_Neuron_vs_Glia", independentFiltering = T, alpha = 0.05) %>% data.frame %>% mutate(PeakName = rownames(.))
DESeqOutCellTypeResultsAnno <- merge(DESeqOutCellTypeResults, BroadPeak, by.x = "PeakName", by.y = "GeneID")

CellTypePeaks <- DESeqOutCellTypeResultsAnno %>% filter(abs(log2FoldChange) > 2, baseMean > 1000)
write.table(CellTypePeaks, "CellTypeH3K27ac/CellTypeH3K27ac.tsv", sep = "\t", row.names = F, col.names = T)



