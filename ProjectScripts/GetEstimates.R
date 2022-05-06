if(!"devtools" %in%  rownames(installed.packages())){
  BiocManager::install("devtools")
}
library(devtools)

source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")
source("ProjectScripts/ProjectFunctions.R")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")
packageF("pheatmap")
if(!"ermineR" %in% rownames(installed.packages())){
  install_github("https://github.com/PavlidisLab/ermineR")
}

library(ermineR)

#Get transcriptome annotations
AnnoLoc = "Annotations"
AssemblyFilename = "gencode.v35.annotation.gff3.gz"

source("ProjectScripts/Annotations.R")

source_url("https://github.com/ltoker/GeneralRscripts/blob/main/Cell_type_PCA.R?raw=T")

packageF("BisqueRNA")

ResultsPath = "EstimatesResults"
if(!ResultsPath %in% list.dirs(full.names = F, recursive = T)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

CountMatrixCombined <- readRDS("Data/txi_combined.Rds") %>% .$counts

rownames(CountMatrixCombined) <- sapply(rownames(CountMatrixCombined), function(x){
  strsplit(x, "\\.")[[1]][1]
})

mitoGenes <- geneNames %>% filter(seqnames == "chrM")

Meta <- readRDS("Data/Metadata.Rds")
NeuronalGenes <- read.table("Data/HumanNeuronVSglia.tsv", header = T, sep = "\t")
Meta$NeuExpRegion <- "Cortex"

region = "Cortex"
CellType_genes <- GetMarkers(region)
CellType_genes$Neurons_Genes <- NeuronalGenes$HugoName

names(CellType_genes) <- paste0("MGP_", names(CellType_genes))
names(CellType_genes) <- sapply(names(CellType_genes), function(x){
  gsub("Oligo", "Oligodendrocyte", x)
})

RNAestimates <- sapply(grep("RNAseq", names(Meta), value = T), function(SampleIDCol){
  if(grepl("1", SampleIDCol)) {
    Prefix = "T1_RNA_"
  } else {
    Prefix = "T2_RNA_"
  }
  
  # Match count matrix names to metadata names
  countMatrix <- CountMatrixCombined[,Meta[[SampleIDCol]]] 
  Metadata = Meta %>%
    mutate(Batch = 1, Group = Neuromics_Classification )
  Metadata$Filename <- Metadata[,SampleIDCol]
  Study = strsplit(SampleIDCol, "_")[[1]][3]
  
  #Remove genes with maximal count < 2 and mitochondrial genes
  Max <- apply(countMatrix, 1, max)
  
  countMatrix <- countMatrix[Max > 10,]
  
  countMito <- countMatrix[rownames(countMatrix) %in% mitoGenes$gene_id2,]
  countMitoSum <- countMito %>% apply(2, sum)
  TotLibSize <- countMatrix  %>% apply(2, sum) 
  ShortGeneSum <-  countMatrix %>% data.frame %>% filter(rownames(.) %in% (geneNames %>% mutate(size = end - start) %>% arrange(size) %>% filter(size < 300) %>% .$gene_id2)) %>% apply(2, sum)
  
  
  MitoCountFiltered <- countMatrix[!rownames(countMatrix) %in% mitoGenes$gene_id2,]
  MitoFiltCountSum = apply(MitoCountFiltered, 2, sum)
  Metadata$MitoCount <- countMitoSum[match(Metadata[[SampleIDCol]], names(countMitoSum))]
  
  #Look at the most highly expressed genes
  TopFiveProportion <- sapply(names(MitoFiltCountSum), function(sbj){
    SubMatrix = data.frame(genes = rownames(MitoCountFiltered), Counts = MitoCountFiltered[,sbj])
    TopFive = SubMatrix %>% arrange(desc(Counts)) %>% head(5)
    TopFive %<>%  mutate(Proportion = Counts/MitoFiltCountSum[sbj])
    Genes <- geneNames[match(TopFive$genes, geneNames$gene_id2),]  %>% select(gene_name, gene_type)
    Genes$Filename = sbj
    temp <- cbind(Genes, TopFive)
    names(temp)[names(temp) == "genes"] <- "ensemblID"
    temp
  }, simplify = FALSE) %>% rbindlist()
  
  
  TopFiveSum <- TopFiveProportion %>% group_by(Filename) %>%
    summarise(TotProp = sum(Proportion)) %>%
    data.frame %>% arrange(TotProp)
  
  
  TopFiveGeneFreq <- TopFiveProportion %>% group_by(ensemblID) %>%
    summarise(n = n()) %>%
    data.frame
  
  TopFiveGeneFreq <- merge(TopFiveGeneFreq, geneNames,
                           by.x = "ensemblID", by.y = "gene_id2", all.x = T, all.y = F, sort = F) 
  
  TopFiveGeneFreq$ensemblID2 <- sapply(TopFiveGeneFreq$ensemblID,  function(x){
    strsplit(x, "\\.")[[1]][1]
  })
  
  TopFiveGeneFreq %<>% mutate(ensemblID2 = paste0(ensemblID2,
                                                  " (", gene_name, ", ", n, ")"))
  TopFiveGeneFreq %<>% arrange(desc(n))
  
  TopFiveProportion$ensemblID2 <- TopFiveGeneFreq$ensemblID2[match(TopFiveProportion$ensemblID,
                                                                   TopFiveGeneFreq$ensemblID)]
  
  #Order subjects based on library size after filtering of mtDNA genes
  TopFiveProportion <- merge(TopFiveProportion, Metadata %>%
                               select(Filename, Batch, Cohort,
                                      DV200, RIN, MitoCount),
                             by = "Filename", sort = F)
  
  TopFiveProportion$Filename <- factor(TopFiveProportion$Filename,
                                       levels = TopFiveSum$Filename)
  TopFiveProportion$ensemblID2 <- factor(TopFiveProportion$ensemblID2,
                                         levels = unique(as.character(TopFiveGeneFreq$ensemblID2)))
  
  
  
  TopFivePlot <- ggplot(TopFiveProportion %>% filter(Proportion > 0.006), aes(Filename, Proportion, fill = ensemblID2)) +
    theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, size = 8)) +
    labs(y = "Proportion of reads", x = "Sample", title = "Top five genes with the highest read count") + 
    scale_fill_manual(values = c(MoviePalettes$MoonRiseKingdomColors[3:9], gray.colors(6)),
                      name = "GeneID (n)") +
    geom_bar(stat = "identity")
  
  TopFivePropPlot <- ggarrange(TopFivePlot,
                               ggplot(TopFiveProportion %>% group_by(Filename) %>%
                                        summarise(TopFiveProp = sum(Proportion),
                                                  DV200 = mean(DV200),
                                                  mtDNAcount = mean(MitoCount)) %>% data.frame(),
                                      aes(DV200, TopFiveProp, color = mtDNAcount)) +
                                 geom_point(),
                               nrow = 2, heights = c(1.5, 1))
  
  #Get the common genes with the highest count in majority of the samples and remove them from the count matrix
  CommonTopGenes <- TopFiveGeneFreq[TopFiveGeneFreq$n > 0.5*ncol(MitoCountFiltered),] %>% filter(!duplicated(ensemblID))
  
  CommonTopGenesSum <- apply(MitoCountFiltered[rownames(MitoCountFiltered )%in% CommonTopGenes$ensemblID,], 2, sum)
  
  countMatrixFiltered <- MitoCountFiltered[!rownames(MitoCountFiltered) %in% as.character(CommonTopGenes$ensemblID),]
  
  # Remove gene with less than 5 counts in  > 80% of the samples
  ZeroCount <- apply(countMatrixFiltered, 1, function(x){
    sum(x < 5)
  })
  
  countMatrixFiltered <- countMatrixFiltered[ZeroCount < 0.8*ncol(countMatrixFiltered),]
  
  
  #Create log2 CPM matrix after removal of mitochondria-encoded genes
  cpmMatrixFiltered <- Count2CPM(countMatrixFiltered) %>% data.frame()
  cpmMatrixFiltered <- apply(cpmMatrixFiltered, c(1,2), function(x) log2(x+1)) %>% data.frame()
  cpmMatrixFiltered <- cbind(rownames(countMatrixFiltered), cpmMatrixFiltered)
  colnames(cpmMatrixFiltered)[1] <- "genes"
  
  
  #Add gene symbols
  GeneSymbolAll <- data.frame(GeneSymbol = geneNames$gene_name[match(rownames(cpmMatrixFiltered), geneNames$gene_id2)],
                              Probe = rownames(cpmMatrixFiltered),
                              ensemblID = rownames(cpmMatrixFiltered))
  
  
  ExpDataCPM <- cbind(GeneSymbolAll, cpmMatrixFiltered[-1])
  
  
  studyFinal <- PreProccessRNAseq(Metadata = Metadata, expData = ExpDataCPM,
                                  SexCol = "Sex", Combat = FALSE, resultsPath = ResultsPath)
  
  #### Calculation of MGPs
  
  #Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
  CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]
  
  #Bootstrap with replacement the samples (90% of the samples)
  SampleNames <- as.character(studyFinal$Metadata$Filename)
  
  PCAresults <- sapply(paste0("Boot_", 1:100), function(boot){
    BootSamples <- sample(SampleNames, 0.9*length(SampleNames), replace = F)
    dataSub <- studyFinal$ExpHigh %>% select(c("GeneSymbol", BootSamples))
    PCA_genes_All_based(dataset_id=Study,
                        dataset=dataSub,
                        CellType_genes=CellType_genes,
                        contName = "SL",SampleReg = "SL",
                        NoiseThershold = studyFinal$NoiseThreshold)
  }, simplify = F)
  
  PCA_resultsMean <- sapply(names(PCAresults[[1]]$modified), function(celltype){
    temp <- data.frame(CommonName = names(PCAresults[[1]]$modified[[celltype]]$x[,1]),
                       Rot = PCAresults[[1]]$modified[[celltype]]$x[,1])
    for(i in 2:length(PCAresults)){
      temp <- merge(temp, data.frame(CommonName = names(PCAresults[[i]]$modified[[celltype]]$x[,1]),
                                     Rot = PCAresults[[i]]$modified[[celltype]]$x[,1]), by = "CommonName", all = TRUE, sort = F)
    }
    names(temp)[2:ncol(temp)] <- paste0("Rot", c(1:c(ncol(temp)-1)))
    temp$MeanRot <- rowMeans(temp[-1], na.rm = T)
    temp
  }, simplify=FALSE)
  
  
  #Add estimation to Metadata 
  AllEstimates <- lapply(PCA_resultsMean, function(x){
    x$MeanRot
  }) %>% do.call(cbind, .) %>% data.frame()
  AllEstimates$CommonName <- PCA_resultsMean[[1]]$CommonName
  
  #Keep only major cell types
  AllEstimates <- AllEstimates %>% select(CommonName, MGP_Astrocyte_Genes, MGP_Endothelial_Genes,
                                          MGP_Microglia_Genes, MGP_Oligodendrocyte_Genes, MGP_Neurons_Genes)
  names(AllEstimates) <- sapply(names(AllEstimates), function(x){
    gsub("_Genes", "", x)
  })
  
  names(AllEstimates)[-1] <- paste0(Prefix, names(AllEstimates)[-1])
  studyFinal$Metadata <- merge(studyFinal$Metadata, AllEstimates, by = "CommonName", sort = F)

  #Calculation of Bisque estimates
  ordered.celltypes <- c("neurons", "oligodendrocytes", "microglia", "endothelial", "astrocytes")
  sc.eset <- readRDS("Data/Bisque//cortex_sc_counts.Rds")
  sc.meta <- read_csv("Data/Bisque/GSE67835_metadata.csv") %>%
    select(Age = AGE, Cell_type, subject_id=experiment_sample_name, sample_id=`Sample Name`) %>% .[grepl("postnatal", .$Age),] %>%
    mutate(Age=as.numeric(gsub(" years", '', gsub("postnatal ", '', Age))))
  sc.meta %<>% filter(Age > 30, Cell_type %in% ordered.celltypes)
  sc.eset <- sc.eset[,sc.meta$sample_id]
  
  markers <- c(read_csv("Data/Bisque/velmeshev_2019_markers.csv")$gene_name,
               read_csv("Data/Bisque/kelley_2018_top100_markers.csv")$gene_name) %>% unique
  
  countMatrixFilteredDF <- countMatrixFiltered %>% data.frame()
  countMatrixFilteredDF$gene_name <- geneNames$gene_name[match(rownames(countMatrixFilteredDF), geneNames$gene_id2)]
  
  SortSample <- names(countMatrixFilteredDF)[1]
  
  countMatrixFilteredDF %<>% select(gene_name, matches("SL")) %>% arrange(desc(across(contains(SortSample)))) %>%
    distinct(gene_name, .keep_all=TRUE)
  rownames(countMatrixFilteredDF) <- countMatrixFilteredDF$gene_name
  countMatrixFilteredDF %<>% select(-gene_name) %>% as.matrix
  
  bulk.eset <- ExpressionSet(assayData=countMatrixFilteredDF)
  estimates <- ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=markers, use.overlap=FALSE)$bulk.prop
  estimates <- t(estimates) %>% as_tibble(rownames="sample_id") %>%
    select(sample_id, Bisque_Astrocytes=astrocytes, Bisque_Endothelial=endothelial,
           Bisque_Microglia=microglia, Bisque_Oligodendrocytes=oligodendrocytes, Bisque_Neurons=neurons)
  names(estimates)[-1] <- paste0(Prefix, names(estimates)[-1])
  
  #MAke all cell type names in singular
  studyFinal$Metadata <- merge(studyFinal$Metadata, estimates, by.x = "CommonName", by.y = "sample_id", sort = F, all.x = T, all.y = T)
  names(studyFinal$Metadata) <- sapply(names(studyFinal$Metadata), function(x){
    gsub("s$", "", x)
  })
  
  estimates2 <- ReferenceBasedDecomposition(bulk.eset, sc.eset, use.overlap=FALSE)$bulk.prop
  estimates2 <- t(estimates2) %>% as_tibble(rownames="sample_id") %>%
    select(sample_id, Bisque_Astrocytes=astrocytes, Bisque_Endothelial=endothelial,
           Bisque_Microglia=microglia, Bisque_Oligodendrocytes=oligodendrocytes, Bisque_Neurons=neurons)
  names(estimates2)[-1] <- paste0(Prefix, names(estimates2)[-1])

  list(TopGenesPlot  = TopFivePropPlot, Estimates = studyFinal$Metadata, cpmMatrix = studyFinal$ExpHigh,
       countMatrixFiltered = countMatrixFiltered, estimates2 = estimates2)
}, simplify = F)

RNAestimates$RNAseq_id_ParkOme1$TopGenesPlot

RNAestimatesCombined <- merge(RNAestimates$RNAseq_id_ParkOme1$Estimates %>% select(-CommonName),
                              RNAestimates$RNAseq_id_ParkOme2$Estimates %>% select(matches("ParkOme1|T._RNA")),
                              by = "RNAseq_id_ParkOme1")

write.table(RNAestimatesCombined, paste0(ResultsPath, "RNAestimatesMGP_Bisque.tsv"),
            row.names = F, col.names = T, sep = "\t")

#This part creates the input for BrainDecode Shiny app (Sutton et al. 2022)
sapply(names(RNAestimates), function(RNAseq){
  temp <- RNAestimates[[RNAseq]]$cpmMatrix %>% select(-GeneSymbol, -Probe)
  temp2 <- RNAestimates[[RNAseq]]$countMatrixFiltered
  temp2 <- temp2[rownames(temp2) %in% temp$ensemblID,]
  temp2CPM <- Count2CPM(temp2)
  temp2DF <- cbind(rownames(temp2CPM), data.frame(temp2CPM))
  RNAseq <- gsub("_id_ParkOme", "", RNAseq)
  names(temp2DF)[1] <- ""
  Samples <- list(A = c(1,2:14),
                  B = c(1, 14:26),
                  C = c(1, 26:38),
                  D = c(1, 38:50))
  lapply(Samples, function(Cols){
    write.table(temp2DF[,Cols], paste0(RNAseq, "_cpm_", Cols[2],
                                    "_", Cols[length(Cols)],  ".csv"),
                row.names = F, col.names = T, sep = ",")
  })
})

