packageF("DESeq2")
packageF("preprocessCore")
packageF("parallel")
packageF("tidyr")
packageF("pheatmap")
packageF("GenomicRanges")
packageF("annotatr")
packageF("readr")

packageF("TxDb.Hsapiens.UCSC.hg19.knownGene")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")

LocateEnhancers <- function(EnhancerFile = "EnhancerData/human_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt"){
  EnhanceExp <- read.table(EnhancerFile, header = T, sep = "\t", quote = "")
  SampleInfo <- read.table("EnhancerData/SampleDescription.txt", header = FALSE, sep = "\t", quote = "")
  BrainSamples <- SampleInfo[grepl("cortex|neuron", SampleInfo$V1),] %>% filter(!grepl("adrena|iPS", .$V1))
  BrainEnhanceExp <- EnhanceExp %>% select(c("Id", as.character(BrainSamples$V2)))
  BrainEnhanceExp$Max <- apply(BrainEnhanceExp %>% select(-Id), 1, max)
  BrainEnhanceExp %<>% filter(Max > 5)
  BrainEnhanceExp$CHR <- sapply(BrainEnhanceExp$Id, function(x) strsplit(as.character(x),  ":")[[1]][1])
  BrainEnhanceExp$START <- sapply(BrainEnhanceExp$Id, function(x){
    x = strsplit(as.character(x),  ":")[[1]][2]
    strsplit(x,"-")[[1]][1]
  })   
  BrainEnhanceExp$END <- sapply(BrainEnhanceExp$Id, function(x){
    x = strsplit(as.character(x),  ":")[[1]][2]
    strsplit(x,"-")[[1]][2]
  })   
  
  BrainEnhanceExp %>% as(.,"GRanges")
}

GetEnhancerData <- function(obj, EnhancerData = NULL, EnhancerFile = "EnhancerData/human_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt"){
  if(is.null(EnhancerData)){
    EnhancerData <- LocateEnhancers(EnhancerFile = EnhancerFile)
  }
  EnhacerLoc <- EnhancerData %>% data.frame() %>% .[1:3] %>% as(.,"GRanges")
  ObjDF <- obj %>% data.frame()
  PeakLoc <- ObjDF %>% .[1:3] %>% as(.,"GRanges")
  OverlapHits <- findOverlaps(query = EnhacerLoc, subject = PeakLoc, maxgap = 5, type = "any", select = "all")
  ObjDF$Enhancer <- "No"
  ObjDF$Enhancer[subjectHits(OverlapHits)] <- "Yes"
  return(ObjDF %>% as(.,"GRanges"))
}

GetGenomeAnno <- function(annots = c('hg19_basicgenes', 'hg19_genes_intergenic','hg19_genes_intronexonboundaries'),
                          genome = "hg19"){
  annoFile <- build_annotations(genome = genome, annotations = annots)
  
  #Remove replicated annotations due to isoform data
  annoFileCollapsed <- annoFile %>% data.frame 
  annoFileCollapsed %<>% mutate(UniqueRegion =  paste0(seqnames, ":", start, "-", end),
                                UniqueRegionGene = paste(UniqueRegion, symbol, sep = "_"),
                                UniqueRegionGeneType = paste(UniqueRegionGene, type, sep = "_"),
                                GeneAnnoType = paste(symbol, type, sep = "_")) %>%
    filter(!duplicated(UniqueRegionGeneType))
  
  
  annoFileCollapsed %<>% as(., "GRanges")
  
  #ADD enhancer information
  #annoFileCollapsed <- GetEnhancerData(annoFileCollapsed)  #this part is removed since for now the enhancer information is not beeing used.
  return(annoFileCollapsed)
}

GetCountMatrixHTseq <- function(countsDF, meta = Metadata, MetaSamleCol = "activemotif_id", countSampleRegEx = "^X", 
                                MetaCol = c("activemotif_id", "SampleID",  "rin", "condition", "sex", "age", "batch", "pm_hours", "library_size", 
                                            grep("Genes", names(Metadata), value = T), 
                                            "h3k27_gapdh", "h3k27_h3", "H3K27gapdh", "H3K27gapdh_Norm"),
                                HouseKeepRegEx = "^GAPDH_|^ACTB_|^UBC_",
                                OtherNormRegEx = NULL){
  names(countsDF)[1] <- "PeakName"
  #Remove counts mapped to contig sequences
  ChrmRM <- countsDF[grepl("GL|hs", countsDF$CHR),] %>% .$CHR %>% as.character %>% unique
  countsDF %<>% filter(!CHR %in% ChrmRM)
  
  ReadsPerChr <- sapply(grep(countSampleRegEx, names(countsDF), value = T), function(Subj){
    subData <- countsDF %>% select(c("CHR", "START", "END", Subj))
    subData %<>% mutate(PeakLength = END-START)
    names(subData)[grepl(Subj, names(subData))] <- "Counts"
    subData %>% group_by(CHR) %>%
      summarise(numPeak = n(), SumReads = sum(Counts),
                MeanLength = mean(PeakLength), MaxLength = max(PeakLength)) %>%
      mutate(GSM = Subj)
  }, simplify = F) %>% rbindlist()
  
  ReadsPerChr <- merge(ReadsPerChr, meta %>% select(MetaCol), by = MetaSamleCol)
  ReadsPerChr$CHR <- factor(ReadsPerChr$CHR, levels = c(as.character(c(1:22)), "X", "Y"))
  
  TotalSampleRead <- ReadsPerChr %>%  group_by(.dots = MetaSamleCol) %>% summarise(TotalCount = sum(SumReads)) %>% data.frame
  TotalSampleRead <- merge(TotalSampleRead, Metadata %>% select(MetaCol), by = MetaSamleCol)
  
  # TotalSampleRead %<>% mutate(FRiP = TotalCount/library_size,
  #                             Background = library_size - TotalCount)  
  
  #Get normalization factor
  # MedianTotalLibrarySize = median(TotalSampleRead$library_size, na.rm = T)
  #MedianRiP = median(TotalSampleRead$TotalCount, na.rm = T)
  #MedianBackground = median(TotalSampleRead$Background, na.rm = T)
  
  # TotalSampleRead %<>% mutate(NormFactAll = library_size/MedianTotalLibrarySize,
  #                             RiP_NormAllCount  = TotalCount/NormFactAll,
  #                             NormFactBackground = Background/MedianBackground,
  #                             RiP_NormBackground = TotalCount/NormFactBackground)
  
  
  
  # Annotate the peaks
  PeakLocation <- countsDF %>% select(CHR, START, END, PeakName) %>% as(., "GRanges")
  seqlevelsStyle(PeakLocation) <- "UCSC"
  PeakLocation <- keepSeqlevels(PeakLocation, paste0("chr",c(1:22, "X", "Y")), pruning.mode = "coarse")
  
  PeakAnnoFile <- mergeByOverlaps(annoFileCollapsed, PeakLocation, maxgap = 0, type = "any", select = "all") %>% data.frame()
  names(PeakAnnoFile)[1:4] <- c("Region.CHR", "Region.START", "Region.END", "Region.width")
  PeakAnnoFile %<>% select(-matches("strand|_id|^id|^PeakName|annoFileCollapsed"))
  
  names(PeakAnnoFile)[grepl("PeakLocation", names(PeakAnnoFile))] <- c("Peak.CHR", "Peak.START", "Peak.END", "Peak.width", "PeakName")
  PeakAnnoFile %<>% mutate(Peak.Location = paste0(Peak.CHR, ":", Peak.START, "-", Peak.END))
  
  
  countsMatrix <- as.matrix(countsDF %>% select(matches(countSampleRegEx)))
  
  rownames(countsMatrix) <- countsDF$PeakName %>% as.character
  
  countsMatrixAnnot <- merge(PeakAnnoFile, countsMatrix, by.x = "PeakName", by.y = "row.names", all.x = F, all.y = T)
  
  # Get the peaks annotated to house keeping genes
  if(is.null(OtherNormRegEx)){
    HouseKeeping <- countsMatrixAnnot %>% .[grepl(HouseKeepRegEx, .$GeneAnnoType),]
  } else {
    HouseKeeping <- countsMatrixAnnot %>% .[grepl(paste0(HouseKeepRegEx, "|", OtherNormRegEx), .$GeneAnnoType),]
  }
  HouseKeeping$Mean <- apply(HouseKeeping %>% select(matches(countSampleRegEx)), 1, mean)
  HouseKeeping %<>% filter(!duplicated(.$PeakName), Mean > 150) 
  
  # Get the ratio to the mean count for the house keeping genes
  HouseKeepingRatio <- sapply(HouseKeeping$PeakName, function(Peak){
    PeakData <- HouseKeeping %>% filter(PeakName == Peak)
    temp <- PeakData %>% select(matches(countSampleRegEx)) %>% unlist
    temp/PeakData$Mean
  }) %>% data.frame
  
  names(HouseKeepingRatio) <- HouseKeeping$PeakName
  names(HouseKeepingRatio) <- paste0(HouseKeeping$symbol,
                                     sapply(names(HouseKeepingRatio), function(x) gsub("(.*)?Peak", "_", x ,ignore.case = T)))
  
  # Get the mean ratio of the house keeping genes
  HouseKeepingRatio$MeanRatioOrg <- apply(HouseKeepingRatio[,grepl(HouseKeepRegEx, names(HouseKeepingRatio))], 1, mean)
  HouseKeepingRatio$MeanRatioAll <- apply(HouseKeepingRatio, 1, mean)
  
  HouseKeepingRatio$SampleName <- rownames(HouseKeepingRatio)
  

  HouseKeepingRatio[[MetaSamleCol]]<- sapply(HouseKeepingRatio$SampleName, function(x){
    gsub("X", "", x)}
  )
  
  SampleInfo = merge(TotalSampleRead, HouseKeepingRatio, by = MetaSamleCol)
  SampleInfo %<>% mutate(RiP_NormMeanRatioOrg = TotalCount/MeanRatioOrg,
                         RiP_NormMeanRatioAll = TotalCount/MeanRatioAll)
  
  HouseKeepingRatioPlot <- merge(HouseKeepingRatio, TotalSampleRead %>% select(c(MetaSamleCol, "TotalCount")), by = MetaSamleCol)
  MeasureCor <- cor(HouseKeepingRatioPlot %>% select(-matches("id|Sample|GSM")), method = "spearman")
  diag(MeasureCor) <- NA
  
  Plot <- pheatmap(MeasureCor, angle_col = 90, na_col = "white", display_numbers = T,
                   filename = paste0(ResultsPath, "HouseKeepingCor", Cohort, ".pdf"), width = 8, height = 8 )
  return(list(countsDF = countsDF,
              countsMatrixAnnot = countsMatrixAnnot,
              SampleInfo = SampleInfo,
              MeasureCor = MeasureCor,
              HeatMap = Plot))
}

GetCollapsedMatrix <- function(countsMatrixAnnot, collapseBy, FilterBy, meta = Metadata, normCol = NULL, title = NULL, samples = "All", CorMethod = "pearson", countSampleRegEx = "^X", MetaSamleCol = "activemotif_id", MetaSamleIDCol = "SampleID", groupCol = "condition"){
  if(is.null(title)){
    title = paste0("Sample correlation (", FilterBy, ")")
  }
  
  Data <- countsMatrixAnnot %>% select(matches(paste0(collapseBy, "|", countSampleRegEx))) %>% group_by(.dots = collapseBy) %>% summarise_if(is.numeric, sum, na.rm = TRUE) %>% data.frame
  names(Data)[grepl(collapseBy, names(Data))] <- "PeakName"
  subData <- Data[grepl(FilterBy, Data$PeakName),]
  rownames(subData) <- subData$PeakName
  subData <- subData[-1]
  
  subData <- apply(subData, c(1, 2), function(x) as.integer(round(x, digits = 0)))
  
  if(!"All" %in% samples){
    subData = subData[,colnames(subData) %in% paste0("X", samples)]
    meta = meta %>% filter(activemotif_id %in% samples)
  }
  if(!is.null(normCol)){
    CPMdata <- sapply(1:ncol(subData), function(i){
      subData[,i]/meta[[normCol]][i]
    })
    colnames(CPMdata) <- colnames(subData)
  } else {
    CPMdata <- Count2CPM(subData)
  }

  SampleCor <- cor(CPMdata, method = CorMethod)
  diag(SampleCor) <- NA
  MedianCor <- apply(SampleCor, 1, function(x) median(x, na.rm = T))
  annoRow = data.frame(Group = meta[[groupCol]],
                       Age = meta$Age,
                       Sex = meta$Sex,
                       row.names = meta[[MetaSamleIDCol]])
  annoCol = data.frame(CETS = meta$CETS,
                       Neuronal_MSP = meta$NeuNall_MSP,
                       Microglia_MSP = meta$Microglia_MSP,
                       Oligo_MSP = meta$Oligo_MSP,
                       row.names = meta[[MetaSamleIDCol]])
  annoColors = list(Group = c(Control = "dodgerblue4" , AD = "chocolate1"),
                    Sex = c(F = "indianred4", M = "cornflowerblue"),
                    Age = c("darkseagreen1", "darkorchid4"),
                    CETS = c("chartreuse4","gray97","maroon"),
                    Neuronal_MSP = c("chartreuse4","gray97","maroon"),
                    Microglia_MSP = c("chartreuse4","gray97","maroon"),
                    Oligo_MSP = c("chartreuse4","gray97","maroon"))
  Plot <- pheatmap(SampleCor, angle_col = 90, na_col = "white",border_color = NA,
                   color = colorRampPalette(c("darkblue", "gold2"))(999),
                   labels_row = meta[[MetaSamleCol]],
                   labels_col = meta[[MetaSamleCol]],show_rownames = F, show_colnames = F,
                   annotation_col = annoCol,
                   annotation_row = annoRow,
                   annotation_colors = annoColors,
                   main = title, filename = paste0(ResultsPath, "SampleCorAllPeaks", Cohort, ".pdf"), useDingbats = F, width = 10, height = 8)
  closeDev()
  return(list(countMatrix = subData,
              Metadata = meta,
              SampleCor = SampleCor[order(MedianCor), order(MedianCor)],
              CPMdata = CPMdata,
              HeatMap = Plot))
}

GetCellularProportions <- function(Metadata, normCol = NULL){
  if("CellTypeH3K27ac.tsv" %in% list.files("CellTypeH3K27ac")){
    CellTypePeaks <- read.table("CellTypeH3K27ac/CellTypeH3K27ac.tsv", header = T, sep = "\t")
  } else {
    source(paste0("ProjectScripts/CellProportions.R"))
  }
  
  #Get relative cell proportion for the samples based on differential NeuN positive and negative cell H3K27ac peaks  
  HTseqCountsSamples <- read.table(CellTypePeakCountLoc, header = T, sep = "\t")
  names(HTseqCountsSamples) <- sapply(names(HTseqCountsSamples), function(x){
    x = gsub(".*bamfiles.0?", "", x)
    strsplit(x,"_")[[1]][1]
  })
  HTseqCountsSamples$NeuronFC <- CellTypePeaks$log2FoldChange[match(HTseqCountsSamples$Geneid, CellTypePeaks$PeakName)]
  
  counMatrixSamples <- HTseqCountsSamples %>% select(as.character(Metadata$GSM)) %>% as.matrix()
  rownames(counMatrixSamples) <- as.character(HTseqCountsSamples$Geneid)
  
  if(!is.null(normCol)){
    samplesCPM <- sapply(1:ncol(counMatrixSamples), function(i){
      counMatrixSamples[,i]/Metadata[[normCol]][i]
    })
    colnames(samplesCPM) <- colnames(counMatrixSamples)
  } else {
    samplesCPM <- Count2CPM(counMatrixSamples)
  }
  counMatrixSamples2 <- log2(samplesCPM +1)
  
  ##### Get expression-based cell type marker genes ###########
  
  if(!"homologene" %in% rownames(installed.packages())){
    install_github("oganm/homologene", force = T)
  }
  library(homologene)
  
  if(!"markerGeneProfile" %in% rownames(installed.packages())){
    install_github("oganm/markerGeneProfile", force = T)
  }
  library(markerGeneProfile)
  data("mouseMarkerGenesCombined")
  
  CellType_genes <- mouseMarkerGenesCombined$Cortex
  for(i in 1:length(CellType_genes)){
    
    CellType_genes[[i]] <- as.vector(mouse2human(CellType_genes[[i]])$humanGene)
    
  }
  
  # Get the cell type specific (neuron vs glia) peaks
  
  if("CellTypeH3K27ac.tsv" %in% list.files("CellTypeH3K27ac")){
    CellTypePeaks <- read.table("CellTypeH3K27ac/CellTypeH3K27ac.tsv", header = T, sep = "\t")
  } else {
    source(paste0("ProjectScripts/CellProportions.R"))
  }
  
  CellTypePeaks %<>% mutate(Peak_Gene = paste(PeakName, symbol, sep = "_"))
  CellTypePeaks %<>% filter(!duplicated(Peak_Gene), !is.na(symbol))
  
  #Identify cell type specific peaks by intersecting between the cell type specific genes and cell type specific peaks
  
  CellTypeSpecificPeaks <-lapply(CellType_genes, function(CellType){
    MarkerPeaks <- CellTypePeaks %>% filter(symbol %in% CellType, !duplicated(Peak_Gene)) %>% select(PeakName, symbol, log2FoldChange, padj) %>% arrange(padj) 
    if(nrow(MarkerPeaks) > 2){
      MarkerPeaks
    } else {
      NULL
    }
  })
  
  #Add general NeuN enriched peaks
  CellTypeSpecificPeaks$NeuNall <- CellTypePeaks %>% filter(log2FoldChange > 3, !duplicated(Peak_Gene)) %>% select(PeakName, symbol, log2FoldChange, padj) %>% arrange(padj)
  
  CellTypeSpecificPeaks <- CellTypeSpecificPeaks[!unlist(lapply(CellTypeSpecificPeaks, function(x) is.null(x)))]
  
  CellTypeSpecificPeakNames <- lapply(CellTypeSpecificPeaks, function(CellType){
    CellType %>% filter(!duplicated(PeakName)) %>% .$PeakName %>% as.character()
  })
  
  ## Run PCA on celltype specific peaks (+ general Neuronal peaks) ##########
  
  PCAcellType <- lapply(CellTypeSpecificPeakNames, function(CellType){
    counMatrixSamples2 <- counMatrixSamples2[rownames(counMatrixSamples2) %in% CellType,]
    temp <- prcomp(t(counMatrixSamples2), scale. = T)
    sumPC1 = sum(temp$rotation[,1])
    if(sumPC1 < 0){
      temp$rotation[,1] <- -1*temp$rotation[,1]
      temp$x[,1] <- -1*temp$x[,1]
    }
    while(sum(temp$rotation[,1] > 0) < nrow(temp$rotation)){
      counMatrixSamples2 <- counMatrixSamples2[rownames(counMatrixSamples2) %in% rownames(temp$rotation)[temp$rotation[,1] > 0],]
      temp <- prcomp(t(counMatrixSamples2), scale. = T)
      sumPC1 = sum(temp$rotation[,1])
      if(sumPC1 < 0){
        temp$rotation[,1] <- -1*temp$rotation[,1]
        temp$x[,1] <- -1*temp$x[,1]
      }
    }
    temp
  })
  
  names(PCAcellType) <- paste0(names(PCAcellType), "_MSP")
  VarianceExplained <- lapply(PCAcellType, function(cell){
    temp <- cell %>% summary
    if(length(temp) > 3){
      temp$importance %>% .[2,1:3] 
    }
  }) %>% do.call(rbind,.) %>% data.frame %>%
    mutate(Peaks = unlist(lapply(PCAcellType, function(x) nrow(x$rotation))))
  rownames(VarianceExplained) <- names(PCAcellType)
  
  MSP_out <- lapply(PCAcellType, function(CellType){
    rescale(CellType$x[,1], to = c(0,1))
  }) %>% do.call(cbind, .) %>% data.frame %>% mutate(Sample = sapply(row.names(.), function(x) gsub("X", "", x)))
  
  #Add general calculation based on ratio between bneuronal and glial peaks, this is for comparing different regions
  counMatrixSamples2 <- counMatrixSamples2[rownames(counMatrixSamples2) %in% (CellTypePeaks %>% filter(abs(log2FoldChange) > 3, !duplicated(PeakName)) %>% .$PeakName),]
  
  CellTypePCA <- prcomp(t(counMatrixSamples2), scale. = T)
  #Make sure that the loadings of genes upregulated in neurons is positive
  DF <- data.frame(Peak = rownames(CellTypePCA$rotation),
                   PC1 = CellTypePCA$rotation[,1],
                   NeuronFC = HTseqCountsSamples$NeuronFC[match(rownames(CellTypePCA$rotation), as.character(HTseqCountsSamples$Geneid))])
  DF$NeuronalPeak <- sapply(DF$NeuronFC, function(x){
    if(x > 0){
      "NeuronalPeak"
    } else {
      "GlialPeak"
    }
  })
  
  ggplot(DF, aes(NeuronalPeak,PC1, color = NeuronalPeak)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.5)
  MedPC1 <- DF %>% group_by(NeuronalPeak) %>% summarise(Median = median(PC1)) %>% data.frame()
  
  if((MedPC1 %>% filter(MedPC1$NeuronalPeak == "NeuronalPeak") %>% .$Median) < (MedPC1 %>% filter(MedPC1$NeuronalPeak == "GlialPeak") %>% .$Median )){
    CellTypePCA$rotation[,1] <- -1*CellTypePCA$rotation[,1]
    CellTypePCA$x[,1] <- -1*CellTypePCA$x[,1]
  }
  
  Metadata$NeuronProp <- rescale(to = c(0,1), x = CellTypePCA$x[,1])
  
  Metadata <- merge(Metadata, MSP_out, by.x = "GSM", by.y = "Sample")
  return(Metadata)
}

RunDESeq <- function(data, meta, normFactor=NULL, sampleToFilter = "none",  FullModel, ReducedModel = "~1", test = "Wald", UseModelMatrix = FALSE, FitType = "parametric", MetaSamleCol = "activemotif_id", SampleNameCol = "SampleName"){
  meta = meta[!grepl(sampleToFilter, meta[[MetaSamleCol]]),] %>% droplevels
  data = data[,as.character(meta[[SampleNameCol]])]
  ModelMatrix <- model.matrix(FullModel, meta)
  
  DESeqDS <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ModelMatrix)
  
  DESeqDS <-  estimateSizeFactors(DESeqDS)
  if(!is.null(normFactor)){
    DESeqDS$sizeFactor <- meta[[normFactor]]
  }
  DESeqDS <- estimateDispersions(DESeqDS, fitType = FitType)
  
  
  if(test == "Wald"){
    if(UseModelMatrix){
      DESeqOut <-  nbinomWaldTest(DESeqDS, modelMatrix = ModelMatrix)
    } else {
      DESeqOut <-  nbinomWaldTest(DESeqDS)
    }
    
  } else if(test == "LRT"){
    ModelMatrixReduced = model.matrix(ReducedModel, meta)
    DESeqOut <-  nbinomLRT(DESeqDS, full = ModelMatrix, reduced = ModelMatrixReduced)
  }
  
  return(DESeqOut)
}

GetDESeqResults <- function(DESeqOut, coef, cutoff = 0.05, indepFilter = TRUE, CountAnnoFile){
  DESeqResults <- results(DESeqOut, name = coef, alpha = cutoff, format = "DataFrame", independentFiltering = indepFilter)
  DESeqResultsDF <- DESeqResults %>% data.frame
  DESeqResultsDF$PeakName = rownames(DESeqResultsDF)
  DESeqResults %>% summary
  
  return(DESeqResultsDF)
}  

AnnotDESeqResult <- function(DESeqResultsDF, CountAnnoFile, by.x, by.y){
  DESeqResultsDF <- merge(DESeqResultsDF, CountAnnoFile,
                          by.x = by.x, by.y = by.y)
  
  names(DESeqResultsDF) <- sapply(names(DESeqResultsDF), function(x) gsub("annoFileCollapsed.", "", x))
  DESeqResultsDF %<>% mutate(Peak_Gene = paste(PeakName, symbol, sep = "_"))
  return(DESeqResultsDF)
}

plotQQ <- function(dist1 = NULL, dist2, xlab = "-log10(uniform)", ylab = "-log10(pvalues)", title){
  dist1 = runif(length(dist2), 0, 1)
  qqplot(-log10(dist1), -log10(dist2), main = title, xlab = xlab, ylab = ylab)
  abline(a=0, b=1)
}

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

GetAdjCountDESeq <- function(dds, Gene,  adjCov){
  GeneRow <- which(rownames(dds) == Gene)
  Mu <- log2(t(t(assays(dds)$mu[GeneRow,])/sizeFactors(dds)))
  Counts <-  log2(t(t(assays(dds)$counts[GeneRow,])/sizeFactors(dds)))
  Resid <- Counts - Mu
  corrFactor <- sizeFactors(dds)
  Coef <- coef(dds)[GeneRow,]
  mod <- attr(dds, "modelMatrix")
  modAdj <-mod
  for(cov in as.character(adjCov$Cov)){
    adjType = adjCov %>% filter(Cov == cov) %>% .$adjType %>% as.character()
    if(adjType == "mean"){
      modAdj[,cov] <- mean(modAdj[,cov], na.rm=T)
      
    } else if (adjType == "base"){
      modAdj[,cov] <- 0
    }
  }
  AdjValue <- (modAdj %*% Coef) + Resid 
  return(AdjValue)
}

ggMA <- function(resultObject, contours=FALSE, minP=0, lims=NULL, trend=NULL, se=FALSE, geneList=NULL, geneColName="GeneSymbol", colour_pvalue=FALSE,
                 logPcolLow = "black", logPcolHigh = "darkorange", contColor = "gray48") {    
  if (!is.null(geneList)){    
    if(!"list" %in% class(geneList)){    
      stop("geneList argument has to be a named list of character arrays, each element corresponding to a pathway")    
    }    
  }                                                                                                                        
  lims <- sort(lims)
  dT <- resultObject %>%    
    mutate(baseMean=log10(baseMean), logPadj=-log10(padj), logFC=log2FoldChange) %>%    
    arrange(logPadj)    
  dT %<>% filter(!is.na(padj))   
  #if (minP > 0) dT <- dT[dT$logP>-log10(minP),]    
  dT$shape <- factor("c", levels=c("c", "tu", "td"))    
  if (!is.null(lims)) {    
    which.low <- which(dT$logFC<lims[1])    
    which.hig <- which(dT$logFC>lims[2])    
    dT[which.low, "logFC"] <- lims[1]    
    dT[which.low,"shape"] <- "td"    
    dT[which.hig,"logFC"] <- lims[2]    
    dT[which.hig,"shape"] <- "tu"    
  }    
  
  gg <- ggplot(dT, aes(x=baseMean, y=logFC)) +
    theme_classic() +
    geom_point(aes(colour=logPadj, shape=shape)) +    
    scale_shape_manual(values=c(c=16, tu=2, td=6), guide=FALSE) +    
    geom_hline(yintercept=0, colour="red") +    
    labs(x="log2(mean of normalised counts)", y="log2(fold change)")    
  if (colour_pvalue) {    
    gg <- gg + scale_colour_gradient(low = logPcolLow, high = logPcolHigh, name = "-log(padj)")    
  } else {    
    gg <- gg + scale_colour_gradientn(colours=c("grey60", "grey60"), guide=FALSE)    
  }    
  if (contours) {    
    gg <- gg + geom_density2d(colour=contColor)    
  }    
  if(!is.null(trend)) {    
    gg <- gg + geom_smooth(method=trend, se=se, linetype="dashed", fill="red", size=1)    
  }    
  if (!is.null(geneList)) {    
    if (is.null(names(geneList))) stop("geneList argument must be a NAMED list")    
    Genes <- do.call(rbind, lapply(names(geneList), function(pathwayName){    
      dT[dT[[geneColName]] %in% geneList[[pathwayName]],] %>%                                         
        mutate(Pathway=pathwayName)
    }))
    gg <- gg + geom_point(data=Genes, aes(x=baseMean, y=logFC, fill=Pathway), colour="black", shape=21, size=2)
  }
  return(gg)
}


manhattan2 <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", "gray60"),
                        pchAll = 20, pointSize = 0.3, pchHighlight = 20, colThresh = "red", SpecificPeaks = NULL,
                        colGenomeWide = "blue", colHighlight = "brown1", Annotate = NULL,
                        chrlabs = NULL, suggestiveline = -log10(1e-05), 
                        genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
                        annotatePval = NULL, annotateTop = TRUE, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = pchAll, cex = pointSize, xlim = c(xmin, xmax), ylim = c(0, 
                                                                                          ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = pchAll, cex = pointSize, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = pchAll, cex = pointSize, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = colThresh)
  if (genomewideline) 
    abline(h = genomewideline, col = colGenomeWide)
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = colHighlight, pch = pchHighlight, 
                             ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625, labs = topHits$SNP, cex = 0.45), 
           ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  
  if(!is.null(SpecificPeaks)){
    SpecificPeaks <- d %>% filter(SNP %in% SpecificPeaks)
    printSNPS <- NULL
    for (i in unique(topHits$CHR)) {
      chrSNPs <- SpecificPeaks[SpecificPeaks$CHR == i, ]
      printSNPS <- rbind(printSNPS, chrSNPs[1, ])
    }
    textxy(printSNPS$pos, -log10(printSNPS$P), offset = 0.625, 
           labs = printSNPS$SNP, cex = 0.5, ...)
  }
  par(xpd = FALSE)
}

GetDESeq2ResultsRNA <- function(DESeqOut, coef, alpha = 0.05, indepFilter = TRUE){
  DEresults <- results(DESeqOut, name = coef, alpha = alpha, format = "DataFrame", independentFiltering = indepFilter)
  DEresults$GeneSymbol <- geneNames$hgnc_symbol[match(rownames(DEresults), geneNames$ensembl_gene_id)]
  DEresults$EnsemblID <- rownames(DEresults)
  DEresults %<>% data.frame %>% filter(GeneSymbol != "")
  return(DEresults)
}
