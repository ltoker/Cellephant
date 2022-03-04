Path = "Data/GSE102538_RAW/"

MarziCounts <- sapply(list.files(path = Path, full.names = T), function(x){
  read.table(x, header = F, sep = "\t")
}, simplify = F)

MarziCountsDF <- lapply(MarziCounts, function(x){
  x$V2
}) %>% do.call(cbind, .)

MarziCountsDF <- cbind(MarziCounts[[1]]["V1"], MarziCountsDF)

names(MarziCountsDF)[1] <- "Geneid"

MarziPeaks <- read.table("Data/GSE102538_H3K27ac_EntorhinalCortex.bed.gz", header = F, sep = " ")
names(MarziPeaks) <- c("CHR",   "START",    "END", "Geneid")

MarziPeaks %<>% mutate(Strand = "+",
                        Length = END - START)

MarziCountsDF <- merge(MarziPeaks, MarziCountsDF, by = "Geneid")

names(MarziCountsDF) <- sapply(names(MarziCountsDF), function(x){
  x = gsub(paste0(Path, "/"), "", x)
  strsplit(x, "_")[[1]][1]
})

write.table(MarziCountsDF, "Data/MarziPaperCounts.tsv", col.names = T, row.names = F, sep = "\t")
