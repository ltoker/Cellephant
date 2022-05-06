if(!AnnoLoc %in% list.dirs(full.names = F, recursive = T)){
  dir.create(AnnoLoc)
}

AnnoLoc = paste0(AnnoLoc, "/")


URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz"
AssemblyFilename = "gencode.v35.annotation.gff3.gz"

if(!file.exists(paste0(AnnoLoc, AssemblyFilename))){
  download.file(URL, destfile = paste0(AnnoLoc, AssemblyFilename))
}

txdbFilename  = paste0(paste0(strsplit(AssemblyFilename, "\\.")[[1]][c(2,4)], collapse = "."), "_DB")
GTF_file = paste0(AnnoLoc, AssemblyFilename)
geneNameFile = paste0(AnnoLoc, strsplit(AssemblyFilename, "\\.")[[1]][2], "_geneNames.Rds")

#Get transcriptome annotations
if(!file.exists(paste0(AnnoLoc, txdbFilename))){
  txdb <- makeTxDbFromGFF(paste0(AnnoLoc, AssemblyFilename))
  saveDb(txdb, paste0(AnnoLoc, txdbFilename))
} else {
  txdb <- loadDb(paste0(AnnoLoc, txdbFilename))
}

if(!file.exists(geneNameFile)){
  geneNames <- rtracklayer::import.gff3(GTF_file) %>%
    data.frame() %>% filter(type == "transcript")  %>%
    select(seqnames, start, end, ID,
           gene_id, gene_name, gene_type, transcript_type)
  saveRDS(geneNames, geneNameFile)
} else {
  geneNames <- readRDS(geneNameFile)
}

geneNames$gene_id2 <- sapply(geneNames$gene_id, function(x){
  strsplit(x, "\\.")[[1]][1]
})
