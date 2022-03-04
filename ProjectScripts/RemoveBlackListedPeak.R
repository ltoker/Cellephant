sapply(list.files("data/Peaks/", pattern = ".narrowPeak.gz"), function(orgFile){
  inFile = paste0("data/Peaks/", orgFile)
  outFile = gsub(".gz", "Clean.gz", inFile)
  Command <- paste0("bedtools intersect -v -a ", inFile,  " -b data/H3K27Ac_black_list.bed | gzip -nc > ", outFile)
  system(Command)
})
  