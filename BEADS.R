
#install bioconductor package
source("http://bioconductor.org/biocLite.R")
biocLite(c('methods','IRanges','BSgenome','digest','Rsamtools','rtracklayer','GenomicRanges','Biostrings'))

#install rbeads directly from github
if (!require("devtools")) install.packages("devtools")
devtools::install_github("przemol/rbeads")

source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Rsamtools")

#run required libraries
library(rbeads)
library(Rsamtools)
library(rtracklayer)
library(AnnotationHub)

#get path to file list
fls <- list.files(path="./alignedBWA/picdup", pattern="*bam$",full.names=T,recursive = F)
fN=4

dir.create("./alignedBWA/normalizedwbeads/")


#run command to normalize the data using the beads algorithm
for (f in seq(1,length(fls),by=2)) {
  print(f)
  #get base of file name from current filename
  fileName1<-strsplit(strsplit(fls[f],"/")[[1]][4],"_")[[1]][1]
  fileName2<-strsplit(strsplit(fls[f+1],"/")[[1]][4],"_")[[1]][1]
  input<-bamFile(paste0("./alignedBWA/picdup/", fileName2) )
  chip<-bamFile(paste0("./alignedBWA/picdup/", fileName1) )
  mappabilitytrack<- import("./ws270.bw")
  ref_fa<- read.fasta("./home/ubelix/izb/bi18k694/genomeversion/ws270/c_elegans.PRJNA13758.WS270.genomic.fa")
  BEADS<-beads(chip, input, mappabilitytrack, ref_fa, uniq = TRUE,mapq_cutoff = 20L, rdata = FALSE, export_er = TRUE, ...)
  export(BEADS,paste0("./alignedBWA/normalizedwbeads/beads_",fileName1,".bw")
}
