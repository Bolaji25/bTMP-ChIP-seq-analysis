library(QuasR)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Celegans.UCSC.ce11)
library(Gviz)


sampleFile<-read.table("./mysamplesmasters.txt")
dir.create("./alignedquasr")
clObj = makeCluster(4)
proj = 
  qAlign(sampleFile= sampleFile,
         genome= "BSgenome.Celegans.UCSC.ce11",
         aligner="Rbowtie", 
         projectName="data",
         cacheDir = "./alignedquasr", clObj = clObj)
qQCReport(proj, pdfFilename="bTMPsupercoilingmastersredo.pdf")

#qExportWig(proj, scaling=F)
qExportWig(proj, binsize=1L, scaling=TRUE, collapseBySample=TRUE)
