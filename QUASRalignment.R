setwd("/")
#to set working directory
library(QuasR)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Celegans.UCSC.ce11)
library(Gviz)
sampleFile<-read.table("samples.txt")
proj = 
  qAlign(sampleFile="samples.txt",
         genome= "BSgenome.Celegans.UCSC.ce11",
         aligner="Rbowtie", 
         projectName="data",
         cacheDir = "/Users/imac/Desktop/Supercoilinganalysis/tocomparefastafiles/cache/")
qQCReport(proj, pdfFilename="tocompare2.pdf")

qExportWig(proj, scaling=F)
qExportWig(proj, binsize=1L, scaling=TRUE, collapseBySample=TRUE)
#to align using QuasR
#to create pdf of the q

#to export genome wig file from alignments
library(GenomicRanges)
#to load genomic ranges library
