library("rtracklayer")
library("BSgenome.Celegans.UCSC.ce11")
library("GenomicAlignments")
library(colorRamps)
library(naturalsort)

# read in list of bam filenames
fls <- list.files(path="../scripts/alignedBWA/subBam/", pattern="*bam$",full.names=T,recursive = F)
fN=5


#make GRanges object from BSgenome (accessed with "Celegans") (i remove the mito chr)
genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6],ranges=IRanges(start=1,
                                                                  end=seqlengths(Celegans)[1:6]), strand="*")


#########################
## generate some bed-graph tracks to look at in IGV.

# read in bam files one at a time, calculate coverage, and create a bedGraph file
#dir.create("./alignedBWA/bedGraph")
#dir.create("./alignedBWA/RDS")
#dir.create("./alignedBWA/bedGraph/normalized")
#dir.create("./alignedBWA/bedGraph/enrichmentchipbwa/")

# normalize counts by dividing by (median coverage + 8)
for (f in fls) {
 print(f)
  #get base of file name from current filename
  fileName<-strsplit(strsplit(f,"/")[[1]][5],"_")[[1]][1]
  #read in bam file
  bamFile<-readGAlignments(f)
  #calculate coverage
  sampleCoverage<-coverage(bamFile)[1:6]
  print(sampleCoverage)
  medCoverage<-median(as.vector(unlist(sampleCoverage)))
  print(medCoverage)
  #add eight to the medianCoverage to deal with samples that have 0 median coverage
  medCoverage <- median(sampleCoverage+8)
  normSampleCoverage<-(sampleCoverage+8)/(medCoverage)
  print(normSampleCoverage)
  #save as RDS to avoid imposing GRanges of random size
  saveRDS(normSampleCoverage,paste0("../scripts/alignedBWA/RDS/norm+8medsub_",fileName,".RDS"))

  #save as a bedgraph file in the bedgraph directory
  export(normSampleCoverage,paste0("../scripts/alignedBWA/bedGraph/normalized/norm+8medsub_",fileName,".bedGraph"))
}

#read files in pairs (input and IP) to calculate enrichment
#for (f in seq(1,length(fls),by=2)) {
  #get base of file name from current filename
  fileName1<-strsplit(strsplit(fls[f],"/")[[1]][5],"_")[[1]][1]
  fileName2<-strsplit(strsplit(fls[f+1],"/")[[1]][5],"_")[[1]][1]
  input<-readRDS(paste0("../scripts/alignedBWA/RDS/norm+8medsub_", fileName2, ".RDS") )
  chip<-readRDS(paste0("../scripts/alignedBWA/RDS/norm+8medsub_", fileName1, ".RDS") )
  #calculate enrichment
  enrichment<- log2(chip/input)
  saveRDS(enrichment,paste0("../scripts//alignedBWA/RDS/enrichmentsub", fileName1, ".RDS"))
  export(enrichment,paste0("../scripts/alignedBWA/bedGraph/enrichmentchipbwa/enrichment+8medsub", fileName1, ".bedGraph"))
}

# get a list of all the files with the enrichment data
normCountsFiles<-list.files(path="./alignedBWA/bedGraph/enrichmentchipbwa/", "bedGraph")
#create GRangesList object with different sized tiles along genome
tile1Mb<-unlist(tile(x=genomeGR,width=1000000))
tile100kb<-unlist(tile(x=genomeGR,width=100000))
tile10kb<-unlist(tile(x=genomeGR,width=10000))
tile1kb<-unlist(tile(x=genomeGR,width=1000))
tile100bp<-unlist(tile(x=genomeGR,width=100))
tileList<-list("tile1Mb"=tile1Mb,"tile100kb"=tile100kb,"tile10kb"=tile10kb,"tile1kb"=tile1kb,"tile100bp"=tile100bp)


# read in the bedgraph files in pairs to calculate enrichment counts at different scales
#get normalised counts
for (f in 1:length(normCountsFiles)) {
normCounts<-import(paste0("../scripts/alignedBWA/bedGraph/enrichmentchipbwa/",normCountsFiles[f]))
sampleName=strsplit(strsplit(normCountsFiles[f],split=".bedGraph",fixed=T)[[1]],split="enrichment+8medsub")[[1]][2]


for (i in 1:length(tileList)) {
  #make a GRangesList from GRanges
  tiles<-split(tileList[[i]],seqnames(tileList[[i]]))

  # convert GRanges back to coverage
  rle<-coverage(normCounts,weight="score")
  names(rle)<-paste0("chr",names(rle))

  #create views corresponding to the tiles on the coverage rle and average counts
  summedCov<-viewMeans(RleViewsList(rangesList=tiles,rleList = rle))

  #now save the values into tiles
  mcols(tileList[[i]])[,sampleName]<-unlist(summedCov)
}
}

## plot box plots by chromosome for all data sets at each tile size
pdf("supercoilingmastersbwa+8medsub.pdf",paper="a4r",width=11,height=8)
par(mfrow=c(2,2))
for (i in 1:length(tileList)) {
  tiles<-tileList[[i]]
  samples<-names(mcols(tiles))
  #plot the raw counts per Mb by chr for each sample (just to look at raw data)
   for (s in samples) {
    boxplot(mcols(tiles)[,s]~as.vector(seqnames(tiles)),main=paste(s,names(tileList)[i]),outline=F,notch=T, ylab="Enrichment (IP$
    abline(h=median(mcols(tiles)[seqnames(tiles)=="chr",s]),lty=2,col="red")
  }
}
dev.off()




## plot scatter plots by chromosome for all data sets at each tile size
pdf("pairwise_smScatterbwa+8medsub.pdf",paper="a4r",width=11,height=8)
par(mfrow=c(4,4))
for (i in 1:length(tileList)) {
  tiles<-tileList[[i]]
  samples<-names(mcols(tiles))
  #plot the raw counts per Mb by chr for each sample (just to look at raw data)
  for (s1 in samples) {
    for (s2 in samples) {
      #if (s1!=s2) {
      myCor<-round(cor(mcols(tiles)[,s1],mcols(tiles)[,s2],use="pairwise.complete.obs"),2)
      smoothScatter(mcols(tiles)[,s1],mcols(tiles)[,s2],colramp=matlab.like,
                    main=paste0(s1," vs ",s2," ",names(tileList)[i]," (R=",myCor,")"),
                    xlab=s2,ylab=s1,cex.main=0.9)
      }
    }
  }
}
dev.off()


########### sliding windows to create smoothed bedgraph file ##########

## function to create list of GRanges with sliding windows of different sizes
makeWinList<-function(GRobj,widths,step=0.1){
 winList<-list()
  for (winSize in widths) {
    # get formatted text for window size
    units<-formatWinSize(winSize)
    # create GRanges of windows of that size and add to list of windows
    winList[paste0("win",units)]<-unlist(slidingWindows(x=GRobj,width=winSize,step=winSize*step))
  }
  return(winList)
}

# function to take a numerical window size and convert to nice text for labels etc
formatWinSize<-function(winSize) {
  if (winSize/1000 < 1) {
    units<-paste0(winSize,"bp")
  } else if (winSize/1e6 >= 1) {
    units<-paste0(winSize/1e6,"Mb")
  } else {
    units<-paste0(winSize/1000,"kb")
  }
  return(units)
}

#sizes of windows you want to create, (e means raised to the power of)
winWidths<-c(10,100,1000,1e4,1e5,5e5)
#create create list of GRanges with sliding windows of different sizes
winList<-makeWinList(genomeGR,winWidths)

# get a list of all the files with the enrichment data
normCountsFiles<-list.files(path="./alignedBWA/bedGraph/enrichmentchipbwa/", "bedGraph")


# read in the bedgraph files in pairs to calculate enrichment counts at different scales
#get normalised counts
  for (f in 1:length(normCountsFiles)) {
    normCounts<-import(paste0("../scripts/alignedBWA/bedGraph/enrichmentchipbwa/",normCountsFiles[f]))
    sampleName=strsplit(strsplit(normCountsFiles[f],split=".bedGraph",fixed=T)[[1]],split="enrichment+8medsub")[[1]][2]

for (i in 1:length(winList)) {
  #make a GRangesList from GRanges
  windows<-split(winList[[i]],seqnames(winList[[i]]))

  # convert GRanges back to coverage
  rle<-coverage(normCounts,weight="score")

  #create views corresponding to the windows on the coverage rle and average counts
  summedCov<-viewMeans(RleViewsList(rangesList=windows,rleList = rle))

  #now save the values into windows
  windows<-unlist(windows)
  mcols(windows)$score<-unlist(summedCov)
  windows<-resize(windows,width=1,fix="center")
  export(windows,paste0("./alignedBWA/bedGraph/smEnrich_",sampleName,"_",names(winList)[i],".bedGraph"))
}
}


################## create GC frequency track for loading to IGV using various window sizes #######

### create GC frequency track for display using different sized views
winSize<-winWidths

# get sequence of CE genome to calculate GC content
CEgenome<-getSeq(Celegans)
#make Granges object with range for every base in the genome
GR<-unlist(slidingWindows(x=genomeGR,width=1,step=1))

# calculate GC frequency over different window sizes and create bedgraph
for (w in winSize){
  gc<-lapply(CEgenome[1:6],function(x) {letterFrequencyInSlidingView(x,view.width=w,letters=c("CG"),as.prob=T)})
  # pad front and back of GC with half w repeats of first and last value so that it has the same
# length as the genome
  padded<-lapply(gc,function(x) {
    newX<-c(rep(x[1],times=w/2),x,rep(x[length(x)],times=(w/2-1)))
    return(newX)
  })
  mcols(GR)$score<-unlist(padded)
  i<-seq(1,length(GR),by=w*0.1)
  export(GR[i],paste0("./alignedBWA/bedGraph/GCfreq_",formatWinSize(w),".bedGraph"))
}

