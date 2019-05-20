library(rtracklayer)
ce<-import("ws260gene.gtf")


ce1<-ce[!(is.na(mcols(ce)$biotype))]
cePC<-ce1[mcols(ce1)$biotype=="protein_coding",]

ceGene<-cePC[mcols(cePC)$type=="gene"]

i<-which(seqnames(ceGene)=="X")

xChr<-ceGene[i]

j<-which(seqnames(ceGene)!="X")
#aChr<-ceGene[j]

export(xChr,"ws260X.gtf","gtf")
export(aChr,"ws260A.gtf","gtf")
#############################################


ce2<-ce[mcols(ce)$type=="five_prime_UTR"]

i<-which(seqnames(ce2)=="X")

xChr<-ce2[i]

j<-which(seqnames(ce2)!="X")
aChr<-ce2[j]

export(xChr,"utr5p_chrX.gtf","gtf")
export(aChr,"utr5p_chrA.gtf","gtf")
#############################################

ce3<-readRDS("low_exp_tss.RDS")
ce3


#i<-which(seqnames(ce3)=="X")

#xChr<-ce3[i]

#j<-which(seqnames(ce3)!="X")
#aChr<-ce3[j]

export(ce3,"lowexpTSS.gtf","gtf")
#export(aChr,"alltss_chrA.gtf","gtf")
#############################################

ce3<-readRDS("/Users/imac/Downloads/Untitled_Message/ChenKreusSaitoTSS_highConf_872.RDS")
ce3


i<-which(seqnames(ce3)=="X")

xChr<-ce3[i]

j<-which(seqnames(ce3)!="X")
aChr<-ce3[j]

export(xChr,"highconftss_chrX.gtf","gtf")
export(aChr,"highconftss_chrA.gtf","gtf")
#############################################

ce4<-import("ws260gene.gtf")


i<-which(seqnames(ce4)=="chrX")

xChr<-ce4[i]

j<-which(seqnames(ce4)!="chrX")
aChr<-ce4[j]

export(xChr,"ws260_chrX.gtf","gtf")
export(aChr,"ws260_chrA.gtf","gtf")


##########################################

epd<-import("/Users/imac/Downloads/hglft_genome_41e5_ad7e40.bed")

ce<-import("./Caenorhabditis_elegans.WBcel235.40.gff3.gz")
ce1<-ce[!(is.na(mcols(ce)$biotype))]
cePC<-ce1[mcols(ce1)$biotype=="protein_coding",]
wbgene<-cePC[mcols(cePC)$type=="gene"]

# want to find overlaps in genes between these two objects
# remove _1 from name of genes in epd
epd$name<-gsub("_1$","",epd$name)

# first subset both vectors to only those genes found in both
inBoth<-epd$name %in% wbgene$Name
epd1<-epd[inBoth]
inBoth<-wbgene$Name %in% epd1$name
wbgene1<-wbgene[inBoth]

# reorder both vectors by gene name
epd2<-epd1[order(epd1$name)]
wbgene2<-wbgene1[order(wbgene1$Name)]
head(epd2$name)
head(wbgene2$Name)

# construct a meta promoter (tss to gene start) depending on strand gene is on
metaPromoter<-epd2
idxPos<-as.vector(strand(epd2)=="+")
# for genes on the positive strand:
start(metaPromoter[idxPos])<-start(epd2[idxPos])
end(metaPromoter[idxPos])<-start(wbgene2[idxPos])
# for genes on the negative strand:
end(metaPromoter)[!idxPos]<-end(epd2[!idxPos])
start(metaPromoter)[!idxPos]<-end(wbgene2[!idxPos])

hist(start(wbgene2[idxPos])-start(epd2[idxPos]),breaks=50)
table((start(wbgene2[idxPos])-start(epd2[idxPos]))>0)

hist(end(epd2[!idxPos])-end(wbgene2[!idxPos]),breaks=50,xlim=c(-2000,2000))
table((end(epd2[!idxPos])-end(wbgene2[!idxPos]))>0)

#############################################

ce3<-readRDS("./GR_ws260gene_L3expnDecile.RDS")
ce3
#to remove NAs(done below)
idx<-is.na(ce3$expnDecile) #all NAs
sum(idx)
ce3<-ce3[!idx] #ce3 without the NAs anymore


i<-which(mcols(ce3)$expnDecile==10)

hegenes<-ce3[i]

j<-which(mcols(ce3)$expnDecile==1)
legenes<-ce3[j]
mcols(hegenes)<-NULL
mcols(legenes)<-NULL
export(legenes,"lowexpggeness.gtf","gtf")
#############################################
library(rtracklayer)
ce<-import("ws260gene.gtf")


ce1<-ce[!(is.na(mcols(ce)$biotype))]
cePC<-ce1[mcols(ce1)$biotype=="protein_coding",]

ceGene<-cePC[mcols(cePC)$type=="gene"]

i<-which(seqnames(ceGene)=="X")

xChr<-ceGene[i]

j<-which(seqnames(ceGene)=="I")
IChr<-ceGene[j]

k<-which(seqnames(ceGene)=="II")
IIChr<-ceGene[k]

l<-which(seqnames(ceGene)=="III")
IIIChr<-ceGene[j]

m<-which(seqnames(ceGene)=="IV")
IVChr<-ceGene[m]

n<-which(seqnames(ceGene)=="V")
VChr<-ceGene[n]
export(xChr,"chrXws260.gtf","gtf")
export(IChr,"chrIws260.gtf","gtf")
export(IIChr,"chrIIws260.gtf","gtf")
export(IIIChr,"chrIIIws260.gtf","gtf")
export(IVChr,"chrIVws260.gtf","gtf")
export(VChr,"chrVws260.gtf","gtf")
#############################################
kreus<-readRDS("KreusTSS_WS235.Rds")
kreus
mcols(kreus)<-NULL
i<-which(seqnames(Kreus)=="chrX")

xChr<-Kreus[i]

j<-which(seqnames(Kreus)!="chrX")
aChr<-Kreus[j]

export(xChr,"kreustsschrX.gtf","gtf")
export(aChr,"kreustsschrA.gtf","gtf")
#############################################
chen<-readRDS("ChenTSS_WS235.Rds")
chen
table(chen$assignmentType)
wormtss<-chen[mcols(chen)$assignmentType=="wormbase_tss"]
i<-which(seqnames(chen)=="chrX")

xChr<-chen[i]

j<-which(seqnames(chen)!="chrX")
aChr<-chen[j]

export(xChr,"chentsschrX.gtf","gtf")
export(aChr,"chentsschrA.gtf","gtf")
#############################################
saito<-readRDS("SaitoTSS_WS235.Rds")
saito
i<-which(seqnames(saito)=="chrX")

xChr<-saito[i]

j<-which(seqnames(saito)!="chrX")
aChr<-saito[j]

export(xChr,"saitotsschrX.gtf","gtf")
export(aChr,"saitotsschrA.gtf","gtf")
############################################
ce3<-readRDS("./GR_ws260gene_L3expnDecile.RDS")
ce3
#to remove NAs(done below)
idx<-is.na(ce3$expnDecile) #all NAs
sum(idx)
ce3<-ce3[!idx] #ce3 without the NAs anymore
ordered<-ce3[order(mcols(ce3)$expnDecile)] #order genes
export(ordered,"orderedgenes.gtf","gtf")
#############################################
nonTSTSSexpordered<-readRDS("nonTS_TSS_expnOrdered.RDS")
nonTSTSSexpordered
i<-which(seqnames(nonTSTSSexpordered)=="X")

xChr<-nonTSTSSexpordered[i]
xChr
j<-which(seqnames(nonTSTSSexpordered)!="X")
aChr<-nonTSTSSexpordered[j]
aChr
export(xChr,"nonTSTSSexporderedchrX.gtf","gtf")
export(aChr,"nonTSTSSxporderedchrA.gtf","gtf")
#export(nonTSTSSexpordered,"nonTSTSSexpordered.gtf","gtf")
#############################################
nonTShiconfTSSexpordered<-readRDS("nonTS_TSShiconf_expnOrdered.RDS")
nonTShiconfTSSexpordered
i<-which(seqnames(nonTShiconfTSSexpordered)=="X")

xChr<-nonTShiconfTSSexpordered[i]
xChr
j<-which(seqnames(nonTShiconfTSSexpordered)!="X")
aChr<-nonTShiconfTSSexpordered[j]
aChr
export(xChr,"nonTShiconfTSSexporderedchrX.gtf","gtf")
export(aChr,"nonTShiconfTSSexporderedchrA.gtf","gtf")
export(nonTShiconfTSSexpordered,"nonTShiconfTSSexpordered.gtf","gtf")
#############################################
ChenTSSshapexpordered<-readRDS("ChenTSS_shape_expnOrdered.RDS")
ChenTSSshapexpordered
i<-which(mcols(ChenTSSshapexpordered)$decileL3expn==10)

hegenes<-ChenTSSshapexpordered[i]
hegenes
j<-which(mcols(ChenTSSshapexpordered)$decileL3expn==1)

legenes<-ChenTSSshapexpordered[j]
legenes
mcols(hegenes)<-NULL
mcols(legenes)<-NULL
export(legenes,"ChenTSSshapelowexpgenes.gtf","gtf")
export(hegenes,"ChenTSSshapehighexpgenes.gtf","gtf")
#i<-which(seqnames(ChenTSSshapexpordered)=="chrX")

#xChr<-ChenTSSshapexpordered[i]

#j<-which(seqnames(ChenTSSshapexpordered)!="chrX")
#aChr<-ChenTSSshapexpordered[j]

#export(xChr,"ChenTSSshapexporderedchrX.gtf","gtf")
#export(aChr,"ChenTSSshapexporderedchrA.gtf","gtf")
#export(ChenTSSshapexpordered, "ChenTSSshapexpordered.gtf", "gtf")
##################################################################
#how to get TES sites from genes list
#Get genomic ranges of genes

genes_A <- import.gff("ws260A.gtf")
genes_X <- import.gff("ws260X.gtf")

#Construct GenomicRanges for the TES
TES_A <- flank(genes_A, width= 1,  start=FALSE)
TES_X <- flank(genes_X, width= 1,  start=FALSE)

#Extend the genomic ranges by 1000 bp on either side, keeping the names
flanking_TES_A <- resize(TES_A, 2000, fix="center", use.names=TRUE)
flanking_TES_X <- resize(TES_X, 2000, fix="center", use.names=TRUE)

#alltss_chrA<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/alltss_chrA.gtf")
#alltss_chrA
#alltss_chrX<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/alltss_chrX.gtf")
#alltss_chrX

#Subsetting TES Genomic ranges based on the gene name
TES_from_alltss_chrA <- flanking_TES_A[flanking_TES_A$Name %in% alltss_chrA$ID]
TES_from_alltss_chrA
TES_from_alltss_chrX <- flanking_TES_X[flanking_TES_X$Name %in% alltss_chrX$ID]
TES_from_alltss_chrX

export(TES_from_alltss_chrA, "TESws260chrA.gtf", "gtf")
export(TES_from_alltss_chrX, "TESws260chrX.gtf", "gtf")

#for epdtss

epdtsschrA<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/epdtss_chrA.gtf")
epdtsschrA
epdtsschrX<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/epdtss_chrX.gtf")
epdtsschrX

TESfromepdtsschrA<-flanking_TES_A[flanking_TES_A$Name %in% epdtsschrA$ID]
TESfromepdtsschrA
TESfromepdtsschrX<-flanking_TES_X[flanking_TES_X$Name %in% epdtsschrX$ID]
TESfromepdtsschrX

export(TESfromepdtsschrA, "TESepdchrA.gtf", "gtf")
export(TESfromepdtsschrX, "TESepdchrX.gtf", "gtf")

#for kruesitss

KruesitsschrA<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/TSSs_from_papers/kreustsschrA.gtf")
KruesitsschrA
KruesitsschrX<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/TSSs_from_papers/kreustsschrX.gtf")
KruesitsschrX

TESfromKruesitsschrA<-flanking_TES_A[flanking_TES_A$Name %in% KruesitsschrA$WBgeneID]
TESfromKruesitsschrA
TESfromKruesitsschrX<-flanking_TES_X[flanking_TES_X$Name %in% KruesitsschrX$WBgeneID]
TESfromKruesitsschrX

export(TESfromKruesitsschrA, "TESKruesichrA.gtf", "gtf")
export(TESfromKruesitsschrX, "TESKruesichrX.gtf", "gtf")

#from chentss
ChentsschrA<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/TSSs_from_papers/chentsschrA.gtf")
ChentsschrA
ChentsschrX<-import("/Users/imac/Desktop/Supercoilinganalysis/featurefiles/TSSs_from_papers/chentsschrX.gtf")
ChentsschrX
table(ChentsschrA$assignmentType) #to select only the wormbase TSSs, exported file with "2" at the end)
ChentsschrAass<-ChentsschrA[mcols(ChentsschrA)$assignmentType=="wormbase_tss"]
ChentsschrAass
table(ChentsschrX$assignmentType) #to select only the wormbase TSSs
ChentsschrXass<-ChentsschrX[mcols(ChentsschrX)$assignmentType=="wormbase_tss"]
ChentsschrXass

TESfromChentsschrAass<-flanking_TES_A[flanking_TES_A$Name %in% ChentsschrAass$WbgeneName]
TESfromChentsschrAass
TESfromChentsschrXass<-flanking_TES_X[flanking_TES_X$Name %in% ChentsschrXass$WbgeneName]
TESfromChentsschrXass

export(TESfromChentsschrAass, "TESfromchenchrA2.gtf", "gtf")
export(TESfromChentsschrXass, "TESfromchenchrX2.gtf", "gtf")
export(ChentsschrAass, "TSSfromChenchrA2.gtf", "gtf")
export(ChentsschrXass, "TSSfromChenchrX2.gtf", "gtf")
########################################################
#making high versus low of TESs
ce3<-readRDS("./GR_ws260gene_L3expnDecile.RDS")
ce3

highlyexpgenes<-which(mcols(ce3)$expnDecile==10)

highlyexpgenes<-ce3[highlyexpgenes]
highlyexpgenes

lowlyexpgenes<-which(mcols(ce3)$expnDecile==1)
lowlyexpgenes<-ce3[lowlyexpgenes]
lowlyexpgenes

#highlyexpgenes<-import("highexpggeness.gtf")
#highlyexpgenes
#lowlyexpgenes<-import("lowexpggeness.gtf")
#lowlyexpgenes

TES_high <- flank(highlyexpgenes, width= 1,  start=FALSE)
TES_low<- flank(lowlyexpgenes, width= 1,  start=FALSE)

flanking_TES_high <- resize(TES_high, 2000, fix="center", use.names=TRUE)
flanking_TES_low <- resize(TES_low, 2000, fix="center", use.names=TRUE)

TES_from_highexpggeness <- flanking_TES_high[flanking_TES_high$Name %in% highlyexpgenes$Name]
TES_from_highexpggeness
TES_from_lowexpggeness <- flanking_TES_low[flanking_TES_low$Name %in% lowlyexpgenes$Name]
TES_from_lowexpggeness

export(TES_from_highexpggeness, "TESfromhighexpgenes.gtf", "gtf")
export(TES_from_lowexpggeness, "TESfromlowexpgenes.gtf", "gtf")
