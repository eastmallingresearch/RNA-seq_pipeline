#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("~/pipelines/RNA_seq/scripts/myfunctions")
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)

#===============================================================================
#       Load features counts data 
#===============================================================================

# load tables into a list of data tables - "." should point to counts directory, e.g. "counts/."
qq <- lapply(list.files(".",".*.txt$",full.names=T,recursive=F),function(x) fread(x)) 

# rename the sample columns (7th column in a feature counts table, saved as the path to the BAM file)
# in the below I'm saving the 8th ([[1]][8]) path depth (which was the informative folder name containg the BAM file)
invisible(lapply(seq(1:length(qq)), function(i) colnames(qq[[i]])[7]<<-strsplit(colnames(qq[[i]])[7],"\\/")[[1]][8]))

# merge the list of data tables into a single data table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

# output "countData"
write.table(m[,c(1,7:(ncol(m))),with=F],"countData",sep="\t",na="",quote=F,row.names=F) 
	    
# output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F) 

#==========================================================================================
#       Read pre-prepared colData,  countData and annotations
##=========================================================================================

# Read sample metadata	    
colData <- read.table("colData",header=T,sep="\t")

# Get counts from m, or load 	    
countData <- data.frame(m[,c(7:length(m)),with=F],row.names=m$Geneid) #or	    
#countData <- fread("countData",sep="\t",header=T,row.names=1) 
	    
countData <- countData[,ccolData$SampleID] # reorder countData columns to same order as colData rows 

#===============================================================================
#       DESeq2 analysis
#		Set alpha to the required significance level. This also effects how
#		DESeq calculated FDR - setting to 0.05 and then extracting results with a
#		significance below 0.01 will give slightly different results form setting
#		alpha to 0.01
#================================================================================

# create DESeq object from count and column data
dds <- 	DESeqDataSetFromMatrix(countData,colData,~1) 
	    
#### Technical replicates only ####
	    
# add grouping factor to identify technical replicates	    
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

# sum replicates (must use same library or library size correction will go wonky)	    
dds <- collapseReplicates(dds,groupby=dds$groupby)

#### end technical replicates ####   
	    	    
# normalise counts for different library size (do after collapsing replicates)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds)) 

# define the DESeq 'GLM' model	    
design=~condition

# add design to DESeq object	    
design(dds) <- design

# Run the DESeq statistical model	    
dds <- DESeq(dds,parallel=T)

# set the significance level for BH adjustment	    
alpha <- 0.05

# calculate the differences - uses the "levels" of the condition factor as the third term for the contrast
# res is a list object containing the DESeq results objects for each contrast
# contrast=c("condition","RH1","RH2") etc. (the below just runs through all of the different sample types (excluding RH1))
res <- results(dds,alpha=alpha)

# get annotations (if you have any)	    

annotations <- read.table("annotations.txt", sep="\t",header=T)	    
	    
# merge results with annotations
res.merged <- left_join(rownames_to_column(as.data.frame(res)),annotations,by=c("rowname"="query_id"))
	
# get significant results
sig.res <- subset(res.merge,padj<=alpha)

# write tables of results
write.table(res.merged,"results.txt",quote=F,na="",row.names=F,sep="\t")
	
#===============================================================================
#       FPKM
#===============================================================================

rowRanges(dds) <- GRangesList(apply(m,1,function(x) GRanges(x[[1]],IRanges(1,as.numeric(x[[6]])),"+")))
myfpkm <- data.table(GeneID=m[,1],length=m[,6],fpkm(dds,robust=T))
write.table(myfpkm,"fpkm.txt",quote=F,na="",sep="\t")
	
#===============================================================================
#       Heirachical clustering
#===============================================================================

clus <- function(X,clusters=10,m=1,name="hclust.pdf") {
	if (m==1) {d <- dist(X, method = "manhattan")}
	else if (m==2) {d <- dist(X, method = "euclidean")}
	else if (m==3) {d <- dist(X, method = "maximum")}
	else if (m==4) {d <- dist(X, method = "canberra")}
	else if (m==5) {d <- dist(X, method = "binary")}
	else d <- {dist(X, method = "minkowski")}
	hc <- hclust(d, method="ward")
	groups <- cutree(hc, k=clusters) # cut tree into n clusters
	pdf(name,height=8,width=8)
	plot(hc)
	rect.hclust(hc,k=clusters)
	dev.off()
	return(list(hc,groups,d))
}

#===============================================================================
#       Graphs
#===============================================================================
	
# PCA 1 vs 2 plot
vst <- varianceStabilizingTransformation(dds,blind=F,fitType="local")
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH1"] <- "02780"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH2"] <- "02793"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH3"] <- "F55"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH4"] <- "10170"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH5"] <- "MWT"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH6"] <- "MOL"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH7"] <- "MKO"
levels(vst@colData$condition)[levels(vst@colData$condition)=="RH8"] <- "TJ"

# calculate PCs				    
mypca <- prcomp(t(assay(vst)))
				    
# calculate variance for each PC
mypca$percentVar <- mypca$sdev^2/sum(mypca$sdev^2)
				    
# create data frame of PCs x variance (sets PCA plot axes to same scale)
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

# set pdf 
pdf("quorn.pca.pdf",height=8,width=8)

# PCA/ordination plotting function 
plotOrd(df,vst@colData,design="condition",xlabel="PC1",ylabel="PC2", pointSize=3,textsize=14)

dev.off()
	
# MA plots	
pdf("MA_plots.pdf")

# plot_ma is an MA plotting function 				    
lapply(res.merged,function(obj) {
	plot_ma(obj[,c(1:5,7]),xlim=c(-8,8))
})
dev.off()
