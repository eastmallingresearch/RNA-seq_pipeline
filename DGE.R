#===============================================================================
#       Load libraries
#===============================================================================

library(DESeq2)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("~/pipelines/RNA-seq/scripts/myfunctions")
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)

#===============================================================================
#       Load data from fetureCounts
#===============================================================================

# load tables into a list of data tables - "." should point to counts directory, e.g. "counts/."
# depending on how your count files are named and stored (in subdirectories or not) the .*txt$ regex and recursive flag may need editing
# the example below assumes all count files have file names ending in .txt and are all in a single folder (featureCounts standard)
qq <- lapply(list.files(".",".*txt$",full.names=T,recursive=F),function(x) fread(x)) 

# rename the sample columns (7th column in a feature counts table, saved as the path to the BAM file)
# in the below I'm saving the 8th ([[1]][8]) path depth (which was the informative folder name containg the BAM file)
invisible(lapply(seq(1:length(qq)), function(i) colnames(qq[[i]])[7]<<-strsplit(colnames(qq[[i]])[7],"\\/")[[1]][8]))

# merge the list of data tables into a single data table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

# output "countData"
# write.table(m[,c(1,7:(ncol(m))),with=F],"countData",sep="\t",na="",quote=F,row.names=F) 
	    
# output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F) 

# Get counts from m
countData <- data.frame(m[,c(7:length(m)),with=F],row.names=m$Geneid)	    
	    

#===============================================================================
#       Load data from SALMON quasi mapping
#===============================================================================

library(tximport)
library(rjson)
library(readr)

# import transcript to gene mapping info
tx2gene <- read.table("trans2gene.txt",header=T,sep="\t")

# import quantification files
# txi.reps <- tximport(list.files(".","quant.sf",full.names=T,recursive=T),type="salmon",tx2gene=tx2gene,txOut=T)
# ...change of plan (it's easier to use list.dirs) though the below is very cool so I've left it here    
# mysamples <- rapply(strsplit(list.files(".","quant.sf",full.names=T,recursive=T),"\\/"), f=`[[`, ...=2, how="unlist")

# import quantification files - will not work if any additional directories are in the specified path
txi.reps <- tximport(paste(list.dirs(".",full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)	    
	    
# get the sample names from the folders	    
mysamples <- list.dirs(".",full.names=F,recursive=F)

# summarise to gene level (this can be done in the tximport step, but is easier to understand in two steps)
txi.genes <- summarizeToGene(txi.reps,tx2gene)

# set the sample names for txi.genes
invisible(sapply(seq(1,3), function(i) colnames(txi.genes[[i]])<<-mysamples))

#==========================================================================================
#       Read pre-prepared sample metadata and annotations
#=========================================================================================

# Read sample metadata	    
colData <- read.table("colData",header=T,sep="\t",row.names=1)
    
# reorder colData for salmon 		 
colData <- colData[mysamples,]
		 
# reorder colData for featureCounts		 
colData <- colData[colnames(countData),]
		 
# get annotations (if you have any)	    
annotations <- read.table("annotations.txt", sep="\t",header=T)	   		 
#===============================================================================
#       DESeq2 analysis
#		Set alpha to the required significance level. This also effects how
#		DESeq calculated FDR - setting to 0.05 and then extracting results with a
#		significance below 0.01 will give slightly different results form setting
#		alpha to 0.01
#================================================================================

# create DESeq object from featureCounts counts and sample metadata
dds <- 	DESeqDataSetFromMatrix(countData,colData,~1) 

# or create from Salmon counts and sample metadata		 
dds <- DESeqDataSetFromTximport(txi.genes,colData,~1)	    
	    
#### Technical replicates only ####
	    
# add grouping factor to identify technical replicates	    
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

# sum replicates (must use same library or library size correction will go wonky)	    
dds <- collapseReplicates(dds,groupby=dds$groupby)

#### end technical replicates ####   
	    	    
# featureCounts only - normalise counts for different library size (do after collapsing replicates)
# NOTE: need to work out what to do if there are technical replicates for salmon workflow
# probably take average of avgTxLength for the summed samples 
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
# contrast=c("condition","S","H") etc - ?results for many more options
res <- results(dds,alpha=alpha)
    
# merge results with annotations
res.merged <- left_join(rownames_to_column(as.data.frame(res)),annotations,by=c("rowname"="query_id"))
	
# get significant results
sig.res <- subset(res.merge,padj<=alpha) 
# sig.res <- res[which(res$padj<=alpha),] # use this is you don't have annotations		 

# write tables of results
write.table(res.merged,"results.txt",quote=F,na="",row.names=F,sep="\t")

# get sequences of significant transcripts - transcripts.fa is the file of transcripts
# this will only work if the fastas are over two lines, i.e. not split every 80 bases (and the names match)
# the below shell script is a method for converting to 2 line fasta
# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' split_transcripts.fa|sed '1d' > transcripts.fa
# seqs <- DNAStringSet(sapply(rownames(sig.res),function(s) {
#	DNAString(system2("grep",c(paste0(s," "),"-A1", "transcripts.fa"),stdout=T)[[2]])
# }))
# alternatively just load the whole of the transcript file (split every 80 bases or not) into a Biostrings object and subset by rowname:
seqs <- readDNAStringSet("transcripts.fa")[rownames(sig.res)]	# o.k this is probably a better method unless transcripts.fa is massive

		 
	
#===============================================================================
#       FPKM
#===============================================================================
# note this is pointless if using salmon or some other pseudo counting aligner
	    
# set GRanges for each gene
# fast method - requires some editing	    
rowRanges(dds) <- GRanges(geneData$Chr,IRanges(geneData$Start,as.numeric(geneData$End)),geneData$Strand)	    	    
# slow method - doesn't require editing	    
rowRanges(dds) <- GRangesList(apply(m,1,function(x) GRanges(x[[1]],IRanges(1,as.numeric(x[[6]])),"+")))

# calculate FPKM values
myfpkm <- data.table(GeneID=m[,1],length=m[,6],fpkm(dds,robust=T))
				    
# write FPKM values				    
write.table(myfpkm,"fpkm.txt",quote=F,na="",sep="\t")
	
#===============================================================================
#       Heirachical clustering
#===============================================================================

# this is out of date - ward.D2 is the prefered method for clustering - actaully the whole function is a bit naff
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

# calculate PCs				    
mypca <- prcomp(t(assay(vst)))
				    
# calculate variance for each PC
mypca$percentVar <- mypca$sdev^2/sum(mypca$sdev^2)
				    
# create data frame of PCs x variance (sets PCA plot axes to same scale)
df <- t(data.frame(t(mypca$x)*mypca$percentVar))

# plotOrd is a PCA/ordination plotting function 
ggsave("pca.pdf",plotOrd(df,vst@colData,design="condition",xlabel="PC1",ylabel="PC2", pointSize=3,textsize=14))

# MA plots	
pdf("MA_plots.pdf")

# plot_ma is an MA plotting function 				    
lapply(res.merged,function(obj) {
	plot_ma(obj[,c(1:5,7]),xlim=c(-8,8))
})
dev.off()
