#===============================================================================
#       Load libraries
#===============================================================================

library(DEXSeq)
library("BiocParallel")
register(MulticoreParam(12))
library(ggplot2)
library(Biostrings)
library(devtools)
load_all("~/pipelines/metabarcoding/myfunctions") # this contains various R scripts for plotting graphs
library(data.table)
library(dplyr)
library(naturalsort)
library(tibble)

#===============================================================================
#       Load featureCounts data 
#===============================================================================

# load tables into a list of data tables - "." should point to counts directory, e.g. "counts/."
qq <- lapply(list.files(".",".*.txt$",full.names=T,recursive=F),function(x) {fread(x,skip=1)}) # don't need the {} round fread(x), but it messes with the git tabbing if you don't 

# rename the sample columns (7th column in a feature counts table, saved as the path to the BAM file)
# in the below this removes everything after the first dot in the 7th column
invisible(lapply(seq(1:length(qq)), function(i) colnames(qq[[i]])[7]<<-sub("\\..*","",colnames(qq[[i]])[7])))

# merge the list of data tables into a single data table
m <- Reduce(function(...) merge(..., all = T,by=c("Geneid","Chr","Start","End","Strand","Length")), qq)

# output "countData"
write.table(m[,c(1,7:(ncol(m))),with=F],"countData",sep="\t",na="",quote=F,row.names=F) 
	    
# output gene details
write.table(m[,1:6,with=F],"genes.txt",sep="\t",quote=F,row.names=F) 

#==========================================================================================
#      Read in data 
##=========================================================================================

colData <- read.table("colData",sep="\t",row.names=1,header=T)

countData <- m[,c(7:length(m)),with=F] # or
# countData <- read.table("countData",sep="\t",header=T,row.names=1) 

# reorder countData columns to same order as colData rows (probably not necessary)
countData <- as.data.frame(countData[,row.names(colData),with=F])

# get the gene data	    
geneData <-  m[,c(1,2:6),with=F]
	    
# add an exon column to geneData by counting occurance of each gene (ordered by position)
geneData[order(Chr,Start,Geneid), Exonid := paste0("exon_",seq_len(.N)), by = Geneid]
    
# design formula
full_design = ~ sample + condition+exon:condition 

# add row ranges    
# featureRanges <- GRanges(m$Geneid,IRanges(m$Start,as.numeric(m$End)),m$Strand)
# the above (which would need a chr name as well) is not necessary  - we're not looking for chimaeras	    
featureRanges <- GRanges(m$Geneid,IRanges(1,as.numeric(m$Length)),m$Strand)

# get the transcript ids from the featureCounts.gtf (or DEXSeq_final.gff - the grep would need to change)	    
transcripts <- fread("grep exon featureCounts.gtf")	  
transcripts <- gsub("\";.*|transcripts \"","",transcripts$V9)	    

#==========================================================================================
#      DEXSeq analysis
##=========================================================================================
	    
# create DEXSeq object	    
dxd <- DEXSeqDataSet(countData,
		     colData,
		     design=full_design,
		     featureID=geneData$Exonid,
		     groupID=geneData$Geneid,
		     featureRanges=featureRanges,
		     transcripts=transcripts)

sizeFactors(dxd) <- sizeFactors(estimateSizeFactors(dxd))
dispersions(dxd) <- dispersions(estimateDispersions(dxd))

reduced_model <- ~ sample + exon
	    
dxd <- testForDEU( dxd,
fullModel = design(dxd),
reducedModel =  reduced_model,BPPARAM=BPPARAM)
