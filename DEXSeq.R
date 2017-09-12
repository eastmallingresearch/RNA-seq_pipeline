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
qq <- lapply(list.files(".",".*.txt$",full.names=T,recursive=F),function(x) fread(x)) 

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
#       Read pre-prepared colData,  countData and annotations
##=========================================================================================

colData <- read.table("colData",header=T,sep="\t")

countData <- m[,c(1,7:length(m)),with=F] # or
# countData <- read.table("countData",sep="\t",header=T,row.names=1) 

# reorder countData columns to same order as colData rows
countData <- countData[,colData$SampleID] 

# design formula
design = ~ 1

flattenedfile = fread("DEXSeq_final.gff",header=F)



dxd <- DEXSeqDataSet(countData,sampleData=colData,design=design,featureID=exons,groupID=genesrle,featureRanges,transcripts)





	    
