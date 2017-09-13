#===============================================================================
#       Load libraries
#===============================================================================

library(DEXSeq)
library("BiocParallel")
BPPARAM=MulticoreParam(12)
register(BPPARAM)
library(data.table)

# library(ggplot2)
# library(Biostrings)
# library(devtools)
# load_all("~/pipelines/metabarcoding/myfunctions") # this contains various R scripts for plotting graphs
# library(dplyr)
# library(naturalsort)
# library(tibble)

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

# Sample data 	    
colData <- read.table("colData",sep="\t",row.names=1,header=T)

# extract counts from m or load them from a file (countData will need to converted to a data frame before calling DEXSeqDataSet)	    
countData <- m[,c(7:length(m)),with=F] # or
# countData <- read.table("countData",sep="\t",header=T,row.names=1) 

# reorder countData columns to same order as colData rows and convert to data.frame 
countData <- as.data.frame(countData[,row.names(colData),with=F])

# get the gene data	    
geneData <-  m[,c(1,2:6),with=F]
	    
# add an exon column to geneData by counting occurance of each gene (ordered by position)
geneData[order(Chr,Start,Geneid), Exonid := paste0("exon_",seq_len(.N)), by = Geneid]

# add row ranges (optional)  
# featureRanges <- GRanges(m$Geneid,IRanges(m$Start,as.numeric(m$End)),m$Strand)
# the above (which would need a chr name as well) might be useful depending on what your're doing 
# but m$Geneid should probably be changed to the chromosome/scaffold/contig name - featureCounts.gtf/DESex_final.gff contain the info	    
featureRanges <- GRanges(m$Geneid,IRanges(1,as.numeric(m$Length)),m$Strand)

# get the transcript ids from featureCounts.gtf	    
transcripts <- fread("grep exon featureCounts.gtf")	  
transcripts <- gsub("\";.*|transcripts \"","",transcripts$V9)	    

#==========================================================================================
#      DEXSeq analysis simple
##=========================================================================================

##### NOTE #####
# estimateDispersions and testForDEU can both be parallelsed by adding 
# BPPARAM=BPPARAM
# BUT parallisation did not work when I was tesing using MulticoreParam 		    
# This is unfortunate as dispersion estimates may take hours to calculate
	    
#### NOTE 2 #####	    
# If you have technical replicates, use the procedure as per DESeq2 to combine (sum) them into single rpelicates detailed below	    

### technical replicates only ###	    

dds <- 	DESeqDataSetFromMatrix(countData,colData,~1)
	    
# add grouping factor to identify technical replicates where sample contains the replicate info   
dds$groupby <- paste(dds$condition,dds$sample,sep="_")

# sum (collapse) replicates (collapseReplicates can be modified if you need mean or some other "summed" values)    
dds <- collapseReplicates(dds,groupby=dds$groupby)

# output summed data  
countData <- assay(dds)	 

# output collapsed sample data	    
colData <- as.data.frame(colData(dds))	 
	    
### end technical replicates ###	 
	    
# design formula
full_design = ~ sample + exon +exon:condition 
	    
# create DEXSeq object     
dxd <- DEXSeqDataSet(countData,
		     colData,
		     design=full_design,
		     featureID=geneData$Exonid,
		     groupID=geneData$Geneid,
		     featureRanges=featureRanges,
		     transcripts=transcripts
)

# calculate size factors 
sizeFactors(dxd) <- sizeFactors(estimateSizeFactors(dxd))

# calculate dispersions    
dxd <- dispersions(estimateDispersions(dxd))
	    
# reduced model	    
reduced_model <- ~ sample + exon

# test for differential exon usage using LRT	    
dxd <- testForDEU( dxd, fullModel = design(dxd), reducedModel =  reduced_model)

# add log fold changes (per sample vs sample 1) to results (I don't find this that useful) 	    
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")	    

# Get the results table
dxr1 <- DEXSeqResults(dxd)	    

    

