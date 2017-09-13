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
countData <- m[order(Chr,Start,Geneid),c(1,7:length(m)),with=F] # or
# countData <- read.table("countData",sep="\t",header=T,row.names=1) 

# reorder countData columns to same order as colData rows and convert to data.frame 
countData <- countData[,c("Geneid",row.names(colData)),with=F]

# get the gene data	    
geneData <-  m[order(Chr,Start,Geneid),c(1,2:6),with=F]
# geneData <- fread(genes.txt) # this MUST be a data.table
	    
# add an exon column to geneData by counting occurance of each gene (ordered by position)
geneData[, Exonid := paste0("exon_",seq_len(.N)), by = Geneid]

# set GRanges for each exon (optional - but set featureRanges to NULL otherwise)	    
featureRanges <- GRanges(geneData$Chr,IRanges(geneData$Start,as.numeric(geneData$End)),geneData$Strand)	    

# retrieve transcipts (optional - but set transcripts to NULL otherwise)	    
features <- fread("grep exon featureCounts.gtf")	  
# get the gene id
features$Geneid <- gsub("\"|.*gene_id \"","",features$V9)	    
# get transcript id
features$Transcriptid <- gsub("\";.*|transcripts \"","",features$V9)	    
# left join geneData and features
geneData[features[,Geneid,Transcriptid],Transcriptid:=i.Transcriptid,on="Geneid"]

transcripts <- geneData$Transcriptid    

#==========================================================================================
#      DEXSeq analysis simple (one factor)
##=========================================================================================

##### NOTE #####
# estimateDispersions and testForDEU can both be parallelsed by adding 
# BPPARAM=BPPARAM
# BUT parallisation did not work when I was tesing using MulticoreParam 		    
# This is unfortunate as dispersion estimates may take hours to calculate
	    
#### NOTE 2 #####	    
# If you have technical replicates, use the procedure as per DESeq2 to combine (sum) them into single rpelicates detailed below	    

### technical replicates only ###	    

dds <- 	DESeqDataSetFromMatrix(countData[,-1],colData,~1)
	    
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
full_design <- ~ sample + exon + condition:exon 
	    
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
reduced_design <- ~ sample + exon

# test for differential exon usage using LRT	    
dxd <- testForDEU( dxd, fullModel = design(dxd), reducedModel =  reduced_design)

# add log fold changes (per sample vs sample 1) to results (I don't find this that useful) 	    
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")	    

# Get the results table
dxr1 <- DEXSeqResults(dxd)	    

#==========================================================================================
#      DEXSeq analysis multiple factors
##=========================================================================================    

# A second factor could be something like batch or block 
# The only change needed is in the models to add a factor:exon interation (which stays in both full and reduced designs)
	    
full_model <- ~ sample + exon + block:exon + condition:exon 	    
reduced_model <- ~ sample + exon + block:exon
	    
# the rest of the analysis is as per one factor
	    
	    
# this is quite nice - 
# compare results from using one factor to multiple factor design     
# then table them out
dxr2 = DEXSeqResults( dxd ) # dxr1 calculated from one factor design
table( before = dxr1$padj < 0.1, now = dxr2$padj < 0.1 )	    
	    
	    
	    
