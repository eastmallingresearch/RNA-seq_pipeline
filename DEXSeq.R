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
#       Read pre-prepared colData,  countData and annotations
##=========================================================================================

colData <- read.table("colData",header=T,sep="\t")

countData <- m[,c(1,7:length(m)),with=F] # or
# countData <- read.table("countData",sep="\t",header=T,row.names=1) 

# reorder countData columns to same order as colData rows
countData <- countData[,colData$SampleID] 

geneData <-  m[,c(1,2:6),with=F]
	    
# design formula
design = ~ 1

flattenedfile = fread("DEXSeq_final.gff",header=F)

dxd <- DEXSeqDataSet(countData,sampleData=colData,design=design,featureID=exons,groupID=genesrle,featureRanges,transcripts)

DEXSeqDataSetFromFC
function (countData, sampleData, design = ~sample + exon + condition:exon,
    flattenedfile = NULL)
{
    if (!all(sapply(countfiles, class) == "character")) {
        stop("The countfiles parameter must be a character vector")
    }
    lf <- lapply(countfiles, function(x) read.table(x, header = FALSE,
        stringsAsFactors = FALSE))
    if (!all(sapply(lf[-1], function(x) all(x$V1 == lf[1]$V1))))
        stop("Count files have differing gene ID column.")
    dcounts <- sapply(lf, `[[`, "V2")
    rownames(dcounts) <- lf[[1]][, 1]
    dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_",
        ]
    rownames(dcounts) <- sub(":", ":E", rownames(dcounts))
    colnames(dcounts) <- countfiles
    splitted <- strsplit(rownames(dcounts), ":")
    exons <- sapply(splitted, "[[", 2)
    genesrle <- sapply(splitted, "[[", 1)
    if (!is.null(flattenedfile)) {
        aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE,
            header = FALSE)
        colnames(aggregates) <- c("chr", "source", "class", "start",
            "end", "ex", "strand", "ex2", "attr")
        aggregates$strand <- gsub("\\.", "*", aggregates$strand)
        aggregates <- aggregates[which(aggregates$class == "exonic_part"),
            ]
        aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
        aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1",
            aggregates$attr)
        transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1",
            aggregates$attr)
        transcripts <- strsplit(transcripts, "\\+")
        exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1",
            aggregates$attr)
        exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start,
            end = aggregates$end), strand = aggregates$strand)
        names(exoninfo) <- paste(aggregates$gene_id, exonids,
            sep = ":E")
        names(transcripts) <- rownames(exoninfo)
        if (!all(rownames(dcounts) %in% names(exoninfo))) {
            stop("Count files do not correspond to the flattened annotation file")
        }
        matching <- match(rownames(dcounts), names(exoninfo))
        stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
        stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
        dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons,
            genesrle, exoninfo[matching], transcripts[matching])
        return(dxd)
    }
    else {
        dxd <- DEXSeqDataSet(dcounts, sampleData, design, exons,
            genesrle)
        return(dxd)
    }
}



	    
