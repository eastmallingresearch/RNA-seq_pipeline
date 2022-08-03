#===============================================================================
#       Load data
#===============================================================================

# read in pfam annotation

annotation <- fread(paste0(home,"/pipelines/common/resources/pfam/names.txt"))
setnames(annotation,"NAME","PFAM_NAME")
pfam_go <- fread(paste0(home,"/pipelines/common/resources/mappings/pfam_go_map"),header=F)

#### bins ####
countData <- fread(paste0(SITE,".countData"))

#setnames(countData,names(countData),sub("(_ND.*_L)([0-9]*)(.*)","_\\2",names(countData)))
setnames(countData,"DOMAINS","PFAM_NAME")

# get number of reads mapped to domains       
print(as.data.frame(colSums(countData[,-1])))

# map bins to pfam and go terms
mapping_pfam <- annotation[countData[,1,with=F],on="PFAM_NAME"]

### end bins ###


### sub bins ###
countData <- fread(paste0(SITE,".bins.countData"))
setnames(countData,names(countData),sub("_L","_",names(countData)))
countData[,PFAM_NAME:=gsub("(k[0-9]+_[0-9]+_)(.*)(_[0-9]+_[0-9]+$)","\\2",SUB_BIN_NAME)]
unsetNA(countData)

mapping_pfam <- annotation[countData[,c(1,ncol(countData)),with=F],on="PFAM_NAME"]

# remove unused columns from countData (BIN_ID can be mapped to mapping_pfam)
colsToDelete <- c("PFAM_NAME") 
countData[, (colsToDelete) := NULL]

### end sub bins###


# annotation
mapping_go   <- copy(mapping_pfam)
mapping_go$ACC <- sub("\\..*","",mapping_go$ACC)
mapping_go <- data.table(left_join(mapping_go,pfam_go,by=(c("ACC"="V1"))))
mapping_go <- mapping_go[complete.cases(mapping_go),]




#===============================================================================
#       Functional analysis (sub bins) 
#===============================================================================
library(topGO)

### bins ###
res_filt <- res_merge[mapping_go,on="PFAM_NAME"]
res_filt <- res_filt[complete.cases(res_filt),]
fwrite(res_filt[,toString(V4),by=list(PFAM_NAME)],"topgo_temp",sep="\t",row.names=F,col.names=F,quote=F)
geneID2GO <- readMappings("topgo_temp")
genes <- unique(res_filt[,c("PFAM_NAME","log2FoldChange","padj")]) 
geneList <- setNames(genes$padj*sign(genes$log2FoldChange),genes$PFAM_NAME)


### sub bins ###
res_filt <- res_merge[mapping_go,on="SUB_BIN_NAME"]
res_filt <- res_filt[complete.cases(res_filt),]
fwrite(res_filt[,toString(V4),by=list(SUB_BIN_NAME)],"topgo_temp",sep="\t",row.names=F,col.names=F,quote=F)
geneID2GO <- readMappings("topgo_temp")
genes <- unique(res_filt[,c("SUB_BIN_NAME","log2FoldChange","padj")]) 
geneList <- setNames(genes$padj*sign(genes$log2FoldChange),genes$SUB_BIN_NAME)


geneSel <- function(X)abs(X)<=0.1

GOdata <- new("topGOdata",ontology = "BP",allGenes = geneList,geneSel = geneSel,annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = 5)

x="Over"
geneSelectionFun(GOdata) <- function(X)abs(X)<=0.1&X>0 # increased expression

x="Under"
geneSelectionFun(GOdata) <- function(X)abs(X)<=0.1&X<0 # decreased expression

# weighted uses the go topology to infer statistical significance
# ks uses strength of signal (i.e. p vlaue) for calculating GO term significance (I'm not convinced this is a good idea)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
#resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultFisher.weight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher, fisherWeighted=resultFisher.weight,orderBy = "fisherWeighted", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))
#allRes <- GenTable(GOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, ksWeighted=resultKS.weight,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = length(GOdata@graph@nodes))
# over_expressed <- allRes[((allRes$Significant)/(allRes$Expected))>1,]

fwrite(allRes,paste(SITE,x,"AOD_GO_RES_SUBBINS_0.5.txt",sep="_"),sep="\t",quote=F)

# get domain list for sig GO terms
qf <- function(GOID,res,GOdata){
  selection<-unlist(genesInTerm(GOdata,GOID))
  res_merge[SUB_BIN_NAME%in%selection,c(1:3,7:10)]
}

test <- lapply(allRes[as.numeric(allRes$fisherWeighted)<=0.1,1],qf,res_merge,GOdata)
names(test) <- allRes[as.numeric(allRes$fisherWeighted)<=0.1,2]
lapply(names(test),function(n) test[n][[1]][,GO:=n])
#test[[1]][,GO:=NA]
OUTPUT <- do.call(rbind,lapply(test,subset,padj<=0.1))

fwrite(OUTPUT,paste(SITE,x,"_AOD_SUBBINS_0.5_OUTPUT.txt",sep="_"),sep="\t",quote=F)

# GO graphs
pdf(paste(SITE,x,"_COD_SUBBINS_0.5_GO_plots.pdf",sep="_"))
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
# showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
# showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
showSigOfNodes(GOdata, score(resultFisher.weight), firstSigNodes = 5, useInfo = 'all')
dev.off()
