## analysis for oryza sativa
# KEGG species names can be found here:
# https://www.genome.jp/kegg/catalog/org_list.html

# NOTE older version of KEGGREST (<1.36) have the incorrect path in .getRootUrl
# either hack up your own version of the function (and link it to KEGGREST) or get the latest version..

# This is copied from : https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html

library(KEGGREST)

# Get the KEGG pathways for rice
pathways.list <- keggList("pathway", "osa")

# get the gene names for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
    function(pwid){
        pw <- keggGet(pwid)
        if (is.null(pw[[1]]$GENE)) return(NA)
        pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
        pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
        return(pw2)
    }
   
)
head(genes.by.pathway)

# get sig genes from test set
test.genes <- res.merge[!is.na(geneID),c(20,2:7)]
test.genes[,geneID:=sub("(GeneID:)([0-9]*)(,.*)","\\2",geneID)]

# Fisher exact test for 2x2 contigency table
#                  sig  |  !sig  | Total
# total genes   |   a   |   b    | a + b 
# pathway genes |   c   |   d    | c + d
# Total         | a + c | b  + d | a + b + c +d

                      

fisher.fun <- function(pathway.genes,exp.genes,padj="padj",cutoff=0.05,geneID="geneID") {
  A=sum(exp.genes[[padj]]<=cutoff,na.rm=T)
  B=nrow(exp.genes)-a
  C=pathway.genes%in%exp.genes[[geneID]][exp.genes[[padj]]<=cutoff]# get these genes
  D=length(pathway.genes)-c
  DF<-data.frame(t=c(A,sum(C)),f=c(B,D))
  f<-fisher.test(DF)
  list(f,DF,pathway.genes[C])   
   # c
}
    
    
wlicox.fun <-function() {
        list.genes.in.pathway <- intersect(test.genes$geneID,pathway.genes)
        list.genes.not.in.pathway <- setdiff(test.genes$geneID, list.genes.in.pathway)
        scores.in.pathway <- geneList[list.genes.in.pathway]
        scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
        if (length(scores.in.pathway) > 0){
            p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
        } else{
            p.value <- NA
        }
        return(c(p.value = p.value, Annotated = length(list.genes.in.pathway)))    
}
pVals.by.pathway <- lapply(genes.by.pathway,fisher.fun,exp.genes)

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
head(outdat)
