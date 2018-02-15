# Functional annotation for transcriptomes

## Interproscan 5

### Installation instructions
Interproscan can be a bit of a pain to get working...

Interproscan 5 requires gcc 4.8 or later (4.7 installed on cluster)
See https://github.com/eastmallingresearch/seq_tools/blob/master/install_progs/instructions.md for instructions for getting gcc upgraded
alternatively use gcc installed in Andy's home directory

These instructions were made using interproscan 5.27-66.0 

#### Download files
```shell
# download interproscan files
mkdir ~/build
cd build
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.27-66.0/interproscan-5.27-66.0-64-bit.tar.gz
tar -zxf interproscan-5.27-66.0-64-bit.tar.gz
```

#### install prereq libraries: 

##### sfld
```shell
cd ~/build/interproscan-5.27-66.0/src/sfld/
sed -i -e 's/^EASEL_DIR=/EASEL_DIR=~\/build\/hmmer-3.1b2-linux-intel-x86_64\/easel/' Makefile # might need to specify full path
make 
mv sfld_preprocess ../../../bin/sfld/.
mv sfld_postprocess ../../../bin/sfld/. 
```
##### cath-resolve-hits
```shell
cd ~/build/interproscan-5.27-66.0/bin/gene3d/4.1.0
wget "https://github.com/UCLOrengoGroup/cath-tools/releases/download/v0.16.2/cath-refine-align.centos6"
cd ../../../
echo -e "\n#cath-resolve settings\ncath.resolve.hits.switches.gene3d=--input-format=hmmer_domtblout" >>interproscan.properties
```

##### rpsbproc
rpsbproc is part of BLAST+. The below uses a prebuilt version which may or may not still be available
```shell
cd ~/build/interproscan-5.27-66.0/bin/blast/ncbi-blast-2.6.0+
wget "https://github.com/ebi-pf-team/interproscan/tree/master/core/jms-implementation/support-mini-x86-32/bin/blast/rpsbproc"
chmod rpsbproc 755

# You also need to create a rpsbproc.ini in the same location:
cat <<EOT>> ~/build/interproscan-5.27-66.0/bin/blast/ncbi-blast-2.6.0+/rpsbproc.ini
cdd.signature.library.release=3.16
cdd.signature.list.path=data/cdd/3.16/data/cddid.tbl
cdd.library.path=data/cdd/3.16/db/Cdd_NCBI
cdd.data.path=data/cdd/3.16/data
EOT
```
The below is the "correct" way to get rpsbproc to work (blast+ install) - requires gcc 4.8+
```shell
mkdir cddblast
cd cddblast
wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-src.tar.gz
tar xvzf ncbi-blast-2.6.0+-src.tar.gz
cd ncbi-blast-2.6.0+-src/c++/src/app/
wget -r --no-parent -l 1 -np -nd -nH -P rpsbproc ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/rpsbproc-src/
#edit Makefile.in and make sure SUB_PROJ is assigned two applications as follows: SUB_PROJ = blast rpsbproc
cd ../../
./compilers/unix/GCC.sh
make && make install 
#after compilation is complete (or update interprocscan.profile to installed file location)
cp ReleaseMT/bin/rpsblast ../../../bin/blast/ncbi-blast-2.6.0+/
cp ReleaseMT/bin/rpsbproc ../../../bin/blast/ncbi-blast-2.6.0+/
```

#### Compiling interproscan (borrowing various other precompiled bits from Andy's home directory)

```shell
cd ~/build/interproscan-5.27-66.0
./configure \
 --prefix=/$HOME/usr/local  \
 --with-gmp="/home/armita/prog/gmp/gmp/6.1.2" \
 --with-mpfr="/home/armita/prog/mpfr/mpfr-3.1.5" \
 --with-mpc="/home/armita/prog/mpc/mpc-0.8.2" \
 --disable-multilib

# change the paths to the location (version 4.8+) of the specified files 
make CC=path_to_gcc/gcc CPP=path_to_g++/g++ CXX=path_to_g++/g++-LD=path_to_g++/g++
```

#### Issues

Gene3D is still not working correctly
The translation option (-t n) is not working correctly (uses EMBOSS getorf which I couldn't get to work).
The included get_all_cds.pl (or get_longest_cds) and translate.pl can convert a nucleotide sequence(s) to protein:
```
get_all_cds.pl <mytranscriptome.fa|translate.pl|sed 's/\*//g' >myproteins.fa
```
It outputs all proteins at least 25 aa long - it's a bit slow mind. The sed gets rid of any leading or trailing stops (\*) as interproscan doesn't like them.   

### Running on cluster
Interproscan has a built in option to run on an SGE cluster - but I couldn't get it to work.  
Therefore I've written my own version which splits a transcriptome into multiple files and runs multiple simultaneous jobs (hardcoded to 10). Compared to the Interproscan inbuilt clustering this will take an additional core for each simultaneous job.  

To run interproscan on the cluster use interproscan.sh. This takes two required + any optional iterproscan arguments.
```shell
interproscan.sh FOLDER FILE <options>
```
FOLDER is the location of the transcriptome file(s) - all output files will be copied here.  
FILE is the name of the transcriptome file (accepts shell expandable wild cards to specify multiple files).  
At least one of the options *must* specify which analyses to perform (-appl) - otherwise the Gene3D analysis will cause the JVM to crash. Below is the list of possible analyses:  
 -appl CDD,COILS,HAMAP,MobiDBLite,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM  
 
 Other interesting options:  
 -iprlookup (provides mappings from matched member database signatures to the InterPro entries that they are integrated into)  
 -goterms  (provides mappings to the Gene Ontology - based on InterPro (IPR) matches)  
 -pa (KEGG, MetaCyc and Reactome pathway lookup - based on InterPro (IPR) matches)  
 -dra (switch off residue level annotations - speeds up analysis)  


The cluster job will create a tempory folder <TEMPFOLDER> in <FOLDER> which will contain the transcriptome split into multiple file containing 1000 sequences, and a file called <FILE>_splitfiles.txt which contains the path to each of the multiple files.
 
Interproscan can take a while to run, if the job needs to be restarted, edit <FILE>_splitfiles.txt removing all files which have already been processed, then run restart_interproscan.sh.  

```shell
restart_interproscan.sh FOLDER FILE TEMPFOLDER <options>
```

### Reshape interproscan results with R
Depending on the analysis options used for interproscan the below script might be useful for combining multiple lines per exon into a single line per gene. Works with tsv output files.

all.txt is a concatenated version of all output files (you could just dread(paste(list.dirs("counts",full.names=T,recursive=F),"/quant.sf",sep="")

```R
library(data.table)
library(dplyr)

# load ip annotations
annotations <- fread("cat *.tsv",header=T,fill=T) # if fread fails could try setting autostart to the line it fails on  (depedent on error)

# add colnames
colnames(annotations) <- c("PROT_ID","MD5","LENGTH","ANALYSIS","MATCH_ID","MATCH_DESC","START","STOP","E_VAL","STATUS","DATE","IPR_ID","IPR_DESC","GO_ID","PATHWAY")
# add gene data (from exon data)
annotations$GENE <- sub("\\..*","",annotations$PROT_ID)
annotations$EXON <- gsub("(.*\\.)([0-9]+)(_.*)","\\2",annotations$PROT_ID)
annotations$FRAME <- gsub("(.*\\.[0-9]+_)([0-9]+)(_.*)","\\2",annotations$PROT_ID)
annotations$ORF <- gsub(".*_","",annotations$PROT_ID)
annotations$DIRECTION <- gsub("(.*\\.[0-9]+_[0-9]+_)([A-Z]+)(_.*)","\\2",annotations$PROT_ID)
annotations$DIRECTION <- sub("F","\\+",annotations$DIRECTION)
annotations$DIRECTION <- sub("RC","\\-",annotations$DIRECTION)
write.table(annotations,"full_annotations.txt",sep="\t",row.names=F,quote=F,na="")

# subset annotations with GeneID an anaylses results
slim_annot <- annotations[,c(16,4,5,6,12,13,14,15)]
# remove MobiDBLite and Coils results (or keep if they're useful to you)
slim_annot <- slim_annot[ANALYSIS!="MobiDBLite"|ANALYSIS!="Coils",]
slim_annot$ANALYSIS <- as.factor(slim_annot$ANALYSIS)
write.table(slim_annot,"slim_annotations_v1.txt",sep="\t",row.names=F,quote=F,na="")

# subset just IPR and GO data
slim_annot2 <- slim_annot[IPR_ID!=""|GO_ID!="",]
slim_annot2 <- slim_annot[,c(1,5,6,7,8)]
slim_annot2
write.table(slim_annot2,"slim_annotations_v2.txt",sep="\t",row.names=F,quote=F,na="")

### reshape data ###

slim_annotations <- slim_annot
slim_annotations$ANALYSIS <- as.factor(slim_annotations$ANALYSIS)

# convert to list of datatables for each analysis
df_annotations <- lapply(levels(slim_annotations$ANALYSIS),function(l) slim_annotations[ANALYSIS==l])

# add names
names(df_annotations) <- levels(slim_annotations$ANALYSIS)

# drop duplicates from all tables
df_annotations <- lapply(df_annotations,function(l) l[!duplicated(l)])

# split "|" seperated strings 
long_annotations <- lapply(df_annotations,function(l) l[, strsplit(as.character(GO_ID), "|", fixed=TRUE), by = .(GENE,MATCH_ID,MATCH_DESC,IPR_ID,IPR_DESC,PATHWAY,GO_ID)][,.(GENE,MATCH_ID,MATCH_DESC,IPR_ID,IPR_DESC,PATHWAY,GO_ID = V1)])


# extract IPR info
IPR <- do.call(rbind, lapply(long_annotations, subset, select=c("GENE", "IPR_ID","IPR_DESC")))
IPR <- IPR[!duplicated(IPR)]

# extract GO info
GO <- do.call(rbind, lapply(long_annotations, subset, select=c("GENE", "GO_ID")))
GO <- GO[!duplicated(GO)]

# extract PATHWAY info
PATHWAY <- do.call(rbind, lapply(long_annotations, subset, select=c("GENE", "PATHWAY")))
PATHWAY <- PATHWAY[!duplicated(PATHWAY)]

# drop none specific analysis stuff from each table
dt_annotations <-   lapply(df_annotations, subset, select=c("GENE", "MATCH_ID","MATCH_DESC"))
# and rename MATCH to name of analysis
lapply(seq_along(dt_annotations),function(i) {
	new_name <- toupper(names(dt_annotations[i]))
	colnames(dt_annotations[[i]])[2] <<- sub("MATCH",new_name,colnames(dt_annotations[[i]])[2])
	colnames(dt_annotations[[i]])[3] <<- sub("MATCH",new_name,colnames(dt_annotations[[i]])[3])
})	

# add go and ipr details to list of datatables
dt_annotations$IPR <- IPR
dt_annotations$GO <- GO
dt_annotations$PATHWAY <- PATHWAY

# collapse genes with more than one annotation for an anlysis
dt_annotations <- lapply(dt_annotations,function(dt) dt[,lapply(.SD, paste,collapse="|"),by = GENE])
# remove "||" errors produced by collapse (not certain why it does this)
dt_annotations <- lapply(dt_annotations,function(dt) 
	as.data.table(lapply(dt, function(x) gsub("^\\||\\|$","",gsub("\\|\\|+", "\\|", x))))
)

# merge list of datatable back together
merged_annotations <- Reduce(function(...) merge(..., all = TRUE), dt_annotations)
# and set NA to ""
merged_annotations[is.na(merged_annotations)] <- ""

annotations <- merged_annotations
write.table(annotations,"ip_annotations.txt",sep="\t",na="",quote=F,row.names=F)
```

