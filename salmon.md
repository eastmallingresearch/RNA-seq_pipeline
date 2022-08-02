# Quantification using Salmon
Salmon can run in two modes, psuedo mapping or alignmnet mode using a pre-aligned SAM/BAM. The BAM <i>must</i> be unsorted and aligned to the transcriptome. 


## Build index for salmon psuedo mapping
Assuming the transcriptome is available...
If you have a genome plus a gff file, I have (somewhere) a scripts which can help make a transcriptome file. gff being a pretty grotty (non) standard, the results need to be checked.

First step is to build a psuedo mapping index file
e.g. transcriptome located at $PROJECT_FOLDER/data/genome/transcriptome.fa
```shell
# Build index from transcript file
salmon index -t $PROJECT_FOLDER/data/genome/transcriptome.fasta -i $PROJECT_FOLDER/data/genome/SALMON_quasi
```

## Mapping and quantification with salmon
is a single process - it accepts gz/bgz compressed files
```shell
for FR in $PROJECT_FOLDER/data/filtered/*_1.filtered.fq.gz; do
 RR=$(echo $FR|sed 's/_1/_2/')
 OUTDIR=$(echo $FR|awk -F"/" '{print $NF}'|sed 's/_.*//')
 $PROJECT_FOLDER/RNA-seq_pipeline/scripts/PIPELINE.sh -c salmon \
 $PROJECT_FOLDER/data/genome/SALMON_quasi \
 $PROJECT_FOLDER/data/counts/$OUTDIR \
 $FR $RR \
 --numBootstraps 1000 --dumpEq --seqBias --gcBias
done
 ```
--writeUnmappedNames --writeMappings=test.sam

From the above the last line includes some of the optional settings for salmon.
see http://salmon.readthedocs.io/en/latest/index.html for description of all options (or; salmon quant --help-reads
). 

### Convert pseudo counts to gene counts
tximport is an R library which can convert salmon transcript "pseudo" counts (or other pseudo mapped counts) to gene counts (looks like later versions of salmon can do this natively if you provide a gff file (with -g flag). I haven't tested this yet).   

First needs a mapping file of transcript to gene.
(Any of the output quant.sf files can be used)
```shell
 awk -F"\t" '{c=$1;sub("\\..*","",$1);print c,$1}' OFS="\t" quant.sf >trans2gene.txt
```

Tximport requires R v 3.3 to install via Bioconductor, otherwise can install from the binary 
Details are in DGE.R

Good for measuring gene level DE derived from transcripts (DTE) - but what about isoforms (different transcript/exon usage - DTU/DEU)?
If you have a well anotated genome with most transcripts already described - or possibly try mapping to exons?
