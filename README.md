# RNA-seq_pipeline

## Index 
 1. Description
 4. Setup pipeline
 7. Quality check
 10. Trim data
 13. Align to ref
 16. Count features
 20. DGE analysis
 

## Description
RNA-seq pipeline for Illumina data. This pipeline is designed to run on a Sun Grid Engine cluster. 

Current implementation uses STAR for aligning, featureCounts for abundance counting and DESeq2 for DGE

Future verions will include k-mer based alignment (Kallisto/Salmon etc.) with their automatic feature abundance estimates and trasnscript quantification/estimation tools (isoEM/eXpress etc.). These tools can be used to analyse differential exon usage. The output is not directly usable by DESeq2 for DGE, but there are methods available to convert estimated abundances to count data.

I might add HISAT2/TOPHAT3 to the mix - once TH3 has been released.

## Setup pipeline
```shell
# set RNSPL variable to pipeline folder
RNSPL=~/RNA-seq_pipeline
# to set permanetly for future (bash) shell sessions (be careful with this, if you have settings in ~/.profile they will no longer load)
echo export RNSPL=~/RNA-seq_pipeline >>~/.bash_profile
```
## Quality check
## Trim data
## Filter data
## Align to ref genome/transcriptome

THIS NEED EDITING - it's wrongheaded

While featureCounts is perfectly capable of counting at the feature/exon level, it has no built-in statistical model for estimating likely isoforms. However the methods which can do this don't produce normal count data - therefore not compatible with DESeq2.
Luckily (as I've just found out), there are tools available which can convert estimated counts to real counts, which can be used in later versions of DESeq (v1.11.23).
https://f1000research.com/articles/4-1521/v2

Update to come with implementation

#### quantification using Salmon
Salmon can run in two modes, psuedo mapping or alignmnet mode using a pre-aligned SAM/BAM. The BAM <i>must</i> be unsorted - aargh, and aligned to a transcriptome not genome -aargh<sup>2</sup>  

O.K. running Salmon with the filtered read files - it's fast...

Test example
```shell
# Build index from transcript file
salmon index -t ../Fven_A3-5_ncbi_final_genes_appended_renamed.cds.fasta -i SALMON_quasi

salmon quant -i ../SALMON_quasi -l A  -1 WTCHG_258645_201.1.fq -2 WTCHG_258645_201.2.fq -o quasi_out_boot \
--numBootstraps 1000 -p 16 --dumpEq --seqBias --gcBias --writeUnmappedNames --writeMappings=test.sam
```
For the above the first line contains required options (for PE reads). The second line includes some of the optional settings
http://salmon.readthedocs.io/en/latest/index.html for all the options to salmon.


tximport is an R library which can convert salmon transcript "pseuod" counts (or other pseudo mapped counts) to gene counts. 
First needs a mapping file of transcript to gene.
```shell
 awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" quant.sf >trans2gene.txt
```

 tximport requires R v 3.3 to install via Bioconductor, otherwise can install from the binary 
```R
library(tximport)
library(rjson)
library(readr)
library(DESeq2)
tx2gene<-read.table("trans2gene.txt",header=T,sep="\t")
txi.reps <- tximport("quant.sf",type="salmon",tx2gene=tx2gene,txOut=T)
txi.genes <- summarizeToGene(txi.reps,tx2gene)
dds<-DESeqDataSetFromTximport(txi.genes,data.frame(S="S1",C="H"),design=~1)
sizeFactors(dds) <- sizeFactors(estimateSizeFactors(dds)
```
Good for measuring gene level DE derived from transcripts (DTE) - but what about isoforms (different transcript/exon usage - DTU/DEU)?
... Got these the wrong way round STAR/featureCounts/DEXSeq of exons for isoform detection and quantification



#### Genome alignment with STAR 
I prefer STAR now - it's performance is not so depedent on choice of input parameters (compared to e.g. Bowtie2).  
An index must first be created
```shell
STAR \
--runMode genomeGenerate \
--genomeDir $QUORN/genome/STAR_illumina \
--genomeFastaFiles $QUORN/genome/Fven_A3-5_ncbi_WT_contigs_unmasked.fa \
--sjdbGTFfile $QUORN/genome/Fven_A3-5_ncb_final_genes_appended_renamed.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene Parent
```

Basic star alignment parameters:
STAR --genomeDir your_out_dir --outFileNamePrefix something --readFilesIn fastq_F fastq_R --outSAMtype SAM --runThreadN 16

Alignment was done using both 2-pass alignment (2-pass alignment and basic alignment) and basic alignment only (though only as I added the two-stage process later, no point in just doing the basic alignment, unless you <i>must</i> have results by tomorrow .  
For two pass, first pass finds extra splice junctions second pass uses these extra annotations for mapping.    
For basic only alignment the below code was modified to comment out "--sjdbFileChrStartEnd $splice_list" (and remove the preceeding line continuation \))

2-pass alignment
```shell
for R1 in $QUORN/filtered/*1.fq; do  
 R2=$(echo $R1|sed -e 's/\.1\./\.2\./');  
 prefix=$(echo $R1|awk -F"/" '{gsub(/\..*/,"",$NF);print $NF}');  
 $QUORN/RNA-seq_pipeline/scripts/PIPELINE.sh -c star \
 $QUORN/genome/STAR_illumina \
 $QUORN/junctions \
 $prefix \
 $R1 \
 $R2 \
 --outSAMmode None
 done
```

basic alignment (with chimera detection)
```shell
splice_list=$(ls $QUORN/junctions/*.tab)
for R1 in $QUORN/filtered/*1.fq; do  
 R2=$(echo $R1|sed -e 's/\.1\./\.2\./');  
 prefix=$(echo $R1|awk -F"/" '{gsub(/\..*/,"",$NF);print $NF}');  
 $QUORN/RNA-seq_pipeline/scripts/PIPELINE.sh -c star \
 $QUORN/genome/STAR_illumina \
 $QUORN/aligned \
 $prefix \
 $R1 \
 $R2 \
 --chimSegmentMin 20 \
 --chimOutType WithinBAM \
 --outSAMtype BAM SortedByCoordinate \
 --sjdbFileChrStartEnd $splice_list; # not used for basic only alignment
 # --outFilterMatchNminOverLread 0.3 # unused parameter - useful for mapping short alignments
 # --outFilterScoreMinOverLread 0.3 # unused parameter - useful for mapping short alignments
done
```
#### Count features
Using featureCounts. 

This can be done in R, but is slow with seqeuntial file import. Outside R each sample is counted seperately.

```shell
featureCounts -o output_file -a gff_file sam_files

```

Can be useful to produce a SAF file from an input gff (as they are not always consistent)
The below will extract exon annotations and output the ninth column stipped of ID= and anything after the first ".", then the first column (chromosome) and etc.
```
grep exon final_genes_appended.gff3|awk -F"\t" '{gsub(/ID=/,"",$NF);gsub(/\..*/,"",$NF);print $NF,$1,$4,$5,$7}' OFS="\t" > $QUORN/counts/exons.SAF
```

The RNA-seq pipeline can be used to run featureCounts:
PIPELINE.sh -c counts annotaions output_dir output_file sam/bam(s) [options]

The below runs with 12 threads (-T 12), counts all multimapping reads (-M) and uses a SAF file for input (-F SAF)
```
for D in $QUORN/aligned/treatment/WTCHG*; do
OUTFILE=$(echo $D|awk -F"/" '{print $(NF)}').counts
$QUORN/RNA-seq_pipeline/scripts/PIPELINE.sh -c counts \
$QUORN/counts/exons.SAF \
$QUORN/counts \
$OUTFILE \
$D/star_aligmentAligned.sortedByCoord.out.bam -T 12 -M -F SAF
done
```
## DGE analysis


