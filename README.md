# EDITING

This is out of date - update required to convert tp slurm and check on current best practices
Don't think it will change dramtically..


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
RNA-seq pipeline for Illumina data. This pipeline is designed to run on a ~~Sun Grid Engine cluster 
~~ Slurm cluster

Current implementation uses Salmon pseudo alignment, with the option of using STAR for full aligning, featureCounts for abundance counting and DESeq2 for DUE

## Setup pipeline
```shell
# set RNSPL variable to pipeline folder
RNSPL=~/RNA-seq_pipeline

# set permanetly for future (bash) shell sessions (be careful with this, if you have settings in ~/.profile they will no longer load)
echo export RNSPL=~/RNA-seq_pipeline >>~/.bash_profile
```

## Common tasks
Create project folder and linked to RNA-seq pipeline 
```shell
PROJECT_FOLDER=~/projects/my_project_folder
mkdir -p $PROJECT_FOLDER

ln -s $RNSPL $PROJECT_FOLDER/RNA-seq_pipeline 
```

Create some sub folders to store files and analysis results
```shell
# folder to store cluster job output
mkdir $PROJECT_FOLDER/cluster

mkdir -p $PROJECT_FOLDER/data/fastq
mkdir $PROJECT_FOLDER/data/trimmed
mkdir $PROJECT_FOLDER/data/filtered
mkdir $PROJECT_FOLDER/data/aligned
mkdir $PROJECT_FOLDER/data/counts

mkdir -p $PROJECT_FOLDER/analysis/quality
mkdir $PROJECT_FOLDER/analysis/DGE
mkdir $PROJECT_FOLDER/analysis/DEU
```

Get hold of the raw data and link to fastq folder.  
```shell
# link
ln -s /data/RNA-seq/this_project_data/*.gz $PROJECT_FOLDER/data/fastq/. 

# Decompress sym linked files (if required)
for FILE in $PROJECT_FOLDER/data/fastq/*.gz; do 
 $PROJECT_FOLDER/RNA-seq_pipeline/scripts/PIPELINE.sh -c unzip $FILE 2
done
```

## Quality check
Qualtiy checking with fastQC (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```shell
for FILE in $PROJECT_FOLDER/data/fastq/*.gz; do 
	$PROJECT_FOLDER/RNA-seq_pipeline/scripts/PIPELINE.sh -c qcheck $FILE $PROJECT_FOLDER/analysis/quality
done
```

## Updates...

Using bbtools for filtering

### Adapter removal/quality filtering and contaminant filtering
BBTools has good options for doing all of this. 

I've merged adapter removal, phix filtering and rRNA removal into a single operation using BBDuk (though it has to run multiple times, rather than a single passthrough). To modify any settings will require editing the mega_duk.sh script. Alternatively the three operations can be run seperately using bbduk (PIPELINE.sc -c bbduk). 


To use the below scripts without editing BBtools will neeed to be installed in ~/pipelines/common/bbtools

#### Adapter/phix/rRNA removal
Runs all three of the options in "Filtering full options" shown at bottom

```shell
for FR in $PROJECT_FOLDER/data/fastq/*_1.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  sbatch --mem=20000 -p medium -c 10 $PROJECT_FOLDER/RNA-seq_pipeline/scripts/mega_duk.sh \
  $PROJECT_FOLDER/RNA-seq_pipeline/common/resources/adapters/truseq.fa \
  $PROJECT_FOLDER/RNA-seq_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/RNA-seq_pipeline/common/resources/contaminants/ribokmers.fa.gz \
  $PROJECT_FOLDER/data/filtered \
  $FR \
  $RR
done 

```
bbduk command line arguments used:  
adapter removal forward; ktrim=l k=23 mink=11 hdist=1 tpe tbo t=10
adapter removal reverse; ktrim=r k=23 mink=11 hdist=1 tpe tbo t=10
phix filtering; k=31 hdist=1 t=4
rRNA filtering; k=31 t=4 

#### Human contaminant removal (BBMap)

Not certain this is relevant for RNAseq...

```shell
for FR in $PROJECT_FOLDER/data/filtered/*_1*.fq.gz; do
  RR=$(sed 's/_1/_2/' <<< $FR)
  sbatch --mem=40000 -p medium -c 20 $PROJECT_FOLDER/metagenomics_pipeline/scripts/slurm/sub_bbmap.sh \
  $PROJECT_FOLDER/metagenomics_pipeline/common/resources/contaminants/bbmap_human \
  $PROJECT_FOLDER/data/cleaned \
  $FR \
  $RR \
  minid=0.95 \
  maxindel=3 \
  bwr=0.16 \
  bw=12 \
  quickmatch \
  fast \
  minhits=2 \
  t=8
done
```

## Align to ref genome/transcriptome

DESeq (v1.11.23).
https://f1000research.com/articles/4-1521/v2

Update to come with implementation

### Quantification using Salmon
Salmon can run in two modes, psuedo mapping or alignmnet mode using a pre-aligned SAM/BAM. The BAM <i>must</i> be unsorted and aligned to the transcriptome. 

O.K. running Salmon with the filtered read files - it's fast...

#### Build index for salmon psuedo mapping
Assuming the transcriptome is available...
If you have a genome plus a gff file, I have (somewhere) a scripts which can help make a transcriptome file. gff being a pretty grotty (non) standard, the results need to be checked.

First step is to build a psuedo mapping index file
e.g. transcriptome located at $PROJECT_FOLDER/data/genome/transcriptome.fa
```shell
# Build index from transcript file
salmon index -t $PROJECT_FOLDER/data/genome/transcriptome.fasta -i $PROJECT_FOLDER/data/genome/SALMON_quasi
```

#### Mapping and quantification with salmon
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

#### Convert pseudo counts to gene counts
tximport is an R library which can convert salmon transcript "pseudo" counts (or other pseudo mapped counts) to gene counts (looks like later versions of salmon can do this natively if you provide a gff file (with -g flag). I haven't tested this yet).   

First needs a mapping file of transcript to gene.
(Any of the output quant.sf files can be used)
```shell
 awk -F"\t" '{c=$1;sub("\..*","",$1);print c,$1}' OFS="\t" quant.sf >trans2gene.txt
```

Tximport requires R v 3.3 to install via Bioconductor, otherwise can install from the binary 
Details are in DGE.R

Good for measuring gene level DE derived from transcripts (DTE) - but what about isoforms (different transcript/exon usage - DTU/DEU)?
If you have a well anotated genome with most transcripts already described - or possibly try mapping to exons?

### Quantification with STAR and featureCounts 

#### Create STAR genome index
A index must first be created
e.g. genome and gff (not essential but useful) located in $PROJECT_FOLDER/data/genome/
```shell
STAR \
--runMode genomeGenerate \
--genomeDir $PROJECT_FOLDER/genome/STAR_illumina \
--genomeFastaFiles $PROJECT_FOLDER/genome/my_genome.fa \
--sjdbGTFfile $PROJECT_FOLDER/genome/my_gff.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene Parent
```

#### Alignment with STAR
Basic star alignment parameters:
STAR --genomeDir your_out_dir --outFileNamePrefix something --readFilesIn fastq_F fastq_R --outSAMtype SAM --runThreadN 16

Alignment using both 2-pass alignment (2-pass alignment and basic alignment) and basic alignment only .     
For two pass, first pass finds extra splice junctions second pass uses these extra annotations for mapping.    
For basic only alignment the below code was modified to comment out "--sjdbFileChrStartEnd $splice_list" (and remove the preceeding line continuation \))

##### 2-pass alignment
```shell
for R1 in $PROJECT_FOLDER/filtered/*1.fq; do  
 R2=$(echo $R1|sed -e 's/\.1\./\.2\./');  
 prefix=$(echo $R1|awk -F"/" '{gsub(/\..*/,"",$NF);print $NF}');  
 $PROJECT_FOLDER/RNA-seq_pipeline/scripts/PIPELINE.sh -c star \
 $PROJECT_FOLDER/genome/STAR_illumina \
 $PROJECT_FOLDER/junctions \
 $prefix \
 $R1 \
 $R2 \
 --outSAMmode None
 done
```

##### single alignment (with chimera detection)
```shell
splice_list=$(ls $PROJECT_FOLDER/junctions/*.tab)
for R1 in $PROJECT_FOLDER/filtered/*1.fq; do  
 R2=$(echo $R1|sed -e 's/\.1\./\.2\./');  
 prefix=$(echo $R1|awk -F"/" '{gsub(/\..*/,"",$NF);print $NF}');  
 $PROJECT_FOLDER/RNA-seq_pipeline/scripts/PIPELINE.sh -c star \
 $PROJECT_FOLDER/genome/STAR_illumina \
 $PROJECT_FOLDER/aligned \
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
This can be done in R, but is slow with seqeuntial file import. Outside R each sample is counted seperately.

```shell
featureCounts -o output_file -a gff_file sam_files
```

Can be useful to produce a SAF file from an input gff (as they are not always consistent)
The below will extract exon annotations and output the ninth column stipped of ID= and anything after the first ".", then the first column (chromosome) and etc.
```
grep exon final_genes_appended.gff3|awk -F"\t" '{gsub(/ID=/,"",$NF);gsub(/\..*/,"",$NF);print $NF,$1,$4,$5,$7}' OFS="\t" > $PROJECT_FOLDER/counts/exons.SAF
```

The RNA-seq pipeline can be used to run featureCounts:
PIPELINE.sh -c counts annotaions output_dir output_file sam/bam(s) [options]

The below runs with 12 threads (-T 12), counts all multimapping reads (-M) and uses a SAF file for input (-F SAF)
```
for D in $PROJECT_FOLDER/aligned/treatment/WTCHG*; do
 OUTFILE=$(echo $D|awk -F"/" '{print $(NF)}').counts
 $PROJECT_FOLDER/RNA-seq_pipeline/scripts/PIPELINE.sh -c counts \
 $PROJECT_FOLDER/counts/exons.SAF \
 $PROJECT_FOLDER/counts \
 $OUTFILE \
 $D/star_aligmentAligned.sortedByCoord.out.bam -T 12 -M -F SAF
done
```

## Differential gene expression with DESeq
Follow script DGE.R

## Differential exon usage with STAR/featureCounts/DEXSeq
NOTE: this section requires editing

DEXSeq requires its own gtf file - they provide a python script to create it, but due to the mess of the gtf/gff "standard" it probably won't work.  

Below is an exampe of a few lines from a typical GFF file  (ignore the header - git requires it)

gff|gff|gff|gff|gff|gff|gff|gff|gff| 
---|---|---|---|---|---|---|---|---
contig_1|AUGUSTUS|gene|1|625|0.61|-|.|ID=g1;
contig_1|AUGUSTUS|mRNA|1|625|0.61|-|.|ID=g1.t1;Parent=g1
contig_1|AUGUSTUS|CDS|1|625|0.61|-|0|ID=g1.t1.CDS1;Parent=g1.t1
contig_1|AUGUSTUS|exon|1|625|.|-|.|ID=g1.t1.exon1;Parent=g1.t1;
contig_1|AUGUSTUS|start_codon|623|625|.|-|0|Parent=g1.t1;
contig_1|AUGUSTUS|gene|2887|5449|0.57|+|.|ID=g2;
contig_1|AUGUSTUS|mRNA|2887|5449|0.57|+|.|ID=g2.t1;Parent=g2
contig_1|AUGUSTUS|start_codon|2887|2889|.|+|0|Parent=g2.t1;
contig_1|AUGUSTUS|CDS|2887|2901|0.57|+|0|ID=g2.t1.CDS1;Parent=g2.t1
contig_1|AUGUSTUS|exon|2887|2901|.|+|.|ID=g2.t1.exon1;Parent=g2.t1;
contig_1|AUGUSTUS|intron|2902|2958|0.57|+|.|Parent=g2.t1;
contig_1|AUGUSTUS|CDS|2959|3208|0.57|+|0|ID=g2.t1.CDS2;Parent=g2.t1
contig_1|AUGUSTUS|exon|2959|3208|.|+|.|ID=g2.t1.exon2;Parent=g2.t1;
contig_1|AUGUSTUS|intron|3209|3260|1|+|.|Parent=g2.t1;
contig_1|AUGUSTUS|CDS|3261|5449|1|+|2|ID=g2.t1.CDS3;Parent=g2.t1
contig_1|AUGUSTUS|exon|3261|5449|.|+|.|ID=g2.t1.exon3;Parent=g2.t1;
contig_1|AUGUSTUS|stop_codon|5447|5449|.|+|0|Parent=g2.t1;

The DEXSeq script requires gene_id and transcript_id in field 9 - for exons only. In the above gff the Parent = transcript_id, but there is no gene_id. But it can be produces by stripping the .t[\d] from the end of Parent. The below script would make a DEXSeq compatible gff:

```shell
grep exon final_genes_appended_renamed.gff3|awk -F"\t" '{gsub(/;$/,"",$9);gsub(/$/,";",$9);x=$9;gsub(/.*=/,"gene_id=",x);gsub(/Parent=/,"transcript_id=",$9);sub(/\..*/,"",x);$9=$9x;print}' OFS="\t"> DEXSeq.gff
```

This can then be fed into the HTSeq script as per DEXSeq manual, or if using featureCounts dexseq_prepare_annotation2.py 
(https://github.com/vivekbhr/Subread_to_DEXSeq)

Youll also need to install the HTSeq python package - if you've got pip installed:
```shell
pip install HTSeq
```
You might need to upgrade numpy as well (numpy-1.9.2 doesn't work)
```shell
dexseq_prepare_annotation2.py -f featureCounts.gtf DEXSeq.gff DEXSeq_final.gff  
```
Then count with featureCounts:
```shell
featureCounts -f -O -p -T 16 \
-F GTF -a featureCounts.gtf \
-o counts.out Test.bam 
```

DEU.r




