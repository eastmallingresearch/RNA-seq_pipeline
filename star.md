
# Quantification with STAR and featureCounts 

## Create STAR genome index
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

## Alignment with STAR
Basic star alignment parameters:
STAR --genomeDir your_out_dir --outFileNamePrefix something --readFilesIn fastq_F fastq_R --outSAMtype SAM --runThreadN 16

Alignment using both 2-pass alignment (2-pass alignment and basic alignment) and basic alignment only .     
For two pass, first pass finds extra splice junctions second pass uses these extra annotations for mapping.    
For basic only alignment the below code was modified to comment out "--sjdbFileChrStartEnd $splice_list" (and remove the preceeding line continuation \))

### 2-pass alignment
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

### single alignment (with chimera detection)
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
## Count features
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
