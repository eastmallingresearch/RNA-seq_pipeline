# EDITING

This is out of date - update required to convert tp slurm and check on current best practices (done

Don't think it will change dramtically - nothing to add


# RNA-seq_pipeline

## Index 
 - Description
 - [Preprocesong](../master/preporcessing.md)
 - [(Pseudo)alignment](../master/preporcessing.md)
   * [Pseudoalignment with Salmon](../master/salmon.md)
   * [Alignment with STAR](../master/star.md)
 - DGE analysis
 - DEU analysis
 

## Description
RNA-seq pipeline for Illumina data. This pipeline is designed to run on a ~~Sun Grid Engine cluster~~ Slurm cluster

Current implementation uses Salmon pseudo alignment + DESeq2, with the option of using STAR for full aligning, featureCounts for abundance counting and DEXSeq for DUE


## Setup pipeline
```shell
# set RNSPL variable to pipeline folder
RNSPL=~/RNA-seq_pipeline

# set permanetly for future (bash) shell sessions (be careful with this, if you have settings in ~/.profile they will no longer load)
echo export RNSPL=~/RNA-seq_pipeline >>~/.bash_profile
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




