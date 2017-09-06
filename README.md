# RNA-seq_pipeline

## Index 
 1. Description
 4. Setup pipeline
 7. Quality check
 10. Trim data
 13. Align to ref
 16. Count features
 20. DGE analysis
 

### Description
RNA-seq pipeline for Illumina data. This pipeline is designed to run on a Sun Grid Engine cluster. 

Current implementation uses STAR for aligning, featureCounts for abundance counting and DESeq2 for DGE

Future verions will include k-mer based alignment (Kallisto/Salmon etc.) with their automatic feature abundance estimates.
Also trasnscript quantification/estimation tools (isoEM/eXpress etc.). If avialable these methods will be analysed with DESeq2 (Note to self to check this).


### Setup pipeline
```shell
# set RNSPL variable to pipeline folder
RNSPL=~/RNA-seq_pipeline
# to set permanetly for future (bash) shell sessions (be careful with this, if you have settings in ~/.profile they will no longer load)
echo export RNSPL=~/RNA-seq_pipeline >>~/.bash_profile
```
### Quality check
### Trim data
### Filter data
### Align to ref genome/transcriptome
Current method uses STAR
### Count features
### DGE analysis

