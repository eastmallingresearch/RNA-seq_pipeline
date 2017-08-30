# RNA-seq_pipeline

## Index 
 1. Description
 4. Setup pipeline
 7. Qualit check
 10. Trim data
 13. Align to ref
 16. Count features
 20. DGE analysis
 

### Description
RNA-seq pipeline for Illumina Mi/HiSeq etc. data. This pipeline is designed to run on a Sun Grid Engine cluster. 

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
### Count Features
### DGE analysis

