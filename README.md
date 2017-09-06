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

Future verions will include k-mer based alignment (Kallisto/Salmon etc.) with their automatic feature abundance estimates and trasnscript quantification/estimation tools (isoEM/eXpress etc.). These tools can be used to analyse differential exon usage. The output is not directly usable by DESeq2 for DGE, but there are methods available to convert estimated abundances to count data.

I might add HISAT2/TOPHAT3 to the mix - once TH3 has been released.

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

