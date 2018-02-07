# Function annotation for transcriptomes

## Interproscan

### Installation instructions
Interproscan is a bit of a pain to get working - install instructions to follow

Gene3D is still not working correctly and Panther has been depricated
### Running on cluster
Interproscan has a built in option to run on an SGE cluster - but I couldn't get it to work.  
Therefore I've written my own version which splits a transcriptome into multiple files and runs multiple simultaneous jobs (hardcoded to 10).    
Compared to the Interproscan inbuilt clustering this will take an additional core for each simultaneous job.  

To run interproscan on the cluster use interproscan.sh. This takes two compulsory + any iterproscan arguments.
```shell
interproscan.sh FOLDER FILE <options>
```
FOLDER is the location of the transcriptome file(s) - all output files will be copied here.  
FILE is the name of the transcriptome file (accepts shell expandable wild cards to specify multiple files)


 -appl CDD,COILS,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM
 -iprlookup
 -goterms
 -pa
 -dra


```shell
# annotation with interproscan
~/pipelines/common/scripts/interproscan.sh \
 FOLDER \
 FILE \
 options
 
# restarting interproscan
 ~/pipelines/common/scripts/restart_interproscan.sh \
  /data/scratch/deakig/blackcurrant/genome \
  ribes.protein.fasta \
  /data/scratch/deakig/blackcurrant/genome/temp /
```  
