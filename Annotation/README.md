# Function annotation for transcriptomes

## Interproscan

### Installation instructions
Interproscan is a bit of a pain to get working - install instructions to follow

Gene3D is still not working correctly and Panther has been depricated
### Running on cluster
Interproscan has a built in option to run on an SGE cluster - but I couldn't get it to work.  
Therefore I've written my own version which splits a transcriptome into multiple files and runs multiple simultaneous jobs (hardcoded to 10). Compared to the Interproscan inbuilt clustering this will take an additional core for each simultaneous job.  

To run interproscan on the cluster use interproscan.sh. This takes two compulsory + any iterproscan arguments.
```shell
interproscan.sh FOLDER FILE <options>
```
FOLDER is the location of the transcriptome file(s) - all output files will be copied here.  
FILE is the name of the transcriptome file (accepts shell expandable wild cards to specify multiple files).  
At least one of the options <i>must<\i> specify which analyses to perform (-appl) - otherwise the Gene3D analysis will cause the JVM to crasdh.Below is the list of possible analyses:  
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
