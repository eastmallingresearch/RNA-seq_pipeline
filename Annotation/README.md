# Function annotation for transcriptomes

## Interproscan

### Installation instructions
Interproscan is a bit of a pain to get working - install instructions to follow

Gene3D is still not working correctly and Panther has been depricated
### Running on cluster



 -appl CDD,COILS,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,ProDom,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM
 -iprlookup
 -goterms
 -pa
 -dra

# annotation with interproscan
~/pipelines/common/scripts/interproscan.sh \
 OUTDIR \
 ribes.protein.fasta \
  /data/scratch/deakig/blackcurrant/genome/temp/

# restarting interproscan
 ~/pipelines/common/scripts/restart_interproscan.sh \
  /data/scratch/deakig/blackcurrant/genome \
  ribes.protein.fasta \
  /data/scratch/deakig/blackcurrant/genome/temp /
  
