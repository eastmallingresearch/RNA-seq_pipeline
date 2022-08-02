# Prepocessing

This describes the steps required to preprocess RNA-seq reads using the RNA-seq pipeline

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
