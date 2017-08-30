#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -pe smp 6

REF=$1; shift
OUTDIR=$1; shift
SUFFIX=$1; shift

cd $TMP

STAR --runThreadN 8 --genomeDir $REF --outFileNamePrefix ./${SUFFIX}. --readFilesIn $@


cp $SUFFIX.* $OUTDIR/.
