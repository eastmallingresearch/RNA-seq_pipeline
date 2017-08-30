#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=2G

ANNOTATIONS=$1
shift
OUTPATH=$1
shift
OUTFILE=$1
shift

cd $TMP

featureCounts -a $ANNOTATIONS -o $OUTFILE.feature_counts.txt $@

cp * $OUTPATH/.