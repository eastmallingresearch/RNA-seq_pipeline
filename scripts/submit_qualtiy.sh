#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G

FILE=$1
shift
OUTDIR=$1
shift

fastqc $FILE -o $OUTDIR $@

