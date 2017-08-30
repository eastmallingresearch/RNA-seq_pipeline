#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

INFILE=$1
shift
NUMREADS=$1
shift
OUTDIR=$1
shift

LINES="$(($NUMREADS * 4))"

FNAME=$(echo $INFILE|awk -F"/" '{print $NF}')

mkdir -p $OUTDIR

TEMPFILE=${FNAME}_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

mkfifo $TEMPFILE
zcat -f -- $INFILE > $TEMPFILE &

split -l $LINES -a 5 --filter='pigz > $FILE.gz' $TEMPFILE $OUTDIR/$FNAME.

rm $TEMPFILE

 

