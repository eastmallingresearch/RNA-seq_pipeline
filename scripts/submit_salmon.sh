#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=4G
#$ -pe smp 6

REF=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=$1; shift

cd $TMP

salmon quant -i $REF -l A -1 $FORWARD -2 $REVERSE -o . -p 6 $@

mkdir -p $OUTDIR

cp -a *  $OUTDIR/
