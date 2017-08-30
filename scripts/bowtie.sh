#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
R2=$2
REF=$3
OUT=$4
I=$5
X=$6
echo "Assemble $REF with Bowtie $R1 \n $R2 \n $REF \n $NAME \n"

echo "R1 is $R1"
echo "R2 is $R2"
echo "REF is $REF"
echo "OUT is $OUT"
qsub $SCRIPT_DIR/submit_bowtie.sh $R1 $R2 $REF $OUT $I $X 