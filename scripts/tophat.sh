#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
R2=$2
SE=$3
REF=$4
OUT=$5
I=$6
X=$7

qsub $SCRIPT_DIR/submit_tophat.sh $R1 $R2 $SE $REF $OUT $I $X 
