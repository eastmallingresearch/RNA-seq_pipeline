#!/bin/bash

#Assemble short reads with Tophat
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

R1=$1
R2=$2
SE=$3
REF=$4
OUT=$5
I=$6
X=$7

echo "Running Tophat with the following in= REF IS '$REF' READ 1 '$R1' READ 2 ' $R2 ' DEST is '$DEST' OUTNAME IS '$OUTNAME' MIN is '$I' MAX is '$X'"

tophat -p 4 -o $OUT $REF $R1 $R2,$SE --b2-I $I --b2-X $X

