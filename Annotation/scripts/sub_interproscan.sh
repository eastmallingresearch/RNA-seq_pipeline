#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=1G
#$ -pe smp 4

SCRIPT_DIR=$1; shift
LOC=$1;shift
DIR=$1;shift
SF=$1;shift

DATA=$(sed -n -e "$SGE_TASK_ID p" $DIR/${SF}_splitfiles.txt)
OUTFILE=$(awk -F"/" '{print $NF}' <<<$DATA)

cd $TMP

$SCRIPT_DIR/../interpro/interproscan-5.27-66.0/interproscan.sh -i $DATA -b ${OUTFILE}.tsv -T $TMP  $@

cp *.* $LOC/.