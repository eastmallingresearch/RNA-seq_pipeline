#!/bin/bash

LOC=$1;shift	
SFPREFIX=$1;shift
TEMPDIR=$1;shift

SCRIPT_DIR=$(readlink -f ${0%/*})

cd $TEMPDIR

TASKS=$(wc -l ${SFPREFIX}_splitfiles.txt|awk -F" " '{print $1}')
qsub -t 1-$TASKS:1 -tc 10  $SCRIPT_DIR/sub_interproscan.sh $SCRIPT_DIR $LOC $TEMPDIR $SFPREFIX $@


