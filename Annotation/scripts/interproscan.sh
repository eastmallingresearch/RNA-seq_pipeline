#!/bin/bash

LOC=$1
shift	
FILES=$1
shift

SCRIPT_DIR=$(readlink -f ${0%/*})

cd $LOC
dir=`mktemp -d -p $LOC`
for f in $FILES
do
	#JOBNAME=OTU_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)

	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}'  <$f > ${f}2
	sed -i -e '1d' ${f}2
	sed -i -e 's/ .*//' ${f}2	
	split -l 1000 ${f}2 -a 4 -d $dir/${f}2.
	cd $dir
	find $PWD -name "${f}2*" >${f}_splitfiles.txt
	TASKS=$(wc -l ${f}_splitfiles.txt|awk -F" " '{print $1}')
    qsub -t 1-$TASKS:1 -tc 10  $SCRIPT_DIR/sub_interproscan.sh $SCRIPT_DIR $LOC $dir $f $@
	cd ..
	#-N ${JOBNAME}_1 
done	


##rm -r $dir


#dir=${LOC}/tmp.LekUBxEJ1e
#cd $dir
#TASKS=$(wc -l ribes11RNA_merge_tophat_g20_0024_N6_newref_cuffmerge_renamed.fasta_splitfiles.txt|awk -F" " '{print $1}')
#TASKS=$(wc -l x_splitfiles.txt|awk -F" " '{print $1}')

#qsub -t 1-$TASKS:1 -tc 10  $SCRIPT_DIR/sub_interproscan.sh $SCRIPT_DIR $LOC $dir $FILES $@

