#!/bin/bash
#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=2G
#$ -pe smp 5

SCRIPT_DIR=$1; shift
REF=$1; shift
OUTDIR=$1; shift
FORWARD=$1; shift
REVERSE=${1:-NOTHING}; shift

cd $TMP

TEMPF=FILTERED_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n1)

bowtie2 -p 4 --no-unal -x $REF -q -U $FORWARD -S /dev/null --un $TEMPF.fq


if [ $REVERSE == "NOTHING" ]; then 
	F=$(echo $FORWARD|awk -F"/" '{print $NF}')
	mv $TEMPF.fq ${F}.f.filtered.fq
	cp ${F}.f.filtered.fq $OUTDIR/.
else
	TEMPR=FILTERED_$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n1)

	bowtie2 -p 4 --no-unal -x $REF -q -U $REVERSE -S /dev/null --un $TEMPR.fq 

	grep -h "^@" $TEMPF.fq $TEMPR.fq| \
	sed -e "s/ 2:/ 1:/"| \
	sed -e 's/^@//'| \
	sort|uniq -d > $TEMPF.fq.l2

	F=$(echo $FORWARD|awk -F"/" '{print $NF}')
	R=$(echo $REVERSE|awk -F"/" '{print $NF}')

	cat $TEMPF.fq.l2|$SCRIPT_DIR/slowx_getseqs.pl $TEMPF.fq > ${F}.f.filtered.fq
	cp ${F}.f.filtered.fq $OUTDIR/.
	
	rm $TEMPF.fq ${F}.f.filtered.fq

	sed -i -e 's/ 1:/ 2:/' $TEMPF.fq.l2
	cat $TEMPF.fq.l2|$SCRIPT_DIR/slowx_getseqs.pl $TEMPR.fq > ${R}.r.filtered.fq

	cp ${R}.r.filtered.fq $OUTDIR/.

#	sed -e 's/^D.*@//'  TEMPF.fq.l1|sort|uniq -d > TEMPF.fq.l2
#	grep "^@" $TEMPR.fq >> $TEMPF.fq.l1
#	sed -i -e "s/ 2:/ 1:/" $TEMPF.fq.l1 

#	sort $TEMPF.fq.l1|uniq -d > $TEMPF.fq.l2
#	sed -i -e 's/^@//' $TEMPF.fq.l2
#	sed -i -e 's/$/\/1/' $TEMPF.fa.l2
#	usearch9  -fastx_getseqs $TEMPF.fq -labels $TEMPF.fq.l2 -fastqout ${F}.f.filtered.fq
#	usearch9  -fastx_getseqs $TEMPR.fq -labels $TEMPF.fq.l2 -fastqout ${R}.r.filtered.fq

fi



#R1=$1
#shift
#R2=$1
#shift
#REF=$1
#shift
#OUTDIR=$1
#shift
#I=$1
#shift
#X=$1
#shift

#bowtie2 -p 8 --no-unal --un-conc $OUTDIR -I $I -X $X -x $REF -1 $R1 -2 $R2 -S /dev/null $@
#cat $OUTDIR/un-conc-mate.1| paste - - - -| sort -k1,1 -S 3G|tr '\t' '\n'|gzip > $OUTDIR/$R1.filtered.gz 
#cat $OUTDIR/un-conc-mate.2| paste - - - -| sort -k1,1 -S 3G|tr '\t' '\n'|gzip > $OUTDIR/$R2.filtered.gz 
#rm $OUTDIR/un-conc-mate.1 $OUTDIR/un-conc-mate.2

