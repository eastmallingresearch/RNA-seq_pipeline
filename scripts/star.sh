#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

read -r -d '' HELP << EOM
#################################################################################
#										#
#	Wrapper script for running STAR on cluster				#
#										#
#	usage: star.sh ref_dir prefix [forward reverse|SE_reads] (options) 	#
#										#
#	Some useful options are:						#
#	--runThreadN <number of threads>					#
#	--readFilesCommand [zcat|pigz -d|bunzip2 -c]				#
#	--outReadsUnmapped Fastx						#
#	--outSAMtype [BAM|SAM|None] [Unsorted|SortedByCoordinate]		#
#										#
#	Use -h for further options						#
#										#
#################################################################################
EOM

function print_help {
	echo;echo "$HELP" >&1;echo;
	STAR 
	exit 1
}

if [ $# -eq 2 ];
then
	echo;echo "$HELP" >&1;echo;
	exit 1   
fi

OPTIND=1

while getopts ":hs:" options; do
	case "$options" in
	s)
	  SCRIPT_DIR=$OPTARG
	  ;;		  
	h)  
	  print_help
	  exit 0
	  ;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [[ -z "$SCRIPT_DIR" ]]; then 
	SCRIPT_DIR=$(readlink -f ${0%/*})
fi

qsub $SCRIPT_DIR/submit_star.sh $@
