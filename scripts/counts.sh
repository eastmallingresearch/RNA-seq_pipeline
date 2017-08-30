#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

read -r -d '' HELP << EOM
#################################################################################
#										#
#	Wrapper script for running featureCounts on cluster				#
#										#
#	usage: counts.sh annotations_file output_file (options) input_file	#
#										#
#	Some useful options are:						#
#	-F SAF	(specifies saf annotations file)				#
#	-T 16	(number of threads)						#
#	-f	(exon rather than gene counting)				#
#										#
#	see featureCounts -h for further options				#
#										#
#################################################################################
EOM

function print_help {
	echo;echo "$HELP" >&1;echo;
	exit 1
}

if [ $# -eq 2 ];
then
   print_help
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

qsub $SCRIPT_DIR/submit_counts.sh $@
