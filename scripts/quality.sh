#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

read -r -d '' HELP << EOM
#########################################################
#								#
#	Wrapper script for checking NGS read quality 	#
#								#
#	usage: quality.sh -p <program> [options]	#
#								#
#	-p ignored, fastqc is only option		#
#								#
#	quality.sh read outdir (options)		#
#								#
#	use -h for options				#
#								#
#########################################################
EOM

function print_help {
	echo;echo "$HELP" >&1;echo;
	fastqc -h
}

if [ $# -eq 2 ];
then
	echo;echo "$HELP" >&1;echo;
	exit 0
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

$SCRIPT_DIR/submit_quality.sh $@
