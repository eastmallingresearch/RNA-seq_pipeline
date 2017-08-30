#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=2G

read -r -d '' HELP << EOM
#################################################################################
#										#
#	Wrapper script for PHIX filtering					#
#										#
#	usage: filter.sh -p <program> [options]					#
#										#
#	-p ignored, bowtie is only option					#
#										#
#	filter.sh Forward Reverse Ref Output I X [options] <SE reads> 		#
#	 									#
#	SE reads should follow the following format:				#
#	-U SE1,SE2,etc. --un-gz Output -S /dev/null [options]  			#							#
#	 									#
#################################################################################
EOM

function print_help {
	echo;echo "$HELP" >&1;echo;
	exit 0
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


# note bowtie will use the same file name for all output files -
# will need to write to seperate directories if submitting multiple pairs
qsub $SCRIPT_DIR/submit_bowtie.sh $SCRIPT_DIR $@
