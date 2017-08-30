#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

#==========================================================
#	Set Help (add message between EOM blocks)
#==========================================================	
read -r -d '' HELP << EOM
#########################################################
#							#
#	RNA-seq pipeline				#
#							#
#	usage: PIPELINE.sh -c <program> [options]	#
#							#
#########################################################

 -c <program>	Program can be any of the defined programs
 -h		display this help and exit	
EOM

function print_help {
	echo;echo "$HELP" >&1;echo;
	exit 1
}

if [ $# -eq 0 ];
then
   print_help
fi

#==========================================================
#	Set command line switch options
#==========================================================

OPTIND=1 

while getopts ":hc:" options; do
	case "$options" in
	c)  
 	    program=$OPTARG
 	    break
	    ;;
	:)
	    echo "Option -$OPTARG requires an argument." >&2
	    exit 1
      	    ;;	
	h)  
	    print_help
	    exit 0
 	    ;;
	?) 
	    echo "Invalid option: -$OPTARG" >&2
	    echo "Call PIPELINE with -h switch to display usage instructions" >&2
	    exit 1
	    ;;
	esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

#==========================================================
#	Set programs (-c)
#==========================================================

case $program in
quality|quality.sh)
	$SCRIPT_DIR/quality.sh -s $SCRIPT_DIR $@
	exit 0
	;;
split|splitfq|splitfq.sh)
	$SCRIPT_DIR/splitfq.sh -s $SCRIPT_DIR $@
	exit 0
	;;
trim|trim.sh)
	$SCRIPT_DIR/trim.sh -s $SCRIPT_DIR $@
	exit 0
;;
filter|filter.sh)
	$SCRIPT_DIR/filter.sh -s $SCRIPT_DIR $@
	exit 0
;;
star|star.sh)
	$SCRIPT_DIR/star.sh -s $SCRIPT_DIR $@ 
	exit 0
;;
align|align.sh)
	$SCRIPT_DIR/align.sh -s $SCRIPT_DIR $@ 
	exit 0
;;
counts|counts.sh)
	$SCRIPT_DIR/counts.sh -s $SCRIPT_DIR $@
	exit 0
;;
anything|anything.sh)
	$SCRIPT_DIR/anything.sh -s $SCRIPT_DIR $@ 
	exit 0
;;
TEST)
	qsub $SCRIPT_DIR/counts.sh -s $SCRIPT_DIR $@ 
	exit 0
;;

*)
	echo "Invalid program: $program" >&2
	exit 1
esac
