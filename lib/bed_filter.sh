#!/bin/bash
## Author S. Monzon
## version v2.0


# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then
   	use_sge=0
else
   	use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u

## Usage

if [ $# != 5 ]; then
   	echo "usage: ............"
   	exit
fi

#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

VCFILE=$1
BEDFILE=$2
DIR=$3
VCFFIL=$4
OUTPUT_DIR=$5

mkdir -p $OUTPUT_DIR/annotation
mkdir -p $OUTPUT_DIR/annotation/bedfilter

bedtools intersect -header -a $DIR/$VCFILE -b $BEDFILE -wa > $OUTPUT_DIR/annotation/bedfilter/$VCFFIL
