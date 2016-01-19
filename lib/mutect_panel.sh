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

 if [ $# != 9 -a "$use_sge" == "1" ]; then
   	echo "usage: ............"
   	exit
 elif [ $# != 10 -a "$use_sge" == "0" ]; then
   	echo "usage: ............"
    	exit
 fi

 #Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
 set -x
 echo `date`

 # Variables

 DIR_BAM=$1
 THREADS=$2
 REF_PATH=$3
 OUTPUT_DIR=$4
 OUTPUT_BAM_NAMES=$5
 CONTROL_NAMES=$6
 OUTPUT_VCF_NAMES=$7
 KNOWN_SNPS=$8
 COSMIC=$9

 if [ "$use_sge" = "1" ]; then
   	sample_number=$SGE_TASK_ID
 else
   	sample_number=${11}
 fi

 CONTROL_NAME=$( echo $CONTROL_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
 CONTROL_BAM=$( echo $OUTPUT_BAM_NAMES | tr ":" "\n" |grep -w "$CONTROL_NAME" )

 mkdir -p $OUTPUT_DIR/panel_control 
 mkdir -p $OUTPUT_DIR/panel_control/$CONTROL_NAME
 
 echo $CONTROL_BAM

 java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
 -T MuTect2 \
 -R $REF_PATH \
 -I:tumor $NORMAL_BAM \
 --dbsnp $KNOWN_SNPS  \
 --cosmic $COSMIC \
 --artifact_detection_mode \
 -o $OUTPUT_DIR/panel_control/$CONTROL_NAME/$OUTPUT_VCF_NAMES \
 -S LENIENT \
 -log $OUTPUT_DIR/OUTPUT_VCF_NAMES-MUTECTPANEL.log
