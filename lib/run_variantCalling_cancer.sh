#!/bin/bash
## Author S. Monzon
## version v2.0

# Help
# usage: run_variantCalling.sh ....
#

##Functions
# function calling_sge  {
#  	CALLING=$( qsub $2 -N $JOBNAME.$3 $1)
#     jobid_calling=$( echo $CALLING | cut -d ' ' -f3 | cut -d '.' -f1 )
#     echo -e "$4:$jobid_calling\n" >> $OUTPUT_DIR/logs/jobids.txt
# }
# 
# function calling_sge_array  {
# 	CALLING=$( qsub $2 -t 1-$3 -N $JOBNAME.$4 $1)
#     jobid_calling=$( echo $CALLING | cut -d ' ' -f3 | cut -d '.' -f1 )
#     echo -e "$5:$jobid_calling\n" >> $OUTPUT_DIR/logs/jobids.txt
# }
# 
# function calling {
#  	for count in `seq 1 $2`
#     do
#      	CALLING=$($1 $count)
#     done
# }

set -e
set -u
set -x

# Variables
USE_SGE=$1
VARIANT_CALLING=$2
VARIANT_CALLER=$3
DUPLICATES=$4
OUTPUT_DIR=$5
REF_PATH=$6
THREADS=$7
SAMPLE_NAMES=$8
OUTPUT_BAM_NAMES=$9
EXOME_ENRICHMENT=${10}
CONTROL_NAMES=${11}
CASE_NAMES=${12}
VCF_NAMES=${13}
STRELKA_CONFIG=${14}
sample_number=${15}

## Folder creation
echo -e "Creating $OUTPUT_DIR/variant_calling"
mkdir -p $OUTPUT_DIR/variant_calling

if [ $DUPLICATES == "YES" ]; then
	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
	CALLING_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid"
else
	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )
	CALLING_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid"
fi

if [ $VARIANT_CALLING == "YES" ];then
	if [ $VARIANT_CALLER == "STRELKA" ];then
		mkdir -p $OUTPUT_DIR/variant_calling/variants_strelka
		CALLING_CMD="$SCRIPTS_DIR/strelka.sh $OUTPUT_DIR/Alignment/BAM $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_strelka $SAMPLE_NAMES $OUTPUT_BAM_NAMES $CONTROL_NAMES $CASE_NAMES $VCF_NAMES $STRELKA_CONFIG"
		if [ "$USE_SGE" = 1 ];then
            CALLING=$( qsub $CALLING_ARGS -t 1-$sample_number -N $JOBNAME $CALLING_CMD )                      
           	jobid_calling=$( echo $CALLING | cut -d ' ' -f3 | cut -d '.' -f1 ) 
           	echo -e "VARIANT_CALLING:$jobid_calling\n" >> $OUTPUT_DIR/logs/jobids.txt       
           	#calling_sge_array "$CALLING_CMD" "$CALLING_ARGS" $sample_number
		#else
           #calling "$CALLING_CMD" $sample_number
		fi

fi


# if [ $VARIANT_CALLER == "MUTECT" ];then
# 	mkdir -p $OUTPUT_DIR/variant_calling/variants_mutect
# 	GATK_PREPROC_CMD="$SCRIPTS_DIR/gatk_preprocessing.sh $OUTPUT_DIR/Alignment/BAM $THREADS $REF_PATH $SAMPLE_NAMES $KNOWN_SNPS $KNOWN_INDELS $OUTPUT_BAM_NAMES $OUTPUT_DIR/variant_calling/variants_mutect $OUTPUT_BAM_REALIGNED_NAMES $OUTPUT_BAM_RECALIBRATED_NAMES $GATK_PATH"
# 	PON_CMD="$SCRIPTS_DIR/mutect_panel.sh $OUTPUT_DIR/variant_calling/variants_mutect/recalibration $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_mutect $OUTPUT_BAM_NAMES $CONTROL_NAMES $VCF_NAMES $KNOWN_SNPS $COSMIC"
# 	PON_COMB_CMD="$SCRIPTS_DIR/mutect_panel_comb.sh $OUTPUT_DIR/variant_calling/variants_mutect/panel_control $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_mutect/panel_control $OUTPUT_PON_NAME"
# 	CALLING_CMD="$SCRIPTS_DIR/mutect.sh $OUTPUT_DIR/Alignment/BAM $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_mutec $SAMPLE_NAMES $OUTPUT_BAM_NAMES $CONTROL_NAMES $CASE_NAMES $VCF_NAMES"
# 	if [ "$USE_SGE" = 1 ];then
# 	   calling_sge_array "$GATK_PREPROC_CMD" "$CALLING_ARGS" $sample_number $SUBJOBNAME "GATK_PREPROC"
#        jobid_gatkpre=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "GATK_PREPROC" | cut -d ':' -f2 )
#        PON_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid_gatkpre" 
#        calling_sge_array "$PON_CMD" "$PON_ARGS" $sample_number $SUBJOBNAME "PON_CREATION"
#        jobid_pon=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "PON_CREATION" | cut -d ':' -f2 ) 
#        PON_COMB_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid_pon"                           
#        calling_sge "$PON_COMB_CMD" "$PON_COMB_ARGS" $SUBJOBNAME "PON_COMB"
#        jobid_poncomb=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "PON_COMB" | cut -d ':' -f2 ) 
#        MUTECT_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid_poncomb"                           
#        calling_sge_array "$CALLING_CMD" "$MUTECT_ARGS" $sample_number $SUBJOBNAME "MUTECT_CALLING"
#        jobid_mutect=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "MUTECT_CALLING" | cut -d ':' -f2 ) 
#     else
#        calling "$PON_CMD" $sample_number
#        CALLING=${$PON_COMB_CMD $sample_number}
#        calling "$CALLING_CMD" $sample_number
#    fi
fi
