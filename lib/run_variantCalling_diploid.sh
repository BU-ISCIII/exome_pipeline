#!/bin/bash
## Author S. Monzon
## version v2.0

# Help
# usage: run_variantCalling.sh ....
#

set -e
set -u
set -x

# Variables
USE_SGE=$1
VARIANT_CALLING=$2
DUPLICATES=$3
OUTPUT_DIR=$4
REF_PATH=$5
THREADS=$6
SAMPLE_NAMES=$7
OUTPUT_BAM_NAMES=$8
EXOME_ENRICHMENT=$9
VCF_NAMES=${10}
sample_number=${11}
OUTPUT_BAM_REALIGNED_NAMES=${12}
OUTPUT_BAM_RECALIBRATED_NAMES=${13}
GATK_PATH=${14}
SNPS_NAME=${15}
SNPS_FIL_NAME=${16}
INDELS_NAME=${17}
INDELS_FIL_NAME=${18}
VCF_FIL_NAME=${19}
VCF_PHASE_NAME=${20}
VCF_BACKED_NAME=${21}
VCF_GTPOS=${22}
VCF_GTPOS_FIL=${23}
VCF_GTPOS_FIL_ANNOT=${24}
KNOWN_SNPS=${25}
KNOWN_INDELS=${26}
SNP_GOLD=${27}
PED_FILE=${28}
VCFFIL=${29}
VCFNOVO=${30}
VCFDOUBLEHIT=${31}
KGGSEQ_PATH=${32}

## Folder creation
echo -e "Creating $OUTPUT_DIR/variant_calling"
mkdir -p $OUTPUT_DIR/variant_calling

if [ "$USE_SGE" = "1" -a $DUPLICATES == "YES" ]; then
  	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
  	PRECALLING_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid"
else
  	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )
  	PRECALLING_ARGS="${SGE_ARGS} -pe openmp $THREADS -hold_jid $jobid"
fi

mkdir -p $OUTPUT_DIR/variant_calling/variants_gatk
PRECALLING_CMD="$SCRIPTS_DIR/gatk_preprocessing.sh $OUTPUT_DIR/Alignment/BAM $THREADS $REF_PATH $SAMPLE_NAMES $KNOWN_SNPS $KNOWN_INDELS $OUTPUT_BAM_NAMES $OUTPUT_DIR/variant_calling/variants_gatk $OUTPUT_BAM_REALIGNED_NAMES  $OUTPUT_BAM_RECALIBRATED_NAMES $GATK_PATH"
CALLING_CMD="$SCRIPTS_DIR/gatk_diploid.sh $OUTPUT_DIR/variant_calling/variants_gatk/recalibration $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_gatk $KNOWN_SNPS $KNOWN_INDELS $SNP_GOLD $OUTPUT_BAM_RECALIBRATED_NAMES $VCF_NAMES $SNPS_NAME $SNPS_FIL_NAME $INDELS_NAME $INDELS_FIL_NAME $VCF_FIL_NAME $VCF_PHASE_NAME $VCF_BACKED_NAME $VCF_GTPOS $VCF_GTPOS_FIL $VCF_GTPOS_FIL_ANNOT $GATK_PATH $PED_FILE"
ANNOTATION_CMD="$SCRIPTS_DIR/annotation.sh $VCF_GTPOS_FIL_ANNOT $EXOME_ENRICHMENT $VCFFIL $VCFNOVO $VCFDOUBLEHIT $OUTPUT_DIR $OUTPUT_DIR/variant_calling/variants_gatk/variants $REF_PATH $PED_FILE $KGGSEQ_PATH $GATK_PATH"

if [ $VARIANT_CALLING == "YES" ]; then
 		if [ "$USE_SGE" = "1" ]; then
            PRECALLING=$( qsub $PRECALLING_ARGS -t 1-$sample_number -N $JOBNAME.CALLING $PRECALLING_CMD)
     		jobid_precalling=$( echo $PRECALLING | cut -d ' ' -f3 | cut -d '.' -f1 )
     		CALLING_ARGS="${SGE_ARGS} -l h_vmem=40g -hold_jid $jobid_precalling"
     		CALLING=$( qsub $CALLING_ARGS -N $JOBNAME.CALLING $CALLING_CMD)
        	jobid_calling=$( echo $CALLING | cut -d ' ' -f3 | cut -d '.' -f1 )
			ANNOTATION_ARGS="${SGE_ARGS} -l h_vmem=40g -hold_jid $jobid_calling"
			ANNOTATION=$( qsub $ANNOTATION_ARGS -N $JOBNAME.ANNOTATION $ANNOTATION_CMD)
        	jobid_annotation=$( echo $ANNOTATION | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "Variant Calling:$jobid_precalling - $jobid_calling \n" >> $OUTPUT_DIR/logs/jobids.txt
            echo -e "Annotation:$jobid_annotation \n" >> $OUTPUT_DIR/logs/jobids.txt
 		else
         	for count in `seq 1 $sample_number`
         	do
         		echo "Running variant calling on sample $count"
         		PRECALLING=$($PRECALLING_CMD $count)
        		done
         	CALLING=$($CALLING_CMD)
			ANNOTATION=$($ANNOTATION_CMD)
       fi
 fi
