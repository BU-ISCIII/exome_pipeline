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

if [ $# != 11 ]; then
   	echo "usage: ............"
   	exit
fi

#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

VCFILE=$1
BEDFILE=$2
VCFFIL=$3
VCFNOVO=$4
VCFDOUBLEHIT=$5
OUTPUT_DIR=$6
DIR=$7
REF_PATH=$8
PEDFILE=$9
KGGSEQ_PATH=${10}
GATK_PATH=${11}

mkdir -p $OUTPUT_DIR/annotation

bedtools intersect -header -a $DIR/$VCFILE -b $BEDFILE -wa > $OUTPUT_DIR/annotation/$VCFFIL

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
	-T VariantsToTable \
	-R $REF_PATH \
	-V $OUTPUT_DIR/annotation/$VCFFIL \
	-F CHROM -F POS -F ID -F QUAL -F AC -F AF -F AN -F BaseQRankSum -F ClippingRankSum -F DB \
	-F DP -F DS -F FS -F MLEAC -F MLEAF -F MQ \
	-F MQRankSum -F PG -F PhasingInconsistent -F hiConfDeNovo -F loConfDeNovo \
	-F set -F SOR -F ReadPosRankSum -F QD -GF GT -GF GQ -GF AD -GF DP \
	-GF FT -GF HP -GF JL -GF PL -GF PP -GF PQ -GF TP \
	-o $OUTPUT_DIR/annotation/vcf_merge.vcf

java -jar $KGGSEQ_PATH/kggseq.jar -Xmx24G \
	-buildver hg19 \
	--vcf-file $OUTPUT_DIR/annotation/$VCFFIL \
	--ped-file $PEDFILE \
	--double-hit-gene-trio-filter \
	--db-gene refgene \
	--db-score dbnsfp \
	--genome-annot \
	--db-filter ESP5400,dbsnp141,1kg201305,exac \
	--rare-allele-freq 0.005 \
	--mendel-causing-predict all \
	--omim-annot \
	--out $OUTPUT_DIR/annotation/$VCFNOVO

java -jar $KGGSEQ_PATH/kggseq.jar -Xmx24G \
	-buildver hg19 \
	--vcf-file $VCFFIL \
	--ped-file $PEDFILE \
	--genotype-filter 4,7 \
	--ignore-homo \
	--db-gene refgene \
	--db-score dbnsfp \
	--genome-annot \
	--db-filter ESP5400,dbsnp141,1kg201305,exac \
	--rare-allele-freq 0.005 \
	--mendel-causing-predict all \
	--omim-annot \
	--out $OUTPUT_DIR/annotation/$VCFDOUBLEHIT
