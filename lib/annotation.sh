#!/bin/bash
## Author S. Monzon
## version v2.0


# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then
   	use_sge=0
   	NSLOTS=1
else
   	use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
#set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
#set -u

## Usage

if [ $# != 10 ]; then
   	echo "usage: ............"
   	exit
fi

#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

VCFILE=$1
VCFNOVO=$2
VCFDOUBLEHIT=$3
VCF_TABLE=$4
OUTPUT_DIR=$5
REF_PATH=$6
PEDFILE=$7
KGGSEQ_PATH=$8
GATK_PATH=$9
DIR=${10}

mkdir -p $OUTPUT_DIR

java -XX:ParallelGCThreads=$NSLOTS -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
	-T VariantsToTable \
	-R $REF_PATH \
	-V $DIR/$VCFILE \
	-F CHROM -F POS -F FILTER -F ID -F QUAL -F AC -F AF -F AN -F BaseQRankSum -F ClippingRankSum -F DB \
	-F DP -F DS -F FS -F MLEAC -F MLEAF -F MQ \
	-F MQRankSum -F PG -F PhasingInconsistent -F hiConfDeNovo -F loConfDeNovo \
	-F set -F SOR -F ReadPosRankSum -F QD -GF GT -GF GQ -GF AD -GF DP \
	-GF FT -GF HP -GF JL -GF PL -GF PP -GF PQ -GF TP \
	--allowMissingData --showFiltered \
	-o $OUTPUT_DIR/vcf_merge.vcf

java -XX:ParallelGCThreads=$NSLOTS -jar $KGGSEQ_PATH/kggseq.jar $JAVA_RAM \
	-buildver hg19 \
	--no-lib-check \
	--disable-resource-update \
	--force-gty-unphased \
	--vcf-file $DIR/$VCFILE \
	--ped-file $PEDFILE \
	--double-hit-gene-trio-filter \
	--db-gene refgene \
	--db-score dbnsfp \
	--genome-annot \
	--db-filter ESP5400,dbsnp141,1kg201305,exac \
	--rare-allele-freq 0.005 \
	--mendel-causing-predict all \
	--omim-annot \
	--out $OUTPUT_DIR/$VCFDOUBLEHIT

java -XX:ParallelGCThreads=$NSLOTS -jar $KGGSEQ_PATH/kggseq.jar $JAVA_RAM \
	-buildver hg19 \
	--no-lib-check \
	--disable-resource-update \
	--force-gty-unphased \
	--vcf-file $DIR/$VCFILE \
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
	--out $OUTPUT_DIR/$VCFNOVO
