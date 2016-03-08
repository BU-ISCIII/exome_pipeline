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
                                                                                                                                                                                                                                     
if [ $# != 5 -a "$use_sge" == "1" ]; then                                                                                                                                                                                         
    	echo "usage: ............"                                                                                                                                                                                                      
    	exit                                                                                                                                                                                                                            
elif [ $# != 6 -a "$use_sge" == "0" ]; then                                                                                                                                                                                       
    	echo "usage: ............"                                                                                                                                                                                                      
     	exit                                                                                                                                                                                                                        
fi                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                     
#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed                     
set -x                                                                                                                                                                                                                             
echo `date`                                                                                                                                                                                                                        
                                                                                                                                                                                                                                     
# Variables                                                                                                                                                                                                                        
                                                                                                                                                                                                                                     
DIR=$1                                                                                                                                                                                                                         
THREADS=$2                                                                                                                                                                                                                         
REF_PATH=$3                                                                                                                                                                                                                        
OUTPUT_DIR=$4                                                                                                                                                                                                                      
OUTPUT_PON_NAME=$5                                                                                                                                                                                                                

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                                      
    	sample_number=$SGE_TASK_ID                                                                                                                                                                                                      
else                                                                                                                                                                                                                               
    	sample_number=${11}                                                                                                                                                                                                             
fi                                                                                                                                                                                                                                 

find $DIR -name "*.vcf" > $DIR/vcf.list

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \                                                                                                                                                       
    -T CombineVariants \
    -R $REF_PATH \
    -V $OUTPUT_DIR/vcf.list \
    -minN 2 \
    --setKey "null" \
    --filteredAreUncalled \
    --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
    -o $OUTPUT_DIR/$OUTPUT_PON_NAME \
	-S LENIENT \                                                                                                                                                                                                                       
  	-log $OUTPUT_DIR/$OUTPUT_PON_NAME-MUTECTPANEL.log                                                                                                                                                                                  
