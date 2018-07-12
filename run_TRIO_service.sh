#!/usr/bin/bash

## Prepare folder structure and necesary files for a TRIO service execution
#
###############################################################################
# USAGE:
# bash run_TRIO_service.sh <Absolute Path for new project>
###############################################################################

# Service to use as reference, files and folder structure will be copied from here
reference_service="/processing_Data/bioinformatics/services_and_colaborations/IIER/genetica_humana/SRVIIER083_20180611_TRIO045_BM_S"
mode=""
echo "Using $reference_service as reference service to copy basic files and folder structure";

# Check that arguments are right and set mode ("create" folder structure where data will be copied, "execute" pipeline where folder and data already exist or "post-processing" for the final steps)
if [ ! -d $1 ]; then
    echo "I need the absolute path to the service folder as argument."
    exit 1
fi
new_service=$1
if [[ $new_service == /processing_Data/bioinformatics/* ]]; then
    if [ -d $new_service ]; then
        if [ -f "$new_service/ANALYSIS/config_diploid.file" ]; then
            echo "Found $new_service/ANALYSIS/config_diploid.file"
            echo "Setting mode post-processing"
            mode="post-processing"
        else
            echo "Found $new_service"
            echo "Setting mode execute"
            mode="execute"
        fi
    else
        echo "Creating new TRIO service in $new_service"
        mkdir $new_service;
        mode="create"
    fi
else
    echo "$new_service is not an absolute path"
    exit 1
fi

# Get services IDs
trio_name=$new_service
trio_name=$( basename "$trio_name" )
trio_id=${trio_name#*_}
trio_id=${trio_id%_*_*}
trio_id=${trio_id}01
reference_name=$reference_service
reference_name=$( basename "$reference_name" )
reference_id=${reference_name#*_}
reference_id=${reference_id%_*_*}
reference_id=${reference_id}01
reference_id="20180220_TRIO04501"

# Create folder structure
if [[ $mode == "create" ]]; then
    mkdir -p "$new_service/ANALYSIS" 
    mkdir -p "$new_service/DOC" 
    mkdir -p "$new_service/RAW" 
    mkdir -p "$new_service/REFERENCES" 
    mkdir -p "$new_service/RESULTS" 
    mkdir -p "$new_service/TMP" 
    
    echo ""
    echo ""
    echo "Project folder has been created in $new_service"
    echo "Copy the RAW data into $new_service/RAW/ and re-execute this script again to run the pipeline"
    echo ""
    echo ""
fi

# Execute pipeline
if [[ $mode == "execute" ]]; then
    
    # Find input inside RAW
    files=$( ls $new_service/RAW/*.fastq.gz )
    echo "Pipeline will execute using $files as input"
    
    # Copy previous service's xgen_exome files into new service's RAW folder
    cp $reference_service/RAW*/xgen-exome-research-panel-* $new_service/RAW/
    echo "Copied xgen-exome files inside RAW";
    
    # Create sample list and symbolic links in 00-reads
    mkdir -p "$new_service/ANALYSIS/00-reads" 
    sample_root=""
    >$new_service/ANALYSIS/samples_id.txt
    for file in $files; do
        link=$( basename "$file")
        strand=""
        if [[ $link == *R1* ]]; then
            strand="R1"
        elif [[ $link == *R2* ]]; then
            strand="R2"
        fi
        link=$( echo $link | cut -f1 -d'_' )
        ln -s $file $new_service/ANALYSIS/00-reads/${link}_${strand}.fastq.gz
        sample_root=${link%?}
        if [[ $strand == "R1" ]]; then
            echo $link >> $new_service/ANALYSIS/samples_id.txt
        fi
    done
    echo "Created symbolinc links for samples inside ANALYSIS/00-reads and listed samples in ANALYSIS/samples_id.txt"
    
    # Copy config_diploid.file and module.sh
    cp $reference_service/ANALYSIS/config_diploid.file $new_service/ANALYSIS/
    cp $reference_service/ANALYSIS/module.sh $new_service/ANALYSIS/
    echo "Copied config_diploid.file and module.sh into $new_service/ANALYSIS/"
    
    # Modify config_diploid.file
    sed -i "s|$reference_service|$new_service|g" $new_service/ANALYSIS/config_diploid.file
    sed -i "s/$reference_id/$trio_id/g" $new_service/ANALYSIS/config_diploid.file
    sed -i "s/ND049/$sample_root/g" $new_service/ANALYSIS/config_diploid.file
    echo "Modified config_diploid.file"
    
    # Copy familypedigri.ped
    cp $reference_service/DOC/familypedigri.ped $new_service/DOC/
    echo "Copied familypedigri.ped into $new_service/DOC/"
    
    # Modify familypedigri.ped
    sed -i "s/ND049/$sample_root/g" $new_service/DOC/familypedigri.ped
    echo "Modified familypedigri.ped"
    
    echo ""
    echo ""
    echo "If your service contains 3 samples, load the environment from \"$new_service/ANALYSIS/module.sh\" and you are good to go."
    echo "If not, adapt the following files manually:"
    echo "$new_service/DOC/familypedigri.ped"
    echo "$new_service/ANALYSIS/config_diploid.file"
    echo ""
    echo "When ready, execute the following command:"
    echo "bash /processing_Data/bioinformatics/pipelines/exome_pipeline/lib/run_exome_germline.sh $new_service/ANALYSIS/config_diploid.file &> $new_service/ANALYSIS/$trio_id.log"
    echo ""
    echo "Wait until the submitted jobs are done and re-execute this script again to run the post-processing steps"
    echo ""
    echo ""
fi

# Post-processing
if [[ $mode == "post-processing" ]]; then
    
    echo "IMPORTANT!!!!"
    echo "Make sure you have loaded the environment from \"$new_service/ANALYSIS/module.sh\" before executing this step"
    echo ""
    
    # Copy merge_parse.R
    cp $reference_service/ANALYSIS/$reference_id/annotation/merge_parse.R $new_service/ANALYSIS/$trio_id/annotation/
    cp $reference_service/ANALYSIS/$reference_id/annotation/merge_parse.R $new_service/ANALYSIS/$trio_id/annotation/bedfilter/
    echo "Copied merfe_parse.R in $new_service/ANALYSIS/$trio_id/annotation/ and $new_service/ANALYSIS/$trio_id/annotation/bedfilter/"
    
    # Execute merge_parse.R
    cd $new_service/ANALYSIS/$trio_id/annotation/
    Rscript merge_parse.R
    cd $new_service/ANALYSIS/$trio_id/annotation/bedfilter/
    Rscript merge_parse.R
    echo "Executer merge_parse.R"
    
    # Unzip FastQC results
    cd $new_service/ANALYSIS/$trio_id/QC/fastqc
    cat ../../../samples_id.txt | xargs -I % echo "cd %;unzip \*.zip;cd .." | bash
    echo "FastQC results unxipped"
    
    # Create stats directories
    mkdir -p "$new_service/ANALYSIS/$trio_id/stats" 
    mkdir -p "$new_service/ANALYSIS/$trio_id/stats/bamstats" 
    mkdir -p "$new_service/ANALYSIS/$trio_id/stats/bedtools" 
    echo "Created stats directories"
    
    # Copy lablog scripts into stats diresctories
    cp $reference_service/ANALYSIS/$reference_id/stats/lablog $new_service/ANALYSIS/$trio_id/stats/
    cp $reference_service/ANALYSIS/$reference_id/stats/bamstats/lablog $new_service/ANALYSIS/$trio_id/stats/bamstats
    cp $reference_service/ANALYSIS/$reference_id/stats/bedtools/lablog $new_service/ANALYSIS/$trio_id/stats/bedtools
    cp $reference_service/ANALYSIS/$reference_id/stats/bedtools/coverage_graphs.R $new_service/ANALYSIS/$trio_id/stats/bedtools
    echo "Copied lablog scripts"
    
    # Execute lablog scripts
    cd $new_service/ANALYSIS/$trio_id/stats/
    bash lablog
    cd $new_service/ANALYSIS/$trio_id/stats/bamstats
    bash lablog
    cd $new_service/ANALYSIS/$trio_id/stats/bedtools
    bash lablog
    echo "Executed lablog scritps"
    
    # Modify coverage_graphs.R
    sample_root=$( $new_service/ANALYSIS/samples_id.txt )
    sample_root=${sample_root%?}
    sed -i "s/ND049/$sample_root/g" $new_service/ANALYSIS/$trio_id/stats/bedtools/coverage_graphs.R
    echo "Modified coverage_graphs.R"
    
    echo ""
    echo ""
    echo "All final scripts successfully created. Check \"len_reads\" if different from 81 and the number of samples if not 3 in $new_service/ANALYSIS/$trio_id/stats/bedtools/coverage_graphs.R "
    echo "Now manually execute in the right order the numbered scripts you will find in:"
    echo "$new_service/ANALYSIS/$trio_id/stats/bamstats"
    echo "$new_service/ANALYSIS/$trio_id/stats/bedtools"
    echo ""
    echo ""
fi
