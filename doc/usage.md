# USAGE

## TRIO SERVICE WRAPPER

Due to the popularity and high demand of TRIO services, we have created a wrapper to automatise most of the analysis, [run_TRIO_service.sh](https://github.com/BU-ISCIII/exome_pipeline/blob/develop/run_TRIO_service.sh). This wrapper needs some extra dependencies apart from the basic scripts inside [lib](https://github.com/BU-ISCIII/exome_pipeline/tree/develop/lib), all of them included in [templates_TRIO_service](https://github.com/BU-ISCIII/exome_pipeline/tree/develop/templates_TRIO_service).

Despite we tried to automatise everything, there are still some points of the pipeline where user supervision is required. For this reason the wrapper is designed to stop at those points and remember the user which things may need to be adapted and which commands should be manually run. Once all the indications of the script have been checked and the user considers it is time to continue executing the pipeline, the wrapper can be invoqued again using the same command as in the previous step and it will auto-detect what step needs to be executed.

The wrapper may be also forced to run a specific step of the pipeline by indicating it as an optional parameter. The only time when specifying the step is requiered is for the last step ("finish"), as it removes intermediate files and prepares the service to be archivated and the results to submitted. For this reason, user manual execution is required.

### SYNTAX

The script needs a mandatory argument, the absolute path where the TRIO service root folder is/will be created. Optionally, a second parameter specifing a mode of execution can be included, in which case the wrapper will be forced to execute that step instead the one it would atomatically detect.

`bash run_TRIO_service.sh <Absolute Path for new project> [mode]`

### MODES

Now we will describe the available modes and which steps of the pipeline are executed in each of them.

#### "create"

This is the first mode of execution. The wrapper will create the right folder structure for the service in the path indicated.

Once the working directories are ready, the wrapper will stop and ask the user to copy the raw fastq.gz files inside the `RAW` folder (or in a subfolder inside RAW). After the input data is succesfully copied, the wrapper can be re-executed to continue the pipeline.

#### "execute"

The script will look for the fastq.gz files in `RAW`, parse their names and prepare the nomenclature for the analysis. Symbolik links of the raw files with the new naming system will be located inside `ANALYSIS/00-reads`. 

`config_diploid.file` and `module.sh` will be generated inside `ANALYSIS`, and `familypedigri.ped` inside `DOC`.

Finally, a `command.sh` will be created in `ANALYSIS` and the wrapper will stop again.

The user will have to check that both `DOC/familypedigri.ped` and `ANALYSIS/config_diploid.file` are right and suit the samples names and number. In case your service have exactly 3 samples with both parents been healthy and the kid affected, you should be goot to go. If there are more than 3 samples you will have to add them in both files, and family group will have to be modified for every group of 3 samples in `DOC/familypedigri.ped`. In case you have less than 3 samples, just remove the samples you do not need from both files.

Once everything is ready, load the modules inside `ANALYSIS/module.sh` and execute `command.sh` and wait for the submited jobs to finish. Then you are ready to re-execute the wrapper for the next step.

#### "post-processing"

IMPORTANT: Make sure you have loaded the environment from `ANALYSIS/module.sh` before executing this step.

This mode will create the scripts for the post-processing steps, and do run of them. For a detailed documentation of each of them, feel free to dive into the code and read each of the scripts involved, it is prety straight forward.

When this step finish, it will ask the user to check and adapt the length of the reads and number of samples inside the script `ANALYSIS/trio_id/stats/bedtools/03_coverage_graphs.R`. Then, the user will have to execute in order the scripts located in `ANALYSIS/trio_id/stats/bamstats` and `ANALYSIS/$trio_id/stats/bedtools`.

#### "finish"

The pipeline is over. Check that everything looks fine and manually execute this mode. It will remove all big intermediate files, rename folders to get the service ready for submition, and create the report.
