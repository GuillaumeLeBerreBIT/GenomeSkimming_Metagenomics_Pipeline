#!/bin/bash
# Want to get the arguments from the command line

## Help function
Help()
{
   # Display Help
   echo "Syntax: scriptTemplate [-h|d|c|s|m|e]"
   echo "options:"
   echo "h     Print this Help."
   echo "d     Give the path to the folder containing the configuration files"
   echo "c     The amount of cores to use"
   echo "s     When used then will use the singularity command"
   echo "m     Need to use the mount togheter with '-s' to save the metadata in a local folder"
   echo "e     Give the conda env with Snakemake installed"
   echo
}

## VARIABLES
# Set the variables as flase, if changed then will be recognized in the IF statement 
conda_env=false
sing=false

# Add a no colon >> Wont need any arguments 
# Add column after to append and use arguments
while getopts "hd:c:sm:e:" flag; do
    case $flag in
        h) Help
            exit ;;                     #Program killed without any errors
        d) configfolder=${OPTARG} ;;
        c) cores=${OPTARG} ;;
        s) sing=true ;;
        m) singmount=${OPTARG} ;;
        e) conda_env=${OPTARG} ;;
        ?) echo "Invalid option: -${OPTARG}."
            exit 1 ;;   # 1 indicates the programmed killed with an error
    esac
done

# If conda env name is passed activate the env containing snakemake
if [ "$conda_env" != false ]
then
    conda activate "$conda_env"
fi


# Check if singularity has been used or not
if [ "$sing" != false ]    # Since the arguments as the options are counted need 7 in total to use the singularity correctly
then
    for file in $configfolder/*; do

        snakemake --use-conda --cores $cores --configfile $file --use-singularity --singularity-args "-B $singmount"
    done
else
    for file in $configfolder/*; do

        snakemake --use-conda --cores $cores --configfile $file
    done
fi 

if [ "$conda_env" != false ]
then
    conda deactivate
fi