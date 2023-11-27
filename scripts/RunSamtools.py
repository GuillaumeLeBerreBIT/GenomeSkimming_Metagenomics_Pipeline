#!/usr/bin/python3
############################# INTRODUCTION #############################
# Author: Guillaume Le Berre
# GitHub: https://github.com/GuillaumeLeBerreBIT
#
# Run the Samtools command of the multiple SAM files after BWA-MEM. Read in folder and process all the SAM files
#
############################# MODULES #############################
import os, subprocess, argparse


parser = argparse.ArgumentParser(description='Generating Genbank files from FASTA & GFF file format')
parser.add_argument('-mf', '--MappingFolder', type = str, required = True, 
                    help ='Give the directory with the mapping file present.')
parser.add_argument('-of', '--OutputFilteredSam', type = str, required = True, 
                    help ='Give the name of the outputfolder for the filtered SAM files.')
args = parser.parse_args()


# If the folder to store the Mapped Indexes is present will skip
if os.path.exists(args.OutputFilteredSam):
    pass
# If the folder is not present it will >>
else:
    # Split up the full path of the folder
    splitted_path = args.OutputFilteredSam.split('/')
    # Cleaning the list containg empty indexes
    splitted_path_cleaned = [d for d in splitted_path if d != '']
    # Start from the root of the path
    path_to_check = '/'    
    # Loop over all the directories
    for directory in splitted_path_cleaned:
        # Iterate over the paths
        path_to_check += f"{str(directory)}/"
        # If it exists do nothing
        if os.path.exists(path_to_check):
            continue
        # If not create the folder
        else: 
            os.mkdir(path_to_check)

# Iterate over the SAM files to remove the unmapped reads from the alignment
for file in os.listdir(args.MappingFolder):
    #COmbining the path + file togheter
    sam_file_path = os.path.join(args.MappingFolder, file)
    basename = os.path.splitext(file)
    # Performing the samtools command
    subprocess.run(f"samtools view -@ 8 -S -F 2308 {sam_file_path} -o {args.OutputFilteredSam}/{basename[0]}_Filtered.sam", shell = True)