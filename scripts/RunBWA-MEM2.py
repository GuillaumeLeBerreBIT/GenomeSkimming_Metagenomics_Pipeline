#!/usr/bin/python3
############################# INTRODUCTION #############################
# Author: Guillaume Le Berre
# GitHub: https://github.com/GuillaumeLeBerreBIT
#
# Run the BWA-MEM of the FASTQ/FASTA file against all the reference organelles which are all placed in the same folder. 
#
############################# MODULES #############################
import os, subprocess, argparse

# Need to provide the [1]FolderReferenceFiles [2]PE_1 [3]PE_2 [4]SampleName

parser = argparse.ArgumentParser(description='Generating Genbank files from FASTA & GFF file format')
parser.add_argument('-f', '--ReferenceFastaFolder', type = str, required = True, 
                    help ='Give the path to the folder containing the reference fasta files.')
parser.add_argument('-1', '--PE_1', type = str, required = True, 
                    help ='The paired-end forward file.')
parser.add_argument('-2', '--PE_2', type = str, required = True, 
                    help ='The paired-end reverse file.')
parser.add_argument('-s', '--sample', type = str, required = True, 
                    help ='Name of the sample')
parser.add_argument('-io', '--IndexOutput', type = str, required = True, 
                    help ='Name of the sample')
parser.add_argument('-mo', '--MappingOutput', type = str, required = True, 
                    help ='Name of the sample')
parser.add_argument('-t', '--Touch', type = str, required = False, 
                    help ='Name of the sample')
args = parser.parse_args()

def check_path(path_to_check): 
    # If the folder to store the Mapped Indexes is present will skip
    if os.path.exists(path_to_check):
        pass
    # If the folder is not present it will >>
    else:
        # Split up the full path of the folder
        splitted_path = path_to_check.split('/')
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

# To see whether the given folder exists or not
check_path(args.IndexOutput)     
check_path(args.MappingOutput)      

# Will read in the folder with the files that need to be indexed
for file in os.listdir(args.ReferenceFastaFolder):
    # Want the index name of the file to be without the extension 
    file_splitted = os.path.splitext(file)
    #Creating the paths required for command line output
    file_index = os.path.join(args.IndexOutput, file_splitted[0])
    file_path = os.path.join(args.ReferenceFastaFolder, file)
    
    # Using the subprocess.run waits till the command is finished upon moving to the next
    # Using Popen would require .communicate() or .wait() to make sure the command is completed
    # It is required to add the argument 'shell = True' to be able to run them
    # Firstly do index the Reference fasta
    subprocess.run(f"bwa-mem2 index -p {file_index} {file_path}", shell = True)
    # If the indexing is done start the alignment of the file. 
    subprocess.run(f"bwa-mem2 mem -t 8 {file_index} {args.PE_1} {args.PE_2} -o {args.MappingOutput}/{args.sample}_{file_splitted[0]}.sam", shell = True)

subprocess.run(f"touch {args.Touch}", shell = True)