#!/usr/bin/python3
############################# INTRODUCTION #############################
# Author: Guillaume Le Berre
# GitHub: https://github.com/GuillaumeLeBerreBIT
# 
# Multiple BLAST Hits concatenated in one fiel >> Want to colelct the unique sequences in a nes file
#
############################# MODULES #############################
import os, argparse
from Bio import SeqIO

############################# COMMAND LINE #############################
parser = argparse.ArgumentParser(description='Combine fasta files to create database of unique sequences.')
parser.add_argument('InputFolder', type=str, 
                    help='Give the path to the folder to read in all files to combine.')
parser.add_argument('OutputFile', type=str, 
                    help='Path/Name of the output file in fasta format.')
args = parser.parse_args()

############################# MAIN SCRIPT #############################
# Create a set/list of unique sequence IDs
unique_records = set()
# Informative to report all sequences processed
total_list = []

# Check wheter the file has the correct file extension.
if args.OutputFile[-6:] != ".fasta":
    filename = args.OutputFile + ".fasta"
else:
    filename = args.OutputFile


# Open the file to write to
file_to_write = open(filename, "w")
# Go over each file in the directory to gather all sequences
for file in os.listdir(args.InputFolder):
    # Open the fasta file to read
    for seq_record in SeqIO.parse(args.InputFolder + file, "fasta"):
        
        # If the sequence ID is not yet in the unique list, write the sequence to a new fasta file
        if seq_record.id not in unique_records:
            
            SeqIO.write(seq_record, file_to_write, "fasta")
        
        # Add the sequence ID to a set (unique list)    
        unique_records.add(seq_record.id)
        # Adding each sequence ID to a list (even duplicates)
        total_list.append(seq_record)
        
# Close the file to write to
file_to_write.close()

print(f"The total amount of sequences: {len(total_list)}\n\
The amount of unique sequences: {len(unique_records)}") 


