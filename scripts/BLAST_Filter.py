#!/usr/bin/python3
############################# INTRODUCTION #############################
# Author: Guillaume Le Berre
# GitHub: https://github.com/GuillaumeLeBerreBIT
# 
# DIAMOND to Assemvly
# 1) Will loop over the folders containing the blast results in csv format, tab seperated
# 2) Gather the matches per gene
# 3) Gather the seq in how many reads found  
# 4) Filter on 20 len, eval 10^-6, 70 percent identity (depending on the input variables)
# 5) Create a final Fasta file for the assembly
#####################################################################
# MODULES
#####################################################################
import os, argparse, csv, re
from Bio import SeqIO   # pip install biopython
import matplotlib.pyplot as plt
import numpy as np

#####################################################################
# COMMAND LINE INPUT
#####################################################################
parser = argparse.ArgumentParser(description='From DIAMOND to Assembly')
parser.add_argument('BlastFile', type=str, 
                    help='Give the (path to and) name of the file containing the BLAST results.')
parser.add_argument('FileSequences', type=str, 
                    help='Give the Fasta file containing the contigs/sequences.')
parser.add_argument('outputAssembly', type=str, 
                    help='The name of the output FASTA file with all sequences.')
parser.add_argument('-i', '--identity', type=int, default = 70, required = False, 
                    help ='Give a number of the percentage identity as the minimum treshold.')
parser.add_argument('-l', '--len', type=int, default = 20, required = False, 
                    help ='Give a number of the min length of residues want to set as minimum treshold.')
# Will later convert it to a float -- > YAML sets it as an string
parser.add_argument('-e', '--evalue', type=str, default = 10e-6, required = False, 
                    help ='Give the minimum value to filter the sequences from DIAMOND.')
args = parser.parse_args()

#####################################################################
# FILE HANDLING
#####################################################################
# VARIABLES & LISTS (set before loop or will reset on each new file)
# Will create a list with all the unique headers validated by the filter parameters
unique_blast_hits = set()

# Empty dictionary for saving the headers
sequence_count = {} 
       
# Can open each csv file
with open(args.BlastFile, "r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    # Each line from the csv file be returned as list
    # 0.  qseqid      query or source (gene) sequence id <-- 
    # 1.  sseqid      subject or target (reference genome) sequence id <-- 
    # 2.  pident      percentage of identical positions <--
    # 3.  length      alignment length (sequence overlap) <--
    # 4.  mismatch    number of mismatches
    # 5.  gapopen     number of gap openings
    # 6.  qstart      start of alignment in query
    # 7.  qend        end of alignment in query
    # 8.  sstart      start of alignment in subject
    # 9.  send        end of alignment in subject
    # 10.  evalue      expect value <--
    # 11.  bitscore    bit score
    for row in csv_reader:
        """
        print(f"Sequence ID: {row[0]}\
            \nTarget Species_Gene: {row[1]}\
            \nPercentage Identical positions: {row[2]}\
            \nLength of sequence overlap: {row[3]}\
            \nE-value: {row[10]}")
        """
        #################
        # FILTER 
        #################
        # First item of each row is the header name 
        header = row[0]
        
        # The float will function in the same way as the integers and can now compare numbers
        # Convert str -- > float/int
        # args.identity == Percentage Identity
        # args.len == Min length of residues 
        # args.evalue == Max e-value a sequence can have as treshold
        if float(row[2]) >= args.identity and\
            int(row[3]) >= args.len and\
            float(row[10]) <= float(args.evalue):

            # Get all the headers that match the set filters
            # The headers that will be used to create FASTA file for assembly 
            # Add a '>' since the fasta files start with a '>' 
            unique_blast_hits.add('>' + header)
        
print(unique_blast_hits)
# Header list with all IDs            
#print(filtered_list)
# Hits per gene
#print(gene_matches)
# The count of how many times a read seen per hit 
#print(sequence_count)
# Can check now the freq of how much a header has been found in a dictionary
#print(header_freq)

#####################################################################
# WRITE TO FASTA
#####################################################################
# HANDLING THE FILE
# Opening a file to write to
with open(args.outputAssembly, "w") as file_to_write:
    # Opening a file to read to
    with open(args.FileSequences, "r") as file_to_read:
        # This will create a list based on the newlines, containing strings
        reading_lines = file_to_read.readlines()
        # Setting a counter
        counter_2 = 0
        counter_3 = 0
        counter_4 = 0
        # Setting a flag
        flag = 0 
        # Iterating over the lines in the list
        for line in reading_lines:
            
            # Check for header lines/ record IDs
            if re.search("^>", line):
                
                line_stripped = line.strip()
                # If the header is in the list set the flag to 1 == to write
                if line_stripped in unique_blast_hits:
                    counter_2 += 1
                    flag = 1
                # If the header is not in the list then set flag to 0 == not write
                elif line_stripped not in unique_blast_hits:
                    flag = 0

            # This flag will allow ot print the lines to a file
            if flag == 1:
                file_to_write.writelines(line)
                counter_3 += 1
            # When flag is 0 it will skip the lines. 
            elif flag == 0:
                counter_4 += 1
                continue
# Closing the files for good practice
    file_to_read.close()
file_to_write.close()

# To get information out the parsed reads.

"""
print(f"File handled for checking matching lines (manually): {counter_2}\
      \nLines printed {counter_3}\
      \nLines skipped {counter_4}")
"""
