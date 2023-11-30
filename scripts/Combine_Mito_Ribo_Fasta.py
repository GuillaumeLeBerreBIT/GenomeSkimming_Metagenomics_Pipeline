#!/usr/bin/python3
############################# INTRODUCTION #############################
# Author: Guillaume Le Berre
# GitHub: https://github.com/GuillaumeLeBerreBIT
#
# Concatenate all the sequences together 
# 
############################# MODULES #############################
#
# Importing necassary modules
from Bio import SeqIO
import os,re, argparse

############################# COMMAND LINE INPUT #############################

parser = argparse.ArgumentParser(description='Generating Genbank files from FASTA & GFF file format')
parser.add_argument('-m', '--MitoFolder', type = str, required = True, 
                    help ='Give the path to the mitochondrion assembly folder produced by GetOrganelle.')
parser.add_argument('-r', '--rDNAFolder', type = str, required = True, 
                    help ='Give the path to the ribosomal assembly folder produced by GetOrganelle.')
parser.add_argument('-o', '--outputFasta', type = str, required = True, 
                    help ='Output name of the Fastafile')
parser.add_argument('-s', '--species', type = str, required = True, 
                    help ='Full name of the species')
args = parser.parse_args()

############################# FUNCTIONS #############################
# Want to use the SAMPLE name to locate the mapping files
# Folder locationsnalways the same -- > 04_GetOrganelle_Mito_Results & 05_GetOrganelle_Ribo_Results
def find_fasta(fasta_folder):
    # Iterate over the FASTA files in a folder 
    for file in os.listdir(fasta_folder):
        
        # Split on the filename extension
        file_ext = os.path.splitext(file)   # This will return [filename, ext]
        
        # only check for the files that have fasta as axtension
        if file_ext[1] == ".fasta":
            
            # Using the break if the wanted file is present will stop the loop and save that as filename
            if re.search('complete.graph1.1', file):
                file_to_read = file
                break
            
            # When a scaffolds file is found will look further for complete file, if nothing else found will use that file
            if re.search('scaffolds.graph1.1', file):
                file_to_read = file
                
    fasta_file_path = os.path.join(fasta_folder, file_to_read)
    return fasta_file_path


def records_to_save(fasta_file, sequence_name, saved_records):
    # Save all the records of a FASTA/ Multi FASTA file
    records_count = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    
    # If only one sequence present 
    if len(records_count) == 1:
    
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # set the sequence name to the sample name >> Have ashortened name
            record.description, record.id, record.name = sequence_name, sequence_name, sequence_name
            # Save all the records
            saved_records.append(record)
    
    else:
        # If multi FASTA iterate over the records and use enumerate to number the contigs
        for i, record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
            # Set a shortened accurate name as seqeunce header
            record.description, record.id, record.name = sequence_name + f'_Contig_{i}', sequence_name + f'_Contig_{i}', sequence_name + f'_Contig_{i}'
                
            saved_records.append(record)
    
    # Do not even need to return > Can just call the function and its list and will append to it   
    #return saved_records


############################# MAIN #############################

mito_fasta = find_fasta(args.MitoFolder)
ribo_fasta = find_fasta(args.rDNAFolder)

# Need [1]InputMito [2]InputRibo [3]OutputMapFile
# File to write the results to
with open(args.outputFasta, 'w') as file_to_write_to:
    content_to_write = []
    # MITO
    records_to_save(mito_fasta,
                    f'{args.species}_Mito',
                    content_to_write)
    
    records_to_save(ribo_fasta,
                    f'{args.species}_Ribo',
                    content_to_write)    

    # Writing the records to a fasta file
    SeqIO.write(content_to_write, file_to_write_to, 'fasta')
    
    