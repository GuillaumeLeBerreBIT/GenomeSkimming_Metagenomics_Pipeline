#!/usr/bin/python3

import os, argparse, re
from Bio import SeqIO

############################# COMMAND LINE INPUT #############################

parser = argparse.ArgumentParser(description='Generating Genbank files from FASTA & GFF file format')
parser.add_argument('-r1', '--R1', type = str, required = True, 
                    help ='The R1 FASTQ file')
parser.add_argument('-r2', '--R2', type = str, required = True, 
                    help ='The R2 FASTQ file')
parser.add_argument('-fs', '--FilterSequences', type = str, required = False, 
                    help ='Give the name to filter the FASTQ sequences on ')
parser.add_argument('-c', '--Counter', action =  "store_true", 
                    help ='To count the amount of sequences in a file')
args = parser.parse_args()


def filter_gamma_delta_reads(inputFastq, outputFastq, FilterSequences):
    
    with open(outputFastq, "w") as file_to_write: 
        
        for record in SeqIO.parse(inputFastq, format = "fastq"):
            
            splitted_header = record.id.split("_", 1)
            
            if re.search(FilterSequences, splitted_header[1]):
                
                SeqIO.write(record, file_to_write, format = "fastq")

def counter(R1, R2):
    
    for file in [R1, R2]:
        count_dict = {}
        
        for record in SeqIO.parse(file, format = "fastq"):
            
            splitted_header = record.id.split("_", 1)
            
            if splitted_header[1] not in count_dict:
                
                count_dict[splitted_header[1]] = 1
                
            else: 
                count_dict[splitted_header[1]] += 1
        
        print(f"{file} -- > {count_dict}")
                
            
                
def main():
    
    if args.Counter:
        counter(args.R1, args.R2)
        
    else: 
        split_ext_1 = os.path.splitext(args.R1)
        split_ext_2 = os.path.splitext(args.R2)
        
        output_1 = split_ext_1[0] + '_' + args.FilterSequences + split_ext_1[1]
        output_2 = split_ext_2[0] + '_' + args.FilterSequences + split_ext_2[1]
        
        filter_gamma_delta_reads(args.R1, output_1, args.FilterSequences)
        filter_gamma_delta_reads(args.R2, output_2, args.FilterSequences)
    
    
main()                


# "/home/genomics/gleberre/01_Research_BAR_ZAND/02_ZAND/02_ZAND_MG/09_Gamma_Delta_Results/MB_MG_ZVL_S42/MB_MG_ZVL_S42_R1_Paired.fastq"
# "/home/genomics/gleberre/01_Research_BAR_ZAND/02_ZAND/02_ZAND_MG/09_Gamma_Delta_Results/MB_MG_ZVL_S42/MB_MG_ZVL_S42_R2_Paired.fastq"



