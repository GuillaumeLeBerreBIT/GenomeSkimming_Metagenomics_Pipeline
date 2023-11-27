#!/usr/bin/python3
############################# INTRODUCTION #############################
#
# Script to add the Taxonomy ID at the start of the sequence header so Kraken can recognize the Taxonomy data
#
#################################### MODULES ####################################
import argparse, re

#################################### COMMAND LINE INPUT ####################################
parser = argparse.ArgumentParser(description='Add the Taxonomy to the sequence headers')                                                         
parser.add_argument('inputFile', type=str, 
                    help='Provide the input file to add the taxid ID to the sequence headers')
parser.add_argument('outputFile', type=str, 
                    help='Provide the (path/)name for the output file')
parser.add_argument('-t', '--taxid', type=str, required = True, 
                    help ='The taxonomy ID for the reads provided.')
args = parser.parse_args()
#################################### FILE HANDLING ####################################

with open(args.outputFile, "w") as file_to_write:
    with open(args.inputFile, "r") as file_to_read:
        file_lines = file_to_read.readlines()

        for line in file_lines:
            if re.search(">", line):
                stripped_line = line.strip()
                seq_taxid = ">kraken:taxid|" + args.taxid + "|" + stripped_line[1:] + "\n"
                file_to_write.writelines(seq_taxid)
            else:
                file_to_write.writelines(line)
    file_to_read.close()
file_to_write.close()