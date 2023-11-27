#!/usr/bin/python3
############################# INTRODUCTION #############################
# Author: Guillaume Le Berre
# GitHub: https://github.com/GuillaumeLeBerreBIT
#
# Using the FASTA file for annotation (ITSx) and sequence (GetOrganelle) togheter with the GFF file to create a Genbank format
#
############################# MODULES #############################
#
# Creating a Genbank file from a GFF file
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import re, argparse, os
import datetime as dt 
#
############################# COMMAND LINE INPUT #############################
#
parser = argparse.ArgumentParser(description='Generating Genbank files from FASTA & GFF file format')
parser.add_argument('-fg', '--fasta_go', type = str, required = True, 
                    help ='Input fasta file sequence to create a GenBank file.')
parser.add_argument('-fi', '--fasta_itsx', type = str, required = True, 
                    help ='Input fasta file sequence to create a GenBank file.')
parser.add_argument('-g', '--gff', type = str, required = True, 
                    help ='Input GFF (folder)file to create a Genbank file.')
parser.add_argument('-gbk', '--genbank', type = str, required = True, 
                    help ='Output name of the Genbank file')
parser.add_argument('-oi', '--organism', type = str, required = True, 
                    help ='Name of the Organism to parse to the Feature table - Source')
parser.add_argument('-l', '--locus', type = str, required = True, 
                    help ='The Locus name to view in the header')
parser.add_argument('-b', '--bioproject', type = str, required = True, 
                    help ='Name of the Biological project')
parser.add_argument('-s', '--sample_type', type = str, required = True, 
                    help ='Type of sample used for the experiment')
parser.add_argument('-t', '--taxonomy', type = str, required = True, 
                    help ='The Taxonomy name of the specie')
parser.add_argument('-oh', '--organism_header', type = str, required = True, 
                    help ='The organism name to view in the Header')
args = parser.parse_args()
# 
############################# FUNCTIONS #############################

def date_create():
        # Using the datetime module to get the current date >> Format the data string to DD-MON-YEAR >> Visualize all caps
        return dt.datetime.now().strftime("%d-%b-%Y").upper()

def determ_strand(split_s):
        # To check if the orientation of the strand is either positive or negative. 
        if split_s == "+": return +1
        else: return -1

# To create the SOURCE of each sequence in that FASTA file
def source_feat(feat_list, seq_rec):
        # Reading in the FASTA file -- > Using Start & End position can extract Sequences from FASTA file
        sequence_str = str(seq_rec)      # Convert the sequence object to a string.
        print(len(sequence_str))
        # Genbank first Feature is the SOURCE. 
        feat_source = SeqFeature(FeatureLocation(0, len(sequence_str), strand = 1), 
                                        type = 'source', 
                                        qualifiers= {'organism': args.organism, # "Clupea harengus"
                                                        'organelle': 'ribosomal', # "mitochondrion"
                                                        'mol_type': 'genomic DNA'
                                                        })
        feat_list.append(feat_source)
                
                # Return a list that has the source Feature & A variable containing the FASTA Sequence 
        return feat_list, sequence_str
    
def features(gene_dict, feat_list):
    # To check if the strand is pos or negative
    for gene in gene_dict.keys():
        
        feat_list.append(SeqFeature(FeatureLocation(gene_dict[gene]["StartP"], gene_dict[gene]["EndP"], strand = gene_dict[gene]["Strand"]), 
                type = "rRNA", 
                qualifiers= {'gene': gene, 
                                'sequence': gene_dict[gene]["Seq"]
                                }))
    return feat_list
        

############################# MAIN SCRIPT ############################# 
      
# CHECK MULTI FASTA OR SINGLE FASTA FILE
seq_id = [record.id for record in SeqIO.parse(args.fasta_go, "fasta")]
# CONSTANT
index = 0

# Iterate over the records in the FASTA file
# Directly calling the SeqIO.parse in the iterator to make it work
for record in SeqIO.parse(args.fasta_go, "fasta"):      
        
        # Fasta file has only One sequence
        if len(seq_id) == 1: 
                
            # Dictionary to store the information of the rDNA genes in
            gene_dict = {}

            # Will iterate over the FASTA file that contains all the records found by ITSX and split >> SSU, ITS1, 5.8S, ITS2, LSU
            for rec_its in SeqIO.parse(args.fasta_itsx, "fasta"): 
                # Split the header line on spaces 
                descr_line = record.description.split()
                
                # Split the START-END positions of ITS regions
                if re.search('ITS1', rec_its.description): its_one_pos = descr_line[-3].split('-')
                if re.search('ITS2', rec_its.description): its_two_pos = descr_line[-3].split('-')
                
                # Creating a Nested dictionary using the Gene name from the header line
                gene_dict[descr_line[2]] = {}
                # Assign the Seq, Start and End position of each gene
                gene_dict[descr_line[2]]["StartP"] = 'StartP'
                gene_dict[descr_line[2]]["EndP"] = 'EndP'
                gene_dict[descr_line[2]]["Seq"] = rec_its.seq

            # Best to have the same and an start positions, it will adjust the positions itself when parsing to a Genbank
            # Encountered this problem when using exact positions >> The start position of following sequences was incorrect
            # set an extra base further
            if 'SSU' in gene_dict:
                gene_dict["SSU"]["EndP"] = int(its_one_pos[0]) - 1 

            if 'ITS1' in gene_dict:
                gene_dict["ITS1"]["StartP"] = int(its_one_pos[0]) - 1
                gene_dict["ITS1"]["EndP"] = int(its_one_pos[1])

            if '5.8S' in gene_dict:
                gene_dict["5.8S"]["StartP"] = int(its_one_pos[1]) 
                gene_dict["5.8S"]["EndP"] = int(its_two_pos[0]) - 1
                
            if 'ITS2' in gene_dict:
                gene_dict["ITS2"]["StartP"] = int(its_two_pos[0]) - 1
                gene_dict["ITS2"]["EndP"] = int(its_two_pos[1]) 

            if 'LSU' in gene_dict:
                gene_dict["LSU"]["StartP"] = int(its_two_pos[1]) 
                
            # Read in the GFF file
            with open(args.gff, "r") as gff_file_to_read:
                gff_lines = gff_file_to_read.readlines()
                
                for feat in gff_lines:
                    
                    splitted_feat = feat.strip().split("\t")
                    
                    # Check if there is circular in the name (can iterate over it does not matter since either if it is circular will always have it in the name)
                    if re.search('circular', splitted_feat[0]):
                            topology = 'circular'
                    else:
                            topology = 'linear'
                    
                    # Set the start position for SSU, Splice the sequence and overwrite it, set the strand. 
                    if re.search("SSU", feat): 
                        # Have to decrement the exact postion minus one, to get the correct position in the Genbank file
                        gene_dict["SSU"]["StartP"] = int(splitted_feat[3]) - 1
                        gene_dict["SSU"]["Seq"] = gene_dict["SSU"]["Seq"][int(splitted_feat[3]) - 1:]   # Because start position 0, decrement it minus 1
                        gene_dict["SSU"]["Strand"] = determ_strand(splitted_feat[6])
                        
                    # Set the start position for LSU, Splice the sequence and overwrite it, set the strand.
                    if re.search("LSU", feat): 
                        gene_dict["LSU"]["EndP"] = int(its_two_pos[1]) + 1 + int(splitted_feat[4])  # Count the end position of ITS2 with the End pos if LSU and one for the index to get the correct end positions
                        gene_dict["LSU"]["Seq"] = gene_dict["LSU"]["Seq"][:int(splitted_feat[4])]
                        gene_dict["LSU"]["Strand"] = determ_strand(splitted_feat[6])

                    if re.search("5.8S", feat): 
                        gene_dict["5.8S"]["Strand"] = determ_strand(splitted_feat[6])
                           
                # Set the strand equal to its neighbouring gene region
                if 'ITS2' in gene_dict: gene_dict["ITS2"]["Strand"] = gene_dict["LSU"]["Strand"]
                if 'ITS1' in gene_dict: gene_dict["ITS1"]["Strand"] = gene_dict["SSU"]["Strand"]
                
                # Need to create a list for the Source features to read in 
                Features = []
                # Get the source feature
                Features, Sequence_str = source_feat(Features, record.seq)
                # Get the rDNA features 
                Features = features(gene_dict, Features)

                        
                # Parsing all information to a GenBank file
                # Header Genbank information  
                locus = args.locus  # Complete_mitochondrial_genome OR Clupea_harengus_mtDNA
                bioproject = args.bioproject # 'Genome_Skimming'  
                taxonomy = [args.taxonomy]   # 'Clupea harengus'

                # Dictionaries with information for SeqFeature qualifiers and SeqRecord annotations  
                header = {'source': args.sample_type, 
                        'organism': args.organism_header, #'Clupea harengus (Atlantic herring)', 
                        'Taxonomy': [args.taxonomy], 
                        'molecule_type': 'DNA', # 'DNA' 
                        'topology': topology, # 'circular'
                        'date': date_create()}  

                # The Introduction to the GenBank file
                record = SeqRecord(Seq(Sequence_str), id = locus, name = locus, \
                        description = bioproject, annotations = header, features = Features)  

                # '/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/scripts/test.gb'  
                with open(args.genbank, 'w') as handle:  
                        SeqIO.write(record, handle, 'genbank')
        
        # If it is a multi FASTA containing multiple contigs >> Need to process each contig in a seperate GENBANK
        else:
            # Increment it each time for new scaffold
            index += 1
            # Need to extract the information of each scaffold per loop >> scaffold_1, Scaffold_2, scaffold_3, ..., scaffold_n
            scaff = 'scaffold_' + str(index)

            # Dictionary to store the information of the rDNA genes in
            gene_dict = {}
            # All the rDNA regions found per scaffold
            rdna_found = []
            
            # Will iterate over the FASTA file that contains all the records found by ITSX and split >> SSU, ITS1, 5.8S, ITS2, LSU
            for rec_its in SeqIO.parse(args.fasta_itsx, "fasta"): 
                # Split the header line on spaces 
                descr_line = rec_its.description.split()
                
                # Since the main loop iterates over the Records can save the information for each scaffold seperatly 
                if re.search(scaff, descr_line[0]):
                    print(descr_line)
                    # Split the START-END positions of ITS regions
                    if re.search('ITS1', rec_its.description): its_one_pos = descr_line[-3].split('-')
                    if re.search('ITS2', rec_its.description): its_two_pos = descr_line[-3].split('-')
                   
                    # Annotated region on the scaffold
                    rdna_found.append(descr_line[2])
                    
                    # Creating a Nested dictionary using the Gene name from the header line
                    gene_dict[descr_line[2]] = {}
                    # Assign the Seq, Start and End position of each gene
                    gene_dict[descr_line[2]]["StartP"] = 'StartP'
                    gene_dict[descr_line[2]]["EndP"] = 'EndP'
                    gene_dict[descr_line[2]]["Seq"] = rec_its.seq
            
            # When the scaffold only has one annotated gene        
            if len(rdna_found) == 1:
                
                # Read in the GFF file
                with open(args.gff, "r") as gff_file_to_read:
                    gff_lines = gff_file_to_read.readlines()
                    # Loop over the lines from the GFF file
                    for feat in gff_lines:
                        # Remove any whitespaces
                        splitted_feat = feat.strip().split("\t")
            
                        # Only process the lines that match a certain scaffold            
                        if re.search(scaff, splitted_feat[0]):
                            
                            # Set the start position for SSU, Splice the sequence and overwrite it, set the strand. 
                            if re.search("SSU", feat): 
                                # Have to decrement the exact postion minus one, to get the correct position in the Genbank file
                                gene_dict["SSU"]["StartP"] = int(splitted_feat[3])
                                gene_dict["SSU"]["EndP"] = int(splitted_feat[4])
                                gene_dict["SSU"]["Seq"] = gene_dict["SSU"]["Seq"][int(splitted_feat[3]) - 1: int(splitted_feat[4])]   # Because start position 0, decrement it minus 1
                                gene_dict["SSU"]["Strand"] = determ_strand(splitted_feat[6])
                                
                            # Set the start position for LSU, Splice the sequence and overwrite it, set the strand.
                            if re.search("LSU", feat): 
                                gene_dict["LSU"]["StartP"] = int(splitted_feat[3])
                                gene_dict["LSU"]["EndP"] = int(splitted_feat[4])
                                gene_dict["LSU"]["Seq"] = gene_dict["LSU"]["Seq"][int(splitted_feat[3]) - 1: int(splitted_feat[4])]   # Because start position 0, decrement it minus 1
                                gene_dict["LSU"]["Strand"] = determ_strand(splitted_feat[6])

                            if re.search("5.8S", feat): 
                                gene_dict["5.8S"]["StartP"] = int(splitted_feat[3])
                                gene_dict["5.8S"]["EndP"] = int(splitted_feat[4])
                                gene_dict["5.8S"]["Seq"] = gene_dict["5.8S"]["Seq"][int(splitted_feat[3]) - 1: int(splitted_feat[4])]   # Because start position 0, decrement it minus 1
                                gene_dict["5.8S"]["Strand"] = determ_strand(splitted_feat[6])
                                   
                # Set the strand equal to its neighbouring gene region
                if 'ITS2' in gene_dict: gene_dict["ITS2"]["Strand"] = gene_dict["LSU"]["Strand"]
                if 'ITS1' in gene_dict: gene_dict["ITS1"]["Strand"] = gene_dict["SSU"]["Strand"]   
                    
                # Need to create a list of features
                Features = []
                # Get the source feature
                Features, Sequence_str = source_feat(Features, record.seq)
                # Get the rDNA features
                Features = features(gene_dict, Features)
                        
                # Parsing all information to a GenBank file
                # Header Genbank information  
                locus = args.locus  # Complete_mitochondrial_genome OR Clupea_harengus_mtDNA
                bioproject = args.bioproject # 'Genome_Skimming'  
                sample_type = args.sample_type # 'ILVO - RMG' 
                taxonomy = [args.taxonomy]   # 'Clupea harengus'
                
                # Dictionaries with information for SeqFeature qualifiers and SeqRecord annotations  
                header = {'source': sample_type, 
                        'organism': args.organism_header, #'Clupea harengus (Atlantic herring)', 
                        'Taxonomy': taxonomy, 
                        'molecule_type': 'DNA', # 'DNA' 
                        'topology': 'linear', # Multiple cotings >> should only be linear 
                        'date': date_create()}  
                
                # The Introduction to the GenBank file
                record = SeqRecord(Seq(Sequence_str), id = locus, name = locus, \
                        description = bioproject, annotations = header, features = Features)  
                # Since provided output name has gb behind >> splittext
                splitted_fileext = os.path.splitext(args.genbank)
                
                # '/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/scripts/test.gb'  
                with open(splitted_fileext[0] + '_' + scaff + splitted_fileext[1], 'w') as handle:  
                        SeqIO.write(record, handle, 'genbank')
            
            # RIBO contig containing multiple genes found
            else:
                
                # ITS1 -- > SSU end, ITS1 start end, 5.8S start
                # ITS2 -- > 5.8S end, ITS2 start end, LSU start
                
                # Check wether the variable for its_two_pos was set or not >> If not use the end position of ITS1 as start for ITS2 
                # If ITS1 was not found THEN using the Start position from ITS2 to set the SSU start and 5.8S end region. 
                try:
                    its_two_pos
                except NameError:
                    its_two_pos = ['Empty']
                else:
                    pass

                # Check wether the variable for its_two_pos was set or not >> If was not found use the ITS2 as end for ITS1
                # If ITS2 was not found THEN they key was not set thus using the end position from ITS1 as start for 5.8S as LSU start region
                try:
                    its_one_pos
                except NameError:
                    its_one_pos = ['Empty']
                else:
                    pass
                
                # If ITS one was found then >> Can Add the lcoations of the following
                if its_one_pos[0] != 'Empty': 
                    
                    if 'SSU' in gene_dict:
                        gene_dict["SSU"]["EndP"] = int(its_one_pos[0]) - 1 

                    if 'ITS1' in gene_dict:
                        gene_dict["ITS1"]["StartP"] = int(its_one_pos[0]) - 1
                        gene_dict["ITS1"]["EndP"] = int(its_one_pos[1])

                    if '5.8S' in gene_dict:
                        gene_dict["5.8S"]["StartP"] = int(its_one_pos[1]) 
                        
                # If ITS2 is found can add the locations of the follwoing rDNA
                if its_two_pos[0] != 'Empty':
                    
                    if '5.8S' in gene_dict:
                        gene_dict["5.8S"]["EndP"] = int(its_two_pos[0]) - 1
                        
                    if 'ITS2' in gene_dict:
                        gene_dict["ITS2"]["StartP"] = int(its_two_pos[0]) - 1
                        gene_dict["ITS2"]["EndP"] = int(its_two_pos[1]) 

                    if 'LSU' in gene_dict:
                        gene_dict["LSU"]["StartP"] = int(its_two_pos[1])   
                
                # Read in again th GFF file >> Go over the scaffold and its rDNA and save it
                with open(args.gff, "r") as gff_file_to_read:
                    gff_lines = gff_file_to_read.readlines()
                    
                    for feat in gff_lines:
                        
                        splitted_feat = feat.strip().split("\t")
                        
                        if re.search(scaff, splitted_feat[0]):
                            
                            # Set the start position for SSU, Splice the sequence and overwrite it, set the strand. 
                            # Set the start position for SSU, Splice the sequence and overwrite it, set the strand. 
                            if re.search("SSU", feat): 
                                # Have to decrement the exact postion minus one, to get the correct position in the Genbank file
                                gene_dict["SSU"]["StartP"] = int(splitted_feat[3]) - 1
                                gene_dict["SSU"]["Seq"] = gene_dict["SSU"]["Seq"][int(splitted_feat[3]) - 1:]   # Because start position 0, decrement it minus 1
                                gene_dict["SSU"]["Strand"] = determ_strand(splitted_feat[6])
                                
                            # Set the start position for LSU, Splice the sequence and overwrite it, set the strand.
                            if re.search("LSU", feat): 
                                gene_dict["LSU"]["EndP"] = int(its_two_pos[1]) + 1 + int(splitted_feat[4])  # Count the end position of ITS2 with the End pos if LSU and one for the index to get the correct end positions
                                gene_dict["LSU"]["Seq"] = gene_dict["LSU"]["Seq"][:int(splitted_feat[4])]
                                gene_dict["LSU"]["Strand"] = determ_strand(splitted_feat[6])

                            if re.search("5.8S", feat): 
                                gene_dict["5.8S"]["Strand"] = determ_strand(splitted_feat[6])
                
                # If th SSU not found >> Then need to set the start position to 0
                # Set the strand to equal to 5.8S
                if  "SSU" in gene_dict and \
                    gene_dict["SSU"]["StartP"] == "StartP":  
                          
                    gene_dict["SSU"]["StartP"] = 0
                    gene_dict["SSU"]["Strand"] = gene_dict["5.8S"]["Strand"] 
                # If the LSU end position is not found set it to len of sequence
                # Set the strand to equal to 5.8S
                if "LSU" in gene_dict and \
                    gene_dict["LSU"]["EndP"] == "EndP":
                        
                    gene_dict["LSU"]["EndP"] == len(record.seq)
                    gene_dict["LSU"]["Strand"] = gene_dict["5.8S"]["Strand"] 
                
                # Set the for non-complexibillity to positve
                if 'ITS2' in gene_dict: gene_dict["ITS2"]["Strand"] = +1
                if 'ITS1' in gene_dict: gene_dict["ITS1"]["Strand"] = +1
                
                # Need to create a list of features
                Features = []
                # Get the source feature
                Features, Sequence_str = source_feat(Features, record.seq)
                # Get the rDNA features
                Features = features(gene_dict, Features)
                        
                # Parsing all information to a GenBank file
                # Header Genbank information  
                locus = args.locus  # Complete_mitochondrial_genome OR Clupea_harengus_mtDNA
                bioproject = args.bioproject # 'Genome_Skimming'  
                sample_type = args.sample_type # 'ILVO - RMG' 
                taxonomy = [args.taxonomy]   # 'Clupea harengus'
                
                # Dictionaries with information for SeqFeature qualifiers and SeqRecord annotations  
                header = {'source': sample_type, 
                        'organism': args.organism_header, #'Clupea harengus (Atlantic herring)', 
                        'Taxonomy': taxonomy, 
                        'molecule_type': 'DNA', # 'DNA' 
                        'topology': 'linear', 
                        'date': date_create()}  
                
                # The Introduction to the GenBank file
                record = SeqRecord(Seq(Sequence_str), id = locus, name = locus, \
                        description = bioproject, annotations = header, features = Features)  
                # Since provided output name has gb behind >> splittext
                splitted_fileext = os.path.splitext(args.genbank)
                
                # '/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/scripts/test.gb'  
                with open(splitted_fileext[0] + f"_{scaff}" + splitted_fileext[1], 'w') as handle:  
                        SeqIO.write(record, handle, 'genbank')
                
    
