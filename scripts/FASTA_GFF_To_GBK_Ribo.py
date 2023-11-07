#!/usr/bin/python3
############################# INTRODUCTION #############################
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


############################# MAIN SCRIPT ############################# 

gene_dict = {}

# The Seq regions from different Subunits   
# "/home/genomics/gleberre/01_Research_BAR_ZAND/02_ZAND/01_ZAND_GS/07_ITSx_Anno_Results/MB_GS_M_johnstoni_S34/MB_GS_M_johnstoni_S34_Concat_All_Regions_ITSx.fasta"
for record in SeqIO.parse(args.fasta_itsx, "fasta"): 
    
    descr_line = record.description.split()
    #print(descr_line)
    
    if re.search('ITS1', record.description): its_one_pos = descr_line[-3].split('-')
    if re.search('ITS2', record.description): its_two_pos = descr_line[-3].split('-')
    
    # Creating a Nested dictionary
    gene_dict[descr_line[2]] = {}
    # Assign the Seq, Start and End position of each gene
    gene_dict[descr_line[2]]["StartP"] = 'StartP'
    gene_dict[descr_line[2]]["EndP"] = 'EndP'
    gene_dict[descr_line[2]]["Seq"] = record.seq

# Check wether the variable for its_two_pos was set or not >> If not use the end position of ITS1 as start for ITS2 
try:
    its_two_pos
except NameError:
    its_two_pos = [its_one_pos[1], its_one_pos[1]]
else:
    pass

# Check wether the variable for its_two_pos was set or not >> If was not found use the other ITS2 as end for ITS1
try:
    its_one_pos
except NameError:
    its_one_pos = [its_two_pos[0], its_two_pos[0]]
else:
    pass

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
    
# Read in the GFF file for the position
# "/home/genomics/gleberre/01_Research_BAR_ZAND/02_ZAND/01_ZAND_GS/08_Barrnap_Anno_Results/MB_GS_M_johnstoni_S34/MB_GS_M_johnstoni_S34_rDNA.gff"
with open(args.gff, "r") as gff_file_to_read:
    gff_lines = gff_file_to_read.readlines()
    
    for feat in gff_lines:
        
        splitted_feat = feat.strip().split("\t")
        
        if re.search('circular', splitted_feat[0]):
                topology = 'circular'
        else:
                topology = 'linear'
        
        if re.search("SSU", feat): 
            gene_dict["SSU"]["StartP"] = int(splitted_feat[3]) - 1
            gene_dict["SSU"]["Seq"] = gene_dict["SSU"]["Seq"][int(splitted_feat[3]) - 1:]
            gene_dict["SSU"]["Strand"] = determ_strand(splitted_feat[6])
        
        if re.search("LSU", feat): 
            gene_dict["LSU"]["EndP"] = int(its_two_pos[1]) + 1 + int(splitted_feat[4])
            gene_dict["LSU"]["Seq"] = gene_dict["LSU"]["Seq"][:int(splitted_feat[4])]
            gene_dict["LSU"]["Strand"] = determ_strand(splitted_feat[6])

        if re.search("5.8S", feat): 
            gene_dict["5.8S"]["Strand"] = determ_strand(splitted_feat[6])
            
    
    if 'ITS2' in gene_dict: gene_dict["ITS2"]["Strand"] = gene_dict["LSU"]["Strand"]
    if 'ITS1' in gene_dict: gene_dict["ITS1"]["Strand"] = gene_dict["SSU"]["Strand"]
       
# CHECK MULTI FASTA OR SINGLE FASTA FILE
seq_id = [record.id for record in SeqIO.parse(args.fasta_go, "fasta")]
# CONSTANT
index = 0

# Iterate over the records in the FASTA file
# Directly calling the SeqIO.parse in the iterator to make it work
for record in SeqIO.parse(args.fasta_go, "fasta"):      
        
        # Fasta file has only One sequence
        if len(seq_id) == 1: 
                
                # Need to create a list for the Source features to read in 
                Features = []

                Features, Sequence_str = source_feat(Features, record.seq)
                
                # To check if the strand is pos or negative
                for gene in gene_dict.keys():
                    
                    Features.append(SeqFeature(FeatureLocation(gene_dict[gene]["StartP"], gene_dict[gene]["EndP"], strand = gene_dict[gene]["Strand"]), 
                            type = "rRNA", 
                            qualifiers= {'gene': gene, 
                                            'sequence': gene_dict[gene]["Seq"]
                                            }))
                        
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
        
        # It is a Multi FASTA file
        else:
                # Need to create a list of features
                Features = []

                Features, Sequence_str = source_feat(Features, record.seq)
                        
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
                        'topology': topology, 
                        'date': date_create()}  

                # The Introduction to the GenBank file
                record = SeqRecord(Seq(Sequence_str), id = locus, name = locus, \
                        description = bioproject, annotations = header, features = Features)  
                # Since provided output name has gb behind >> splittext
                splitted_fileext = os.path.splitext(args.genbank)
                
                # '/home/genomics/gleberre/01_Research_BAR_ZAND/03_Pipeline_development/scripts/test.gb'  
                with open(splitted_fileext[0] + f"_{index}" + splitted_fileext[1], 'w') as handle:  
                        SeqIO.write(record, handle, 'genbank')
        
        # Increment the indx if there would be multiple GFF files present >> Will be in subfolders created by MITOS
        index += 1
