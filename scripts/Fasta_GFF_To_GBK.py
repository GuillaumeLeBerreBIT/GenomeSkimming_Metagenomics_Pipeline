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

parser = argparse.ArgumentParser(description='Generating Genbank files from FASTA & GFF file format')
parser.add_argument('-f', '--fasta', type = str, required = True, 
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
                                                        'organelle': 'mitochondrion', # "mitochondrion"
                                                        'mol_type': 'genomic DNA'
                                                        })
        feat_list.append(feat_source)
                
                # Return a list that has the source Feature & A variable containing the FASTA Sequence 
        return feat_list, sequence_str

# Add the different gene features to the list
def gene_feat(feat_list, gff_folder, Sequence_str, index):
              
        # If it is an FASTA file with only one sequence
        if os.path.exists(gff_folder + f"/result.gff"): 
                
                file_path =  gff_folder + f"/result.gff"
                
        # If it is an Multi FASTA file        
        elif os.path.exists(gff_folder + f"/{index}/result.gff"):
                        
                file_path =  gff_folder + f"/{index}/result.gff"
                    
        topology = ''
        
        # Reading in the GFF file -- > Getting the Features
        with open(file_path, "r") as gff_file_read:
                # Create a list of each GFF line from the file
                gff_lines = gff_file_read.readlines()
                # Read in each line from the file/list seperately
                for gff in gff_lines:
                        # Split the GFF file on tabular spaces
                        splitted_gff_feat = gff.strip().split("\t")
                        
                        # Each line is in a seperate list
                        # [0] = SeqName, [1] = Source, [2] = Feature, [3] = Start, [4] = End, [5] = Score, [6] = Strand, [7] = Frame, [8] = Attribute
                        # Check if the header line either contains circular or not (can iterate over the lines since they either all or none have circular in the name)
                        if re.search('circular', splitted_gff_feat[0]):
                                topology = 'circular'
                        else:
                                topology = 'linear'
                        
                        
                        if re.search("^gene", splitted_gff_feat[2]) or\
                        re.search("tRNA", splitted_gff_feat[2]) or \
                        re.search("rRNA", splitted_gff_feat[2]) or \
                        re.search("origin_of_replication", splitted_gff_feat[2]):
                                        
                                # Extracting the start and end positions. Need to substract one base postion since index starts from 0.
                                # The End can be the same due to when splicing it will not take the last base
                                start, end = (int(splitted_gff_feat[3]) - 1), int(splitted_gff_feat[4])
                                                                
                                # Extract the gene name
                                attribute = splitted_gff_feat[8].split(';')
                                
                                # To get the key containing the gene name 
                                gene_name = [gene for gene in attribute if re.search('^gene_id', gene)]
                                
                                # Need to split the gene name on =
                                splitted_gene_name = gene_name[0].split("=")
                                
                                # Need to change the Origin_Of_Replication name 
                                if re.search("origin_of_replication", splitted_gff_feat[2]): gene_type = "OOR"
                                else: gene_type = splitted_gff_feat[2]
                                
                                # When circular and the annotation loops through the end tot the start again
                                if start > end:
                                        # To check if the strand is pos or negative
                                        strand = determ_strand(splitted_gff_feat[6])
                                        # Take the spliced sequence from both sides and paste toghter
                                        seq_fromatted = Sequence_str[start:] + Sequence_str[0:end]
                                        
                                        # Get the features from both side of the sequence
                                        FeatureLoc1 = FeatureLocation(start, len(Sequence_str), strand = strand)
                                        FeatureLoc2 = FeatureLocation(0, end, strand = strand)
                                        # Parse the joined FeatureLocations to the Genbank file to CompoundLocation using a list. THis will combine them both
                                        feat_list.append(SeqFeature(CompoundLocation([FeatureLoc1, FeatureLoc2]), 
                                                type = gene_type, 
                                                qualifiers= {'gene':splitted_gene_name[1], 
                                                                'sequence': seq_fromatted
                                                                }))
                                        
                                else:   
                                        # To check if the strand is pos or negative
                                        strand = determ_strand(splitted_gff_feat[6])
                                        
                                        feat_list.append(SeqFeature(FeatureLocation(start, end, strand = strand), 
                                                type = gene_type, 
                                                qualifiers= {'gene':splitted_gene_name[1], 
                                                                'sequence': Sequence_str[start:end]
                                                                }))            
                return feat_list, topology

############################# MAIN SCRIPT #############################        

# CHECK MULTI FASTA OR SINGLE FASTA FILE
seq_id = [record.id for record in SeqIO.parse(args.fasta, "fasta")]
# CONSTANT
index = 0

# Iterate over the records in the FASTA file
# Directly calling the SeqIO.parse in the iterator to make it work
for record in SeqIO.parse(args.fasta, "fasta"):      
        
        # Fasta file has only One sequence
        if len(seq_id) == 1: 
                
                # Need to create a list for the Source features to read in 
                Features = []

                # Get the source feature and whole FASTA sequence
                Features, Sequence_str = source_feat(Features, record.seq)
                
                # Returns a feature list and Topology
                Features, topology = gene_feat(Features, args.gff, Sequence_str, index)
                        
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
                
                # Get the source feature and whole FASTA sequence         
                Features, Sequence_str = source_feat(Features, record.seq)
                
                # Returns a feature list and Topology
                Features, topology = gene_feat(Features, args.gff, Sequence_str, index)
                        
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
