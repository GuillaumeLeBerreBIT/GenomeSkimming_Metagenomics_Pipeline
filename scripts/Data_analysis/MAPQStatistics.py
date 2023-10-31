import pandas as pd
#from plotnine import aes, geom_histogram, ggplot, scale_x_continuous
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
import os, re

#################################### EXTRACT MAPQ VALUES ####################################
def extract_MAPQ(file_path):
    # Empty list
    list_for_mapq = []
    # Open the GFF-3 file with the MAPQ values
    with open(file_path, "r") as gff_to_read:
        # Convert the file to readlines list
        gff_lines = gff_to_read.readlines()#[0:5]
        
        # Print each item from the list
        for line in gff_lines:
            # Split on tabular spaces -- > New list with sepreate values per header
            splitted_line = line.split("\t")
            # Add the MAPQ values to a new list
            list_for_mapq.append(int(splitted_line[4]))

    gff_to_read.close()
    
    return list_for_mapq    # Return a list to save the values in a variable
    
#################################### P platessa bucket on P platessa Reference Mito ####################################
P_platessa_P_platessa_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/P_platessa_P_platessa_Ref_Mapped_Only.sam")
    
#################################### P platessa bucket on Sscombrus Reference Mito ####################################
P_platessa_S_scombrus_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/P_platessa_S_scombrus_Ref_Mapped_Only.sam") 

C_harengus_bucket_P_platessa_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/C_harengus_P_platessa_Ref_Mapped_Only.sam")
C_harengus_bucket_S_scombrus_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/C_harengus_S_scombrus_Ref_Mapped_Only.sam")

P_platessa_bucket_C_harengus_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/P_platessa_C_harengus_Ref_Mapped_Only.sam")
S_scombrus_bucket_C_harengus_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/S_scombrus_C_harengus_Ref_Mapped_Only.sam")
S_scombrus_bucket_P_platessa_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/S_scombrus_P_platessa_Ref_Mapped_Only.sam")
S_scombrus_bucket_S_scombrus_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/S_scombrus_S_scombrus_Ref_Mapped_Only.sam")
C_harengus_bucket_C_harengus_Ref_Mapped_Only_Mito = extract_MAPQ("Mito/C_harengus_C_harengus_Ref_Mapped_Only.sam")
Field_1_C_harengus_Ref_Mapped_Only = extract_MAPQ("Mito/Field_1_C_harengus_Ref_Mapped_Only.sam")
Field_1_P_platessa_Ref_Mapped_Only = extract_MAPQ("Mito/Field_1_P_platessa_Ref_Mapped_Only.sam")
Field_1_S_scombrus_Ref_Mapped_Only = extract_MAPQ("Mito/Field_1_S_scombrus_Ref_Mapped_Only.sam")
Field_2_C_harengus_Ref_Mapped_Only = extract_MAPQ("Mito/Field_2_C_harengus_Ref_Mapped_Only.sam")
Field_2_P_platessa_Ref_Mapped_Only = extract_MAPQ("Mito/Field_2_P_platessa_Ref_Mapped_Only.sam")
Field_2_S_scombrus_Ref_Mapped_Only = extract_MAPQ("Mito/Field_2_S_scombrus_Ref_Mapped_Only.sam")
Field_3_C_harengus_Ref_Mapped_Only = extract_MAPQ("Mito/Field_3_C_harengus_Ref_Mapped_Only.sam")
Field_3_P_platessa_Ref_Mapped_Only = extract_MAPQ("Mito/Field_3_P_platessa_Ref_Mapped_Only.sam")
Field_3_S_scombrus_Ref_Mapped_Only = extract_MAPQ("Mito/Field_3_S_scombrus_Ref_Mapped_Only.sam")

#################################### P platessa bucket on Sscombrus Reference Ribo ####################################
P_platessa_P_platessa_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/P_platessa_P_platessa_Ref_Mapped_Only.sam") 

#################################### P platessa bucket on Sscombrus Reference Ribo ####################################
P_platessa_S_scombrus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/P_platessa_S_scombrus_Ref_Mapped_Only.sam") 

C_harengus_bucket_P_platessa_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/C_harengus_P_platessa_Ref_Mapped_Only.sam")
C_harengus_bucket_S_scombrus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/C_harengus_S_scombrus_Ref_Mapped_Only.sam")
P_platessa_bucket_C_harengus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/P_platessa_C_harengus_Ref_Mapped_Only.sam")
S_scombrus_bucket_C_harengus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/S_scombrus_C_harengus_Ref_Mapped_Only.sam")
S_scombrus_bucket_P_platessa_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/S_scombrus_P_platessa_Ref_Mapped_Only.sam")
Field_1_C_harengus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_1_C_harengus_Ref_Mapped_Only.sam")
Field_1_P_platessa_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_1_P_platessa_Ref_Mapped_Only.sam")
Field_1_S_scombrus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_1_S_scombrus_Ref_Mapped_Only.sam")
Field_2_C_harengus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_2_C_harengus_Ref_Mapped_Only.sam")
Field_2_P_platessa_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_2_P_platessa_Ref_Mapped_Only.sam")
Field_2_S_scombrus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_2_S_scombrus_Ref_Mapped_Only.sam")
Field_3_C_harengus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_3_C_harengus_Ref_Mapped_Only.sam")
Field_3_P_platessa_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_3_P_platessa_Ref_Mapped_Only.sam")
Field_3_S_scombrus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/Field_3_S_scombrus_Ref_Mapped_Only.sam")

S_scombrus_bucket_S_scombrus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/S_scombrus_S_scombrus_Ref_Mapped_Only.sam")
C_harengus_bucket_C_harengus_Ref_Mapped_Only_Ribo = extract_MAPQ("rDNA/C_harengus_C_harengus_Ref_Mapped_Only.sam")


#################################### STATISTICAL VISUALIZATION ####################################


############# HISTOGRAM: 
#Here do not need return since the figure is saved, do not need a variable to store.
def plot_histogram(mapq_list, name):
    
    mapq_arry = np.array(mapq_list)
    style.use('ggplot')
    # Set a predefined picture size to save in format 
    plt.figure(figsize=(7,5))
    # Plot a histogram of sequence lengths after DIAMOND & Filtering BLAST matches
    # Setting amount of bins & range of the graph. 
    plt.hist(mapq_arry, bins = 10, range = [min(mapq_arry), 70], color = 'darkorchid', edgecolor='black')
    # Setting title, x and y labels. 
    plt.title(f"MAPQ scores of {name}")
    plt.xlabel('MAPQ score')
    plt.ylabel('Amount of reads')
    # Determining to show the interval of x-axis ticks. 
    plt.xticks(np.arange(0, 70, 10))
    plt.savefig(f"{name}.png", dpi=200, bbox_inches='tight')
    # Savefig does not close the plot. 
    plt.clf()
    plt.close()

plot_histogram(P_platessa_P_platessa_Ref_Mapped_Only_Mito, "P. platessa bucket reads against P. platessa mitochondrial genome")
plot_histogram(P_platessa_S_scombrus_Ref_Mapped_Only_Mito, "P. platessa bucket reads against S. scombrus mitochondrial genome")
plot_histogram(P_platessa_P_platessa_Ref_Mapped_Only_Ribo, "P. platessa bucket reads against P. platessa ribosomal repeat")
plot_histogram(P_platessa_S_scombrus_Ref_Mapped_Only_Ribo, "P. platessa bucket reads against S. scombrus ribosomal repeat")
plot_histogram(C_harengus_bucket_P_platessa_Ref_Mapped_Only_Mito, "C. harengus bucket reads against P. platessa mitochondrial genome")
plot_histogram(C_harengus_bucket_S_scombrus_Ref_Mapped_Only_Mito, "C. harengus bucket reads against S. scombrus mitochondrial genome")
plot_histogram(P_platessa_bucket_C_harengus_Ref_Mapped_Only_Mito, "P. platessa bucket reads against C. harengus mitochondrial genome")
plot_histogram(S_scombrus_bucket_C_harengus_Ref_Mapped_Only_Mito, "S. scombrus bucket reads against C. harengus mitochondrial genome")
plot_histogram(S_scombrus_bucket_P_platessa_Ref_Mapped_Only_Mito, "S. scombrus bucket reads against P. platessa mitochondrial genome")
plot_histogram(C_harengus_bucket_P_platessa_Ref_Mapped_Only_Ribo, "C. harengus bucket reads against P. platessa ribosomal repeat")
plot_histogram(C_harengus_bucket_S_scombrus_Ref_Mapped_Only_Ribo, "C. harengus bucket reads against S. scombrus ribosomal repeat")
plot_histogram(P_platessa_bucket_C_harengus_Ref_Mapped_Only_Ribo, "P. platessa bucket reads against C. harengus ribosomal repeat")
plot_histogram(S_scombrus_bucket_C_harengus_Ref_Mapped_Only_Ribo, "S. scombrus bucket reads against C. harengus ribosomal repeat")
plot_histogram(S_scombrus_bucket_P_platessa_Ref_Mapped_Only_Ribo, "S. scombrus bucket reads against P. platessa ribosomal repeat")
plot_histogram(Field_1_C_harengus_Ref_Mapped_Only_Ribo, "Field 1 reads against C. harengus ribosomal repeat")
plot_histogram(Field_1_P_platessa_Ref_Mapped_Only_Ribo, "Field 1 reads against P. platessa ribosomal repeat")
plot_histogram(Field_1_S_scombrus_Ref_Mapped_Only_Ribo, "Field 1 reads against S. scombrus ribosomal repeat")
plot_histogram(Field_2_C_harengus_Ref_Mapped_Only_Ribo, "Field 2 reads against C. harengus ribosomal repeat")
plot_histogram(Field_2_P_platessa_Ref_Mapped_Only_Ribo, "Field 2 reads against P. platessa ribosomal repeat")
plot_histogram(Field_2_S_scombrus_Ref_Mapped_Only_Ribo, "Field 2 reads against S. scombrus ribosomal repeat")
plot_histogram(Field_3_C_harengus_Ref_Mapped_Only_Ribo, "Field 3 reads against C. harengus ribosomal repeat")
plot_histogram(Field_3_P_platessa_Ref_Mapped_Only_Ribo, "Field 3 reads against P. platessa ribosomal repeat")
plot_histogram(Field_3_S_scombrus_Ref_Mapped_Only_Ribo, "Field 3 reads against S. scombrus ribosomal repeat")


def plot_combined_histogram(mapq_list_1, mapq_list_2, mapq_list_3, name):
    
    mapq_arry_1 = np.array(mapq_list_1)
    mapq_arry_2 = np.array(mapq_list_2)
    mapq_arry_3 = np.array(mapq_list_3)
    
    style.use('ggplot')
    # Set a predefined picture size to save in format 
    # Once a figure set can draw multiple plots on the same figure
    plt.figure(figsize=(7,5))
    
    # Setting amount of bins & range of the graph, set to 61 since last value not taken in.
    # Define for both histograms a different color to seperate them. Alpha value to make them more transparent  
    plt.hist(mapq_arry_1, bins = 10, alpha = 0.5, range = [min(mapq_arry_1), 70], color = 'Red', label = 'P. platessa Ref')
    plt.hist(mapq_arry_2, bins = 10, alpha = 0.5, range = [min(mapq_arry_2), 70], color = 'Green', label = 'S. scombrus Ref')
    plt.hist(mapq_arry_3, bins = 10, alpha = 0.5, range = [min(mapq_arry_2), 70], color = 'Purple', label = 'C. harengus Ref')
    # Setting title, x and y labels. 
    plt.title(f"MAPQ scores of {name}")
    plt.xlabel('MAPQ score')
    plt.ylabel('Amount of reads')
    # Determining to show the interval of x-axis ticks. 
    plt.xticks(np.arange(0, 70, 10))
    plt.legend(loc = 'upper left')
    plt.savefig(f"{name}.png", dpi=200, bbox_inches='tight')
    # Savefig does not close the plot. 
    plt.clf()
    plt.close()

plot_combined_histogram(P_platessa_P_platessa_Ref_Mapped_Only_Mito, P_platessa_S_scombrus_Ref_Mapped_Only_Mito, P_platessa_bucket_C_harengus_Ref_Mapped_Only_Mito, "P platessa reads against Mitochondrial genome")
plot_combined_histogram(P_platessa_P_platessa_Ref_Mapped_Only_Ribo, P_platessa_S_scombrus_Ref_Mapped_Only_Ribo, P_platessa_bucket_C_harengus_Ref_Mapped_Only_Ribo, "P platessa reads against Ribosomal repeat")
plot_combined_histogram(S_scombrus_bucket_P_platessa_Ref_Mapped_Only_Mito, S_scombrus_bucket_S_scombrus_Ref_Mapped_Only_Mito, S_scombrus_bucket_C_harengus_Ref_Mapped_Only_Mito, "S. scombrus reads against Mitochondrial genome")
plot_combined_histogram(S_scombrus_bucket_P_platessa_Ref_Mapped_Only_Ribo, S_scombrus_bucket_S_scombrus_Ref_Mapped_Only_Ribo, S_scombrus_bucket_C_harengus_Ref_Mapped_Only_Ribo, "S. scombrus reads against Ribosomal repeat")
plot_combined_histogram(C_harengus_bucket_P_platessa_Ref_Mapped_Only_Mito, C_harengus_bucket_S_scombrus_Ref_Mapped_Only_Mito, C_harengus_bucket_C_harengus_Ref_Mapped_Only_Mito, "C. harengus reads against Mitochondrial genome")
plot_combined_histogram(C_harengus_bucket_P_platessa_Ref_Mapped_Only_Ribo, C_harengus_bucket_S_scombrus_Ref_Mapped_Only_Ribo, C_harengus_bucket_C_harengus_Ref_Mapped_Only_Ribo, "C. harengus reads against Ribosomal repeat")



############################### BARPLOT ###############################
# Pleuronectes_platessa_Mito mapped sequences: 23342
# Pleuronectes_platessa_Ribo mapped sequences: 64153
# --------------------------------------------------------------------------
# Pleuronectes_platessa_Bucket on Scomber_scombrus_Mito_Reference: 3443
# Pleuronectes_platessa_Bucket on Scomber_scombrus_Ribo_Reference: 37254
# --------------------------------------------------------------------------
# P_platessa_bucket_C_harengus_Ref Mito mapped sequences: 2129
# P_platessa_bucket_C_harengus_Ref Ribo mapped sequences: 23899


# Define the species and their corresponding data
organelle = ("Mitochondrion", "Ribosomal")
reference = {
    'P. platessa Ref': (23342, 64153),
    'S. scombrus Ref': (3443, 37254),
    'C. harengus Ref': (2129, 23899)
}

# Create an array of x-axis locations for the bars
x = np.arange(len(organelle))

# Set the width of each bar
width = 0.25

# Multiplier to position the x-axis ticks in the middle of the grouped bars
multiplier = 0

# Create a figure and axis for the plot
fig, ax = plt.subplots()

# Changing the fiure size
fig.set_size_inches(8, 5)

# Define a list of colors for each set of bars
colors = ['coral', 'mediumaquamarine', 'skyblue']

# Iterate through the penguin_means dictionary to plot the bars
for i, (attribute, measurement) in enumerate(reference.items()):
    offset = width * multiplier  # Calculate the offset for each set of bars
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=colors[i], edgecolor = 'Black')
    ax.bar_label(rects, padding=3)  # Add labels to the bars
    multiplier += 1  # Increment the multiplier to separate the next set of bars

# Add labels and titles to the plot
ax.set_ylabel('Mapped reads')
ax.set_title('Amount of mapped P. platessa reads from bucket samples')
ax.set_xticks(x + width, organelle)  # Set custom x-axis tick positions and labels
ax.legend(loc='upper left', ncols=2)  # Add a legend to the plot
ax.set_ylim(0, 170000)  # Set the y-axis limit
ax.set_facecolor("floralwhite") # The color of the background

# Display the plot
#plt.show()

plt.savefig(f"TotalMappedReads_Pp_Bucket.png", dpi=200, bbox_inches='tight')
plt.clf()

############################### BARPLOT ###############################
# C_harengus_bucket_P_platessa_Ref Mito mapped sequences: 73
# C_harengus_bucket_S_scombrus_Ref Mito mapped sequences: 50

# C_harengus_bucket_P_platessa_Ref Ribo mapped sequences: 334
# C_harengus_bucket_S_scombrus_Ref Ribo mapped sequences: 687

# Clupea_harengus_Mito Mito mapped sequences: 582
# Clupea_harengus_Ribo Ribo mapped sequences: 371

# Define the species and their corresponding data
organelle = ("Mitochondrion", "Ribosomal")
reference = {
    'C. harengus Ref': (582, 371),
    'P. platessa Ref': (73, 334),
    'S. scombrus Ref': (50, 687)
}

# Create an array of x-axis locations for the bars
x = np.arange(len(organelle))

# Set the width of each bar
width = 0.25

# Multiplier to position the x-axis ticks in the middle of the grouped bars
multiplier = 0

# Create a figure and axis for the plot
fig, ax = plt.subplots()

# Changing the fiure size
fig.set_size_inches(8, 5)

# Define a list of colors for each set of bars
colors = ['coral', 'mediumaquamarine', 'skyblue']

# Iterate through the penguin_means dictionary to plot the bars
for i, (attribute, measurement) in enumerate(reference.items()):
    offset = width * multiplier  # Calculate the offset for each set of bars
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=colors[i], edgecolor = 'Black')
    ax.bar_label(rects, padding=3)  # Add labels to the bars
    multiplier += 1  # Increment the multiplier to separate the next set of bars

# Add labels and titles to the plot
ax.set_ylabel('Mapped reads')
ax.set_title('Amount of C. harengus mapped reads from bucket samples')
ax.set_xticks(x + width, organelle)  # Set custom x-axis tick positions and labels
ax.legend(loc='upper left', ncols=2)  # Add a legend to the plot
ax.set_ylim(0, 2000)  # Set the y-axis limit
ax.set_facecolor("floralwhite") # The color of the background

# Display the plot
#plt.show()

plt.savefig(f"TotalMappedReads_Ch_Bucket.png", dpi=200, bbox_inches='tight')
plt.clf()

############################### BARPLOT ###############################
# S_scombrus_bucket_C_harengus_Ref Mito mapped sequences: 6850
# S_scombrus_bucket_P_platessa_Ref Mito mapped sequences: 14410

# Scomber_scombrus_Mito Mito mapped sequences: 60767
# Scomber_scombrus_Ribo Ribo mapped sequences: 40111

# S_scombrus_bucket_C_harengus_Ref Ribo mapped sequences: 12493
# S_scombrus_bucket_P_platessa_Ref Ribo mapped sequences: 13532

# Define the species and their corresponding data
organelle = ("Mitochondrion", "Ribosomal")
reference = {
    'S. scombrus Ref': (60767, 40111),
    'P. platessa Ref': (14410, 13532),
    'C. harengus Ref': (6850, 12493)
}

# Create an array of x-axis locations for the bars
x = np.arange(len(organelle))

# Set the width of each bar
width = 0.25

# Multiplier to position the x-axis ticks in the middle of the grouped bars
multiplier = 0

# Create a figure and axis for the plot
fig, ax = plt.subplots()

# Changing the fiure size
fig.set_size_inches(8, 5)

# Define a list of colors for each set of bars
colors = ['coral', 'mediumaquamarine', 'skyblue']

# Iterate through the penguin_means dictionary to plot the bars
for i, (attribute, measurement) in enumerate(reference.items()):
    offset = width * multiplier  # Calculate the offset for each set of bars
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=colors[i], edgecolor = 'Black')
    ax.bar_label(rects, padding=3)  # Add labels to the bars
    multiplier += 1  # Increment the multiplier to separate the next set of bars

# Add labels and titles to the plot
ax.set_ylabel('Mapped reads')
ax.set_title('Amount of S. scombrus mapped reads from bucket samples')
ax.set_xticks(x + width, organelle)  # Set custom x-axis tick positions and labels
ax.legend(loc='upper right', ncols=2)  # Add a legend to the plot
ax.set_ylim(0, 70000)  # Set the y-axis limit
ax.set_facecolor("floralwhite") # The color of the background

# Display the plot
#plt.show()

plt.savefig(f"TotalMappedReads_Ss_Bucket.png", dpi=200, bbox_inches='tight')
plt.clf()

################################ BARPLOT ################################
# For each Field sample need to how many reads with a MAPQ of 60 map on 1 or more fishes.
#A dictionary where the read  

def count_mapping(reg_expr):
    # Empty list
    mapping_dict = {}
    # The folder to read in
    for file in os.listdir("rDNA"):
        # Only want certain files to open
        if re.match(reg_expr, file):
            # Open the GFF-3 file with the MAPQ values
            file_path = os.getcwd() + "/rDNA/" + file
            #Open the file to read in
            with open(file_path, "r") as gff_to_read:
                # Convert the file to readlines list
                gff_lines = gff_to_read.readlines()#[0:5]
                
                # Print each item from the list
                for line in gff_lines:
                    # Split on tabular spaces -- > New list with sepreate values per header
                    splitted_line = line.split("\t")
                    
                    # If the MAPQ is 60 then will append the hit where the read mapped to the dictionary
                    if int(splitted_line[4]) == 60:
                        
                        if splitted_line[0] not in mapping_dict:
                            mapping_dict[splitted_line[0]] = [splitted_line[2]]
                        else:
                            mapping_dict[splitted_line[0]].append(splitted_line[2])
                    
            gff_to_read.close()
            
    return mapping_dict

Field_1 = count_mapping("^Field_1")
Field_2 = count_mapping("^Field_2")
Field_3 = count_mapping("^Field_3")

# Example on how to print a dicitonary
# print(Field_1['LH00260:17:22CG72LT3:5:1178:16284:17839_1:N:0:CTTGTACACC+AAGCGCGCTT'][2])

def Map_per_species(field_list):
    # The count per species -- > Graph
    Pp_count = 0
    Ss_count = 0
    Ch_count = 0 
    PpSs_count = 0
    PpCh_count = 0
    SsCh_count = 0
    PpChSs_count = 0 

    # Need to iterate over the dictionary
    for hits in field_list.values():
        #print(f"{read}: {hits}")
        
        # For now want to have only the results of a read mapping to differnt species, possible same reads maps well (MAPQ == 60) in the same species but difference in CIGAR string
        unique_hits = set()
        for hit in hits:
            
            # Split the item in the list on underscore > Take the 2 first items > Concatenate them back togheter
            unique_hits.add("_".join(hit.split("_")[:2]))
        # Create a code in form of a string of unique matches 
        code = ""
        for match in unique_hits:
            code += match
        
        #print(code)
        
        score = 0
        # Set up the possible flags what to do.
        if re.search("Pleuronectes_platessa", code):
            score += 1
        if re.search("Scomber_scombrus", code):
            score += 2
        if re.search("Clupea_harengus", code):
            score += 4
        
        if score == 1: Pp_count += 1
        elif score == 2: Ss_count += 1
        elif score == 3: PpSs_count += 1
        elif score == 4: Ch_count += 1
        elif score == 5: PpCh_count += 1
        elif score == 6: SsCh_count += 1
        elif score == 7: PpChSs_count += 1

    return {'P platessa': Pp_count, 'S scombrus':Ss_count, 'Clupea harengus':Ch_count, 'P platessa &\nS scombrus':PpSs_count, 'P platessa &\nC harengus':PpCh_count, 'S scombrus &\nC harengus':SsCh_count, 'P platessa &\nS scombrus &\n C harengus':PpChSs_count}

Map_per_species_1 = Map_per_species(Field_1)
Map_per_species_2 = Map_per_species(Field_2)
Map_per_species_3 = Map_per_species(Field_3)


def Field_Mapping_Species(sample, dictionary):
    style.use('ggplot')
    plt.figure(figsize = (18,10))
    plt.bar(dictionary.keys(), dictionary.values(), 
            edgecolor = 'black',
            color = 'mediumaquamarine')
    plt.title(f'Mapping of a read in {sample} sample',
              fontsize = 22)
    plt.xlabel('On which species a read maps',
               fontsize = 18,
               weight = 'bold')
    plt.xticks(fontsize = 16)
    plt.ylabel('Amount of mapped reads',
               fontsize = 18,
               weight = 'bold')
    plt.yticks(fontsize = 16)
    
    plt.tight_layout()
    plt.savefig(f"MappingPerSpecies{sample}.png", dpi=200, bbox_inches='tight')
    
Field_Mapping_Species("Field 1", Map_per_species_1)
Field_Mapping_Species("Field 2", Map_per_species_2)
Field_Mapping_Species("Field 3", Map_per_species_3)


################################ GROUPED BARPLOT ################################
def uniqe_mult_map(Dict):
    count_unique_mapping = 0
    count_mult_mapping = 0

    for v in Dict.values():
        if v == 1: count_unique_mapping += 1
        else: count_mult_mapping += 1
    
    return count_unique_mapping, count_mult_mapping


count_unique_mapping = 0
count_mult_mapping = 0

for v in Field_1.values():
    if v == 1: count_unique_mapping += 1
    else: count_mult_mapping += 1
    
Field_1_Unique_Map, Field_1_Mult_Map = uniqe_mult_map(Field_1)
Field_2_Unique_Map, Field_2_Mult_Map = uniqe_mult_map(Field_2)
Field_3_Unique_Map, Field_3_Mult_Map = uniqe_mult_map(Field_3)


# Define the species and their corresponding data
organelle = ("Field 1", "Field 2", "Field 3")
reference = {
    'Maps Once': (Field_1_Unique_Map, Field_2_Unique_Map, Field_3_Unique_Map),
    'Maps Multiple times': (Field_1_Mult_Map, Field_2_Mult_Map, Field_3_Mult_Map)
}

style.use('ggplot')

# Create an array of x-axis locations for the bars
x = np.arange(len(organelle))

# Set the width of each bar
width = 0.25

# Multiplier to position the x-axis ticks in the middle of the grouped bars
multiplier = 0.5

# Create a figure and axis for the plot
fig, ax = plt.subplots()

# Changing the fiure size
fig.set_size_inches(8, 5)

# Define a list of colors for each set of bars
colors = ['coral', 'mediumaquamarine']

# Iterate through the penguin_means dictionary to plot the bars
for i, (attribute, measurement) in enumerate(reference.items()):
    offset = width * multiplier  # Calculate the offset for each set of bars
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=colors[i], edgecolor = 'Black')
    ax.bar_label(rects, padding=3)  # Add labels to the bars
    multiplier += 1  # Increment the multiplier to separate the next set of bars

# Add labels and titles to the plot
ax.set_ylabel('Mapped reads')
ax.set_title('Amount of Mapped reads with MAPQ of 60 per Field sample')
ax.set_xticks(x + width, organelle)  # Set custom x-axis tick positions and labels
ax.legend(loc='upper right', ncols=2)  # Add a legend to the plot
ax.set_ylim(0, 400)  # Set the y-axis limit
ax.set_facecolor("floralwhite") # The color of the background

# Display the plot
#plt.show()

plt.savefig(f"MappingField.png", dpi=200, bbox_inches='tight')
plt.clf()





