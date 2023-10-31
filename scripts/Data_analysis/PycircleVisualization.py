# Importing the correct modules
from pygenomeviz import GenomeViz


# ITSx: SSU: 1-1586, ITS1: 1587-2294, 5.8S: 2295-2452, ITS2: 2453-2930, LSU: 2931-10293
# Barrnap: SSU: 8709 - 10289, 5.8S: 7846 - 7998, LSU: 3358 - 8006

genome_list = (
    {"name": "P. platessa Ribosomal - Barrnap", "size": 10293, "rdna_list": ((8709, 10289, -1), (7846, 7998, -1), (3358, 8006, -1))},
    {"name": "P. platessa Ribosomal - ITSx", "size": 10293, "rdna_list": ((1587, 2294, 1), (2453, 2930, 1))},
    {"name": "P. platessa Reverse Comp. Ribosomal - Barrnap", "size": 10293, "rdna_list": ((3004, 4841, 1), (5544, 10192, 1), (5552, 5704, 1))},
)

label_1 = ["SSU","5.8S","LSU"]
label_2 = ["ITS1", "ITS2"]
label_3 = ["SSU", "LSU", "5.8S"]

color_1 = "orange"
color_2 = "tomato"

count = 1

gv = GenomeViz(tick_style="axis")

for genome in genome_list:
    name, size, rdna_list = genome["name"], genome["size"], genome["rdna_list"]
    track = gv.add_feature_track(name, size)
    
    # Using the correct labels for the annotation
    if count == 1: label, color = label_1, color_1
    elif count == 2: label, color = label_2, color_2
    elif count == 3: label, color = label_3, color_1
    #Incrementing the count
    count += 1
    
    for idx, rrna in enumerate(rdna_list, 1):
    
        start, end, strand = rrna
        track.add_feature(start, end, strand, label = label[idx - 1], linewidth=1, labelrotation=0, labelvpos="top", labelhpos="center", labelha="center", facecolor = color)
        
gv.savefig("P_platessa_Ribosomal.png")



# ITSx: SSU: 1-2925, ITS1: 2926-3574, 5.8S: 3575-3732, ITS2: 3733-4202, LSU: 4203-10818
# Barrnap: SSU: 1089 - 2924, 5.8S: 3576- 3728, LSU: 3571- 8115

genome_list = (
    {"name": "S. scombrus Ribosomal - Barrnap", "size": 10818, "rdna_list": ((1089, 2925, 1), (3575, 3732, 1), (4203, 8115, 1))},
    {"name": "S. scombrus Ribosomal - ITSx", "size": 10818, "rdna_list": ((2926, 3574, 1), (3733, 4202, 1))},
)

label_1 = ["SSU","5.8S","LSU"]
label_2 = ["ITS1", "ITS2"]

color_1 = "orange"
color_2 = "tomato"

count = 1

gv = GenomeViz(tick_style="axis")

for genome in genome_list:
    name, size, rdna_list = genome["name"], genome["size"], genome["rdna_list"]
    track = gv.add_feature_track(name, size)
    
    # Using the correct labels for the annotation
    if count == 1: label, color = label_1, color_1
    elif count == 2: label, color = label_2, color_2
    #Incrementing the count
    count += 1
    
    for idx, rrna in enumerate(rdna_list, 1):
    
        start, end, strand = rrna
        track.add_feature(start, end, strand, label = label[idx - 1], linewidth=1, labelrotation=0, labelvpos="top", labelhpos="center", labelha="center", facecolor = color)
 
gv.savefig("S_scombrus_Ribosomal.png")



# Setting the name and genome size
# ITSx 2715 bp: SSU: 1-1859, ITS1: 1860-2244, 5.8S: 2245-2401, ITS2: 2402-2715, LSU: Not found
# ITSx 4936 bp: SSU: Not found, ITS1: Not found, 5.8S: Not found, ITS2: 1-69, LSU: 70-4936, Broken or partial sequence, no 5.8S!
# Barrnap:
# Scaffold 1: SSU: 858 - 2709, 5.8S: 319 - 470 
# Scaffold 2:LSU: 69 - 4028


genome_list = (
    {"name": "Contig 1 - Barrnap", "size": 2715, "rdna_list": ((858, 2709, 1), (319, 470, 1))},
    {"name": "Contig 1 - ITSx", "size": 2715, "rdna_list": ((1860, 2244, 1), (2402, 2715, 1))},
    {"name": "Contig 2 - Barrnap", "size": 4936, "rdna_list": ((69, 4028, 1),)},
    {"name": "Contig 2 - ITSx", "size": 4936, "rdna_list": ((1, 69, 1),)},
)

label_1 = ["SSU","5.8S",]
label_2 = ["ITS1", "ITS2"]
label_3 = ["LSU"]
label_4 = ["ITS2"]

color_1 = "orange"
color_2 = "tomato"

count = 1

gv = GenomeViz(tick_style="axis")
for genome in genome_list:
    name, size, rdna_list = genome["name"], genome["size"], genome["rdna_list"]
    track = gv.add_feature_track(name, size)
    
    # Using the correct labels for the annotation
    if count == 1: label, color = label_1, color_1
    elif count == 2: label, color = label_2, color_2
    elif count == 3: label, color = label_3, color_1
    elif count == 4: label, color = label_4, color_2
    #Incrementing the count
    count += 1
    
    for idx, rrna in enumerate(rdna_list, 1):
        #print(idx, rrna)
        start, end, strand = rrna
        track.add_feature(start, end, strand, label = label[idx - 1], linewidth=1, labelrotation=0, labelvpos="top", labelhpos="center", labelha="center", facecolor = color)


gv.savefig("C_harengus_Ribosomal.png")