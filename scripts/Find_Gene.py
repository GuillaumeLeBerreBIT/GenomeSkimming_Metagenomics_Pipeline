from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

gb_file = open("test.gb", "r")

for gb_record in SeqIO.parse(gb_file, "genbank"):
    # now do something with the record
    for feat in gb_record.features:
        if feat.type == "gene" \
            and 'cox1' in feat.qualifiers["gene"]:
                print(feat.qualifiers["sequence"], feat.qualifiers["gene"])     


var = ["file", ".fastq"]

print(var[1][1:])       
            