############################# WORKFLOW #############################
# 
# The AllAssemblies_workflow will perform the 2 possible assembly methods of Genome Skimming pipeline, without building the KRAKEN Custom DB. 
# 1) Assembling mitochondrial reference genomes, ribosomal repeats generating a Genbank file with the results summerized. 
# 2) Performing the Assembly free method to save all reads after contamination filtering of Species
#
###################################################################
#
############################# RULES - PROGRAMS #############################
include:
    '../rules/PreprocessingStrict.smk'
include:
    '../rules/MitoAssembly.smk'
include:
    '../rules/RibosomalAssembly.smk'
include:
    '../rules/PreAssemblies.smk'
include:
    '../rules/MitoAnnotation.smk'
include:
    '../rules/RibosomalAnnotation.smk'
include:
    '../rules/MergePE.smk'
include:
    '../rules/ContaminationFiltering.smk'
include:
    '../rules/ExtractKrakenReads.smk'
include:
    '../rules/CreateGenbank.smk'
include: 
    '../rules/CreateReferenceGenomeMapping.smk'
    
############################# RULE - RESULTING ANALYSIS #############################
# This functions as an "rule all" where the input is expected to be the finalized output. 
# By using touch can create an empty file to fake the output. 
rule Analysis:
    input:
        expand(
            [
                os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_paired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_for_unpaired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_paired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/02_FastQC_Results/{SAMPLE}_back_unpaired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/06_MITOS_Results/{SAMPLE}/AnnotationMitosDone.txt"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/08_Barrnap_Anno_Results/{SAMPLE}/{SAMPLE}_rDNA.gff"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/14_Kraken2_Library_Reads/{TAXID}_{SAMPLE}.fasta"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/16_Genbank/{TAXID}_{SAMPLE}.done"),
                os.path.join(DATA_DIR_GS, "{PROJECT}/17_Reference_Mapping_MG/{TAXID}_{SPECIES}.fasta"),
                
            ],
            PROJECT = config["genome_skimming"]["project"],
            SAMPLE = config["genome_skimming"]["sample"],
            TAXID = config['genome_skimming']['taxid'],
            SPECIES = config["genome_skimming"]["species"]
        )
    output:
        touch('AllAssembly_analysis.done')


