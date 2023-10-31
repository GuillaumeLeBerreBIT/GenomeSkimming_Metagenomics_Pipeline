############################# INTRODUCTION #############################
# 
# This file will contain all the rules necessary to perform the analysis of the Assembly. 
# The rules are split up in different files for efficient combining in different pipelines. 
#
###################################################################

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
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_paired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_for_unpaired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_paired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/02_FastQC_Results/{sample}_back_unpaired_fastqc.html"),
                os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
                os.path.join(DATA_DIR_GS, "{project}/08_Barrnap_Anno_Results/{sample}/{sample}_rDNA.gff"),
                os.path.join(DATA_DIR_GS, "{project}/16_Genbank/{TAXID}_{sample}.done"),
                os.path.join(DATA_DIR_GS, "{project}/17_Reference_Mapping_MG/{TAXID}_{species}.fasta"),
            ],
            project = config["genome_skimming"]["project"],
            sample = config["genome_skimming"]["sample"],
            TAXID = config['genome_skimming']['taxid'],
            species = config["genome_skimming"]["species"]
        )
    output:
        touch('Assembly_analysis.done')
