# #!/usr/bin/python3
############################# INTRODUCTION #############################
# 
# This is to perform the annotation of mitochondrial genome by MITOS. 
# Here MITOS is pulled from a doker container to perform the analysis. 
#
###################################################################
# 
############################# MODULES #############################
import os, glob

############################# RULES #############################
# The mitochondrial annotation doen by MITOS
# Since the input fasta file has a variable name and only the last part is consistent: *.path_sequence.fasta
# Due to wanting to annotate only the first path can use >> *.1.1.path_sequence.fasta 
# MITOS has a variable output, if the input fasta file has only 1 sequence then it will create the output directly in {sample}/
# BUT when the input is a multifasta it will create subfolders 0/, 1/, ... n/ so the output is not consistent
# Will use a dummy file as output file which can be done with 'touch' command from bash

rule MitosMito:
    input:
        expand(os.path.join(DATA_DIR_GS, "{project}/04_GetOrganelle_Mito_Results/{sample}/get_org.log.txt"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            )
    output:
    # Use a dummy file as target output to avoid the complications of variable output
        expand(os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            )

    #conda:
    #    "../envs/mitos.yaml"
    singularity:
        "docker://nanozoo/mitos:2.1.0--4cc80aa"

    params:
        # Define folders under params and not as output otherwise will not work. 
        MitosFolder = expand(os.path.join(DATA_DIR_GS, "{project}/06_MITOS_Results/{sample}/"),
            project = config["genome_skimming"]["project"], sample = config["genome_skimming"]["sample"]
            ),
        RefFolder = expand(os.path.join(DB_DIR, "databases/MITOS_DB/"),
            DB_DIR = config['mitos']['database']),
        MitRefSeq = "refseq89m/",
        gencode = config['mitos']['gencode'],
        # This is to only take one FASTA file as input for the annotation
        #AnnoFile = os.path.join(DATA_DIR, "{project}/04_GetOrganelle_Mito_Results/{sample}/*1.1.path_sequence.fasta"),
        # Want to only take the first matching file from GetOrganelle >> Alphabetically "Complete" comes before "Contigs" thus will preffered take "Complete" fasta file

        ## DOCKER PATHS
        InputSing = expand(os.path.join(DK_DIR, "04_GetOrganelle_Mito_Results/{sample}/*1.1.path_sequence.fasta"),
            sample = config["genome_skimming"]["sample"]
            ),
        OutputSing = expand(os.path.join(DK_DIR, "06_MITOS_Results/{sample}/AnnotationMitosDone.txt"),
            sample = config["genome_skimming"]["sample"]
            ),
        FolderSing = expand(os.path.join(DK_DIR, "06_MITOS_Results/{sample}/"),
            sample = config["genome_skimming"]["sample"]
            )
        # The output file will be saved in those folders, thus being accesable on your local system
    
    # unset _JAVA_OPTIONS

    shell:
        """
        first_matching_file=$(ls -1 {params.InputSing} | head -n 1)
        runmitos.py -i $first_matching_file -c {params.gencode} -r {params.MitRefSeq} -o {params.FolderSing} --refdir {params.RefFolder}
        touch {params.OutputSing}
        """
# When encountering the following error >> DELETE THE EXISTING FOLDER 
#
#ERROR:root:mitfi exception
#Exception in thread "main" java.lang.StringIndexOutOfBoundsException: begin 0, end 10, length 0
#	at java.base/java.lang.String.checkBoundsBeginEnd(String.java:3319)
#	at java.base/java.lang.String.substring(String.java:1874)
#	at mitfi.Main.main(Main.java:347)

