import os, glob


rule BWA_Index:
    input:
        PAIRED_1 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        #UNPAIRED_1 = expand(
        #    os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R1_Unpaired.fastq"),
        #    project = config["metagenomics"]["project"], sample=config["metagenomics"]["sample"]
        #    ),
        PAIRED_2 = expand(
            os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Paired.fastq"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            ),
        #UNPAIRED_2 = expand(
        #    os.path.join(DATA_DIR_MG, "{project}/01_Trimmomatic_Results/{sample}_R2_Unpaired.fastq"),
        #    project = config["metagenomics"]["project"], sample=config["metagenomics"]["sample"]
        #    )        
    output:
        expand(
            os.path.join(DATA_DIR_MG, "{project}/07_BWA_Mapped_Sequences/{sample}_Mapping.done"),
            project = config["metagenomics"]["project"], sample = config["metagenomics"]["sample"]
            )
        
    conda:
        "../envs/bwa-mem2.yaml"

    params:
        FolderReferenceFasta =  config["metagenomics"]["FolderReferenceFasta"],
        OutputFolderIndexes = expand(os.path.join(DATA_DIR_MG, "{project}/06_Mapped_Indexes/"),
            project = config["metagenomics"]["project"]
            ),
        Sample = config["metagenomics"]["sample"],
        IndexOutput = "/home/genomics/gleberre/01_Research_BAR_ZAND/01_BAR/02_BAR_MG/06_Mapped_Indexes/"
    

    shell:
        """
        python3 scripts/RunBWA-MEM2.py -f {params.FolderReferenceFasta} -1 {input.PAIRED_1} -2 {input.PAIRED_2} -s {params.Sample} -io {params.IndexOutput}
        touch {output}
        """

#    for file in "{params.FolderReferenceFasta}/*"
#        do
#            echo $file
#            filename=$(basename "$file" | cut -d. -f1)
#            bwa-mem2 index -p "{params.OutputFolderIndexes}/$filename" $file
#            bwa-mem2 mem "{params.OutputFolderIndexes}/$filename" "{input.PAIRED_1}" "{input.PAIRED_2}" -o {params.Sample}_"$filename".sam
#        done
#    touch "{output}"