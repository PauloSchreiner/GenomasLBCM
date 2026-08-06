# ==============================================================================
# assembly.smk 
# contains the rules necessary for annotation. gets called by snakefile
# ==============================================================================

"""
anotação funcional: 
eggNOG-mapper para atribuição de GO
InterProScan para identificação de domínios conservados
antiSMASH para buscar metabólitos secundários (BGCs)
dbCAN/HMMER - para buscar CAZymes 

"""

# placeholder
rule placeholder_annotation:
    input:
        "results/02_assembly/{sample}/final_assembly.fasta"
    output:
        "results/03_annotation/{sample}/placeholder.txt"
    shell:
        """
        echo "annotation for {wildcards.sample} !" > {output}
        """