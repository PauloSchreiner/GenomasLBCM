"""

anotação funcional: 
eggNOG-mapper para atribuição de GO
InterProScan para identificação de domínios conservados
antiSMASH para buscar metabólitos secundários (BGCs)
dbCAN/HMMER - para buscar CAZymes 

"""

# Regra temporária para não quebrar a execução. Removeremos assim que criarmos a do eggNOG.
rule placeholder_annotation:
    input:
        "results/02_assembly/{sample}/final_assembly.fasta"
    output:
        "results/03_annotation/{sample}/placeholder.txt"
    shell:
        """
        echo "Anotação da amostra {wildcards.sample} virá aqui!" > {output}
        """