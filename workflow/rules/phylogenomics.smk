# ==============================================================================
# assembly.smk 
# contains the rules necessary for phylogenomics. gets called by snakefile
# ==============================================================================

"""
Input:
Genoma montado + GCFs de referência para o clado (FASTAs) + GCF do outgroup do clado

Métodos: 
Concatenar tudo e rodar uma árvore só
Rodar uma árvore para cada gene e juntar

Output:
Árvore filogenética (gerar png e .treefile)

Dúvidas:
Como fazer a amostragem de táxons??? Devo pegar o máximo de genomas próximos de boa qualidade? 
Já começar a prospectar genomas para os clados 


Passos:
1) criar uma pasta dentro de phylogenomics para o clado que vc quer gerar
2) criar um arquivo .yaml contendo:

Samples (from pipeline - names should be formatted as ..):
- mrakii

OTUs/taxons: (pensar no nome)
- GCF123123
- GCA123124

Outgroup:
- GCF123123

"""

