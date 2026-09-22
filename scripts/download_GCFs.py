'''
input: 
1) tsv with list of GCF codes
2) directory
'''

import pandas as pd

df = pd.read_csv(sep="\t", )

# abrir tsv
# para cada linha, rodar:
f"datasets download genome accession {gcf} --include genome --filename {taxon}_genome.gz"