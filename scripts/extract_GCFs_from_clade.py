import subprocess
import json
import pandas as pd
from datetime import datetime

tag = datetime.now().strftime("%Y-%m-%d_%H:%M")

caminho_tsv = f"./GCFs_Tremellales_{tag}.tsv"

df = pd.DataFrame(columns=["taxon", "GCF"])

grupos_taxonomicos = ["Tremellales"]

for grupo in grupos_taxonomicos:
    print("\nExtraindo espécies de " + grupo + "...")

    result = subprocess.run(["datasets", "summary", "genome", "taxon", grupo,
                "--reference", 
                "--assembly-source", "refseq",
                "--as-json-lines"],
                capture_output=True, text=True)

    if result.stdout != "":

        data_list = result.stdout.split("\n")

        for species_data in data_list:
            species_dict = json.loads(species_data)

            print(json.dumps(species_dict, sort_keys=True, indent=4))


            gcf = species_dict["current_accession"] 
            species = species_dict["organism"]["organism_name"] ## O erro ta aqui

            print(species + " possui genoma de referência de código " + gcf)

            df.loc[len(df)] = [species, gcf]


print(df)

# print(f"Salvando resultados em {caminho_tsv}...")
# df.to_csv(caminho_tsv, index=False, sep="\t")

# print(f"\nDe {len(df)}, {(df['GCF'] != "Not found").sum()} têm genoma de referência!")
