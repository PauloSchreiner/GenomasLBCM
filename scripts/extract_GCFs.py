import subprocess
import json
import pandas as pd

caminho_tsv = "./especies_filogenomica.tsv"

df = pd.read_csv(caminho_tsv, sep="\t")

# inicializa a coluna GCF com "" como padrão se ela nao existir
if "GCF" not in df.columns:
    df["GCF"] = None

# puxa as especies
especies = list(df["especies"])


for e in especies:
    # pula se ja contem
    if (df.loc[df["especies"] == e, "GCF"] != None).any():
        print("\nProcessando " + e + "...")

        result = subprocess.run(["datasets", "summary", "genome", "taxon", e,
                    "--reference", 
                    "--assembly-source", "refseq",
                    "--as-json-lines"],
                    capture_output=True, text=True)

        if result.stdout != "":
            data_dict = json.loads(result.stdout)
            gcf = data_dict["current_accession"] 

            print(e + " POSSUI genoma de refereência de código " + gcf)
            df.loc[df["especies"] == e, "GCF"] = gcf

            # print(json.dumps(data_dict, sort_keys=True, indent=4))
        else: 
            df.loc[df["especies"] == e, "GCF"] = "Not found"
            print(e + " não possui genoma de referencia")

print(df)

print(f"Salvando resultados em {caminho_tsv}...")
df.to_csv(caminho_tsv, index=False, sep="\t")

print(f"\nDe {len(df)}, {(df['GCF'] != "Not found").sum()} têm genoma de referência!")