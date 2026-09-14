import subprocess
import json
import pandas as pd

caminho_tsv = "./especies_filogenomica.tsv"

df = pd.read_csv(caminho_tsv)

# inicializa a coluna GCF com "" como padrão se ela nao existir
if "GCF" not in df.columns:
    df["GCF"] = ""


# puxa as especies
especies = list(df["especies"])

for e in especies:
    # pula se
    if df.loc[[df["especies"] == e], "GCF"] != "":    
        print("\nProcessando " + e + "...")

        result = subprocess.run(["datasets", "summary", "genome", "taxon", e,
                    "--reference", 
                    "--assembly-source", "refseq",
                    "--as-json-lines"],
                    capture_output=True, text=True)

        if result.stdout != "":
            data_dict = json.loads(result.stdout)
            gcf = data_dict["current_accession"] 

            # data_pretty = json.dumps(data_dict, sort_keys=True, indent=4)
            # print(data_pretty)

            print(e + " is " + gcf)

            df.loc[df["especies"] == e, "GCF"] = gcf
        else: 
            df.loc[df["especies"] == e, "GCF"] = "None"
            print(e + " não possui genoma de referencia")

print(df)

print(f"\nDe {len(df)}, {(df['GCF'] != "None").sum()} têm genoma de referência!")