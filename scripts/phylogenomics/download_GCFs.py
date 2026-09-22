import os
import shutil
import subprocess
import zipfile
from pathlib import Path
import pandas as pd
import gzip

diretorio_destino = Path("data/phylogenomics/clade1/ref_genomes/")
caminho_planilha = "GCFs_Tremellales.tsv"

# Cria o diretório de destino caso não exista
diretorio_destino.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(caminho_planilha, sep="\t")

for index, row in df.iterrows():
    gcf = str(row["GCF"]).strip()
    taxon = str(row["taxon"]).replace(" ", "_").strip()

    print(f"[{index+1}/{len(df)}] Baixando: {taxon} ({gcf})...")

    zip_temp = diretorio_destino / f"{taxon}_temp.zip"
    fasta_gz_final = diretorio_destino / f"{taxon}_genome.fa.gz"

    # Se já existir o .fa.gz, pula
    if fasta_gz_final.exists():
        print(f"  -> Arquivo {fasta_gz_final.name} já existe. Pulando...")
        continue

    cmd_download = [
        "datasets",
        "download",
        "genome",
        "accession",
        gcf,
        "--include",
        "genome",
        "--filename",
        str(zip_temp),
    ]

    try:
        subprocess.run(cmd_download, check=True)

        with zipfile.ZipFile(zip_temp, "r") as zip_ref:
            fasta_in_zip = [
                f for f in zip_ref.namelist() if f.endswith(".fna")
            ]

            if fasta_in_zip:
                # Extrai do ZIP e comprime diretamente para .fa.gz
                with (
                    zip_ref.open(fasta_in_zip[0]) as f_in,
                    gzip.open(fasta_gz_final, "wb") as f_out,
                ):
                    shutil.copyfileobj(f_in, f_out)
                print(f"  -> Salvo e compactado: {fasta_gz_final.name}")
            else:
                print(f"  [AVISO] Nenhum FASTA encontrado no ZIP de {gcf}")

    except subprocess.CalledProcessError as e:
        print(f"  [ERRO] Falha ao baixar o acesso {gcf}: {e}")
    finally:
        if zip_temp.exists():
            os.remove(zip_temp)


