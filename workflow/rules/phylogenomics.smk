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


PASSOS PIPELINE
- Rodar BUSCO para todos genomas
- Extrair os ortólogos comuns para todos os genomas
- Ver quais genomas tem baixa completude do BUSCO e não incluir eles na análise para evitar baixar muito a quantidade de genes
- Proteoma? faa?
- Rodar árvore concatenada E rodar árvore por coalescência

"""


import pandas as pd
from pathlib import Path
from snakemake.io import glob_wildcards

configfile: "config/config.yaml"
include: "rules/assembly.smk"

DATA_ROOT = Path("data/phylogenomics")

# ---------- Descoberta ----------
def list_clades():
    return sorted(p.name for p in DATA_ROOT.iterdir() if p.is_dir())

def list_genomes(clade):
    """Todos os FASTAs em data/phylogenomics/<clade>/genomes/."""
    pattern = str(DATA_ROOT / clade / "genomes" / "{sample}.fa.gz")
    _, samples = glob_wildcards(pattern)
    return sorted(samples)

def list_by_role(clade, role):
    """role in {'SAMPLE', 'REF', 'OUTGROUP'}"""
    return [s for s in list_genomes(clade) if s.startswith(f"{role}_")]

def list_outgroups(clade):
    return list_by_role(clade, "OUTGROUP")


# ---------- SANITY CHECK (roda primeiro) ----------
rule validate_layout:
    """Verifica que todo genoma tem um prefixo reconhecido e que há >=1 outgroup."""
    input:
        lambda w: [str(DATA_ROOT / w.clade / "genomes" / f"{s}.fa.gz")
                   for s in list_genomes(w.clade)],
    output:
        "results/phylogenomics/{clade}/layout.ok"
    run:
        valid_prefixes = ("SAMPLE_", "REF_", "OUTGROUP_")
        bad = [f for f in input if not Path(f).name.startswith(valid_prefixes)]
        if bad:
            raise ValueError(
                f"Genomas sem prefixo reconhecido em {wildcards.clade}: {bad}\n"
                f"Use SAMPLE_, REF_ ou OUTGROUP_."
            )
        if not list_outgroups(wildcards.clade):
            raise ValueError(
                f"Nenhum OUTGROUP_* em {wildcards.clade}. "
                f"A árvore não pode ser enraizada."
            )
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).touch()


# ---------- BUSCO ----------
rule busco:
    input:
        assembly   = lambda w: str(DATA_ROOT / w.clade / "genomes" / f"{w.sample}.fa.gz"),
    output:
        outdir = directory("results/phylogenomics/{clade}/busco/{sample}"),
    params:
        out_path = "results/phylogenomics/{clade}/busco/",
        lineage  = lambda w: config["clades"][w.clade]["busco_lineage"],
    container:
        config["containers"]["busco"],
    threads: 8,
    log:
        "results/logs/busco/{clade}_{sample}.log",
    shell:
        """
        export PYTHONUNBUFFERED=1
        export OMP_NUM_THREADS={threads}

        busco -i {input.assembly} \
              -l {params.lineage} \
              -o {wildcards.sample} \
              --out_path {params.out_path} \
              --download_path data/dbs/busco \
              --offline --opt-out-run-stats \
              -m genome -c {threads} -f \
              > {log} 2>&1
        """


# ---------- INVENTÁRIO AUTOMÁTICO ----------
rule inventory:
    """Gera um TSV com metadados a partir do layout + resultados BUSCO.
       Este TSV é derivado — NÃO é fonte de verdade."""
    input:
        buscos  = lambda w: expand(
            "results/phylogenomics/{clade}/busco/{sample}",
            clade=[w.clade],
            sample=list_genomes(w.clade),
        ),
        layout  = "results/phylogenomics/{clade}/layout.ok",
    output:
        tsv = "results/phylogenomics/{clade}/inventory.tsv",
    run:
        rows = []
        for s in list_genomes(wildcards.clade):
            role = next(r for r in ("SAMPLE", "REF", "OUTGROUP") if s.startswith(r + "_"))
            gcf  = s.split("_", 1)[1] if role != "SAMPLE" else None
            rows.append({"sample": s, "role": role, "gcf": gcf})
        df = pd.DataFrame(rows)
        df.to_csv(output.tsv, sep="\t", index=False)


# ---------- rule all ----------
rule all:
    input:
        # layout validado para todos os clados
        expand("results/phylogenomics/{clade}/layout.ok", clade=list_clades()),
        # BUSCO para todos os genomas de todos os clados
        [f"results/phylogenomics/{clade}/busco/{s}"
         for clade in list_clades()
         for s in list_genomes(clade)],
        # inventário
        expand("results/phylogenomics/{clade}/inventory.tsv", clade=list_clades()),
        # resultados finais
        expand("results/phylogenomics/{clade}/trees/astral.treefile",
               clade=list_clades()),