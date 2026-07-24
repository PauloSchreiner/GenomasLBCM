# YeastAssemble: Pipeline de Montagem e Anotação Genômica de Leveduras

Uma pipeline bioinformática automatizada, modular e reprodutível construída em **Snakemake**. Desenvolvida para processar genomas de leveduras selvagens (Ascomicetos e Basidiomicetos), lidando com os desafios de montagens heterozigotas a partir de dados Illumina (paired-end).

---

## 🧬 Arquitetura da Pipeline

A pipeline é dividida em dois grandes módulos lógicos:

### 1. Módulo de Montagem e QC (`assembly.smk`)
* **Controle de Qualidade (QC) e Trimming:** Limpeza de dados brutos e remoção de adaptadores utilizando o `fastp`.
* **Montagem *de novo*:** Montagem do genoma utilizando o `SPAdes`. O modo de montagem pode ser customizado por amostra (ex: `--isolate`, `--meta`).
* **Filtragem de Contigs:** Remoção de fragmentos pequenos e ruidosos (< 500 pb) utilizando o `seqtk`.
* **Avaliação da Montagem:** 
  * Métricas de contiguidade (N50, L50, tamanho total) via `QUAST`.
  * Avaliação de completude gênica via `BUSCO` (banco de dados customizável por amostra).
* **Consolidação de Relatórios:** Geração de um dashboard interativo com `MultiQC` e um sumário executivo em texto puro via script Python customizado.

### 2. Módulo de Anotação Funcional (`annotation.smk`) - *[Em Desenvolvimento]*
* **eggNOG-mapper:** Atribuição de ontologia gênica (GO) e ortologia.
* **InterProScan:** Identificação de domínios proteicos conservados.
* **antiSMASH:** Mineração de clusters de biossíntese de metabólitos secundários (BGCs).
* **dbCAN/HMMER:** Busca por enzimas ativas em carboidratos (CAZymes) - *Foco: Metabolismo de Xilose.*
* **RagTag:** Ordenação de scaffolds baseada em genoma de referência (opcional).

---

## 📂 Estrutura do Projeto

```text
YeastAssemble/
├── raw_data/                  # Coloque os arquivos brutos aqui (ex: amostra_1.fastq.gz, amostra_2.fastq.gz)
├── results/                   # Pasta gerada automaticamente pela pipeline
│   ├── 01_clean_data/
│   ├── 02_assembly/
│   ├── 03_annotation/
│   ├── qc_reports/
│   └── logs/
├── scripts/
│   └── generate_summary.py    # Script gerador do boletim médico da montagem
├── workflow/
│   ├── Snakefile              # Orquestrador principal da pipeline
│   └── rules/
│       ├── assembly.smk       # Regras da montagem
│       └── annotation.smk     # Regras da anotação
├── config.yaml                # Arquivo de configuração global
└── samples.tsv                # Banco de dados de amostras e parâmetros específicos
```

---

## ⚙️ Configuração

Antes de rodar a pipeline, configure os dois arquivos principais na raiz do projeto:

### 1. `config.yaml`
Controla os parâmetros globais das ferramentas:
```yaml
fastp:
  qual_threshold: 20
  min_length: 50

spades:
  mode: "--isolate"  # Default, mas é sobrescrito pelo samples.tsv
# mode: "--careful" 
# mode: "--meta" 

filtering:
  min_contig_length: 500

busco:
  lineage: "saccharomycetes_odb10"  # Default, mas é sobrescrito pelo samples.tsv
# lineage: "ascomycota_odb10"
# lineage: "basidiomycota_odb10"
```

### 2. `samples.tsv`
Define quais amostras serão processadas e permite ajuste fino para genomas com biologia distinta (separação por *Tab*).
```tsv
sample	spades_mode	busco_lineage
pseudolambica	--isolate	saccharomycetes_odb10
levedura_red_1	--meta	basidiomycota_odb10
```
*(Nota: Os arquivos na pasta `raw_data/` devem ter exatamente o mesmo prefixo definido na coluna `sample`)*

---

## 🚀 Como Executar

A pipeline foi modularizada para permitir execuções parciais. 

**Para rodar a pipeline inteira (Montagem + Anotação):**
```bash
snakemake -c 8
```
*(Substitua `8` pelo número de núcleos de processamento disponíveis na sua máquina)*

**Para rodar APENAS a Montagem e o QC:**
Ideal para a etapa de iteração e refinamento de parâmetros.
```bash
snakemake assemble -c 8
```

**Para rodar APENAS a Anotação Funcional:**
Ideal para quando os genomas finais já estão montados e filtrados.
```bash
snakemake annotate -c 8
```

---

## 📊 Resultados

Os relatórios finais da montagem estarão disponíveis em:
* **`results/qc_reports/multiqc/multiqc_report.html`**: Dashboard interativo completo.
* **`results/qc_reports/{sample}_assembly_summary.txt`**: Resumo limpo e direto da amostra, ideal para anexar em cadernos de laboratório.