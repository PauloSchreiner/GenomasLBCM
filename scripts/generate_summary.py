import sys
import json
import re
import datetime
import os

# 1. Argumentos da linha de comando
quast_file = sys.argv[1]
busco_file = sys.argv[2]
fastp_file = sys.argv[3]
kraken_file = sys.argv[4]
out_txt = sys.argv[5]
out_tsv = sys.argv[6]
sample = sys.argv[7]
spades_mode = sys.argv[8]
kmers = sys.argv[9]
min_len = sys.argv[10]
host_flag = sys.argv[11]
downsample = sys.argv[12]
pypolca_file = sys.argv[13]

timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# 2. Parsing do PyPOLCA
polca_snps = "0"
polca_indels = "0"
if os.path.exists(pypolca_file):
    with open(pypolca_file, 'r') as f:
        for line in f:
            if "Substitution errors found" in line:
                polca_snps = line.strip().split(":")[-1].strip()
            elif "Insertion/Deletion errors found" in line:
                polca_indels = line.strip().split(":")[-1].strip()

# 3. Parsing do Hostile (Trata Lista e Dicionário)
hostile_removed = "N/A"
hostile_json = f"results/qc_reports/hostile/{sample}_hostile.json"
if str(host_flag).strip().lower() in ["human", "yes", "true"] and os.path.exists(hostile_json):
    try:
        with open(hostile_json, 'r') as f:
            h_data = json.load(f)
            h_dict = h_data[0] if isinstance(h_data, list) and len(h_data) > 0 else (h_data if isinstance(h_data, dict) else {})
            hostile_removed = str(h_dict.get("reads_removed", "0"))
    except Exception:
        hostile_removed = "Error"

# 4. Parsing do Fastp
with open(fastp_file, 'r') as f:
    fastp_data = json.load(f)

reads_before = fastp_data['summary']['before_filtering']['total_reads']
reads_after = fastp_data['summary']['after_filtering']['total_reads']
dup_rate = fastp_data['duplication']['rate'] * 100
q30_after = fastp_data['summary']['after_filtering']['q30_rate'] * 100
insert_size = fastp_data.get('insert_size', {}).get('peak', 'N/A')

# 5. Parsing do QUAST
quast_metrics = {}
with open(quast_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            quast_metrics[parts[0]] = parts[1]

total_len_1000 = quast_metrics.get('Total length (>= 1000 bp)', 'N/A')
total_contigs = quast_metrics.get('# contigs', 'N/A')
total_contigs_1000 = quast_metrics.get('# contigs (>= 1000 bp)', 'N/A')
largest_contig = quast_metrics.get('Largest contig', 'N/A')
n50 = quast_metrics.get('N50', 'N/A')
l50 = quast_metrics.get('L50', 'N/A')
n90 = quast_metrics.get('N90', 'N/A')
l90 = quast_metrics.get('L90', 'N/A')
ns_per_100kbp = quast_metrics.get("N's per 100 kbp", 'N/A')
gc = quast_metrics.get('GC (%)', 'N/A')

# 6. Parsing do Kraken2
unclassified_pct = "0.0%"
top_hits = []
with open(kraken_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 6: continue
        pct = float(parts[0].strip())
        rank = parts[3].strip()
        name = parts[5].strip()
        if rank == 'U':
            unclassified_pct = f"{pct}%"
        elif rank in ['S', 'G']:
            top_hits.append((pct, name))

top_hits = sorted(top_hits, key=lambda x: x[0], reverse=True)
unique_hits = []
seen_names = set()
for pct, name in top_hits:
    if name not in seen_names:
        unique_hits.append(f"{name} ({pct}%)")
        seen_names.add(name)
    if len(unique_hits) == 3: break
while len(unique_hits) < 3:
    unique_hits.append("N/A")

# 7. Parsing do BUSCO
busco_metrics = {}
with open(busco_file, 'r') as f:
    for line in f:
        if "Total BUSCO groups searched" in line:
            busco_metrics['Total'] = line.strip().split()[0]
        elif "C:" in line and "%" in line:
            summary_str = line.strip()
            busco_metrics['C_pct'] = re.search(r'C:(.*?)%', summary_str).group(1)
            busco_metrics['S_pct'] = re.search(r'S:(.*?)%', summary_str).group(1)
            busco_metrics['D_pct'] = re.search(r'D:(.*?)%', summary_str).group(1)
            busco_metrics['F_pct'] = re.search(r'F:(.*?)%', summary_str).group(1)
            busco_metrics['M_pct'] = re.search(r'M:(.*?)%', summary_str).group(1)

# 8. Geração do Relatório TXT
with open(out_txt, 'w') as out:
    out.write(f"======================================================\n")
    out.write(f" GENOMIC ASSEMBLY REPORT - SAMPLE: {sample.upper()}\n")
    out.write(f"======================================================\n\n")
    out.write(f"[ ASSEMBLY PARAMETERS ]\n")
    out.write(f"Decontaminate Human       : {host_flag}\n")
    out.write(f"Host Reads Removed        : {hostile_removed}\n")
    out.write(f"Reads Downsampled to      : {downsample}\n")
    out.write(f"SPAdes Mode               : {spades_mode}\n")
    out.write(f"SPAdes k-mers             : {kmers}\n")
    out.write(f"Min Contig Length Filter  : {min_len} bp\n\n")
    out.write(f"[ READ FILTERING & POLISHING ]\n")
    out.write(f"Total Reads (Raw)         : {reads_before}\n")
    out.write(f"Total Reads (Clean)       : {reads_after}\n")
    out.write(f"Q30 Rate (Clean)          : {q30_after:.2f}%\n")
    out.write(f"Estimated Insert Size     : {insert_size} bp\n")
    out.write(f"Duplication Rate          : {dup_rate:.2f}%\n")
    out.write(f"POLCA Base Corrections    : {polca_snps} SNPs | {polca_indels} Indels\n\n")
    out.write(f"[ CONTAMINATION CHECK (KRAKEN2) ]\n")
    out.write(f"Unclassified Fragments    : {unclassified_pct}\n")
    out.write(f"Top Hit 1                 : {unique_hits[0]}\n")
    out.write(f"Top Hit 2                 : {unique_hits[1]}\n")
    out.write(f"Top Hit 3                 : {unique_hits[2]}\n\n")
    out.write(f"[ CONTIGUITY STATISTICS (QUAST) ]\n")
    out.write(f"Total Length (>= 1000 bp) : {total_len_1000} bp\n")
    out.write(f"Total Contigs (All sizes) : {total_contigs}\n")
    out.write(f"Total Contigs (>= 1000 bp): {total_contigs_1000}\n")
    out.write(f"Largest Contig            : {largest_contig} bp\n")
    out.write(f"N50                       : {n50} bp\n")
    out.write(f"L50                       : {l50}\n")
    out.write(f"N90                       : {n90} bp\n")
    out.write(f"L90                       : {l90}\n")
    out.write(f"GC Content                : {gc}%\n\n")
    out.write(f"[ GENE COMPLETENESS (BUSCO) ]\n")
    out.write(f"Total BUSCOs searched     : {busco_metrics.get('Total', 'N/A')}\n")
    out.write(f"Complete (C)              : {busco_metrics.get('C_pct', 'N/A')}%\n")
    out.write(f"  - Single-copy (S)       : {busco_metrics.get('S_pct', 'N/A')}%\n")
    out.write(f"  - Duplicated (D)        : {busco_metrics.get('D_pct', 'N/A')}%\n")
    out.write(f"Fragmented (F)            : {busco_metrics.get('F_pct', 'N/A')}%\n")
    out.write(f"Missing (M)               : {busco_metrics.get('M_pct', 'N/A')}%\n")
    out.write(f"======================================================\n")

# 9. Geração da Linha TSV Estruturada (28 Colunas Alinhadas)
tsv_data = [
    sample,
    timestamp,
    host_flag,
    hostile_removed,
    downsample,
    spades_mode,
    kmers,
    min_len,
    reads_before,
    reads_after,
    f"{q30_after:.2f}",
    insert_size,
    polca_snps,
    polca_indels,
    unclassified_pct,
    unique_hits[0],
    unique_hits[1],
    total_len_1000,
    total_contigs_1000,
    largest_contig,
    n50,
    l50,
    gc,
    busco_metrics.get('C_pct', 'N/A'),
    busco_metrics.get('S_pct', 'N/A'),
    busco_metrics.get('D_pct', 'N/A'),
    busco_metrics.get('F_pct', 'N/A'),
    busco_metrics.get('M_pct', 'N/A')
]

with open(out_tsv, 'w') as f:
    f.write("\t".join(map(str, tsv_data)) + "\n")