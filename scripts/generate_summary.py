import sys
import json
import re
import datetime

# Command-line arguments mapping.
# The script is designed to act as a bridge between the raw outputs of various
# bioinformatics tools and the structured data needs of the project's tracking system.
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
host_remove = sys.argv[11]
downsample_to = sys.argv[12]

# Standardized timestamp to link this report to its specific execution environment capsule.
timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


# --- 1. PRE-PROCESSING METRICS (FASTP) ---
# Evaluates library complexity and sequencing machine accuracy. 
# High duplication rates or low Q30 scores often explain downstream assembly fragmentation.
with open(fastp_file, 'r') as f:
    fastp_data = json.load(f)

reads_before = fastp_data['summary']['before_filtering']['total_reads']
reads_after = fastp_data['summary']['after_filtering']['total_reads']
dup_rate = fastp_data['duplication']['rate'] * 100
q30_after = fastp_data['summary']['after_filtering']['q30_rate'] * 100
# Insert size is critical for evaluating whether paired-end reads overlap effectively.
insert_size = fastp_data.get('insert_size', {}).get('peak', 'N/A')


# --- 2. STRUCTURAL METRICS (QUAST) ---
# Extracts contiguity metrics. Filtering out contigs under 1000bp removes noise 
# (like residual adapters or unresolved repetitive elements) to provide a clearer 
# picture of the biologically relevant assembly size.
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


# --- 3. TAXONOMIC METRICS (KRAKEN2) ---
# Identifies potential environmental or cross-sample contamination.
# We restrict parsing to Species ('S') or Genus ('G') ranks to avoid flooding the 
# report with high-level domains (e.g., "Bacteria") that lack actionable specificity.
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

# Enforce uniqueness to prevent the report from displaying multiple overlapping strains
# of the exact same species, which clutters the overview.
top_hits = sorted(top_hits, key=lambda x: x[0], reverse=True)
unique_hits = []
seen_names = set()

for pct, name in top_hits:
    if name not in seen_names:
        unique_hits.append(f"{name} ({pct}%)")
        seen_names.add(name)
    if len(unique_hits) == 3: break
        
# Pad data structure to maintain TSV column consistency if the sample is highly pure.
while len(unique_hits) < 3:
    unique_hits.append("N/A")


# --- 4. BIOLOGICAL COMPLETENESS (BUSCO) ---
# Assesses the presence of single-copy orthologs. A highly contiguous assembly (good N50) 
# can still be biologically poor if it collapses essential genes or artificially duplicates them.
busco_metrics = {}
with open(busco_file, 'r') as f:
    for line in f:
        if "Total BUSCO groups searched" in line:
            busco_metrics['Total'] = line.strip().split()[0]
        # BUSCO outputs are unstructured text files. Regex provides a resilient method 
        # to extract the C/S/D/F/M percentages regardless of minor version formatting changes.
        elif "C:" in line and "%" in line:
            summary_str = line.strip()
            busco_metrics['C_pct'] = re.search(r'C:(.*?)%', summary_str).group(1)
            busco_metrics['S_pct'] = re.search(r'S:(.*?)%', summary_str).group(1)
            busco_metrics['D_pct'] = re.search(r'D:(.*?)%', summary_str).group(1)
            busco_metrics['F_pct'] = re.search(r'F:(.*?)%', summary_str).group(1)
            busco_metrics['M_pct'] = re.search(r'M:(.*?)%', summary_str).group(1)


# --- 5. HUMAN-READABLE OUTPUT GENERATION ---
# Formats the parsed data into a quick-reference text file for immediate visual inspection.
with open(out_txt, 'w') as out:
    out.write(f"======================================================\n")
    out.write(f" GENOMIC ASSEMBLY REPORT - SAMPLE: {sample.upper()}\n")
    out.write(f"======================================================\n\n")
    
    out.write(f"[ ASSEMBLY PARAMETERS ]\n")
    out.write(f"Host Contamination Removed: {host_remove}\n")
    out.write(f"Reads Downsampled to      : {downsample_to}\n")
    out.write(f"SPAdes Mode               : {spades_mode}\n")
    out.write(f"SPAdes k-mers             : {kmers}\n")
    out.write(f"Min Contig Length Filter  : {min_len} bp\n\n")
    
    out.write(f"[ RAW DATA & FILTERING (FASTP) ]\n")
    out.write(f"Total Reads (Before)      : {reads_before}\n")
    out.write(f"Total Reads (After)       : {reads_after}\n")
    out.write(f"Q30 Rate (After)          : {q30_after:.2f}%\n")
    out.write(f"Estimated Insert Size     : {insert_size} bp\n")
    out.write(f"Duplication Rate          : {dup_rate:.2f}%\n\n")
    
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
    out.write(f"N's per 100 kbp           : {ns_per_100kbp}\n")
    out.write(f"GC Content                : {gc}%\n\n")
    
    out.write(f"[ GENE COMPLETENESS (BUSCO) ]\n")
    out.write(f"Total BUSCOs searched     : {busco_metrics.get('Total', 'N/A')}\n")
    out.write(f"Complete (C)              : {busco_metrics.get('C_pct', 'N/A')}%\n")
    out.write(f"  - Single-copy (S)       : {busco_metrics.get('S_pct', 'N/A')}%\n")
    out.write(f"  - Duplicated (D)        : {busco_metrics.get('D_pct', 'N/A')}%\n")
    out.write(f"Fragmented (F)            : {busco_metrics.get('F_pct', 'N/A')}%\n")
    out.write(f"Missing (M)               : {busco_metrics.get('M_pct', 'N/A')}%\n")
    out.write(f"======================================================\n")


# --- 6. PROGRAMMATIC DATABASE PAYLOAD ---
# Maps the structured variables strictly to the headers defined in the Bash shell rule.
# This ensures a flat architecture where data can be seamlessly appended to the persistent 
# database across multiple runs without complex merging scripts.
tsv_data = [
    sample,
    timestamp,
    host_remove,
    downsample_to,
    spades_mode,
    kmers,
    min_len,
    reads_before,
    reads_after,
    f"{q30_after:.2f}",
    insert_size,
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