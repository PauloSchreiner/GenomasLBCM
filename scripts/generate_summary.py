import sys
import re

# Fetch arguments passed by Snakemake/terminal
quast_file = sys.argv[1]
busco_file = sys.argv[2]
output_file = sys.argv[3]
sample_name = sys.argv[4]

# Initial logic and parameters are validated. Standard processing applied.
quast_metrics = {}
with open(quast_file, 'r') as q:
    for line in q:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            quast_metrics[parts[0]] = parts[1]

busco_stats = {"C": "N/A", "S": "N/A", "D": "N/A", "F": "N/A", "M": "N/A", "n": "N/A"}
with open(busco_file, 'r') as b:
    for line in b:
        if "C:" in line and "F:" in line and "M:" in line:
            # Jumping directly to regex pattern matching for BUSCO parsing
            match = re.search(r"C:(.*?)%\[S:(.*?)%,D:(.*?)%\],F:(.*?)%,M:(.*?)%,n:(\d+)", line)
            if match:
                busco_stats["C"] = match.group(1)
                busco_stats["S"] = match.group(2)
                busco_stats["D"] = match.group(3)
                busco_stats["F"] = match.group(4)
                busco_stats["M"] = match.group(5)
                busco_stats["n"] = match.group(6)
            break

# Final transformation and report formatting
report = f"""
======================================================
 GENOMIC ASSEMBLY REPORT - SAMPLE: {sample_name.upper()}
======================================================

[ CONTIGUITY STATISTICS (QUAST) ]
Total Length (>= 1000 bp) : {quast_metrics.get('Total length (>= 1000 bp)', 'N/A')} bp
Total Contigs (>= 1000 bp): {quast_metrics.get('# contigs (>= 1000 bp)', 'N/A')}
Largest Contig            : {quast_metrics.get('Largest contig', 'N/A')} bp
N50                       : {quast_metrics.get('N50', 'N/A')} bp
L50                       : {quast_metrics.get('L50', 'N/A')}
GC Content                : {quast_metrics.get('GC (%)', 'N/A')}%

[ GENE COMPLETENESS (BUSCO) ]
Total BUSCOs searched     : {busco_stats['n']}
Complete (C)              : {busco_stats['C']}%
  - Single-copy (S)       : {busco_stats['S']}%
  - Duplicated (D)        : {busco_stats['D']}%
Fragmented (F)            : {busco_stats['F']}%
Missing (M)               : {busco_stats['M']}%

[ QUICK DIAGNOSTICS ]
- An N50 below 10,000 bp indicates a highly fragmented assembly.
- A GC content deviating significantly from the expected range for the genus may suggest contamination.
- A high percentage of Fragmented (F) BUSCOs is usually a direct consequence of a low N50.
======================================================
"""

with open(output_file, 'w') as out:
    out.write(report.strip() + "\n")