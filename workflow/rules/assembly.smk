rule fastp:
    input:
        r1 = "raw_data/{sample}_1.fastq.gz",
        r2 = "raw_data/{sample}_2.fastq.gz"
    output:
        r1_clean = "results/01_clean_data/{sample}_1.fastq.gz",
        r2_clean = "results/01_clean_data/{sample}_2.fastq.gz",
        html_report = "results/qc_reports/fastp/{sample}_fastp.html",
        json_report = "results/qc_reports/fastp/{sample}_fastp.json"
    params:
        qual = config["fastp"]["qual_threshold"],
        length = config["fastp"]["min_length"]
    threads: 4
    log:
        "results/logs/01_fastp_{sample}.log"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1_clean} -O {output.r2_clean} \
            -h {output.html_report} -j {output.json_report} \
            -q {params.qual} -l {params.length} \
            --detect_adapter_for_pe \
            -c \
            --thread {threads} \
            2> {log}
        """

rule spades:
    input:
        r1 = "results/01_clean_data/{sample}_1.fastq.gz",
        r2 = "results/01_clean_data/{sample}_2.fastq.gz"
    output:
        contigs = "results/02_assembly/{sample}/contigs.fasta",
        scaffolds = "results/02_assembly/{sample}/scaffolds.fasta"
    params:
        outdir = "results/02_assembly/{sample}",
        # A função lambda busca o modo do SPAdes específico para esta amostra na tabela
        mode = lambda wildcards: samples_df.loc[wildcards.sample, "spades_mode"]
    threads: 8
    log:
        "results/logs/02_spades_{sample}.log"
    shell:
        """
        spades.py \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            -t {threads} \
            {params.mode} \
            > {log} 2>&1
        """

# REMOVE CONTIGS PEQUENAS (MENORES DE 500pb)
rule filter_contigs:
    input:
        scaffolds = "results/02_assembly/{sample}/scaffolds.fasta"
    output:
        filtered = "results/02_assembly/{sample}/scaffolds_filtered.fasta"
    params:
        min_len = config["filtering"]["min_contig_length"]
    log:
        "results/logs/02b_filter_{sample}.log"
    shell:
        """
        seqtk seq -L {params.min_len} {input.scaffolds} > {output.filtered} 2> {log}
        """

rule quast:
    input:
        assembly = "results/02_assembly/{sample}/scaffolds_filtered.fasta"
    output:
        report = "results/qc_reports/quast/{sample}/report.tsv"
    params:
        outdir = "results/qc_reports/quast/{sample}"
    threads: 4
    log:
        "results/logs/03_quast_{sample}.log"
    shell:
        """
        quast.py {input.assembly} \
                 -o {params.outdir} \
                 -t {threads} \
                 > {log} 2>&1
        """

rule busco:
    input:
        assembly = "results/02_assembly/{sample}/scaffolds_filtered.fasta"
    output:
        # Usamos uma função para definir o output exato baseado na linhagem da tabela
        summary = lambda wildcards: f"results/qc_reports/busco/{wildcards.sample}/short_summary.specific.{samples_df.loc[wildcards.sample, 'busco_lineage']}.{wildcards.sample}.txt"
    params:
        out_path = "results/qc_reports/busco",
        lineage = lambda wildcards: samples_df.loc[wildcards.sample, "busco_lineage"]
    threads: 8
    log:
        "results/logs/03_busco_{sample}.log"
    shell:
        """
        busco -i {input.assembly} \
              -l {params.lineage} \
              -o {wildcards.sample} \
              --out_path {params.out_path} \
              -m genome \
              -c {threads} \
              -f \
              > {log} 2>&1
        """

rule multiqc:
    input:
        fastp = expand("results/qc_reports/fastp/{sample}_fastp.json", sample=SAMPLES),
        quast = expand("results/qc_reports/quast/{sample}/report.tsv", sample=SAMPLES),
        busco = BUSCO_OUTPUTS
    output:
        report = "results/qc_reports/multiqc/multiqc_report.html"
    params:
        search_dir = "results/qc_reports",
        outdir = "results/qc_reports/multiqc"
    log:
        "results/logs/05_multiqc.log"
    shell:
        """
        multiqc {params.search_dir} \
            -o {params.outdir} \
            -f \
            > {log} 2>&1
        """

rule custom_report:
    input:
        quast = "results/qc_reports/quast/{sample}/report.tsv",
        busco = lambda wildcards: f"results/qc_reports/busco/{wildcards.sample}/short_summary.specific.{samples_df.loc[wildcards.sample, 'busco_lineage']}.{wildcards.sample}.txt"
    output:
        report = "results/qc_reports/{sample}_assembly_summary.txt"
    shell:
        """
        python scripts/generate_summary.py {input.quast} {input.busco} {output.report} {wildcards.sample}
        """

# Incluir uma rule RagTag para amostras que possuam genoma de referencia (para alinhar as contigs)