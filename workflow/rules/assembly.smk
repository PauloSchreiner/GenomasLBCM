# ==============================================================================
# assembly.smk 
# contains the rules necessary for QC and assembly. gets called by snakefile
# ==============================================================================

# trimming and pre-processing
rule fastp:
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz"
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

# uses minimap2 and samtools to remove all reads that correspond to human contamination.
# to activate this, add "human" as the value of the column "host_remove" in samples.tsv.
# as of now, it only works for human genomes. in the future it may be necessary to expand 
rule decontaminate_host:
    input:
        r1="results/01_clean_data/{sample}_1.fastq.gz",
        r2="results/01_clean_data/{sample}_2.fastq.gz"
    output:
        r1="results/01b_decontam/{sample}_1.fastq.gz",
        r2="results/01b_decontam/{sample}_2.fastq.gz"
    params:
        host=lambda wildcards: samples.loc[wildcards.sample, "Host_Remove"],
        tool=config.get("decontaminator", "minimap2"),
        fasta=config.get("db_human_fasta", ""),
        bt2_idx=config.get("db_human_bowtie2", "")
    log:
        "results/logs/01b_decontam_{sample}.log"
    threads: 8
    shell:
        """
        set -euo pipefail
        
        # 1. Verifica se a amostra precisa de descontaminação
        if [ "{params.host}" == "human" ]; then
            
            # 2a. Trilho do Servidor (Minimap2)
            if [ "{params.tool}" == "minimap2" ]; then
                echo "[INFO] Iniciando Minimap2 para {wildcards.sample}..." > {log}
                minimap2 -a -x sr -t {threads} -I 8g -K 500M {params.fasta} {input.r1} {input.r2} > results/01b_decontam/{wildcards.sample}_tmp.sam 2>> {log}
                
                samtools fastq -f 12 -@ {threads} -1 {output.r1} -2 {output.r2} -s /dev/null -0 /dev/null results/01b_decontam/{wildcards.sample}_tmp.sam >> {log} 2>&1
                
                rm results/01b_decontam/{wildcards.sample}_tmp.sam
            
            # 2b. Trilho do Laptop (Bowtie2)
            elif [ "{params.tool}" == "bowtie2" ]; then
                echo "[INFO] Iniciando Bowtie2 para {wildcards.sample}..." > {log}
                # O parâmetro --un-conc-gz já cospe os FASTQs limpos magicamente!
                bowtie2 -x {params.bt2_idx} -1 {input.r1} -2 {input.r2} -p {threads} --un-conc-gz results/01b_decontam/{wildcards.sample}_%.fastq.gz > /dev/null 2>> {log}
            
            else
                echo "[ERROR] Ferramenta {params.tool} nao reconhecida!" >> {log}
                exit 1
            fi

        # 3. Se nao houver contaminacao, apenas cria o atalho (symlink)
        else
            echo "[INFO] Nenhuma descontaminacao solicitada. Criando symlinks..." > {log}
            ln -sf $(realpath {input.r1}) {output.r1}
            ln -sf $(realpath {input.r2}) {output.r2}
        fi
        """


# uses seqtk to downsample reads, which is necessary when there are too many. 
# to activate this, add the desired amount of reads in the column "downsample_to"
# of samples.tsv. Recommended amount is 4000000
rule downsample_reads:
    input:
        r1 = "results/01b_decontam/{sample}_1.fastq.gz",
        r2 = "results/01b_decontam/{sample}_2.fastq.gz"
    output:
        r1 = "results/01c_downsample/{sample}_1.fastq.gz",
        r2 = "results/01c_downsample/{sample}_2.fastq.gz"
    params:
        # selects corresponding value from samples.tsv
        target = lambda wildcards: samples_df.loc[wildcards.sample, "downsample_to"]
    threads: 2
    log:
        "results/logs/01c_downsample_{sample}.log"
    shell:
        """
        if [ "{params.target}" == "none" ]; then
            ln -sf $(realpath {input.r1}) {output.r1}
            ln -sf $(realpath {input.r2}) {output.r2}
        else
            # Selects randomly a defined amount of reads (seed 42 for reproducibility)
            seqtk sample -s 42 {input.r1} {params.target} | gzip > {output.r1}
            seqtk sample -s 42 {input.r2} {params.target} | gzip > {output.r2}
        fi
        """

# assembly
rule spades:
    input:
        r1 = "results/01c_downsample/{sample}_1.fastq.gz",
        r2 = "results/01c_downsample/{sample}_2.fastq.gz"
    output:
        contigs = "results/02_assembly/{sample}/contigs.fasta",
        scaffolds = "results/02_assembly/{sample}/scaffolds.fasta"
    params:
        outdir = "results/02_assembly/{sample}",
        # locates mode from samples.tsv
        mode = lambda wildcards: samples_df.loc[wildcards.sample, "spades_mode"],
        # selects k-mer sizes from samples.tsv. works for both auto and specified k values
        kmers = lambda wildcards: "" if str(samples_df.loc[wildcards.sample, "spades_k"]).strip() == "auto" else f"-k {samples_df.loc[wildcards.sample, 'spades_k']}",
        # Lógica para corte de cobertura ()
        cov = lambda wildcards: "" if str(samples_df.loc[wildcards.sample, "cov_cutoff"]).strip() == "auto" else f"--cov-cutoff {samples_df.loc[wildcards.sample, 'cov_cutoff']}"
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
            {params.kmers} \
            {params.cov} \
            > {log} 2>&1
        """

# REMOVE CONTIGS PEQUENAS (MENORES DE 500pb) 
rule filter_contigs:
    input:
        scaffolds = "results/02_assembly/{sample}/scaffolds.fasta"
    output:
        filtered = "results/02_assembly/{sample}/scaffolds_filtered.fasta"
    params:
        min_len = lambda wildcards: samples_df.loc[wildcards.sample, "filter_len"]
    log:
        "results/logs/02b_filter_{sample}.log"
    shell:
        """
        seqtk seq -L {params.min_len} {input.scaffolds} > {output.filtered} 2> {log}
        """


# Chamada so se a amostra tiver algo diferente de "none" na coluna "reference_genome" em samples.tsv
# RagTag usa um genoma de referência pra montar os scaffolds
rule ragtag_scaffold:
    input:
        assembly = "results/02_assembly/{sample}/scaffolds_filtered.fasta",
        reference = lambda wildcards: samples_df.loc[wildcards.sample, "reference_genome"]
    output:
        scaffolded = "results/02_assembly/{sample}/ragtag_output/ragtag.scaffold.fasta"
    params:
        outdir = "results/02_assembly/{sample}/ragtag_output"
    threads: 4
    log:
        "results/logs/02c_ragtag_{sample}.log"
    shell:
        """
        ragtag.py scaffold {input.reference} {input.assembly} -o {params.outdir} -t {threads} -u > {log} 2>&1
        """

# Padroniza o nome do arquivo final pra aceitar tanto 
# o resultado do RagTag quanto do seqtk 
rule unify_assembly:
    input:
        # se for "none", pega do seqtk. Se não, pega do RagTag
        lambda wildcards: (
            "results/02_assembly/{sample}/ragtag_output/ragtag.scaffold.fasta"
            if str(samples_df.loc[wildcards.sample, "reference_genome"]).strip().lower() != "none"
            else "results/02_assembly/{sample}/scaffolds_filtered.fasta"
        )
    output:
        final_assembly = "results/02_assembly/{sample}/final_assembly.fasta"
    shell:
        """
        # Cria um simlink com um nome unificado
        ln -sf $(realpath {input}) {output.final_assembly}
        """


# Checks for contamination 
rule kraken2:
    input:
        assembly = "results/02_assembly/{sample}/final_assembly.fasta"
    output:
        report = "results/qc_reports/kraken2/{sample}_kraken2_report.txt",
        output_txt = "results/qc_reports/kraken2/{sample}_kraken2_output.txt"
    params:
        db = config["kraken2"]["db"],
        mem = config["kraken2"]["mem"]
    threads: 4
    log:
        "results/logs/03b_kraken2_{sample}.log"
    shell:
        """
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            {params.mem} \
            --report {output.report} \
            --output {output.output_txt} \
            {input.assembly} \
            > {log} 2>&1
        """

rule quast:
    input:
        assembly = "results/02_assembly/{sample}/final_assembly.fasta"
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
        assembly = "results/02_assembly/{sample}/final_assembly.fasta"
    output:
        summary = "results/qc_reports/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt"
    params:
        out_path = "results/qc_reports/busco",
        lineage = lambda wildcards: samples_df.loc[wildcards.sample, "busco_lineage"]
    threads: 8
    log:
        "results/logs/03_busco_{sample}_{lineage}.log"
    shell:
        """
        busco -i {input.assembly} \
              -l {wildcards.lineage} \
              -o {wildcards.sample} \
              --out_path {params.out_path} \
              --download_path data/dbs/busco \
              -m genome \
              -c {threads} \
              -f \
              > {log} 2>&1
        """

        
rule multiqc:
    input:
        fastp = expand("results/qc_reports/fastp/{sample}_fastp.json", sample=SAMPLES),
        quast = expand("results/qc_reports/quast/{sample}/report.tsv", sample=SAMPLES),
        busco = BUSCO_OUTPUTS,
        kraken = expand("results/qc_reports/kraken2/{sample}_kraken2_report.txt", sample=SAMPLES)
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


# Apaga os arquivos intermediarios gerados pelo spades
# Eles ficam em pastas como K21/ e K33/
rule clean_spades_intermediates:
    input:
        "results/02_assembly/{sample}/scaffolds.fasta"
    output:
        # Cria um arquivo sinal
        touch("results/02_assembly/{sample}/.cleaned")
    shell:
        """
        rm -rf results/02_assembly/{wildcards.sample}/K[0-9]*
        rm -rf results/02_assembly/{wildcards.sample}/tmp
        """



# Generates both a human-readable text report and a programmatic TSV line for each sample.
# Aggregating metrics across multiple tools (Fastp, Quast, Kraken, BUSCO) into a single 
# step ensures that downstream analysis doesn't require parsing raw tool outputs manually.
rule custom_report:
    input:
        quast = "results/qc_reports/quast/{sample}/report.tsv",
        busco = lambda wildcards: f"results/qc_reports/busco/{wildcards.sample}/short_summary.specific.{samples_df.loc[wildcards.sample, 'busco_lineage']}.{wildcards.sample}.txt",
        fastp = "results/qc_reports/fastp/{sample}_fastp.json",
        kraken = "results/qc_reports/kraken2/{sample}_kraken2_report.txt"
    output:
        report = "results/qc_reports/{sample}_assembly_summary.txt",
        tsv_line = "results/qc_reports/{sample}_summary_line.tsv"
    params:
        mode = lambda wildcards: samples_df.loc[wildcards.sample, "spades_mode"],
        kmers = lambda wildcards: samples_df.loc[wildcards.sample, "spades_k"],
        filter_len = lambda wildcards: samples_df.loc[wildcards.sample, "filter_len"],
        host = lambda wildcards: samples_df.loc[wildcards.sample, "host_remove"],
        downsample = lambda wildcards: samples_df.loc[wildcards.sample, "downsample_to"]
    shell:
        """
        python scripts/generate_summary.py \
            {input.quast} \
            {input.busco} \
            {input.fastp} \
            {input.kraken} \
            {output.report} \
            {output.tsv_line} \
            {wildcards.sample} \
            {params.mode} \
            "{params.kmers}" \
            {params.filter_len} \
            {params.host} \
            {params.downsample}
        """

# Implements a persistent data-logging architecture. 
# Bioinformatics working directories (like results/) are ephemeral and often deleted 
# to save storage. This rule captures a permanent snapshot of the run's metadata, 
# configurations, and metrics in the project root, enabling longitudinal performance 
# tracking across different pipeline executions without risking data loss.
rule compile_run_history:
    input:
        txts = expand("results/qc_reports/{sample}_assembly_summary.txt", sample=SAMPLES),
        tsvs = expand("results/qc_reports/{sample}_summary_line.tsv", sample=SAMPLES)
    output:
        "run_history/run_latest.stamp"
    shell:
        """
        # Isolate the current execution metadata into a unique timestamped capsule
        mkdir -p run_history/capsules
        DATE_TAG=$(date +"%Y-%m-%d_%H-%M-%S")
        CAPSULE_DIR="run_history/capsules/run_${{DATE_TAG}}"
        mkdir -p $CAPSULE_DIR

        # Initialize the persistent tracking database if it doesn't exist
        MASTER_TSV="run_history/master_assembly_log.tsv"
        if [ ! -f "$MASTER_TSV" ]; then
            echo -e "Sample\tRun_Timestamp\tHost_Remove\tDownsample_To\tSpades_Mode\tKmers\tMin_Contig\tReads_Raw\tReads_Clean\tQ30_Pct\tInsert_Size\tUnclassified_Pct\tTop_Hit_1\tTop_Hit_2\tTotal_Length_bp\tContigs_>1k\tLargest_Contig\tN50\tL50\tGC_Pct\tBUSCO_C\tBUSCO_S\tBUSCO_D\tBUSCO_F\tBUSCO_M" > $MASTER_TSV
        fi

        # Append standardized metrics from all samples to the persistent database
        cat {input.tsvs} >> $MASTER_TSV

        # Snapshot the workflow configuration to guarantee full reproducibility of these metrics
        cp samples.tsv $CAPSULE_DIR/
        cp config.yaml $CAPSULE_DIR/
        cp {input.txts} $CAPSULE_DIR/
        
        touch {output}
        """

