# ==============================================================================
# assembly.smk 
# contains the rules necessary for QC and assembly. gets called by snakefile
# ==============================================================================


# downloads busco lineage db for lineages specified in samples.tsv
rule download_busco_lineage:
    output:
        lineage_dir = directory("data/dbs/busco/lineages/{lineage}")
    log:
        "results/logs/dbs/busco_download_{lineage}.log"
    shell:
        """
        mkdir -p data/dbs/busco/information data/dbs/busco/lineages
        
        # Download version manifest if absent
        if [ ! -f data/dbs/busco/information/file_versions.tsv ]; then
            wget -q -O data/dbs/busco/information/file_versions.tsv "https://busco-data.ezlab.org/v5/data/file_versions.tsv" || \
            wget -q -O data/dbs/busco/information/file_versions.tsv "https://busco-data2.ezlab.org/v5/data/file_versions.tsv"
        fi

        # Stream download and uncompress directly into lineages directory
        wget -qO- "https://busco-data.ezlab.org/v5/data/lineages/{wildcards.lineage}.2024-01-08.tar.gz" 2> {log} | \
        tar -xz -C data/dbs/busco/lineages/
        """


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
    container:
        config["containers"]["fastp"]
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


# DECONTAMINATION ROUTING
# Evaluates 'decontaminate_human' column in samples.tsv to determine input path.
def get_reads_for_downsampling(wildcards):
    """
    Dynamically routes clean FASTQs depending on whether host decontamination
    was requested for the given sample.
    """
    host_status = str(samples_df.loc[wildcards.sample, "decontaminate_human"]).strip().lower()
    
    if host_status in ["human", "yes", "true"]:
        # Hostile appends '.clean_1' and '.clean_2' to the original input basenames
        return {
            "r1": f"results/01b_decontam/{wildcards.sample}_1.clean_1.fastq.gz",
            "r2": f"results/01b_decontam/{wildcards.sample}_2.clean_2.fastq.gz"
        }
    else:
        # Bypass decontamination and consume fastp clean reads directly
        return {
            "r1": f"results/01_clean_data/{wildcards.sample}_1.fastq.gz",
            "r2": f"results/01_clean_data/{wildcards.sample}_2.fastq.gz"
        }


# HOSTILE DECONTAMINATION
# Removes human host reads using masked T2T-CHM13v2.0 index and Bowtie2.
rule hostile_decontam:
    input:
        r1 = "results/01_clean_data/{sample}_1.fastq.gz",
        r2 = "results/01_clean_data/{sample}_2.fastq.gz"
    output:
        r1 = "results/01b_decontam/{sample}_1.clean_1.fastq.gz",
        r2 = "results/01b_decontam/{sample}_2.clean_2.fastq.gz",
        json_log = "results/qc_reports/hostile/{sample}_hostile.json"
    container:
        config["containers"]["hostile"]
    params:
        cache_dir = "data/dbs/hostile",
        index = "human-t2t-hla",
        out_dir = "results/01b_decontam"
    threads: 8
    log:
        "results/logs/01b_hostile_{sample}.log"
    shell:
        """
        # Override default user home cache directory to project workspace
        export HOSTILE_CACHE_DIR="{params.cache_dir}"

        hostile clean \
            --fastq1 {input.r1} \
            --fastq2 {input.r2} \
            --aligner bowtie2 \
            --index {params.index} \
            --output {params.out_dir} \
            --threads {threads} \
            --airplane \
            --force \
            > {output.json_log} 2> {log}
        """


# DOWNSAMPLING
# Normalizes sequencing depth across samples prior to de novo assembly.
rule downsample_reads:
    input:
        unpack(get_reads_for_downsampling)
    output:
        r1 = "results/01c_downsample/{sample}_1.fastq.gz",
        r2 = "results/01c_downsample/{sample}_2.fastq.gz"
    container:
        config["containers"]["seqtk"]
    params:
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
            # Sample reads randomly with fixed seed (42) for reproducibility
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
    container:
        config["containers"]["spades"]
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
    container:
        config["containers"]["seqtk"]
    params:
        min_len = lambda wildcards: samples_df.loc[wildcards.sample, "filter_len"]
    log:
        "results/logs/02b_filter_{sample}.log"
    shell:
        """
        seqtk seq -L {params.min_len} {input.scaffolds} > {output.filtered} 2> {log}
        """


# POLISHING
# Refines the assembly by correcting base errors and small indels using the original short reads.
rule pypolca:
    input:
        assembly = "results/02_assembly/{sample}/scaffolds_filtered.fasta",
        r1 = "results/01c_downsample/{sample}_1.fastq.gz",
        r2 = "results/01c_downsample/{sample}_2.fastq.gz"
    output:
        corrected = "results/02_assembly/{sample}/pypolca/{sample}_corrected.fasta",
        report = "results/02_assembly/{sample}/pypolca/{sample}.report"
    container:
        config["containers"]["pypolca"]
    params:
        outdir = "results/02_assembly/{sample}/pypolca",
        prefix = "{sample}"
    threads: 8
    log:
        "results/logs/02c_pypolca_{sample}.log"
    shell:
        """
        pypolca run \
            -a {input.assembly} \
            -1 {input.r1} \
            -2 {input.r2} \
            -t {threads} \
            -o {params.outdir} \
            -p {params.prefix} \
            --careful \
            -f \
            > {log} 2>&1
        """


# Chamada so se a amostra tiver algo diferente de "none" na coluna "reference_genome" em samples.tsv
# RagTag usa um genoma de referência pra montar os scaffolds
rule ragtag_scaffold:
    input:
        assembly = "results/02_assembly/{sample}/pypolca/{sample}_corrected.fasta", # << AQUI
        reference = lambda wildcards: samples_df.loc[wildcards.sample, "reference_genome"]
    output:
        scaffolded = "results/02_assembly/{sample}/ragtag_output/ragtag.scaffold.fasta"
    container:
        config["containers"]["ragtag"]
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
            else "results/02_assembly/{sample}/pypolca/{sample}_corrected.fasta" # << AQUI
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
    container:
        config["containers"]["kraken2"]
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
    container:
        config["containers"]["quast"]
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
        assembly = "results/02_assembly/{sample}/final_assembly.fasta",
        lineage_db = lambda wildcards: f"data/dbs/busco/lineages/{samples_df.loc[wildcards.sample, 'busco_lineage']}"
    output:
        summary = "results/qc_reports/busco/{sample}/short_summary.specific.{lineage}.{sample}.txt"
    container:
        config["containers"]["busco"]
    params:
        out_path = "results/qc_reports/busco",
        lineage = lambda wildcards: samples_df.loc[wildcards.sample, "busco_lineage"]
    threads: 8
    log:
        "results/logs/03_busco_{sample}_{lineage}.log"
    shell:
        """
        # Force Python unbuffered logging and limit OpenMP threads to prevent deadlocks
        export PYTHONUNBUFFERED=1
        export OMP_NUM_THREADS={threads}

        busco -i {input.assembly} \
              -l {input.lineage_db} \
              -o {wildcards.sample} \
              --out_path {params.out_path} \
              --download_path data/dbs/busco \
              --offline \
              --opt-out-run-stats \
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
    container:
        config["containers"]["multiqc"]
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
        kraken = "results/qc_reports/kraken2/{sample}_kraken2_report.txt",
        pypolca = "results/02_assembly/{sample}/pypolca/{sample}.report" # NOVO INPUT
    output:
        report = "results/qc_reports/{sample}_assembly_summary.txt",
        tsv_line = "results/qc_reports/{sample}_summary_line.tsv"
    params:
        mode = lambda wildcards: samples_df.loc[wildcards.sample, "spades_mode"],
        kmers = lambda wildcards: samples_df.loc[wildcards.sample, "spades_k"],
        filter_len = lambda wildcards: samples_df.loc[wildcards.sample, "filter_len"],
        host = lambda wildcards: samples_df.loc[wildcards.sample, "decontaminate_human"],
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
            {params.downsample} \
            {input.pypolca} 
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
        stamp = "run_history/run_latest.stamp",
        dashboard = "results/full_report.html"
    shell:
        """
        mkdir -p run_history/capsules
        DATE_TAG=$(date +"%Y-%m-%d_%H-%M-%S")
        CAPSULE_DIR="run_history/capsules/run_${{DATE_TAG}}"
        mkdir -p $CAPSULE_DIR

        MASTER_TSV="run_history/master_assembly_log.tsv"
        if [ ! -f "$MASTER_TSV" ]; then
            echo -e "Sample\tRun_Timestamp\tHost_Remove\tHost_Removed_Reads\tDownsample_To\tSpades_Mode\tKmers\tMin_Contig\tReads_Raw\tReads_Clean\tQ30_Pct\tInsert_Size\tPolca_SNPs\tPolca_Indels\tUnclassified_Pct\tTop_Hit_1\tTop_Hit_2\tTotal_Length_bp\tContigs_>1k\tLargest_Contig\tN50\tL50\tGC_Pct\tBUSCO_C\tBUSCO_S\tBUSCO_D\tBUSCO_F\tBUSCO_M" > $MASTER_TSV
        fi

        cat {input.tsvs} >> $MASTER_TSV

        cp samples.tsv $CAPSULE_DIR/
        cp config.yaml $CAPSULE_DIR/
        cp {input.txts} $CAPSULE_DIR/
        
        # Gera o Dashboard HTML interativo
        python scripts/build_html_report.py $MASTER_TSV {output.dashboard} {input.txts}
        
        # Salva uma cópia do dashboard na cápsula do tempo
        cp {output.dashboard} $CAPSULE_DIR/

        touch {output.stamp}
        """

