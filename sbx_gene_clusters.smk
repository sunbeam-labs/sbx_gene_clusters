# -*- mode: Snakemake -*-
#
# Chunyu Zhao 2018-07-23
# Reads mapping to gene clusters/pathways/gene families of interest:
#   Rules for Diamond or BLASTx reads against protein databases

try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")

GENES_DIR = Cfg["sbx_gene_clusters"]["genes_fp"]
GENES_KEY = [PurePath(f.name).stem for f in GENES_DIR.glob("*.fasta")]
GENES_VAL = [str(GENES_DIR) + "/" + g + ".fasta" for g in GENES_KEY]
GENES_DICT = dict(zip(GENES_KEY, GENES_VAL))
print(f"sbx_gene_clusters::INFO Found these genes dbs: {str(GENES_DICT)}")

TARGET_GENES = expand(
    str(MAPPING_FP / "sbx_gene_family" / "{gene}" / "{sample}_1.txt"),
    gene=GENES_DICT.keys(),
    sample=Samples.keys(),
)


rule all_gene_family:
    input:
        TARGET_GENES,


rule merge_pairs:
    input:
        r1=str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
        r2=str(QC_FP / "decontam" / "{sample}_2.fastq.gz"),
    output:
        r1=str(MAPPING_FP / "merged" / "{sample}.fastq"),
    threads: Cfg["sbx_gene_clusters"]["threads"]
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        vsearch \
        --fastq_mergepairs {input.r1} --reverse {input.r2} \
        --fastqout {output.reads} --threads {threads} \
        --fastq_allowmergestagger --fastq_maxdiffs 5 \
        --fastq_minovlen 10 --fastq_minmergelen 100
        """


rule build_gene_clusters_diamond_db:
    input:
        lambda wildcards: GENES_DICT[wildcards.gene],
    output:
        expand(str(GENES_DIR / "{{gene}}.fasta.{index}"), index=["dmnd"]),
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        diamond makedb --in {input} -d {input} 
        """


rule build_blast_db:
    input:
        lambda wildcards: GENES_DICT[wildcards.gene],
    output:
        expand(str(GENES_DIR / "{{gene}}.fasta.{index}"), index=["psq", "pin", "phr"]),
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        makeblastdb -in {input} -dbtype prot
        """


rule fq_2_fa:
    input:
        str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
    output:
        str(MAPPING_FP / "R1" / "{sample}_1.fasta"),
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        seqtk seq -a < <(gzip -cd {input}) > {output}
        """


rule diamond_reads:
    input:
        read=str(MAPPING_FP / "R1" / "{sample}_1.fasta"),
        db=expand(str(GENES_DIR / "{{gene}}.fasta.{index}"), index=["dmnd"]),
    output:
        str(MAPPING_FP / "sbx_gene_family" / "{gene}" / "{sample}_1.m8"),
    threads: Cfg["sbx_gene_clusters"]["threads"]
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        diamond blastx \
            --db {input.db} --query {input.read} \
            --threads {threads} --evalue 1e-6 \
            --max-target-seqs 0 \
            --out {output} \
            --outfmt 6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore
        """


rule gene_hits:
    input:
        aln_fp=str(MAPPING_FP / "sbx_gene_family" / "{gene}" / "{sample}_1.m8"),
        db_annot_fp=expand(str(GENES_DIR / "{{gene}}.{index}"), index=["tsv"]),
    output:
        str(MAPPING_FP / "sbx_gene_family" / "{gene}" / "{sample}_1.txt"),
    benchmark:
        BENCHMARK_FP / "gene_hits_{gene}_{sample}.tsv"
    log:
        LOG_FP / "gene_hits_{gene}_{sample}.log",
    params:
        evalue=float(Cfg["sbx_gene_clusters"]["evalue"]),
        alnLen=Cfg["sbx_gene_clusters"]["alnLen"],
        mismatch=Cfg["sbx_gene_clusters"]["mismatch"],
    script:
        "scripts/gene_hits.py"


rule blastx_reads:
    input:
        read=str(MAPPING_FP / "R1" / "{sample}_1.fasta"),
        db=expand(
            str(GENES_DIR / "{{gene}}.fasta.{index}"), index=["psq", "pin", "phr"]
        ),
    output:
        str(MAPPING_FP / "sbx_gene_family" / "{gene}" / "{sample}_1.blastx"),
    params:
        db=lambda wildcard: GENES_DICT[wildcard.gene],
    threads: Cfg["sbx_gene_clusters"]["threads"]
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        blastx -query {input.read} -db {params.db} \
               -num_threads {threads} -evalue 1e-6 \
               -max_target_seqs 5000 \
               -out {output} \
               -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore"
        """


rule uniref50_download:
    output:
        str(MAPPING_FP / "sbx_gene_family" / "databases" / "uniref50.fasta"),
    params:
        str(MAPPING_FP / "sbx_gene_family" / "databases"),
    shell:
        """
        set +o pipefail
        mkdir -p {params}
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P {params}
        """
