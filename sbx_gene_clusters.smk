from pathlib import PurePath


try:
    SBX_GENE_CLUSTERS_VERSION = get_ext_version("sbx_gene_clusters")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_GENE_CLUSTERS_VERSION = "0.0.0"

GENES_DIR = Cfg["sbx_gene_clusters"]["genes_fp"]
GENES_KEY = [PurePath(f.name).stem for f in GENES_DIR.glob("*.fasta")]
GENES_VAL = [str(GENES_DIR) + "/" + g + ".fasta" for g in GENES_KEY]
GENES_DICT = dict(zip(GENES_KEY, GENES_VAL))
print(f"sbx_gene_clusters::INFO Found these genes dbs: {str(GENES_DICT)}")


rule all_gene_clusters:
    input:
        expand(
            MAPPING_FP / "sbx_gene_clusters" / "{gene}" / "{sample}_1.txt",
            gene=GENES_DICT.keys(),
            sample=Samples.keys(),
        ),


rule merge_pairs:
    input:
        r1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        r2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        reads=MAPPING_FP / "merged" / "{sample}.fastq",
    benchmark:
        BENCHMARK_FP / "merge_pairs_{sample}.tsv"
    log:
        LOG_FP / "merge_pairs_{sample}.log",
    threads: Cfg["sbx_gene_clusters"]["threads"]
    conda:
        "envs/sbx_gene_clusters_env.yml"
    shell:
        """
        vsearch \
        --fastq_mergepairs {input.r1} \
        --reverse {input.r2} \
        --fastqout {output.reads} \
        --threads {threads} \
        --fastq_allowmergestagger \
        --fastq_maxdiffs 5 \
        --fastq_minovlen 10 \
        --fastq_minmergelen 100 \
        > {log} 2>&1
        """


rule build_gene_clusters_diamond_db:
    input:
        lambda wildcards: GENES_DICT[wildcards.gene],
    output:
        expand(GENES_DIR / "{{gene}}.fasta.{index}", index=["dmnd"]),
    benchmark:
        BENCHMARK_FP / "build_gene_clusters_diamond_db_{gene}.tsv"
    log:
        LOG_FP / "build_gene_clusters_diamond_db_{gene}.log",
    conda:
        "envs/sbx_gene_clusters_env.yml"
    shell:
        """
        diamond makedb --in {input} -d {input} > {log} 2>&1
        """


rule build_gene_clusters_blast_db:
    input:
        lambda wildcards: GENES_DICT[wildcards.gene],
    output:
        expand(GENES_DIR / "{{gene}}.fasta.{index}", index=["psq", "pin", "phr"]),
    benchmark:
        BENCHMARK_FP / "build_gene_clusters_blast_db_{gene}.tsv"
    log:
        LOG_FP / "build_gene_clusters_blast_db_{gene}.log",
    conda:
        "envs/sbx_gene_clusters_env.yml"
    shell:
        """
        makeblastdb -in {input} -dbtype prot > {log} 2>&1
        """


rule fq_2_fa:
    input:
        QC_FP / "decontam" / "{sample}_1.fastq.gz",
    output:
        MAPPING_FP / "R1" / "{sample}_1.fasta",
    benchmark:
        BENCHMARK_FP / "fq_2_fa_{sample}.tsv"
    log:
        LOG_FP / "fq_2_fa_{sample}.log",
    conda:
        "envs/sbx_gene_clusters_env.yml"
    shell:
        """
        (seqtk seq -a < <(gzip -cd {input}) > {output}) > {log} 2>&1
        """


rule gene_hits:
    input:
        read=MAPPING_FP / "R1" / "{sample}_1.fasta",
        db=expand(GENES_DIR / "{{gene}}.fasta.{index}", index=["dmnd"]),
        db_annot_fp=expand(GENES_DIR / "{{gene}}.{index}", index=["txt"]),
    output:
        MAPPING_FP / "sbx_gene_clusters" / "{gene}" / "{sample}_1.txt",
    benchmark:
        BENCHMARK_FP / "gene_hits_{gene}_{sample}.tsv"
    log:
        diamond_log=LOG_FP / "gene_hits_diamond_{gene}_{sample}.log",
        script_log=LOG_FP / "gene_hits_{gene}_{sample}.log",
    params:
        evalue=float(Cfg["sbx_gene_clusters"]["evalue"]),
        alnLen=Cfg["sbx_gene_clusters"]["alnLen"],
        mismatch=Cfg["sbx_gene_clusters"]["mismatch"],
        m8=str(MAPPING_FP / "sbx_gene_clusters" / "{gene}" / "{sample}_1.m8"),
    threads: Cfg["sbx_gene_clusters"]["threads"]
    conda:
        "envs/sbx_gene_clusters_env.yml"
    script:
        "scripts/gene_hits.py"


rule blastx_reads:
    input:
        read=MAPPING_FP / "R1" / "{sample}_1.fasta",
        db=expand(GENES_DIR / "{{gene}}.fasta.{index}", index=["psq", "pin", "phr"]),
    output:
        MAPPING_FP / "sbx_gene_clusters" / "{gene}" / "{sample}_1.blastx",
    benchmark:
        BENCHMARK_FP / "blastx_reads_{gene}_{sample}.tsv"
    log:
        LOG_FP / "blastx_reads_{gene}_{sample}.log",
    params:
        db=lambda wildcard: GENES_DICT[wildcard.gene],
    threads: Cfg["sbx_gene_clusters"]["threads"]
    conda:
        "envs/sbx_gene_clusters_env.yml"
    shell:
        """
        blastx -query {input.read} -db {params.db} \
               -num_threads {threads} -evalue 1e-6 \
               -max_target_seqs 5000 \
               -out {output} \
               -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore" \
                > {log} 2>&1
        """


rule uniref50_download:
    output:
        MAPPING_FP / "sbx_gene_clusters" / "databases" / "uniref50.fasta",
    shell:
        """
        set +o pipefail
        mkdir -p $(dirname {output})
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P $(dirname {output})
        """
