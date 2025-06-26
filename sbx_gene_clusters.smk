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
        LOG_FP / "gene_hits_{gene}_{sample}.log",
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
