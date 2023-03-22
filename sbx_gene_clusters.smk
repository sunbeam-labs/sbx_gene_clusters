# -*- mode: Snakemake -*-
#
# Chunyu Zhao 2018-07-23
# Reads mapping to gene clusters/pathways/gene families of interest:
#   Rules for Diamond or BLASTx reads against protein databases

import csv
import os
from collections import Counter, OrderedDict

GENES_DIR = Cfg['sbx_gene_clusters']['genes_fp']
GENES_KEY = [PurePath(f.name).stem for f in GENES_DIR.glob('*.fasta')]
GENES_VAL = [str(GENES_DIR) + '/' + g+'.fasta' for g in GENES_KEY]
GENES_DICT = dict(zip(GENES_KEY, GENES_VAL))

TARGET_GENES = expand(str(MAPPING_FP/'sbx_gene_family'/'{gene}'/'{sample}_1.txt'), 
                     gene=GENES_DICT.keys(), sample=Samples.keys())

rule all_gene_family:
    input:
        TARGET_GENES


def read_alignment(reader, aln_chunk):
    #Read a chunk of alignment for a single query.
    for line in reader:
        if aln_chunk[0][0] != line[0]: #found complete chunk
            return aln_chunk, [line]
        else:
            aln_chunk.append(line) #grow the chunk
    else:
        return aln_chunk, [] #final chunk

def filter_evalue(alignment_chunk, e_cutoff):
    return [aln for aln in alignment_chunk if float(aln[12]) < e_cutoff]

def filter_mismatch(alignment_chunk, mismatch_cutoff):
    return [aln for aln in alignment_chunk if ( int(aln[6]) - int(alignment_chunk[0][6]) ) < mismatch_cutoff]

def filter_aln_length(alignment_chunk, aln_len_cutoff):
    return [aln for aln in alignment_chunk if int(aln[5]) > aln_len_cutoff]

def filter_chunk(alignment_chunk, e_cutoff, aln_len_cutoff, mismatch_cutoff):
    alignment_chunk = filter_evalue(alignment_chunk, e_cutoff)
    alignment_chunk = filter_mismatch(alignment_chunk, mismatch_cutoff)
    alignment_chunk = filter_aln_length(alignment_chunk, aln_len_cutoff)
    return alignment_chunk

def get_best_hit(alignment_chunk):
    if alignment_chunk:
        return alignment_chunk.pop(0)

def get_gene_function(aln, db):
    if aln:
        entries = db.get(aln[1])
        if entries:
            return [aln + entry for entry in entries]

def write_gene_hits(in_fp, out_fp, db_annot_fp, evalue, alnLen, mismatch):

    ## Organize the gene information    
    db_organized = OrderedDict()
    with open(db_annot_fp) as db_in:
        db = csv.DictReader(db_in, delimiter='\t')
        for row in db:
            proteinID = row.get('proteinID')
            entry = db_organized.get(proteinID)
            if entry:
                db_organized.get(proteinID).append([row.get('geneID'), float(row.get('weight')), row.get('taxon')])
            else: # one protein can map to multiple classifications
                db_organized[proteinID] = [[row.get('geneID'), float(row.get('weight')), row.get('taxon')]]
    
    ## Read in the alignment file to count
    counter_genes = Counter()
    if os.path.getsize(in_fp) > 0:    
        with open(in_fp) as in_file:
            reader = csv.reader(in_file, delimiter='\t') # initialize the csv reader
            aln_chunk = [next(reader)] #initialize the seed chunk
        
            while aln_chunk:
                current_chunk, aln_chunk = read_alignment(reader, aln_chunk)
                current_chunk = filter_chunk(current_chunk, evalue, alnLen, mismatch)
                current_chunk = get_best_hit(current_chunk) # is there a way to modify the current chunk and "pop" all but the first item?
                current_chunk = get_gene_function(current_chunk, db_organized)
                if current_chunk:
                    for chunk in current_chunk:    
                        counter_genes[(chunk[14], chunk[16])] += chunk[15]
    
    ## Write the counts to file
    with open(out_fp, 'w') as out_file:
        writer=csv.writer(out_file, delimiter='\t')
        writer.writerow(["geneID", "taxon", "count"])
        for key, value in counter_genes.items():
            writer.writerow(list(key) + [value])

rule gene_hits:
    input:
        aln_fp = str(MAPPING_FP/'sbx_gene_family'/'{gene}'/'{sample}_1.m8'),
        db_annot_fp = expand(str(GENES_DIR/'{{gene}}.{index}'), index=['txt'])
    output:
        str(MAPPING_FP/'sbx_gene_family'/'{gene}'/'{sample}_1.txt')
    params:
        evalue = float(Cfg['sbx_gene_clusters']['evalue']),
        alnLen = Cfg['sbx_gene_clusters']['alnLen'],
        mismatch = Cfg['sbx_gene_clusters']['mismatch']
    run:
        write_gene_hits(input.aln_fp, output[0], input.db_annot_fp[0], params.evalue, params.alnLen, params.mismatch)

rule merge_pairs:
    input:
        r1 = str(QC_FP/'decontam'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'decontam'/'{sample}_2.fastq.gz')
    output:
        r1 = str(MAPPING_FP/'merged'/'{sample}.fastq')
    threads:
        Cfg['sbx_gene_clusters']['threads']
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
        lambda wildcards: GENES_DICT[wildcards.gene]
    output:
        expand(str(GENES_DIR/'{{gene}}.fasta.{index}'),index=['dmnd'])
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        diamond makedb --in {input} -d {input} 
        """

rule build_blast_db:
    input:
        lambda wildcards: GENES_DICT[wildcards.gene]
    output:
        expand(str(GENES_DIR/'{{gene}}.fasta.{index}'),index=['psq','pin','phr'])
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        makeblastdb -in {input} -dbtype prot
        """

rule fq_2_fa:
    input:
        str(QC_FP/'decontam'/'{sample}_1.fastq.gz')
    output:
        str(MAPPING_FP/'R1'/'{sample}_1.fasta')
    conda:
        "sbx_gene_clusters_env.yml"
    shell:
        """
        seqtk seq -a < <(gzip -cd {input}) > {output}
        """

rule diamond_reads:
    input:
        read = str(MAPPING_FP/'R1'/'{sample}_1.fasta'),
        db = expand(str(GENES_DIR/'{{gene}}.fasta.{index}'), index=['dmnd'])
    output:
        str(MAPPING_FP/'sbx_gene_family'/'{gene}'/'{sample}_1.m8')
    threads:
        Cfg['sbx_gene_clusters']['threads']
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

rule blastx_reads:
    input:
        read = str(MAPPING_FP/'R1'/'{sample}_1.fasta'),
        db = expand(str(GENES_DIR/'{{gene}}.fasta.{index}'), index=['psq','pin','phr'])
    output:
        str(MAPPING_FP/'sbx_gene_family'/'{gene}'/'{sample}_1.blastx')
    params:
        db=lambda wildcard: GENES_DICT[wildcard.gene]
    threads:
        Cfg['sbx_gene_clusters']['threads']
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
        str(MAPPING_FP/'sbx_gene_family'/'databases'/'uniref50.fasta')
    params:
        str(MAPPING_FP/'sbx_gene_family'/'databases')
    shell:
        """
        set +o pipefail
        mkdir -p {params}
        wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P {params}
        """
