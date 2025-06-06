<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_gene_clusters

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_gene_clusters/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_gene_clusters/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_gene_clusters)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_gene_clusters/)
<!-- Badges end -->

Reads-level based alignment to gene clusters of interest, e.g. bai operon or butyrate producing genes. Please refer to [sunbeam_database](https://github.com/zhaoc1/sunbeam_databases.git) for details. Make a diamond database from your proteins of interest fasta file and provide a text annotation file with the following columns: geneID, proteinID, ARO, taxon, weight.

## Configuration

 - threads: Is the number of threads to run parallel processes with
 - genes_fp: Is the path to the downloaded database
 - evalue: 
 - alnLen: 
 - mismatch: 

Take [**UniRef50** database](https://www.uniprot.org/downloads) as an example. Download it and point `sunbeam_config.yml` to it:

 ```bash
 mkdir -p /path/to/uniref50/
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P /path/to/uniref50/
 ```

## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).
