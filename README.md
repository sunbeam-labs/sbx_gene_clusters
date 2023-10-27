<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_gene_clusters

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_gene_clusters/actions/workflows/main.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_gene_clusters/actions/workflows/main.yml)
[![Super-Linter](https://github.com/sunbeam-labs/sbx_gene_clusters/actions/workflows/linter.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_gene_clusters/actions/workflows/linter.yml)
<!-- Badges end -->

Reads-level based alignment to gene clusters of interest, e.g. bai operon or butyrate producing genes. Please refer to [sunbeam_database](https://github.com/zhaoc1/sunbeam_databases.git) for details. 

## Installation

To install, activate your conda environment (using the name of your environment) and use `sunbeam extend`:

    conda activate <i>sunbeamX.X.X</i>
    sunbeam extend https://github.com/sunbeam-labs/sbx_gene_clusters.git

Now take [**UniRef50** database](https://www.uniprot.org/downloads) as an example. First download the uniref50.fasta into your current `sunbeam_output/mapping/sbx_gene_family/databases/`.

 ```bash
 mkdir -p sunbeam_output/mapping/sbx_gene_family/database/
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P sunbeam_output/mapping/sbx_gene_family/database/
 ```
Be sure to update the `config.yml` with the proper path.

## Usage
 
To generate alignments,

 ```bash
 sunbeam run --profile /path/to/project all_gene_clusters
 ```

## Configuration

 - threads: Is the number of threads to run parallel processes with
 - genes_fp: Is the path to the downloaded database
 - evalue: 
 - alnLen: 
 - mismatch: 

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_gene_clusters.git extensions/sbx_gene_clusters

and then include it in the config for any given project with:

    cat extensions/sbx_gene_clusters/config.yml >> /path/to/project/sunbeam_config.yml
 
