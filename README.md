# sbx_gene_family

This extension can be used to perform functional mapping, i.e. mapping metagenomic reads to proteins. The database to map against could be [UniRef50](https://www.uniprot.org/downloads), all prokaryotic proteins from [KEGG](https://www.kegg.jp/kegg/download/), or more targeted [databases](https://github.com/zhaoc1/sunbeam_databases.git), e.g. bai operon or butyrate producing genes.

Take [**UniRef50** database](https://www.uniprot.org/downloads) as an example. First download the uniref50.fasta into your current `sunbeam_output/mapping/sbx_gene_family/databases/`.

 ```bash
 mkdir -p sunbeam_output/mapping/sbx_gene_family/database/
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P sunbeam_output/mapping/sbx_gene_family/database/
 ```
 Second, update the `config.yml` with the proper path.

## Usage

 With you [sunbeam](https://github.com/sunbeam-labs/sunbeam) conda environemnt activated, 
 
 1. Clone into your Sunbeam extensions directory:
 
  ```bash
  git clone https://github.com/sunbeam-labs/sbx_gene_clusters
  ```
  
 2. Add the new config options to your config file
 
  ```bash
  cat sunbeam/extensions/sbx_gene_clusters/config.yml >> sunbeam_config.yml
  ```
 
 3. Install the requirements:
 
  ```bash
  conda install --file extensions/sbx_gene_clusters/requirements.txt
  ```
  
 4. Run time

By default, mapping uses DIAMOND, but this extension also supports using BLAST.
 
  ```bash
  sunbeam run -- --configfile sunbeam_config.yml all_gene_family
  ```
 
