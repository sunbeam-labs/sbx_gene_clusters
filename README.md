# sbx_gene_clusters

Reads-level based alignment to gene clusters of interest, e.g. bai operon or butyrate producing genes. Please refer to [sunbeam_database](https://github.com/zhaoc1/sunbeam_databases.git) for details. 

Take [**UniRef50** database](https://www.uniprot.org/downloads) as an example. First download the uniref50.fasta into your current `sunbeam_output/mapping/sbx_gene_family/databases/`.

 ```bash
 mkdir -p sunbeam_output/mapping/sbx_gene_family/database/
 wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz -P sunbeam_output/mapping/sbx_gene_family/database/
 ```
 Second, update the `config.yml` with the proper path.

## Usage

With your [sunbeam](https://github.com/sunbeam-labs/sunbeam) conda environemnt activated, 

1. Install the extension:

 ```bash
 sunbeam extend https://github.com/sunbeam-labs/sbx_gene_clusters
 ```

2. (Optional) update the config for an existing project (future projects will have these config options by default)

 ```bash
 sunbeam config update -i /path/to/config.yaml
 ```
 
3. Run time

 ```bash
 sunbeam run --configfile sunbeam_config.yml all_gene_family
 ```

## For sunbeam <3.0

 With your [sunbeam](https://github.com/sunbeam-labs/sunbeam) conda environemnt activated, 
 
 1. Clone into your Sunbeam extensions directory:
 
  ```bash
  cd $SUNBEAM_DIR
  git clone https://github.com/sunbeam-labs/sbx_gene_clusters
  ```
  
 2. Add the new config options to your config file
 
  ```bash
  cat sunbeam/extensions/sbx_gene_clusters/config.yml >> sunbeam_config.yml
  ```
 
 3. Install the requirements (optional, alternatively run `sunbeam` with the `--use-conda` and `--conda-prefix` flags as shown below):
 
  ```bash
  conda install --file extensions/sbx_gene_clusters/gene_clusters_env.yaml
  ```
  
 4. Run time

 - Use diamond
 
  ```bash
  sunbeam run -- --configfile sunbeam_config.yml --use-conda --conda-prefix $SUNBEAM_DIR/.snakemake all_gene_family
  ```
 
