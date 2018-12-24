# sbx_gene_families

Reads-level based alignment to gene clusters of interest, e.g. bai operon or butyrate producing genes. Please refer to [sunbeam_database](https://github.com/zhaoc1/sunbeam_databases.git) for details.

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

 - Use diamond
 
  ```bash
  sunbeam run -- --configfile sunbeam_config.yml _all_genes
  ```
 
 - Use blast
 
  ```bash
  sunbeam run -- --configfile sunbeam_config.yml _all_genes_blast
  ```
