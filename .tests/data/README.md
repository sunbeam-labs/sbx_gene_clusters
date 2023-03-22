To create `.tests/dbs/`:

```bash
    cat raw/*.fna > .tests/dbs/bacteria.fa
    makeblastdb -dbtype nucl -in .tests/dbs/bacteria.fa
    makeblastdb -dbtype prot -in .tests/dbs/prot.fa
```

where `raw/*.fna` includes genomes of bacteria you want to include in the nucl database and `prot.fa` is a file containing proteins, some of which should be identified by tests.