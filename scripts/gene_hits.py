import csv
import os
from collections import Counter, OrderedDict


def read_alignment(reader, aln_chunk):
    # Read a chunk of alignment for a single query.
    for line in reader:
        if aln_chunk[0][0] != line[0]:  # found complete chunk
            return aln_chunk, [line]
        else:
            aln_chunk.append(line)  # grow the chunk
    else:
        return aln_chunk, []  # final chunk


def filter_evalue(alignment_chunk, e_cutoff):
    return [aln for aln in alignment_chunk if float(aln[12]) < e_cutoff]


def filter_mismatch(alignment_chunk, mismatch_cutoff):
    return [
        aln
        for aln in alignment_chunk
        if (int(aln[6]) - int(alignment_chunk[0][6])) < mismatch_cutoff
    ]


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


def write_gene_hits(in_fp, out_fp, db_annot_fp, evalue, alnLen, mismatch, log):
    ## Organize the gene information
    db_organized = OrderedDict()
    with open(db_annot_fp) as db_in:
        db = csv.DictReader(db_in, dialect="excel-tab")
        for row in db:
            log.write(f"{str(row)}\n")
            proteinID = row.get("proteinID")
            entry = db_organized.get(proteinID)
            if entry:
                db_organized.get(proteinID).append(
                    [row.get("geneID"), float(row.get("weight")), row.get("taxon")]
                )
            else:  # one protein can map to multiple classifications
                db_organized[proteinID] = [
                    [row.get("geneID"), float(row.get("weight")), row.get("taxon")]
                ]

    ## Read in the alignment file to count
    counter_genes = Counter()
    if os.path.getsize(in_fp) > 0:
        with open(in_fp) as in_file:
            reader = csv.reader(in_file, delimiter="\t")  # initialize the csv reader
            aln_chunk = [next(reader)]  # initialize the seed chunk

            while aln_chunk:
                current_chunk, aln_chunk = read_alignment(reader, aln_chunk)
                current_chunk = filter_chunk(current_chunk, evalue, alnLen, mismatch)
                current_chunk = get_best_hit(
                    current_chunk
                )  # is there a way to modify the current chunk and "pop" all but the first item?
                current_chunk = get_gene_function(current_chunk, db_organized)
                if current_chunk:
                    for chunk in current_chunk:
                        counter_genes[(chunk[14], chunk[16])] += chunk[15]

    ## Write the counts to file
    with open(out_fp, "w") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        writer.writerow(["geneID", "taxon", "count"])
        for key, value in counter_genes.items():
            writer.writerow(list(key) + [value])

    ## Remove the bulky .m8 file
    os.remove(in_fp)


with open(snakemake.log[0], "w") as log:
    write_gene_hits(
        snakemake.input.aln_fp,
        snakemake.output[0],
        snakemake.input.db_annot_fp[0],
        snakemake.params.evalue,
        snakemake.params.alnLen,
        snakemake.params.mismatch,
        log,
    )
