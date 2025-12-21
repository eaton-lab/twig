#!/usr/bin/env python

"""Parse genome files to write table w/ chromosome info.

This takes the genome and annotation files as input and will check
that the gene names match between the two. It can often be useful to
shorten gene names before running downstream analyses. See the tool
... for relabeling gene names in gff and sequence file headers which
you would want to run before this step.

Example
-------
$ python get_chromosome_table.py --genome A.fa --gff A.gff

Options
-------
- handle GTF format
- filter overlapping transcripts to keep longest?

Output
------
Chr01   100001  500
Chr02   90002   300
...

TODO
-----
Do we always parse on gene or are there cases where we should parse
on mRNA? Some annotations are marked by mRNAs but not as genes, as
a result of different annotation sources...
"""


import sys
import gzip
from pathlib import Path
from loguru import logger
import pandas as pd


def get_chrom_names_and_lengths(genome: Path) -> dict[str, int]:
    """Return dict mapping scaffold names to their sequence lengths"""
    seq_lengths = {}
    header = None
    current_length = 0
    if genome.suffix == ".gz":
        xopen = gzip.open
    else:
        xopen = open
    with xopen(genome, mode='rt') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seq_lengths[header] = current_length
                # store with ">" removed and reset counter
                header = line[1:]
                current_length = 0
            else:
                current_length += len(line)
        # capture last sequence length
        if header:
            seq_lengths[header] = current_length
    return seq_lengths


def get_genes_per_chrom(gff: Path) -> dict[str, int]:
    """Return a dict mapping scaffold names to their number of annotated genes"""
    ngenes = {}
    if gff.suffix == ".gz":
        xopen = gzip.open
    else:
        xopen = open
    with xopen(gff, mode='rt') as file:
        for line in file:
            line = line.split()
            if line[2].upper() == "MRNA":
                name = line[0]
                if name not in ngenes:
                    ngenes[name] = 1
                else:
                    ngenes[name] += 1
    return ngenes


def run_genome_table(args):
    """..."""
    # parse data
    clens = get_chrom_names_and_lengths(args.genome)
    if not args.annotation:
        ngenes = {i: 0 for i in clens}
    else:
        ngenes = get_genes_per_chrom(args.annotation)

        # validate scaffold names match, and set ngenes to 0 if no genes.
        bad_gene_names = set(ngenes) - set(clens)
        if bad_gene_names:
            raise ValueError(f"gene names found in annotation but not genome fasta:\n{bad_gene_names}")
        for name in clens:
            if not ngenes.get(name):
                ngenes[name] = 0

    # organize table
    table = pd.DataFrame({
        "name": list(clens),
        "length": [clens[i] for i in clens],
        "ngenes": [ngenes[i] for i in clens],
    })
    if args.sort_alpha:
        table.sort_values(by="name", inplace=True)
        logger.debug(f"sorting scaffolds by name. First={table.iloc[0].tolist()[:2]}")        
    if args.sort_len:
        table.sort_values(by="length", inplace=True, ascending=False)
        logger.debug(f"sorting scaffolds by length. Longest={table.iloc[0].tolist()[:2]}")
    if args.subset:
        table = table.iloc[:args.subset]
        logger.debug(f"subselecting scaffolds in first N ({args.subset}) while retaining current order")
    if args.subset_idx:
        # entered by user as 1-indexed
        subset_idx = [i - 1 for i in table.index if i in args.subset_idx]
        table = table.iloc[subset_idx]
        logger.debug("subselecting scaffolds in subset indices while retaining current order")        
    if args.subset_names:
        table = table.loc[[i in args.subset_names for i in table.name]]
        logger.debug("subselecting scaffolds in subset names while retaining current order")                
    table.to_csv(sys.stdout, sep="\t", index=False)


def main():
    parser = get_parser_genome_table()
    args = parser.parse_args()
    run_genome_table(args)


def test():
    """Example run"""


if __name__ == "__main__":
    main()
