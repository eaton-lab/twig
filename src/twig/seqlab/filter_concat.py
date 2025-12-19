#!/usr/bin/env python

"""Concatenate a collection of fasta files into one large file (fa, phy, nex)

# concatenate genes from many files subsampled by header names
$ twig seq-filter -i AA/*.aa -c --subsample-genes GNAMES >

# concatenate alignments
$ twig seq-filter -i MSAs/*.fa -c -F phy > CONCAT.phy

# concatenate alignments with tips delim relabeled
$ twig concatenate -i A.fa B.fa C.fa -d "|" -di 0 -F phy > out.phy

# concatenate and write to nex
$ twig concatenate -i {A,B,C}.fa -F nex > out.phy

"""

import sys
from loguru import logger
from twig.seqlab.macse_refine import parse_fasta_to_dict
from twig.utils.logger_setup import set_log_level


def check_args(args):
    # if subsampling by gene name then not subsampling by IMAP
    if args.subsample_genes and (args.imap or args.delim):
        raise ValueError("--subsample-genes and (--imap or --delim) are mutually exclusive methods to subsample genes")

    # require sequences to be aligned (check a few) if concatenating
    if args.concatenate:
        pass


def stream(args):
    """Stream filtered loci to out"""

    # parse set of gene headers to subsample [this can be large >1M strs]
    keep_set = set()
    if args.subsample_genes:
        with args.subsample_genes.open() as hin:
            keep_set = set(i.strip() for i in hin)
    logger.warning(f"{list(keep_set)[:10]}")

    # iterate over the input fastas in order
    loci = []
    for fasta in args.input:
        seqs = parse_fasta_to_dict(fasta)
        logger.warning(f"{list(seqs.keys())}")

        # filter genes from fasta
        seqs = {i: j for (i, j) in seqs.items() if i in keep_set}

        # filter this locus?
        if len(seqs) >= args.min_samples:
            loci.append(seqs)

        # yield as fasta
        if len(loci) > 1000:
            fasta = "\n".join([f">{i}\n{j}" for (i, j) in seqs.items()] for seqs in loci)
            yield "\n".join(fasta).encode("utf-8")
            loci = []
    if loci:
        fasta = "\n".join([f">{i}\n{j}" for (i, j) in seqs.items()] for seqs in loci)
        yield "\n".join(fasta).encode("utf-8")


def run_filter_concat(args):
    """..."""
    set_log_level(args.log_level)
    check_args(args)
    logger.warning("RUN")

    # concatenate pipeline
    if args.concatenate:
        concatenate(args)
    else:
        if args.out:
            out = args.out.open('wb')
        else:
            out = sys.stdout
        for chunk in stream(args):
            sys.stdout.write(chunk)
        if args.out:
            out.close()



def main():
    from twig.cli.subcommands import get_parser_filter_concat
    parser = get_parser_filter_concat()
    args = parser.parse_args()
    run_filter_concat(args)


if __name__ == "__main__":
    pass
