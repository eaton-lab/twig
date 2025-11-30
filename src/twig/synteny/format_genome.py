#!/usr/bin/env python

# DEPRECATED FOR format_fasta.py
"""Format genome fasta file.

This takes a genome fasta file and returns a modified copy.

Example
-------
$ python format_genome.py GENOME ...
$ twig format-genome GENOME ...

Output
------
Writes to stdout by default unless a -o option is provided.
"""

import sys
import gzip
import textwrap
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
from loguru import logger


def get_parser_format_genome(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="format-genome",
        usage="%(prog)s fasta [options]",
        help="format a genome fasta to sort, subselect, or relabel headers",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=textwrap.dedent("""
            -------------------------------------------------------------------
            | format-genome: sort/subselect/relabel headers & change masking  |
            -------------------------------------------------------------------
            | This tool is used to modify a genome fasta file by (1) sorting; |
            | (2) subselecting; and/or (3) relabeling headers; and (4)        |
            | optionally modifying sequences, with operations performed in    |
            | that order. The modified genome is written to STDOUT. A subset  |
            | of scaffolds can be returned in any order using --subset-idx or |
            | --subset-names, but note that --sort-alpha or --sort-len will   |
            | override the ordering.                                          |
            -------------------------------------------------------------------
        """),
        epilog=textwrap.dedent("""
            Examples
            --------
            $ format-genome GENOME --subset 3 --relabel-list A B C > genome.fa
            $ format-genome GENOME --subset-idx 0 1 2 3 --relabel-prefix-idx A_ > subgenomeA.fa
            $ format-genome GENOME --subset-names Chr1 Chr2 --relabel-list chromosome_1 chromosome_2 > genome.fa
            $ format-genome GENOME --sort-len --subset 8 --max-width 80 > chromosomes.fa
        """)
    )

    # create parser or connect as subparser to cli parser
    if parser:
        kwargs['name'] = kwargs.pop("prog")
        parser = parser.add_parser(**kwargs)
    else:
        kwargs.pop("help")
        parser = ArgumentParser(**kwargs)

    # add arguments
    parser.add_argument("genome", type=Path, help="a fasta genome sequence (can be .gz)")
    sort_options = parser.add_mutually_exclusive_group()
    sort_options.add_argument("--sort-alpha", action="store_true", help="sort chromosomes by name alphanumerically")
    sort_options.add_argument("--sort-len", action="store_true", help="sort chromosomes by length (longest to shortest)")

    relabel_options = parser.add_mutually_exclusive_group()
    relabel_options.add_argument("--relabel-prefix", type=str, metavar="str", default="", help="append prefix label to start of existing scaffold labels (default='')")
    relabel_options.add_argument("--relabel-prefix-idx", type=str, metavar="str", default="", help="overwrite labels with {{prefix}}{{counter}} on scaffolds in order (default='')")
    relabel_options.add_argument("--relabel-list", type=str, metavar="str", nargs="+", help="overwrite labels with new labels provided as a list")

    sub_options = parser.add_mutually_exclusive_group()
    sub_options.add_argument("--subset", type=int, metavar="int", default=None, help="subselect the first N scaffolds after optional sorting")
    sub_options.add_argument("--subset-idx", type=int, metavar="int", nargs="+", help="subselect one or more scaffolds by 1-based index after optional sorting")
    sub_options.add_argument("--subset-names", type=str, metavar="str", nargs="+", help="subselect one or more scaffolds by name after optional relabeling")

    parser.add_argument("--max-width", type=int, metavar="int", default=None, help="write sequences with max width textwrapping (e.g., 80)")
    parser.add_argument("--unmask", action="store_true", help="convert any soft-masked sequences to upper case")
    return parser


def get_headers_to_seqs_dict(genome: Path) -> dict[str, str]:
    """Return dict mapping scaffold names to their sequences"""
    # keep track of sequence lengths
    headers_to_seqs = {}

    # select opener func
    if genome.suffix == ".gz":
        xopen = gzip.open
    else:
        xopen = open

    # read through genome file
    with xopen(genome, mode='rt') as file:
        for line in file:
            line = line.strip()
            # skip empty lines
            if not line:
                continue
            # start new scaffold
            if line.startswith(">"):
                header = line[1:]
                headers_to_seqs[header] = []
            # append to existing scaffold
            else:
                headers_to_seqs[header].append(line)
    return headers_to_seqs


def sort_headers_to_seqs_dict(seqs: dict[str, str], sort_len: bool, sort_alpha: bool) -> dict[str, str]:
    """..."""
    if sort_len:
        return {i: seqs[i] for i in sorted(seqs, key=lambda x: len("".join(seqs[x])), reverse=True)}
    if sort_alpha:
        return {i: seqs[i] for i in sorted(seqs)}
    return seqs


def subselect(seqs: dict[str, str], subset: int | None, subset_idx: list[int], subset_names: list[str]) -> dict[str, str]:
    """Return a subset of {names: seqs} while retaining sorted order."""
    if subset:
        return {i: seqs[i] for i in list(seqs)[:subset]}

    # index and filter...
    if subset_idx:
        if 0 in subset_idx:
            raise ValueError("subset-idx args should be 1-indexed; 0 is an invalid option.")
        names = list(seqs)
        return {i: j for (i, j) in seqs.items() if names.index(i) - 1 in subset_idx}
        # return {names[i]: seqs[names[i]] for i in sorted(subset_idx)}

    # ...
    if subset_names:
        bad_names = set(subset_names).difference(set(seqs))
        if bad_names:
            raise ValueError(f"scaffold names {bad_names} not found in the genome file")
        logger.debug("subselecting scaffolds with headers matching subset-names. Example: ")
        return {i: seqs[i] for i in seqs if i in subset_names}
    return seqs


def relabel(seqs: dict[str, str], relabel_prefix: str, relabel_prefix_idx: str, relabel_list: list[str]) -> dict[str, str]:
    """..."""
    # relabel numerically with prefix
    if relabel_prefix:
        logger.debug(f"relabeling scaffolds with prefix '{relabel_prefix}{{existing}}'")
        seqs = {f"{relabel_prefix}{i + 1}": j for i, j in enumerate(seqs.values())}

    # relabel numerically with prefix
    if relabel_prefix_idx:
        logger.debug(f"relabeling scaffolds with prefix '{relabel_prefix_idx}{{count}}'")
        seqs = {f"{relabel_prefix_idx}{i + 1}": j for i, j in enumerate(seqs.values())}

    # relabel individually with list of names
    if relabel_list:
        assert len(relabel_list) == len(seqs), "relabel list must be same lengths as the number of selected scaffolds"
        logger.debug(f"relabeling scaffolds with provided labels, e.g., {list(seqs.values())[0]} -> {relabel_list[0]}")
        seqs = dict(zip(relabel_list, seqs.values()))
    return seqs


def run_genome_format(args):
    """Run the formatting options and write the result to stdout"""
    sdict = get_headers_to_seqs_dict(args.genome)
    sdict = sort_headers_to_seqs_dict(sdict, args.sort_len, args.sort_alpha)
    sdict = subselect(sdict, args.subset, args.subset_idx, args.subset_names)
    sdict = relabel(sdict, args.relabel_prefix, args.relabel_prefix_idx, args.relabel_list)

    # print genome to stdout
    for header, seq in sdict.items():
        sys.stdout.write(f">{header}\n")
        if args.unmask:
            seq = seq.upper()
        chunk_size = 1024 * 1024 * 2  # 2 Mb
        bigseq = "".join(seq)
        if args.max_width:
            for i in range(0, len(bigseq), chunk_size):
                chunk = bigseq[i: i + chunk_size]
                wrapped_chunk = textwrap.fill(chunk, width=args.max_width)
                sys.stdout.write(wrapped_chunk + "\n")
        else:
            for i in range(0, len(bigseq), chunk_size):
                chunk = bigseq[i: i + chunk_size]
                sys.stdout.write(chunk + "\n")


def main():
    parser = get_parser_format_genome()
    args = parser.parse_args()
    # allow graceful exit if output is piped but pipe is broken
    try:
        run_genome_format(args)
    except BrokenPipeError:
        sys.exit(0)



if __name__ == "__main__":
    main()
