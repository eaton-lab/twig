#!/usr/bin/env python

"""Format multi-fasta file.

This takes a fasta file and returns a modified copy.

Example
-------
$ python format_genome.py GENOME ...
$ twig format-genome GENOME ...

Output
------
Writes fasta to stdout, but can also write additional metadata files
to a prefix path.
"""

import sys
import gzip
from io import StringIO
from pathlib import Path
from loguru import logger


def get_headers_to_seqs_dict_subset(genome: Path, min_len: int, subset: int, subset_idx: list[int], subset_names: list[str]) -> dict[str, str]:
    """Return dict mapping scaffold names to their sequences.

    This method is faster than the one below b/c it does not read the
    entire fasta file, but instead will stop once the subset of selected
    scaffolds has been found. This can make a big difference in speed
    for very large files, especially when selecting the first few loci.
    """
    # keep track of sequence lengths
    headers_to_seqs = {}
    
    # select opener func
    if genome.suffix == ".gz":
        xopen = gzip.open
    else:
        xopen = open

    # read through genome file
    idx = 0
    paused = False
    with xopen(genome, mode='rt') as file:
        for line in file:
            line = line.strip()
            # skip empty lines
            if not line:
                continue
            # start new scaffold
            if line.startswith(">"):
                idx += 1
                header = line[1:]

                # breaking conditions for faster stopping
                if subset:
                    if len(headers_to_seqs) == subset:
                        break
                elif subset_idx is not None:
                    if not subset_idx:
                        break
                    if idx in subset_idx:
                        subset_idx.remove(idx)
                        paused = False
                    else:
                        paused = True
                elif subset_names is not None:
                    if not subset_names:
                        break
                    if header in subset_names:
                        subset_names.remove(header)
                        paused = False
                    else:
                        paused = True

                # store header
                if not paused:
                    headers_to_seqs[header] = []
            # append to existing scaffold
            else:
                if not paused:
                    headers_to_seqs[header].append(line)

    # remove scaffs shorter than min_len
    if min_len:
        start = len(headers_to_seqs)
        headers_to_seqs = {i: j for (i, j) in headers_to_seqs.items() if sum(len(s) for s in j) > min_len}
        removed = start - len(headers_to_seqs)
        logger.info(f"format-fasta removed {removed} scaffolds shorter than {min_len} bp")
    return headers_to_seqs


def get_headers_to_seqs_dict(genome: Path, min_len: int) -> dict[str, str]:
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

    # remove scaffs shorter than min_len
    if min_len:
        start = len(headers_to_seqs)
        headers_to_seqs = {i: j for (i, j) in headers_to_seqs.items() if sum(len(s) for s in j) > min_len}
        removed = start - len(headers_to_seqs)
        logger.info(f"format-fasta removed {removed} scaffolds shorter than {min_len} bp")
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


def reverse_complement(seq: str, dna: bool) -> str:
    """Return reverse complement of DNA, or reverse sequence of AA."""
    if not dna:
        return seq[::-1]
    else:
        # define the complement translation map for DNA
        complement_map = str.maketrans("ATCGatcg", "TAGCtagc")
    
        # translate sequence to get complement, then reverse it
        return seq.translate(complement_map)[::-1]


def contains_invalid_dna_letters(sequence: str) -> bool:
    """Return True if not DNA sequence"""
    valid_nucleotides = {"A", "C", "G", "T"}
    return any(letter not in valid_nucleotides for letter in sequence.upper())


def run_format_fasta(args):
    """Run the formatting options and write the result to stdout"""
    try:
        # get subset dict, no sort, faster option
        if any([args.subset, args.subset_idx, args.subset_names]) and not any([args.sort_len, args.sort_alpha]):
            sdict = get_headers_to_seqs_dict_subset(args.fasta, args.min_len, args.subset, args.subset_idx, args.subset_names)
            rdict = relabel(sdict, args.relabel_prefix, args.relabel_prefix_idx, args.relabel_list)

        # get full dict, sort, subset, and modify
        else:            
            sdict = get_headers_to_seqs_dict(args.fasta, args.min_len)
            sdict = sort_headers_to_seqs_dict(sdict, args.sort_len, args.sort_alpha)
            sdict = subselect(sdict, args.subset, args.subset_idx, args.subset_names)
            rdict = relabel(sdict, args.relabel_prefix, args.relabel_prefix_idx, args.relabel_list)

        # optional write translation table
        if args.map:
            logger.debug(f"fasta labels translation .tsv written to {args.map}")
            with open(args.map, 'w') as file:
                file.write("\n".join("\t".join([i, j]) for (i, j) in zip(rdict, sdict)))

        # if doing revcomp check if it is DNA or AA sequences
        dna = True
        if args.revcomp:
            for i in list(rdict)[:100]:
                seq = "".join(rdict[i])
                if contains_invalid_dna_letters(seq):
                    dna = False

        # print fasta buffered to stdout
        buffer = StringIO()
        flush_size = 1024 * 1024 * 20  # 2 Mb
        for idx, (header, seqlist) in enumerate(rdict.items()):
            # sys.stdout.write(f">{header}\n")
            buffer.write(f">{header}\n")

            # get full locus sequence
            seq = "".join(seqlist)

            # optionally remove soft masking
            if args.unmask:
                seq = seq.upper()

            # optionally revcomp (1-indexed selector)
            if args.revcomp:
                seq = reverse_complement(seq, dna)

            # 
            if args.max_width:
                for i in range(0, len(seq), args.max_width):
                    buffer.write(seq[i:i + args.max_width] + "\n")
            else:
                buffer.write(seq + "\n")
                # for i in range(0, len(seq), chunk_size):
                #     chunk = seq[i: i + chunk_size]
                #     sys.stdout.write(chunk + "\n")

            # dump buffer to stdout
            if buffer.tell() >= flush_size:
                sys.stdout.write(buffer.getvalue())
                buffer.seek(0)  # Reset buffer position to the start
                buffer.truncate(0)  # Clear buffer contents

        # dump any remaining buffer to stdout
        sys.stdout.write(buffer.getvalue())
        buffer.close()

    # allow graceful exit when pipe is broken (e.g., ... | head -n 10)
    except BrokenPipeError:
        sys.exit(0)
    except KeyboardInterrupt:
        logger.error("KeyboardInterrupt by user")
        sys.exit(0)


def main():
    """module-level main cli"""
    from ..cli.subcommands import get_parser_format_fasta
    parser = get_parser_format_fasta()
    args = parser.parse_args()
    # allow graceful exit if output is piped but pipe is broken
    run_format_fasta(args)


if __name__ == "__main__":
    main()
