#!/usr/bin/env python

"""...


Example
-------
$ twig diamond-blast-all -d A.faa B.faa C.faa -o ./data

Output
------
./data
|_ A_B.diamond.tsv
|_ A_C.diamond.tsv
|_ B_C.diamond.tsv

Format
------
...

"""

import sys
import textwrap
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
import itertools
from tempfile import gettempdir
from loguru import logger
# import numpy as np
# from twig.utils.path_utils import expand_multiple_paths
# from concurrent.futures import ProcessPoolExecutor, as_completed

DIAMOND_BIN = Path(sys.prefix) / "bin" / "diamond"
TMPDIR = gettempdir()


def get_parser_diamond_pw(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="diamond-pw",
        usage="%(prog)s paths [args]",
        help="write .tsv blast hits for all pairs of fasta sequence files",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=textwrap.dedent("""
            -------------------------------------------------------------------
            | diamond-pw: write .tsv's of all-by-all pairwise diamond blastp  |
            -------------------------------------------------------------------
            | Perform all-by-all blastp on a set of fasta protein sequence
            | files. There are many more options available from the `diamond
            | blast` tool itself. This is intended only as a simple wrapper.  |
            -------------------------------------------------------------------
        """),
        epilog=textwrap.dedent("""
            Examples
            --------
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5 --no-self
            $ diamond-pw -d d1.fa d2.fa d3.fa -o blast/ -e 1e-5 -k 1 -j 10 -t 4

            # Example result dir (blast/) contains
            # d1_d1.tsv  d2_d2.tsv  d3_d3.tsv
            # d1_d2.tsv  d1_d3.tsv  d2_d3.tsv
            # d2_d1.tsv  d3_d1.tsv  d3_d2.tsv
            # sample_map.tsv
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
    parser.add_argument("-i", "--input", type=Path, metavar="path", help="fasta sequence files for all-by-all blastp")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", default="./blast", help="directory to write results (created if necessary) [%(default)s]")
    parser.add_argument("-e", "--evalue", type=float, metavar="float", default=1e-5, help="max evalue to report an alignment [%(default)s]")
    parser.add_argument("-m", "--min-bitscore", type=float, metavar="float", default=None, help="min bit-score to report an alignment (overrides -e) [%(default)s]")
    parser.add_argument("-k", "--max-target-seqs", type=int, metavar="int", default=25, help="max n target seqs to report hits for [%(default)s]")
    parser.add_argument("-d", "--dbdir", type=Path, metavar="path", default=TMPDIR, help="directory for temp .db files [%(default)s]")
    parser.add_argument("-t", "--threads", type=int, default=4, help="number of threads per diamond job [%(default)s]")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="number of diamond jobs to run in parallel [%(default)s]")
    parser.add_argument("-n", "--no-self", action="store_true", help="do not perform self-self search")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite .db files if they exist in dbdir")
    # parser.add_argument("-r", "--relabel-headers", type=int, default=1, help="...")
    # parser.add_argument("-r", "--relabel-taxa", type=int, default=1, help="...")
    # parser.add_argument("-I", "--imap", type=Path, help="map file to translate sample labels.")
    return parser


def run_diamond_pw(args):
    """..."""
    logger.info(args.data)



def main():
    parser = get_parser_diamond_pw()
    args = parser.parse_args()
    run_diamond_pw(args)


if __name__ == "__main__":
    pass
