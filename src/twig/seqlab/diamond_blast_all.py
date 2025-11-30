#!/usr/bin/env python

"""...


Example
-------
$ twig diamond-blast-all A.faa B.faa C.faa -o ./data

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
from subprocess import Popen, STDOUT, PIPE
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from pathlib import Path
from tempfile import gettempdir
from loguru import logger
import numpy as np
from twig.utils.make_wide import make_wide
from twig.utils.path_utils import expand_multiple_paths
from concurrent.futures import ProcessPoolExecutor, as_completed


logger = logger.bind(name="twig")
DIAMOND_BIN = Path(sys.prefix) / "bin" / "diamond"
TMPDIR = gettempdir()


def get_parser(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    kwargs = dict(
        prog="diamond-blast-all",
        usage="%(prog)s paths [args]",
        help="write .tsv blast hits for all pairs of fasta sequence files",
        formatter_class=make_wide(RawDescriptionHelpFormatter),
        description=textwrap.dedent("""
            -------------------------------------------------------------------
            | diamond-blast-all: write .tsv's of all-by-all diamond blast hits |
            -------------------------------------------------------------------
            | ...  
            | There are many more options available from the `diamond blast`  |
            | tool itself. This is intended only as a simple wrapper.         |
            -------------------------------------------------------------------
        """),
        epilog=textwrap.dedent("""
            Examples
            --------
            $ diamond-blast-all d1.faa d2.faa d3.faa -o blast/ -e 1e-5
            $ diamond-blast-all d1.faa d2.faa d3.faa -o blast/ -e 1e-5 --no-self
            $ diamond-blast-all d1.faa d2.faa d3.faa -o blast/ -e 1e-5 --k 1 
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
    parser.add_argument("query", type=Path, help="a fasta sequence to align to the database")
    parser.add_argument("database", type=Path, help="a fasta sequence or diamond db file as the database")    
    parser.add_argument("-e", "--evalue", type=float, default=1e-5, help="max evalue to report an alignment")
    parser.add_argument("-m", "--min-bitscore", type=float, default=None, help="min bit-score to report an alignment (overrides -e)")
    parser.add_argument("-k", "--max-target-seqs", type=int, default=25, help="max number of target sequences to report alignments for")
    parser.add_argument("-d", "--dbdir", type=Path, default=TMPDIR, help="directory for .db files")
    parser.add_argument("-f", "--force", type=bool, help="overwrite .db files if they exist in dbdir")
    parser.add_argument("-t", "--threads", type=int, default=4, help="number of threads per diamond job")
    parser.add_argument("-j", "--jobs", type=int, default=1, help="number of diamond jobs to run in parallel")
    # parser.add_argument("-m", "--map", type=Path, help="map file to translate sample labels.")
    return parser


