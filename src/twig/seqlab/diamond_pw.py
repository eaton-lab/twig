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
from pathlib import Path
import itertools
from tempfile import gettempdir
from loguru import logger
# import numpy as np
# from twig.utils.path_utils import expand_multiple_paths
# from concurrent.futures import ProcessPoolExecutor, as_completed

DIAMOND_BIN = Path(sys.prefix) / "bin" / "diamond"
TMPDIR = gettempdir()


def run_diamond_pw(args):
    """..."""
    logger.info(args.data)


def main():
    from ..cli.subcommands import get_parser_diamond_pw
    parser = get_parser_diamond_pw()
    args = parser.parse_args()
    run_diamond_pw(args)


if __name__ == "__main__":
    pass
