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



def call_diamond_blastp_parallel(proj: Project, pool: ProcessPoolExecutor) -> None:
    """

    """
    # default threading if not provided as a pool attr
    threads = getattr(pool, 'threads', 1)

    # get all pairwise searches
    searches = []
    for name1 in proj.samples:
        for name2 in proj.samples:
            searches.append((proj.samples[name1], proj.samples[name2]))
    total = len(searches)

    # remove any previously unfinished searches
    to_search = []
    for (samp1, samp2) in searches:
        tsv = proj.blastdir / f"{samp1.name}_to_{samp2.name}.tsv.gz"
        if not tsv.exists():
            to_search.append((samp1, samp2))
        else:
            logger.debug(f"{tsv} {tsv.stat().st_size}")
            if not tsv.stat().st_size:
                to_search.append((samp1, samp2))
                logger.debug(f"removed previously unfinished file {tsv}")
                tsv.unlink()
    n_to_search = len(to_search)
    finished = total - n_to_search

    # if any new searches will be done then remove the graphdir
    if n_to_search:
        logger.info(f"running all-by-all blastp searches (total={total}; completed={finished})")
        old_dir = proj.rundir / "graph_dir"
        if old_dir.exists():
            rmtree(old_dir)

    # submit enough jobs to fill the pool
    futures = []
    for i in range(pool.cores):
        if to_search:
            samp1, samp2 = to_search.pop()
            args = (proj, samp1, samp2, threads)
            fut = pool.submit(call_diamond_blastp, *args)
            fut.names = (samp1.name, samp2.name)
            futures.append(fut)

    # loop until all unaligned files have been aligned
    checkpoint = max(5, n_to_search // 10)
    while finished < total:
        # block until a future finishes, check errors.
        future = next(as_completed(futures))
        path = future.result()

        # store result path
        s1, s2 = future.names
        proj.samples[s1].files.blast[s2] = path

        # cleanup; incr counter; save
        futures.remove(futures[futures.index(future)])
        del future
        finished += 1
        proj.save_json()

        # submit a new job to queue
        if to_search:
            samp1, samp2 = to_search.pop()
            args = (proj, samp1, samp2, threads)
            fut = pool.submit(call_diamond_blastp, *args)
            fut.names = (samp1.name, samp2.name)
            futures.append(fut)

        # log progress occasionally
        if not finished % checkpoint:
            logger.info(f"finished blastp searches {finished}/{total}")
    logger.info(f"all pairwise blastp searches completed ({total})")
    proj.save_json()


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
