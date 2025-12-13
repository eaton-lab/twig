#!/usr/bin/env python

"""Command-line interface for twig.

Examples
--------
twig -h
twig [subcommand]
"""

from typing import Optional
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from textwrap import dedent
from loguru import logger
from ..synteny.format_fasta import run_format_fasta, get_parser_format_fasta
from ..synteny.format_gff import run_format_gff, get_parser_format_gff
from ..synteny.genome_table import run_genome_table, get_parser_genome_table
# from ..synteny.format_data import run_format_data, get_parser_format_data
# from ..synteny.dot_draw import run_dot_draw, get_parser_dot_draw
# from ..synteny.dot_coords import run_dot_coords, get_parser_dot_coords
from ..seqlab.diamond_bl import run_diamond_bl, get_parser_diamond_bl
from ..seqlab.diamond_pw import run_diamond_pw, get_parser_diamond_pw
from ..seqlab.macse_prep import run_macse_prep, get_parser_macse_prep
from ..seqlab.macse_align import run_macse_align, get_parser_macse_align
from ..seqlab.macse_refine import run_macse_refine, get_parser_macse_refine
from ..seqlab.csubst import run_csubst, get_parser_csubst
from ..treelab.tree_filter import run_tree_filter, get_parser_tree_filter
from ..treelab.tree_rooter import run_tree_rooter, get_parser_tree_rooter
#from ..utils.logger_setup import set_log_level
# from ..utils.make_wide import make_wide
from twig import __version__ as VERSION

THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def setup_parsers() -> ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = ArgumentParser(
        "twig",
        usage="twig [subcommand] --help",
        formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        description=dedent("""
            -----------------------------------------------------
            |  %(prog)s: tree-based workflows for integrative genomics |
            -----------------------------------------------------
            """),
        epilog=dedent(r"""
            Tutorials
            ---------
            # genome file formatting
            https://eaton-lab.org/twig/formatting

            # gene orthogroup analysis
            https://eaton-lab.org/twig/orthology

            # gene synteny analysis
            https://eaton-lab.org/twig/synteny

            # sequence rate/convergence/etc analysis
            https://eaton-lab.org/twig/evolution
        """)
    )
    parser.add_argument("-v", "--version    ", action='version', version=f"twig {VERSION}")
    subparsers = parser.add_subparsers(
        prog="%(prog)s", required=True, 
        title="subcommands", 
        dest="subcommand", 
        metavar="--------------", 
        help="-----------------------------------------------------",
    )

    get_parser_csubst(subparsers)
    get_parser_genome_table(subparsers)
    get_parser_format_fasta(subparsers)
    get_parser_format_gff(subparsers)
    # get_parser_format_data(subparsers)
    # get_parser_dot_coords(subparsers)
    # get_parser_dot_draw(subparsers)
    get_parser_diamond_bl(subparsers)
    get_parser_diamond_pw(subparsers)
    # get_parser_lastal_align(subparsers)
    # get_parser_ortholog_graph(subparsers)
    # get_parser_ortholog_...(subparsers)    
    get_parser_macse_prep(subparsers)
    get_parser_macse_align(subparsers)
    get_parser_macse_refine(subparsers)
    get_parser_tree_filter(subparsers)
    get_parser_tree_rooter(subparsers)    

    # OLD
    # _add_parser_seqlab_run(subparsers)
    # _add_parser_branch(subparsers)
    # _add_parser_continue(subparsers)
    # _add_parser_info(subparsers)
    # _add_parser_seqlab_prep(subparsers)
    # _add_parser_seqlab_score(subparsers)
    # _add_parser_seqlab_graph(subparsers)

    return parser


def main(cmd: Optional[str] = None) -> int:
    """Command line tool.

    """
    parser = setup_parsers()
    args = parser.parse_args(cmd.split() if cmd else None)

    dispatch = {
        # "dot-coords": run_dot_coords,
        # "dot-draw": run_dot_draw,
        "csubst": run_csubst,
        "diamond-bl": run_diamond_bl,
        "diamond-pw": run_diamond_pw,
        "format-fasta": run_format_fasta,
        "format-gff": run_format_gff,
        # "format-data": run_format_data,
        "genome-table": run_genome_table,
        "macse-prep": run_macse_prep,
        "macse-align": run_macse_align,
        "macse-refine": run_macse_refine,
        "tree-filter": run_tree_filter,
        "tree-rooter": run_tree_rooter,
    }

    # in cli pin numpy to single-threading
    for k in THREAD_ENV_VARS:
        os.environ.setdefault(k, '1')

    # run function and handle errors
    if args.subcommand:
        try:
            run_func = dispatch[args.subcommand]
            run_func(args)
            return 0
        except KeyboardInterrupt:
            logger.warning("interrupted by user")
        except Exception as exc:
            logger.error(exc)
            raise
        return 1

    # unreachable
    parser.print_help()
    return 1


if __name__ == "__main__":
    
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except RuntimeError as exc:
        logger.error(exc)
    except Exception as exc:
        logger.error(exc)
