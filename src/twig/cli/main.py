#!/usr/bin/env python

"""Command-line interface for twig.

Examples
--------
twig -h
twig [subcommand]
"""

from typing import Optional
import argparse
import textwrap
from loguru import logger
from ..synteny.format_fasta import run_format_fasta, get_parser_format_fasta
from ..synteny.format_gff import run_format_gff, get_parser_format_gff
# from ..synteny.format_data import run_format_data, get_parser_format_data
# from ..synteny.genome_table import run_genome_table, get_parser_genome_table
# from ..synteny.dot_draw import run_dot_draw, get_parser_dot_draw
# from ..synteny.dot_coords import run_dot_coords, get_parser_dot_coords
from ..seqlab.diamond_blastp import run_diamond_blastp, get_parser_diamond_blastp
from ..seqlab.macse import run_macse, get_parser_macse
from ..treelab.tree_filter import run_tree_filter, get_parser_tree_filter
from ..treelab.tree_rooter import run_tree_rooter, get_parser_tree_rooter
#from ..utils.logger_setup import set_log_level
# from ..utils.make_wide import make_wide
from twig import __version__ as VERSION

CLI_NAME = "twig"
# VERSION = "0.0.1"


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = argparse.ArgumentParser(
        CLI_NAME,
        usage=f"{CLI_NAME} [subcommand] --help",
        formatter_class=lambda prog: argparse.RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
        # formatter_class=CustomHelpFormatter,
        description=textwrap.dedent("""
            -----------------------------------------------------
            |  %(prog)s: comparative genomics on trees              |
            -----------------------------------------------------
            """),
        epilog=textwrap.dedent(r"""
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

    # get_parser_genome_table(subparsers)
    # get_parser_format_data(subparsers)
    get_parser_format_fasta(subparsers)
    get_parser_format_gff(subparsers)
    # get_parser_dot_coords(subparsers)
    # get_parser_dot_draw(subparsers)
    get_parser_diamond_blastp(subparsers)
    # get_parser_diamond_blastp-all(subparsers)    
    # get_parser_lastal_align(subparsers)
    # get_parser_ortholog_graph(subparsers)
    # get_parser_ortholog_...(subparsers)    
    get_parser_macse(subparsers)
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
        "diamond-blastp": run_diamond_blastp,
        # "format-data": run_format_data,
        "format-fasta": run_format_fasta,
        "format-gff": run_format_gff,
        # "genome-table": run_genome_table,
        "macse": run_macse,
        "tree-filter": run_tree_filter,
        "tree-rooter": run_tree_rooter,
    }

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
            logger.exception(exc)            
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
