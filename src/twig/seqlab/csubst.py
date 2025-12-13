#!/usr/bin/env python

"""csubst wrapper

conda create -n csubst
conda install csubst 'iqtree<3' -c bioconda
"""

import subprocess
import sys
from textwrap import dedent
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
from twig.utils.logger_setup import set_log_level

INSTALL_MSG = """\
A 'csubst' binary was not found in a conda env named '{env}'.

Due to conflicting dependencies, I recommend creating a separate conda
env to install csubst into, and which can be called from twig by using
the -e parameter to specify the env name, or using the name 'csubst',
as in the installation command below. Note, to call csubst from twig,
do not activate the csubst env, just specify the name with -e.

$ conda create -n csubst csubst 'iqtree<3' Python=3.10 --strict-channel -c conda-forge -c bioconda
$ twig csubst -a MSA -t NWK -g FG -o OUT/PRE -e csubst

"""

KWARGS = dict(
    prog="csubst",
    usage="twig csubst -i MSA -t TREE -g FG [options]",
    help="run csubst in wrapper",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=dedent("""
        -------------------------------------------------------------------
        |
        -------------------------------------------------------------------
        |
        |
        |
        -------------------------------------------------------------------
    """),
    epilog=dedent("""
        Examples
        --------
        $ twig csubst -a MSA -t NWK -g FG -o OUT/PRE   # OUT/PRE.cb...

        $ twig csubst -a MSA -t NWK -g FG -o OUT/PRE   # OUT/PRE.cb...
    """)
)


def get_parser_csubst(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path args
    parser.add_argument("-a", "--alignment", type=Path, metavar="path", required=True, help="input CDS alignment")
    parser.add_argument("-t", "--tree", type=Path, metavar="path", required=True, help="input rooted tree file")
    parser.add_argument("-g", "--foreground", type=Path, metavar="path", help="foreground file")
    # parser.add_argument("-o", "--out", type=Path, metavar="path", help="out prefix; default is input path [{input}]")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", help="output directory. Created if it doesn't exist")
    parser.add_argument("-p", "--prefix", type=str, metavar="str", help="optional outfile prefix. If None the cds filename is used")

    parser.add_argument("-m", "--max-arity", type=int, metavar="int", default=2, help="max combinatorial number of branches (K)")
    parser.add_argument("-u", "--exhaustive-until", type=int, metavar="int", default=1, help="perform exhaustive (non-heuristic) search up N branch combs")
    parser.add_argument("-c", "--cutoff-stat", type=str, metavar="str", default="OCNany2spe,2.0|omegaCany2spe,5.0", help="Cutoff stats for searching higher-order branch combs [%(default)s]")
    parser.add_argument("-F", "--foreground-table", action="store_true", help="foreground file is a table (fg_format=2)")

    # others
    parser.add_argument("-e", "--env", type=Path, metavar="path", help="conda env name where 'csubst' in installed [csubst]")
    parser.add_argument("-j", "--threads", type=int, metavar="int", default=1, help="number of threads")
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def run_csubst(args):
    """..."""
    set_log_level(args.log_level)

    BIN = Path(sys.prefix).parent / f"{args.env}" / "bin"
    BIN_CSUBST = str(BIN / "csubst")
    # BIN_IQTREE = str(BIN / "iqtree")

    # check that macse is in PATH
    if not Path(BIN_CSUBST).exists():
        logger.error(INSTALL_MSG)
        sys.exit(1)

    # ensure outdir exists
    args.outdir.mkdir(exist_ok=True)
    args.prefix = args.prefix if args.prefix is not None else args.alignment.name

    # check to create workdir
    args.workdir = args.outdir / args.prefix
    if args.workdir.exists() and args.workdir.is_dir():
        if not args.force:
            raise ValueError(f"path exists at {args.workdir}/{args.prefix}. Use --force to overwrite.")
        else:
            # check if it contains only csubst files
            for item in args.workdir.iterdir():
                if item.is_file() and (item.name.startswith("csubst") or item.name.startswith("tmp.csubst")):
                    item.unlink()
            args.workdir.rmdir()
            logger.info(f"removed existing results in {args.workdir}")
    args.workdir.mkdir(exist_ok=True)

    # ...
    args.alignment = args.alignment.expanduser().absolute()
    args.tree = args.tree.expanduser().absolute()
    args.foreground = args.foreground.expanduser().absolute()
    for argpath in [args.alignment, args.tree, args.foreground]:
        assert argpath.exists(), f"path '{argpath}' does not exist"

    # bail out if final file exists
    # result = args.outdir / (args.prefix + ".msa.nt.fa")
    # if result.exists() and not args.force:
    #     logger.info(f"[{args.prefix}] [skipping] {result} already exists")
    #     return 0

    # run it
    call_csubst(args)
    logger.info(f"[{args.prefix}] csubst result written to {args.workdir}/")


def call_csubst(args):
    """..."""
    """Run Alignment step with default settings"""
    cmd = [
        BIN_CSUBST, "analyze",
        "--alignment_file", str(args.alignment),
        "--rooted_tree_file", str(args.tree),
        "--foreground", str(args.foreground),
        "--fg_format", "2" if args.foreground_table else "1",
        "--threads", str(args.threads),
        "--max_arity", str(args.max_arity),
        "--exhaustive_until", str(args.exhaustive_until),
        "--iqtree_exe", str(BIN_IQTREE),
    ]
    logger.debug(f"[{args.prefix}] " + " ".join(cmd))
    if args.verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True, cwd=args.workdir)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=args.workdir)
    if proc.returncode:
        raise Exception(proc.stderr)
        # raise subprocess.CalledProcessError(cmd, proc.stderr)
    # (outdir / f"{prefix}.tmp.msa.aa.fa").unlink()
    # return proc.returncode


def main():
    parser = get_parser_csubst()
    args = parser.parse_args()
    run_csubst(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
        raise
