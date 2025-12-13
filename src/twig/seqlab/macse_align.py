#!/usr/bin/env python

"""Prep align sequences

"""

import subprocess
import sys
import textwrap
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
# from twig.utils.parallel import run_pipeline  # , run_with_pool
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")
ISOFORM_REGEX_DEFAULT = r"^([^|]+)\|.*?__(.+?)_i\d+"
# group 1 (shared part up until |)
# group 2 (after __ and up until first _i)


KWARGS = dict(
    prog="macse-align",
    usage="macse-align -i CDS -o OUTDIR [options]",
    help="run macse alignment (CDS/AA)",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | macse-align: ...
        -------------------------------------------------------------------
        | Macse ...
        | This implements `macse -prog AlignSequences` and
        | additional steps to ...
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        $ twig macse-align -i CDS -o OUT -p TEST

        # run parallel jobs on many cds files
        $ parallel -j 10 'twig macse-prep -i {} -o OUT -p {/.}' ::: CDS/*.fa
    """)
)


def get_parser_macse_align(parser: ArgumentParser | None = None) -> ArgumentParser:
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
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input CDS sequence (aligned or unaligned)")
    parser.add_argument("-o", "--out", type=Path, metavar="path", required=True, help="out prefix; parent dirs created if necessary [{input}]")
    # others
    parser.add_argument("-m", "--max-refine-iter", type=int, metavar="int", default=-1, help="max refinement iterations during optimizing [default -1 = no limit]")
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    # parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def run_macse_align(args):
    """..."""
    set_log_level(args.log_level)#, args.log_file)

    # check that macse is in PATH
    assert Path(BIN_MACSE).exists(), f"macse binary not found. Checked: {BIN_MACSE}"

    # check in files
    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")

    # ensure outpath and pdir exists
    args.outprefix = args.out
    if args.outprefix is None:
        args.outprefix = args.input
    args.outprefix.parent.mkdir(exist_ok=True)

    # bail out if final file exists
    result = args.outprefix.with_suffix(args.outprefix.suffix + ".msa.nt.fa")
    if result.exists() and not args.force:
        logger.warning(f"[{args.input.name}] [skipping] {result} already exists. Using --force to overwrite")
        return 0
    call_macse_align(args.input, args.outprefix, args.max_refine_iter, args.verbose)
    logger.info(f"[{args.input.name}] alignment written to {args.outprefix}.nt.fa")


def call_macse_align(cds_fasta: Path, outdir: Path, prefix: str, max_iter: int, verbose: bool):
    """Run Alignment step with default settings"""
    cmd = [
        BIN_MACSE, "-prog", "alignSequences",
        "-seq", str(cds_fasta),
        "-out_NT", str(outdir / f"{prefix}.msa.nt.fa"),
        "-out_AA", str(outdir / f"{prefix}.tmp.msa.aa.fa"),
        "-max_refine_iter", str(max_iter),
    ]
    logger.debug(f"[{prefix}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    (outdir / f"{prefix}.tmp.msa.aa.fa").unlink()
    return proc.returncode


def main():
    parser = get_parser_macse_align()
    args = parser.parse_args()
    run_macse_align(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
