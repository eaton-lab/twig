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
    parser.add_argument("-i", "--cds", type=Path, metavar="path", required=True, help="input CDS sequence (aligned or unaligned)")
    parser.add_argument("-o", "--outdir", type=Path, metavar="path", required=True, help="output directory, created if it doesn't exist")
    parser.add_argument("-p", "--prefix", type=str, metavar="str", help="optional outfile prefix. If None the cds filename is used")
    # parser.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="optional names or glob to exclude one or more sequences")
    # parser.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="optional names or glob to include only a subset sequences")
    # options
    # parser.add_argument("-ac", "--aln-trim-ends-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment edges to where a min percent of samples have data [%(default)s]")
    # parser.add_argument("-as", "--aln-trim-window-size", type=int, metavar="int", default=5, help="trim alignment using a sliding 'half_window_size' defined as ... [%(default)s]")
    # parser.add_argument("-r", "--refine", action="store_true", help="run refineAlignment instead of AlignSequences")

    # others
    # parser.add_argument("-B", "--binary", type=Path, metavar="path", help="path to macse binary if not in $PATH")
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

    # ensure outdir exists
    args.outdir.mkdir(exist_ok=True)
    args.prefix = args.prefix if args.prefix is not None else args.cds.name

    # bail out if final file exists
    result = args.outdir / (args.prefix + ".msa.nt.fa")
    if result.exists() and not args.force:
        logger.info(f"[{args.prefix}] [skipping] {result} already exists")
        return 0
    call_macse_align(args.cds, args.outdir, args.prefix, args.max_refine_iter, args.verbose)
    logger.info(f"[{args.prefix}] alignment written to {args.outdir}/{args.prefix}.nt.fa")


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
