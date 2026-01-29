#!/usr/bin/env python

"""Prep align sequences

If you call:
$ twig macse-align -i CDS -o OUT/ID.msa.nt.fa

It will produce:
- OUT/ID.msa.nt.fa
"""

import subprocess
import sys
from pathlib import Path
from loguru import logger
# from twig.utils.parallel import run_pipeline  # , run_with_pool
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")


def run_macse_align(args):
    """..."""
    set_log_level(args.log_level)

    # check that macse is in PATH
    assert Path(BIN_MACSE).exists(), f"macse binary not found. Checked: {BIN_MACSE}"

    # check in files
    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")

    # ensure outpath and pdir exists
    if args.outpath.is_dir():
        raise IOError("outpath must be a file path not dir")
    args.outpath.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    if args.outpath.exists() and not args.force:
        logger.warning(f"[{args.outpath.name}] [skipping] {args.outpath} already exists. Use --force to overwrite")
        return 0
    call_macse_align(args.input, args.outpath, args.max_refine_iter, args.verbose)
    logger.info(f"[{args.outpath.name}] alignment written to {args.outpath}")


def call_macse_align(cds_fasta: Path, outpath: str, max_iter: int, verbose: bool):
    """Run Alignment step with default settings"""
    cmd = [
        BIN_MACSE, "-prog", "alignSequences",
        "-seq", str(cds_fasta),
        "-out_NT", f"{outpath}",
        "-out_AA", f"{outpath}.tmp.aa.fa",
        "-max_refine_iter", str(max_iter),
    ]
    logger.debug(f"[{outpath.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise Exception(proc.stderr + proc.stdout)
    (outpath.with_suffix(outpath.suffix + ".tmp.aa.fa")).unlink()
    return proc.returncode


def main():
    from ..cli.subcommands import get_parser_macse_align
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
