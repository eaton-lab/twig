#!/usr/bin/env python

"""csubst wrapper

conda create -n csubst
conda install csubst 'iqtree<3' -c bioconda
"""

import subprocess
import sys
from pathlib import Path
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


def run_csubst(args):
    """..."""
    set_log_level(args.log_level)

    BIN = Path(sys.prefix).parent / f"{args.env}" / "bin"
    global BIN_CSUBST, BIN_IQTREE
    BIN_CSUBST = str(BIN / "csubst")
    BIN_IQTREE = str(BIN / "iqtree")

    # check that macse is in PATH
    if not Path(BIN_CSUBST).exists():
        logger.error(INSTALL_MSG.format(args.env))
        sys.exit(1)

    # ensure outdir exists
    # args.outdir.mkdir(exist_ok=True)
    # args.prefix = args.prefix if args.prefix is not None else args.alignment.name

    # check to create workdir
    args.workdir = args.outdir
    if args.workdir.exists() and args.workdir.is_dir():
        if not args.force:
            raise ValueError(f"path exists at {args.workdir}/. Use --force to overwrite.")
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
    logger.info(f"[{args.alignment}] csubst result written to {args.workdir}/")


def call_csubst(args):
    """..."""
    """Run Alignment step with default settings"""
    cmd = [
        BIN_CSUBST, "analyze",
        "--alignment_file", str(args.alignment),
        "--rooted_tree_file", str(args.tree),
        "--foreground", str(args.foreground),
        "--fg_format", "2" if args.foreground_table else "1",
        "--fg_exclude_wg", "yes" if args.fg_exclude_wg else "no",
        "--fg_stem_only", "yes" if args.fg_stem_only else "no",
        "--threads", str(args.threads),
        "--max_arity", str(args.max_arity),
        "--exhaustive_until", str(args.exhaustive_until),
        "--iqtree_exe", str(BIN_IQTREE),
    ]
    logger.debug(f"[{args.alignment}] " + " ".join(cmd))
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
    from ..cli.subcommands import get_parser_csubst
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
