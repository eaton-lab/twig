#!/usr/bin/env python

"""Subselect, refine, trim, and export a CDS alignment to CDS/AA.

"""

from typing import List
import re
import sys
import textwrap
import subprocess
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
import toytree
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")
# ISOFORM_REGEX_DEFAULT = r'^(.+_i)(\d+)(\..*)?$'
# group 1 (shared part up until _i)
# group 2 (isoform index)
# group 3 (optional suffix)
ISOFORM_REGEX_DEFAULT = r"^([^|]+)\|.*?__(.+?)_i\d+"
# group 1 (shared part up until |)
# group 2 (after __ and up until first _i)


KWARGS = dict(
    prog="macse-refine",
    usage="macse-refine -i CDS [options]",
    help="...",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | macse-refine: CDS/AA jointly and filter low homology seqs
        -------------------------------------------------------------------
        | Macse codon-aware alignment of CDS and AA sequences. This method
        | also includes options to filter and trim sequences to remove low
        | homology fragments or sequences; sequences that are too short;
        | and extra isoforms. If a result with -p PREFIX name exists in -o
        | OUTDIR it will be skipped, allowing for easy checkpointing to
        | run in parallel on a large set of sequences and restart if
        | interrupted.
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        $ twig macse-refine -i CDS -o OUT -p TEST
        $ twig macse-refine -i CDS -o OUT -mh 0.1 -mi 0.5 -ti 50 -te 50 -mc 15 -c 20
        $ twig macse-refine -i CDS -o OUT -mh 0.3 -mi 0.8 -ti 25 -te 25 -mc 15 -c 20
        $ twig macse-refine -i CDS -o OUT -mh 0.5 -mc 10 -k -xa -ml 200 -e '^sppA.*'
        $ twig macse-refine -i CDS -o OUT -s '^sppA.*'

        # full pipeline
        $ twig macse-prep -i CDS -o OUT -p PRE
        $ twig macse-align -i OUT/PRE.nt.fa -o OUT -p PRE
        $ twig macse-refine -i OUT/PRE.msa.nt.fa -o OUT
    """)
)


def get_parser_macse_refine(parser: ArgumentParser | None = None) -> ArgumentParser:
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
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="input CDS alignment")
    parser.add_argument("-o", "--out", type=Path, metavar="path", required=True, help="out prefix; parent dirs created if necessary [{input}]")
    parser.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="optional names or glob to exclude one or more sequences")
    parser.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="optional names or glob to include only a subset sequences")
    parser.add_argument("-t", "--tree", type=Path, metavar="path", help="optional newick file to subsample genes present in tree")
    # options
    parser.add_argument("-ml", "--min-length", type=int, metavar="int", default=0, help="min length of non-missing sequence in a sample [%(default)s]")
    parser.add_argument("-ac", "--aln-trim-ends-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment edges to where a min percent of samples have data [%(default)s]")
    parser.add_argument("-as", "--aln-trim-window-size", type=int, metavar="int", default=5, help="trim alignment using a sliding 'half_window_size' defined as ... [%(default)s]")
    parser.add_argument("-if", "--codon-int-fs", type=str, metavar="str", default="NNN", help="codon to sub for internal frame shift [NNN]")
    parser.add_argument("-ef", "--codon-ext-fs", type=str, metavar="str", default="NNN", help="codon to sub for external frame shift [NNN]")
    parser.add_argument("-fs", "--codon-final-stop", type=str, metavar="str", default="NNN", help="codon to sub for final stop [NNN]")
    parser.add_argument("-is", "--codon-int-stop", type=str, metavar="str", default="NNN", help="codon to sub for internal stop [NNN]")

    parser.add_argument("-r", "--refine-alignment", action="store_true", help="refine alignment")
    parser.add_argument("-R", "--refine-alignment-if", action="store_true", help="refine alignment only if >=1 sequences are filtered out")
    parser.add_argument("-ri", "--max-iter-refine-alignment", type=int, metavar="int", default=-1, help="max iterations in refine alignment [%(default)s]")

    # others
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def filter_sequences(cds_fasta: Path, outprefix: Path, exclude: List[str], subsample: List[str], subsample_tree: Path, min_length: int, force: bool) -> None:
    """Write fasta with only one isoform per gene. When multple are present
    the one with greatest homology to other sequences is retained, with
    ties broken by length, and then order.
    """
    # logger.info("filtering alignment")
    pre = cds_fasta.name
    out = outprefix.with_suffix(outprefix.suffix + ".tmp.msa.nt.fa")

    # parse trimmed fasta file
    seqs = {}
    with open(cds_fasta, 'r') as datain:
        for line in datain.readlines():
            if line:
                if line.startswith(">"):
                    uname = line.strip()[1:]
                    seqs[uname] = ""
                else:
                    seqs[uname] += line.strip()

    # raise an exception if not aligned
    lengths = [len(i) for i in seqs.values()]
    if not len(set(lengths)) == 1:
        raise ValueError("input contains sequences of variable lengths (not aligned)")

    # expand exclude list of names/globs to select all name matches
    matched = []
    names = list(seqs)
    if exclude:
        patterns = [re.compile(p) for p in exclude]
        matched = [s for s in names if any(r.search(s) for r in patterns)]
    if subsample:
        patterns = [re.compile(p) for p in subsample]
        matched = [s for s in names if any(r.search(s) for r in patterns)]
        matched = list(set(names) - set(matched))
    if subsample_tree:
        matched = toytree.tree(subsample_tree).get_tip_labels()
        matched = list(set(names) - set(matched))
        for t in matched:
            if t not in seqs:
                logger.warning(t)
            else:
                logger.info(t)

    # group sequences by isoform regex
    f = {"min_length": 0, "user": 0}
    keep = {}
    for name, seq in seqs.items():

        # exclude if user-excluded
        if name in matched:
            logger.debug(f"[{pre}] {name} excluded by user args")
            f["user"] += 1
            continue

        # exclude if too short
        seq = seqs[name]
        nbases = sum(1 for i in seq if i != "-")
        if nbases < min_length:
            logger.debug(f"[{pre}] {name} excluded by min_length ({len(seq)})")
            f["min_length"] += 1
            continue
        keep[name] = seq

    # report
    logger.info(f"[{pre}] {len(seqs)} seqs -> {len(keep)} seqs, filtered by [min_length={f['min_length']}, user={f['user']}])")

    # write output
    if sum(f.values()):
        with open(out, 'w') as hout:
            for uname in keep:
                hout.write(f">{uname}\n{keep[uname]}\n")
        return out, True
    return cds_fasta, False


def call_macse_refine_alignment(data: Path, outprefix: str, force: bool, max_iter: int, verbose: bool):
    """Run Alignment step with default settings"""
    cmd = [
        BIN_MACSE, "-prog", "refineAlignment",
        "-align", str(data),
        "-out_NT", f"{outprefix}.rmsa.nt.fa",
        "-out_AA", f"{outprefix}.rmsa.aa.fa",
        "-max_refine_iter", str(max_iter),
    ]
    logger.info("refining alignment")
    logger.debug(f"[{data.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    return outprefix.with_suffix(outprefix.suffix + ".rmsa.nt.fa")


def call_macse_trim_alignment(data: Path, outprefix: str, half_window_size: int, min_percent_at_ends: float, verbose: bool, force: bool):
    """..."""
    dataout = outprefix.with_suffix(outprefix.suffix + ".tmp.msa.trimmed.nt.fa")
    datainfo = outprefix.with_suffix(outprefix.suffix + ".msa.trimmed.info")
    cmd = [
        BIN_MACSE, "-prog", "trimAlignment",
        "-align", str(data),
        "-respect_first_RF_ON",
        "-half_window_size", str(half_window_size),
        "-min_percent_NT_at_ends", str(min_percent_at_ends),
        "-out_trim_info", str(datainfo),
        "-out_NT", str(dataout),
    ]
    logger.info("trimming alignment")
    logger.debug(f"[{data.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    return dataout


def call_macse_export_alignment(data: Path, outprefix: str, codon_efs, codon_ifs, codon_fst, codon_ist, verbose, force):
    """..."""
    out_nt = outprefix.with_suffix(".msa.refined.nt.fa")
    out_aa = outprefix.with_suffix(".msa.refined.aa.fa")
    cmd = [
        BIN_MACSE, "-prog", "exportAlignment",
        "-align", str(data),
        "-codonForExternalFS", str(codon_efs),
        "-codonForInternalFS", str(codon_ifs),
        "-codonForFinalStop", str(codon_fst),
        "-codonForInternalStop", str(codon_ist),
        "-out_NT", str(out_nt),
        "-out_AA", str(out_aa),
    ]
    logger.debug(f"[{outprefix.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    return out_nt


def run_macse_refine(args):
    """..."""
    set_log_level(args.log_level)#, args.log_file)

    # check that macse is in PATH
    assert Path(BIN_MACSE).exists(), f"macse binary not found. Checked: {BIN_MACSE}"

    # check infiles
    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")
    if args.tree and not (args.tree.exists() and args.tree.is_file()):
        raise IOError(f"{args.tree} not found")

    # only one or the other allowed
    nargs = len([i for i in [args.exclude, args.subsample, args.tree] if i])
    if nargs > 1:
        raise ValueError("choose one of --exclude, --subample, or --tree")

    # ensure outpath and pdir exists
    args.outprefix = args.out
    if args.outprefix is None:
        args.outprefix = args.input
    args.outprefix.parent.mkdir(exist_ok=True)

    # bail out if final file exists
    result = args.outprefix.with_suffix(args.outprefix.suffix + ".msa.refined.nt.fa")
    if result.exists() and not args.force:
        logger.warning(f"[{args.input.name}] [skipping] {result} already exists. Using --force to overwrite")
        return

    # filter by minimum length
    data, filtered = filter_sequences(args.input, args.outprefix, args.exclude, args.subsample, args.tree, args.min_length, args.force)
    if args.refine_alignment:
        data = call_macse_refine_alignment(data, args.outprefix, args.force, args.max_iter_refine_alignment, args.verbose)
    if args.refine_alignment_if and filtered:
        data = call_macse_refine_alignment(data, args.outprefix, args.force, args.max_iter_refine_alignment, args.verbose)
    data = call_macse_trim_alignment(data, args.outprefix, args.aln_trim_window_size, args.aln_trim_ends_min_coverage, args.verbose, args.force)
    data = call_macse_export_alignment(data, args.outprefix, args.codon_int_fs, args.codon_ext_fs, args.codon_final_stop, args.codon_int_stop, args.verbose, args.force)

    # clean up tmp files
    suffices = [
        ".rmsa.nt.fa",
        ".rmsa.aa.fa",
        ".tmp.msa.trimmed.nt.fa",
    ]
    if not args.keep:
        for suffix in suffices:
            path = args.outprefix.with_suffix(args.outprefix.suffix + suffix)
            if path.exists():
                path.unlink()
    logger.info(f"[{args.input.name}] alignment written to {data}")


def main():
    parser = get_parser_macse_refine()
    args = parser.parse_args()
    run_macse_refine(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
