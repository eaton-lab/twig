#!/usr/bin/env python

"""Prep sequences for macse alignment (trim,filter,iso-collapse)

"""

from typing import List
import re
import sys
import textwrap
import subprocess
from collections import defaultdict
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
import pandas as pd
# from twig.utils.parallel import run_pipeline  # , run_with_pool
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")
ISOFORM_REGEX_DEFAULT = r"^([^|]+)\|.*?__(.+?)_i\d+"
# group 1 (shared part up until |)
# group 2 (after __ and up until first _i)


KWARGS = dict(
    prog="macse-prep",
    usage="macse-prep -i CDS -o OUTDIR [options]",
    help="prepare CDS for macse alignment (trim,filter,iso-collapse)",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=120, max_help_position=120),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | macse-prep: prepare CDS for alignment (trim, filter, iso-collapse)
        -------------------------------------------------------------------
        | Macse ...
        | This implements `macse -prog trimNonHomologousFragments` and
        | additional steps to ...
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        $ twig macse-prep -i CDS -o OUT -p TEST
        $ twig macse-prep -i CDS -o OUT -mh 0.1 -mi 0.5 -ti 50 -te 50 -mc 15 -c 20
        $ twig macse-prep -i CDS -o OUT -mh 0.3 -mi 0.8 -ti 25 -te 25 -mc 15 -c 20
        $ twig macse-prep -i CDS -o OUT -mh 0.5 -mc 10 -k -xa -ml 200 -e '^sppA.*'
        $ twig macse-prep -i CDS -o OUT -s '^sppA.*'

        # run parallel jobs on many cds files
        $ parallel -j 10 'twig macse-prep -i {} -o OUT -p {/.}' ::: CDS/*.fa
    """)
)


def get_parser_macse_prep(parser: ArgumentParser | None = None) -> ArgumentParser:
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
    parser.add_argument("-p", "--prefix", type=str, metavar="str", help="optional outfile prefix. If None the cds filename is used")
    parser.add_argument("-e", "--exclude", type=str, metavar="str", nargs="*", help="optional names or glob to exclude one or more sequences")
    parser.add_argument("-s", "--subsample", type=str, metavar="str", nargs="*", help="optional names or glob to include only a subset sequences")
    # options
    parser.add_argument("-ml", "--min-length", type=int, metavar="int", default=0, help="min nt sequence length after trimming [%(default)s]")
    parser.add_argument("-hf", "--min-homology-full", type=float, metavar="float", default=0.1, help="min homology required w/ >=mc others across full sequence [%(default)s]")
    parser.add_argument("-hi", "--min-homology-internal", type=float, metavar="float", default=0.5, help="min homology required w/ >=mc others in the internal sequence [%(default)s]")
    parser.add_argument("-hc", "--min-homology-coverage", type=int, metavar="int", default=3, help="min samples a seq must share homology with at >= mh and mi [%(default)s]")
    parser.add_argument("-ti", "--min-trim-length-homology-internal", type=int, metavar="int", default=50, help="trim fragments w/ len <ti and homology <mi w/ <mc sequences [%(default)s]")
    parser.add_argument("-te", "--min-trim-length-homology-external", type=int, metavar="int", default=50, help="trim fragments w/ len <tx and homology <mh w/ <mc sequences [%(default)s]")
    parser.add_argument("-mm", "--mem-length", type=int, metavar="int", default=6, help="homology is the prop of aa Maximum Exact Matches of this length [%(default)s]")
    parser.add_argument("-is", "--isoform-regex", type=re.compile, metavar="str", default=ISOFORM_REGEX_DEFAULT, help="regex used to group isoform sequences ['%(default)s']")
    # parser.add_argument("-ac", "--aln-trim-ends-min-coverage", type=float, metavar="float", default=0.4, help="trim alignment edges to where a min percent of samples have data [%(default)s]")
    # parser.add_argument("-as", "--aln-trim-window-size", type=int, metavar="int", default=5, help="trim alignment using a sliding 'half_window_size' defined as ... [%(default)s]")
    # parser.add_argument("-xa", "--skip-alignment", action="store_true", help="skip alignment step and only trim/filter sequences")

    # choose one or more
    # parser.add_argument("-xt", "--skip-trim-and-filter", action="store_true", help="skip trim and filter step")
    parser.add_argument("-xi", "--skip-isoform-collapse", action="store_true", help="skip isoform collapse step")

    # others
    parser.add_argument("-v", "--verbose", action="store_true", help="print macse progress info to stderr")
    parser.add_argument("-f", "--force", action="store_true", help="overwrite existing result files in outdir")
    parser.add_argument("-k", "--keep", action="store_true", help="keep tmp files (for debugging)")
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def call_macse_trim_non_homologous_fragments(
    cds_fasta: Path,
    min_homology_to_keep_seq: float,
    min_internal_homology_to_keep_seq: float,
    min_cov: int,
    min_trim_ext: int,
    min_trim_in: int,
    min_mem_length: int,
    outdir: Path,
    prefix: str,
    # kwargs: ...,
    verbose: bool,
    force: bool,
):
    """Run macse 'trimNonHomologousFragments' on a CDS fasta.

    Parameters
    ----------
    ...

    Command
    -------
    >>> macse -prog trim... -seq CDS.fa ... -out_NT CDS.trim.fa
    """
    cmd = [
        BIN_MACSE, "-prog", "trimNonHomologousFragments",
        "-seq", str(cds_fasta),
        "-min_homology_to_keep_seq", str(min_homology_to_keep_seq),
        "-min_internal_homology_to_keep_seq", str(min_internal_homology_to_keep_seq),
        "-min_cov", str(min_cov),
        "-min_trim_ext", str(min_trim_ext),
        "-min_trim_in", str(min_trim_in),
        "-min_MEM_length", str(min_mem_length),
        "-out_NT", str(outdir / f"{prefix}.trim"),
        "-out_AA", str(outdir / f"{prefix}.tmp.aa.fa"),
        "-out_trim_info", str(outdir / f"{prefix}.trim_info"),
        "-out_mask_detail", str(outdir / f"{prefix}.tmp.trim_mask"),
    ]
    logger.debug(f"[{prefix}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    return proc.returncode


def filter_sequences(
    outprefix: str,
    isoform_regex: re.compile,
    exclude: List[str],
    subsample: List[str],
    min_length: int,
    force: bool,
) -> None:
    """Write fasta with only one isoform per gene. When multple are present
    the one with greatest homology to other sequences is retained, with
    ties broken by length, and then order.
    """
    # use trim file if present
    trim = outprefix.with_suffix(outprefix.suffix + ".trim")
    tnfo = outprefix.with_suffix(outprefix.suffix + ".trim_info")
    tout = outprefix.with_suffix(outprefix.suffix + ".nt.fa")

    # keep track of filtered-by
    f = {"homology": 0, "min_length": 0, "isoform": 0, "user": 0}

    # get table with lengths and homology scores
    info = {}
    with open(tnfo, 'r') as indata:
        _header = indata.readline()
        for line in indata.readlines():
            data = line.strip().split(";")
            name, length, kept, trimmed, itrimmed, phomology_internal, phomology_total, kept_seq = data
            info[name] = {
                "name": name,
                "length": int(length),
                "bp_kept": int(kept),
                "bp_trim": int(trimmed),
                "bp_trim_i": int(itrimmed),
                "homology_internal": float(phomology_internal),
                "homology_total": float(phomology_total),
                "discarded": bool(kept_seq == "false")
            }

    # report stats on trimmed
    data = pd.DataFrame(info).T.set_index("name")
    data['discarded'] = data['discarded'].astype(bool)
    f['homology'] = data.loc[data['discarded']].shape[0]
    low_homology = data.loc[data['discarded']]
    for i in low_homology.index:
        logger.debug(f"[{trim.name}] {i} excluded by low_homology ({low_homology.loc[i, 'homology_internal']:.3f},{low_homology.loc[i, 'homology_total']:.3f})")

    # parse trimmed fasta file
    seqs = {}
    with open(trim, 'r') as datain:
        for line in datain.readlines():
            if line:
                if line.startswith(">"):
                    uname = line.strip()[1:]
                    seqs[uname] = ""
                else:
                    seqs[uname] += line.strip()

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

    # group sequences by isoform regex
    pat = isoform_regex
    groups = defaultdict(list)
    for name, seq in seqs.items():

        # exclude if user-excluded
        if name in matched:
            logger.debug(f"[{trim.name}] {name} excluded by user args")
            f["user"] += 1
            continue

        # exclude if too short
        seq = seqs[name]
        if len(seq) < min_length:
            logger.debug(f"[{trim.name}] {name} excluded by min_length ({len(seq)})")
            f["min_length"] += 1
            continue

        # assign to a group
        m = pat.match(name)
        if m:
            group_key = m.group(1)  # which capture group from regex pattern
        else:
            # no match, treat as its own group
            group_key = name

        # store grouped data
        idata = (name, seq, info[name]['homology_total'])
        groups[group_key].append(idata)

    # ...
    collapsed = {}
    for group_key, members in groups.items():
        # members: list of (name, idict)
        # keep only one per group by [highest homology_total, longest_length, or first]
        best_name, best_seq, best_score = max(
            members,
            key=lambda x: (x[2], len(x[1]))  # (score, length)
        )
        collapsed[best_name] = (best_seq, best_score)
    f["isoform"] = len(seqs) - f["min_length"] - f["user"] - len(collapsed)

    # report isoform collapsed
    for group in groups:
        for name, *_ in groups[group]:
            if name not in collapsed:
                logger.debug(f"[{trim.name}] {name} excluded by isoform collapse")

    # write output
    with open(tout, 'w') as hout:
        for uname in sorted(collapsed):
            hout.write(f">{uname}\n{seqs[uname]}\n")

    # report filtering stats
    keys = list(collapsed)
    mean_length = data.loc[keys, "bp_kept"].mean()
    mean_trimmed = data.loc[keys, "bp_trim"].mean()
    mean_homology = data.loc[keys, "homology_total"].mean()
    logger.info(f"[{trim.name}] {len(info)} seqs -> {len(collapsed)} seqs, filtered by [min_homology={f['homology']}, min_length={f['min_length']}, user={f['user']}, isoform={f['isoform']}])")
    logger.info(f"[{trim.name}] stats of retained sequences: mean_nt_length={mean_length:.2f}; mean_nt_trimmed={mean_trimmed:.2f}; mean_homology={mean_homology:.2f}")



def run_macse_prep(args):
    """..."""
    set_log_level(args.log_level)#, args.log_file)

    # check that macse is in PATH
    assert Path(BIN_MACSE).exists(), f"macse binary not found. Checked: {BIN_MACSE}"

    # check in files
    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")

    # only one or the other allowed
    if args.exclude and args.subsample:
        raise ValueError("choose one of --exclude or --subample, but not both")

    # ensure outpath and pdir exists
    args.outprefix = args.out
    if args.outprefix is None:
        args.outprefix = args.input
    args.outprefix.parent.mkdir(exist_ok=True)

    # bail out if final file exists
    result = args.outprefix.with_suffix(args.outprefix.suffix + ".nt.fa")
    if result.exists() and not args.force:
        logger.warning(f"[{args.input.name}] [skipping] {result} already exists. Using --force to overwrite")
        return

    # if skip isoform collapse then set isoform grouper to arbitrary str
    if args.skip_isoform_collapse:
        args.isoform_regex = re.compile("@@@@@")

    # trim sequences
    call_macse_trim_non_homologous_fragments(
        args.input,
        args.min_homology_full,
        args.min_homology_internal,
        args.min_homology_coverage,
        args.min_trim_length_homology_external,
        args.min_trim_length_homology_internal,
        args.mem_length,
        args.outprefix,
        args.verbose,
        args.force,
    )

    # filter by minimum length
    filter_sequences(
        args.outprefix,
        args.isoform_regex,
        args.exclude,
        args.subsample,
        args.min_length,
        args.force,
    )

    # clean up tmp files
    suffices = [
        ".trim",
        ".trim_info",
        ".tmp.aa.fa",
        ".tmp.trim_mask",
    ]
    if not args.keep:
        for suffix in suffices:
            path = args.outdir / f"{args.prefix}{suffix}"
            if path.exists():
                path.unlink()
    logger.info(f"[{args.prefix}] trimmed/filtered sequences written to {args.outdir}/{args.prefix}.nt.fa")


def main():
    parser = get_parser_macse_prep()
    args = parser.parse_args()
    run_macse_prep(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
