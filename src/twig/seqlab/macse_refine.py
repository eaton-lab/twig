#!/usr/bin/env python

"""Subselect, refine, trim, and export a CDS alignment to CDS/AA.

"""

from typing import List
import re
import sys
import subprocess
from pathlib import Path
from loguru import logger
import toytree
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")


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
    out_nt = outprefix.with_suffix(outprefix.suffix + ".msa.refined.nt.fa")
    out_aa = outprefix.with_suffix(outprefix.suffix + ".msa.refined.aa.fa")
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
    from ..cli.subcommands import get_parser_macse_refine
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
