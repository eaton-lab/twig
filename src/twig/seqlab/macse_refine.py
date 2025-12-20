#!/usr/bin/env python

"""Subselect, refine, trim, and export a CDS alignment to CDS/AA.

1. subsample by user selection and re-align
2. trim edges
3. filter to remove short or terraced samples and re-align
4. export

If you call:
$ twig macse-refine -i CDS -o OUT/ID.refined

It will produce:
- OUT/ID.refined.aa.fa
- OUT/ID.refined.nt.fa
"""

from __future__ import annotations
from typing import Dict, List, Tuple
import re
import sys
import subprocess
from pathlib import Path
from loguru import logger
import toytree
import numpy as np
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")
MISSING = set(["-", "N", "n", "?", "."])


def parse_fasta_to_dict(path: str) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError("FASTA starts with sequence before any header.")
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def present_matrix(seqs: Dict[str, str], missing=MISSING) -> Tuple[List[str], np.ndarray]:
    names = list(seqs.keys())
    L = len(seqs[names[0]])
    for n in names[1:]:
        if len(seqs[n]) != L:
            raise ValueError("All sequences must have the same alignment length.")
    arr = np.empty((len(names), L), dtype=np.bool_)
    for i, n in enumerate(names):
        s = seqs[n]
        arr[i] = np.fromiter((c not in missing for c in s), count=L, dtype=np.bool_)
    return names, arr


def overlap_matrix(present: np.ndarray) -> np.ndarray:
    x = present.astype(np.uint16, copy=False)
    return x @ x.T  # (n,n), diagonal = present counts


def prune_to_pairwise_min_overlap(
    names: List[str],
    overlap: np.ndarray,
    min_ov: int,
    present_counts: np.ndarray | None = None,
    min_keep: int = 2,
) -> Tuple[List[str], List[Tuple[str, int, str]], bool]:
    """
    Keep a set S such that overlap(i,j) >= min_ov for EVERY i!=j in S.
    Drop one sample at a time until satisfied, or fail (kept < min_keep).

    Returns:
      kept_names,
      removed_log: (removed_name, min_overlap_at_removal, partner_name_at_that_min),
      success: bool
    """
    n = len(names)
    if overlap.shape != (n, n):
        raise ValueError("overlap must be (n,n) matching names length.")

    active = np.ones(n, dtype=bool)
    removed: List[Tuple[str, int, str]] = []

    while active.sum() >= min_keep:
        idx = np.flatnonzero(active)
        sub = overlap[np.ix_(idx, idx)].astype(int, copy=False)

        k = sub.shape[0]
        if k < 2:
            break

        offdiag = ~np.eye(k, dtype=bool)

        # global minimum overlap among remaining pairs (excluding diagonal)
        min_val = sub[offdiag].min()
        if min_val >= min_ov:
            return [names[i] for i in idx], removed, True

        # one worst pair (i,j) achieving min_val
        wi, wj = np.argwhere((sub == min_val) & offdiag)[0]
        gi, gj = int(idx[wi]), int(idx[wj])

        # decide which endpoint to drop:
        # drop the one with smaller mean overlap to others (tie-break: fewer present sites)
        row_sums = sub.sum(axis=1) - np.diag(sub)         # exclude self
        means = row_sums / (k - 1)

        drop_local = wi if means[wi] < means[wj] else wj

        if means[wi] == means[wj] and present_counts is not None:
            pci, pcj = present_counts[gi], present_counts[gj]
            drop_local = wi if pci <= pcj else wj

        drop_global = int(idx[drop_local])
        partner_global = gj if drop_global == gi else gi

        removed.append((names[drop_global], int(min_val), names[partner_global]))
        active[drop_global] = False

    kept = [names[i] for i in np.flatnonzero(active)]
    return kept, removed, False


def filter_by_selection(fasta: Path, outprefix: Path, exclude: List[str], subsample: List[str], subsample_tree: Path) -> None:
    """Write fasta with only one isoform per gene. When multple are present
    the one with greatest homology to other sequences is retained, with
    ties broken by length, and then order.
    """
    pre = fasta.name
    seqs = parse_fasta_to_dict(fasta)

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
    f = {"user": 0}
    keep = {}
    for name, seq in seqs.items():
        if name in matched:
            logger.debug(f"[{pre}] {name} excluded")
            f["user"] += 1
        else:
            keep[name] = seq

    # write output
    out = outprefix.with_suffix(outprefix.suffix + ".tmp.msa.nt.fa")
    if sum(f.values()):
        with open(out, 'w') as hout:
            for uname in keep:
                hout.write(f">{uname}\n{keep[uname]}\n")
        return out, True
    # nothing was filtered, return orig file and False
    return fasta, False


def filter_by_min_overlap(fasta: Path, min_ov: int, min_samples: int):
    seqs = parse_fasta_to_dict(fasta)
    names, arr = present_matrix(seqs)
    ovarr = overlap_matrix(arr)
    kept, removed, filtered = prune_to_pairwise_min_overlap(names, ovarr, min_ov, arr.sum(axis=1), min_samples)
    return kept, removed, filtered


def call_macse_refine_alignment(data: Path, outprefix: str, max_iter: int, verbose: bool):
    """Run Alignment step with default settings"""
    cmd = [
        BIN_MACSE, "-prog", "refineAlignment",
        "-align", str(data),
        "-out_NT", f"{outprefix}.tmp.rmsa.nt.fa",
        "-out_AA", f"{outprefix}.tmp.rmsa.aa.fa",
        "-max_refine_iter", str(max_iter),
    ]
    logger.info("refining alignment")
    logger.debug(f"[{data.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise Exception(proc.stderr)
    return outprefix.with_suffix(outprefix.suffix + ".tmp.rmsa.nt.fa")


def call_macse_trim_alignment(data: Path, outprefix: str, half_window_size: int, min_percent_at_ends: float, verbose: bool, ):
    """..."""
    cmd = [
        BIN_MACSE, "-prog", "trimAlignment",
        "-align", str(data),
        "-respect_first_RF_ON",
        "-half_window_size", str(half_window_size),
        "-min_percent_NT_at_ends", str(min_percent_at_ends),
        "-out_trim_info", f"{outprefix}.tmp.trimaln.info",
        "-out_NT", f"{outprefix}.tmp.trimaln.nt.fa",
    ]
    logger.info("trimming alignment")
    logger.debug(f"[{data.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise Exception(proc.stderr)
    return outprefix.with_suffix(outprefix.suffix + ".tmp.trimaln.nt.fa")


def call_macse_export_alignment(data: Path, outprefix: str, codon_efs, codon_ifs, codon_fst, codon_ist, verbose):
    """..."""
    out_nt = outprefix.with_suffix(outprefix.suffix + ".nt.fa")
    out_aa = outprefix.with_suffix(outprefix.suffix + ".aa.fa")
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
        raise Exception(proc.stderr)
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
    if args.outprefix.is_dir():
        raise IOError("outprefix must be a file path prefix, not a directory")
    args.outprefix.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    result = args.outprefix.with_suffix(args.outprefix.suffix + ".nt.fa")
    if result.exists() and not args.force:
        logger.warning(f"[{args.input.name}] [skipping] {result} already exists. Using --force to overwrite")
        return
    ####################################################################

    # REFINE ALIGNMENT ON SUBSAMPLE

    # ITERATIVELY REFINE
    data = args.input
    while 1:
        # FILTER TO OPTIONALLY SUBSAMPLE -> tmp.msa.nt.fa
        data, filtered = filter_by_selection(data, args.outprefix, args.exclude, args.subsample, args.tree)
        # mask orig input args for these
        args.exclude = args.subsample = args.tree = None

        # REFINE ALIGNMENT IF ANYTHING CHANGED
        if args.refine_alignment:
            data = call_macse_refine_alignment(data, args.outprefix, args.max_iter_refine_alignment, args.verbose)
        elif args.refine_alignment_if and filtered:
            data = call_macse_refine_alignment(data, args.outprefix, args.max_iter_refine_alignment, args.verbose)
        # TRIM EDGES
        data = call_macse_trim_alignment(data, args.outprefix, args.aln_trim_window_size, args.aln_trim_ends_min_coverage, args.verbose)
        # DROP MIN_OV SAMPLES
        kept, removed, success = filter_by_min_overlap(data, args.min_overlap, args.min_samples)
        args.exclude = [i[0] for i in removed]
        for name, minov, _ in removed:
            logger.info(f"[{args.input.name}] removed {name} by min-ov ({minov}) < min overlap")
        if not success:
            raise Exception("locus filtered")
        if not removed:
            break

    # WRITE TRANSLATED ALIGNMENT
    data = call_macse_export_alignment(data, args.outprefix, args.codon_int_fs, args.codon_ext_fs, args.codon_final_stop, args.codon_int_stop, args.verbose)

    # clean up tmp files
    suffices = [
        ".tmp.msa.nt.fa",
        ".tmp.rmsa.nt.fa",
        ".tmp.rmsa.aa.fa",
        ".tmp.trimaln.nt.fa",
        ".tmp.trimaln.info",
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




# # ---- example usage ----
# fasta = "../CSUBST/HOG_20018/msa.fa"
# seqs = parse_fasta_to_dict(fasta)
# names, pres = present_matrix(seqs)
# ov = overlap_matrix(pres)
# kept, removed, success = prune_to_pairwise_min_overlap(
#     names, ov, min_ov=50, present_counts=pres.sum(axis=1), min_keep=2
# )

# print("success:", success)
# print("kept:", len(kept), kept[:10], "..." if len(kept) > 10 else "")
# print("removed:")
# for n, mn, partner in removed:
#     print(f"  - {n} (min overlap={mn} with {partner})")

