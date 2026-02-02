#!/usr/bin/env python

"""Trim and filter sequences for macse alignment.

Command
-------
$ twig macse-prep -i CDS -o OUT/ID.nt.fa

Output
------
- OUT/ID.trim.nt.fa
- OUT/ID.trim.nt.fa.info.tsv
"""

from typing import List
import re
import sys
import subprocess
from pathlib import Path
from loguru import logger
import pandas as pd
from twig.utils.logger_setup import set_log_level

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")


def call_macse_trim_non_homologous_fragments(
    cds_fasta: Path,
    min_homology_to_keep_seq: float,
    min_internal_homology_to_keep_seq: float,
    min_cov: int,
    min_trim_ext: int,
    min_trim_in: int,
    min_mem_length: int,
    outpath: str,
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
        "-out_NT", f"{outpath}",
        "-out_AA", f"{outpath}.tmp.aa.fa",
        "-out_trim_info", f"{outpath}.trim_info",
        "-out_mask_detail", f"{outpath}.tmp.trim.mask",
    ]
    logger.debug(f"[{outpath.name}] " + " ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    return proc.returncode


def filter_sequences(
    outpath: Path,
    exclude: List[str],
    subsample: List[str],
    min_length: int,
    min_count: int,
    # min_cov: int,
    force: bool,
) -> None:
    """Write fasta with filtered sequences and write info tsv.
    """
    # use trim file if present
    info_orig = outpath.with_suffix(outpath.suffix + ".trim_info")
    info_new = outpath.with_suffix(outpath.suffix + ".info.tsv")

    # keep track of sequences filtered-by
    f = {"homology": 0, "min_length": 0, "user": 0}

    # get table with lengths and homology scores
    info = {}
    with open(info_orig, 'r') as indata:
        _header = indata.readline()
        for line in indata.readlines():
            data = line.strip().split(";")
            name, length, kept, trimmed, itrimmed, phomology_internal, phomology_total, kept_seq = data
            info[name] = {
                "name": name,
                "homology_total": float(phomology_total),
                "homology_internal": float(phomology_internal),
                "length": int(length),
                "bp_kept": int(kept),
                "bp_trim": int(trimmed),
                "bp_trim_i": int(itrimmed),
                "discarded": bool(kept_seq == "false")
            }

    # report stats on trimmed
    data = pd.DataFrame(info).T.set_index("name")
    data['discarded'] = data['discarded'].astype(bool)
    f['homology'] = data.loc[data['discarded']].shape[0]
    low_homology = data.loc[data['discarded']]
    for i in low_homology.index:
        logger.debug(f"[{outpath.name}] {i} excluded by low_homology ({low_homology.loc[i, 'homology_internal']:.3f},{low_homology.loc[i, 'homology_total']:.3f})")

    # write this as a TSV and remove the awful semicolon-delimited macse info file
    data.to_csv(info_new, sep="\t")
    info_orig.unlink()

    # parse trimmed fasta file
    seqs = {}
    with open(outpath, 'r') as datain:
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
    passed = []
    for name, seq in seqs.items():

        # exclude if user-excluded
        if name in matched:
            logger.debug(f"[{outpath.name}] {name} excluded by user args")
            f["user"] += 1
            continue

        # exclude if too short
        seq = seqs[name]
        if len(seq) < min_length:
            logger.debug(f"[{outpath.name}] {name} excluded by min_length ({len(seq)})")
            f["min_length"] += 1
            continue

        passed.append((name, seq))

    # filter the locus?
    filtered = len(passed) < min_count

    # write output
    if not filtered:
        with open(outpath, 'w') as hout:
            for pname, pseq in sorted(passed, key=lambda x: x[0]):
                hout.write(f">{pname}\n{pseq}\n")

    # report filtering stats
    keys = sorted([i[0] for i in passed])
    mean_length = data.loc[keys, "bp_kept"].mean()
    mean_trimmed = data.loc[keys, "bp_trim"].mean()
    mean_homology = data.loc[keys, "homology_total"].mean()
    logger.info(f"[{outpath.name}] {len(info)} seqs -> {len(passed)} seqs, filtered by [min_homology={f['homology']}, min_length={f['min_length']}, user={f['user']}])")
    logger.info(f"[{outpath.name}] stats of retained sequences: mean_nt_length={mean_length:.2f}; mean_nt_trimmed={mean_trimmed:.2f}; mean_homology={mean_homology:.2f}")
    if filtered:
        logger.info(f"[{outpath.name}] locus did not pass 'min-count' filter ({len(passed)} < {min_count})")
        outpath.unlink()


def run_macse_trim(args):
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
    if args.outpath.is_dir():
        raise IOError("outpath should be a file path not dir")
    args.outpath.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    if args.outpath.exists() and not args.force:
        logger.warning(f"[{args.input.name}] [skipping] {args.outpath} already exists. Using --force to overwrite")
        return

    # trim sequences
    call_macse_trim_non_homologous_fragments(
        args.input,
        args.min_homology_full,
        args.min_homology_internal,
        args.min_homology_coverage,
        args.min_trim_length_homology_external,
        args.min_trim_length_homology_internal,
        args.mem_length,
        args.outpath,
        args.verbose,
        args.force,
    )

    # filter by minimum length
    filter_sequences(
        args.outpath,
        args.exclude,
        args.subsample,
        args.min_length,
        args.min_count,
        args.force,
    )

    # clean up tmp files
    suffices = [
        # ".trim_info",
        ".tmp.aa.fa",
        ".tmp.trim.mask",
    ]
    if not args.keep:
        for suffix in suffices:
            path = args.outpath.parent / (args.outpath.name + f"{suffix}")
            if path.exists():
                path.unlink()
    logger.info(f"[{args.outpath.name}] trimmed/filtered sequences written to {args.outpath}")


def main():
    from twig.cli.subcommands import get_parser_macse_prep
    parser = get_parser_macse_prep()
    args = parser.parse_args()
    run_macse_trim(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
