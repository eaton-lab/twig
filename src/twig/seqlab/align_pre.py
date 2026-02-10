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
from dataclasses import dataclass
from loguru import logger
import numpy as np
from twig.utils.logger_setup import set_log_level
from twig.utils.exceptions import TwigError

BIN = Path(sys.prefix) / "bin"
BIN_MACSE = str(BIN / "macse")


def parse_bool(value: str) -> bool:
    """Parse 'true'/'false' tokens from MACSE trim_info output."""
    norm = str(value).strip().lower()
    if norm == "true":
        return True
    if norm == "false":
        return False
    raise ValueError(f"expected 'true' or 'false', got: {value!r}")


@dataclass
class TrimInfo:
    name: str
    initial_length: int
    bp_kept: int
    bp_trim: int
    bp_trim_informative: int
    homology_total: float
    homology_internal: float
    passed: bool

    @property
    def informative_bp(self) -> int:
        return max(0, self.initial_length - self.bp_trim_informative)

    @property
    def effective_informative_homolog_bp(self) -> float:
        # primary score (bigger is better)
        return self.informative_bp * float(self.homology_internal)


def get_scores(
    fasta_path: Path,
    use_trim_info: bool = True,
    warn_on_existing_trim_info: bool = False,
):
    """Return dict mapping {sequence: dict} with trim info data."""
    info_file = None
    if use_trim_info:
        info_file = fasta_path.with_suffix(fasta_path.suffix + ".trim_info")
        # MACSE trim writes {outpath}.trim_info while sequences are read from
        # {outpath}.tmp.trimmed.nt.fa. Detect this naming pattern and resolve it.
        if (
            not info_file.exists()
            and fasta_path.name.endswith(".tmp.trimmed.nt.fa")
        ):
            base_name = fasta_path.name[:-len(".tmp.trimmed.nt.fa")]
            candidate = fasta_path.parent / f"{base_name}.trim_info"
            if candidate.exists():
                info_file = candidate
    scores = {}
    if info_file and info_file.exists():
        if warn_on_existing_trim_info:
            logger.warning(
                f"using existing trim info file without running trimming: {info_file}"
            )
        with open(info_file, 'r') as datain:
            for line in datain.readlines():
                if line:
                    # skip header
                    if line.startswith("seqName"):
                        continue
                    score = line.strip().split(";")
                    info = TrimInfo(
                        name=str(score[0]),
                        initial_length=int(score[1]),
                        bp_kept=int(score[2]),
                        bp_trim=int(score[3]),
                        bp_trim_informative=int(score[4]),
                        homology_internal=float(score[5]),
                        homology_total=float(score[6]),
                        passed=parse_bool(score[7]),
                    )
                    scores[str(info.name)] = info
    else:
        logger.warning("no info file detected, falling back to selecting the longest length isoforms")
        seqs = parse_fasta(fasta_path)
        for name, seq in seqs.items():
            length = len(seq.replace("-", "")) # in case it is aligned
            scores[name] = TrimInfo(
                name=name,
                initial_length=length,
                bp_kept=length,
                bp_trim=0,
                bp_trim_informative=0,
                homology_total=1.0,
                homology_internal=1.0,
                passed=True,
            )
    # logger.debug(f"scores={scores}")
    return scores


def get_patterns_dict_from_split(split: str, split_file: Path, imap: dict[str, list[str]]):
    """Return dict mapping {group-name: regex-pattern} for each imap group."""
    patterns = {}
    if split_file:
        if split_file.is_file() and split_file.exists():
            with open(split_file, 'r') as datain:
                for line in datain.readlines():
                    if line:
                        # it should be fine if header is present, just a key with no matches...

                        # parse regex file line
                        name, delim, *x = line.strip().split()
                        delim = delim.strip("'").strip('"')
                        idx = int(x[0]) if x else 1

                        # store the {name: (regex, key)}
                        patterns[str(name)] = (delim, idx)

            # if using split_file then a pattern must be present for every sequence
            if len(imap) > 1:
                diff = set(imap) - set(patterns)
                if any(diff):
                    raise TwigError(f"split-file tsv does not contain the following (-d) delimited names from the input: [{diff}]")
            else:
                diff = set(imap["ALL-SEQUENCES"]) - set(patterns)
                if any(diff):
                    raise TwigError(f"split-file tsv does not contain the following sequence names from the input: [{diff}]")
        else:
            raise TwigError(f"split file is malformed: {split_file}")
    else:
        if split:
            patterns = {i: (str(split[0]), int(split[1])) for i in imap}
        else:
            patterns = {i: ("", 0) for i in imap}
    return patterns


def get_groups_from_split(split: list[str, int], split_file: Path, imap: dict[str, list[str]]):

    # get patterns dict {str: re.compile} of {group: regex}
    patterns = get_patterns_dict_from_split(split, split_file, imap)

    # iterate over (taxon-group, regular expression) tuples
    observed = set()
    groups = {}
    for taxon_group, (delim, idx) in patterns.items():

        # select only the sequences in this taxon_group
        subgroup = imap.get(taxon_group)
        logger.warning((taxon_group, delim, idx, subgroup))

        # skip if this locus does not contain this subgroup
        if not subgroup:
            continue

        # apply this groups specific regular expression to find isoforms
        for sname in subgroup:

            # sanity check
            if sname in observed:
                raise ValueError(f"{sname} matched to more than one group. Check your delim pattern.")
            observed.add(sname)

            # get the group key (up one group)
            group_key = None
            if delim:
                parts = sname.split(delim)
                parts_left_of_nth = parts[:idx]
                group_key = delim.join(parts_left_of_nth)

            # if no key split from name then use full name (no collapse)
            if not group_key:
                group_key = sname

            # store group results
            if group_key in groups:
                groups[group_key].append(sname)
            else:
                groups[group_key] = [sname]
    return groups


def parse_fasta(inpath: Path) -> dict[str, str]:
    seqs = {}
    with open(inpath, 'r') as datain:
        for line in datain.readlines():
            if line:
                if line.startswith(">"):
                    uname = line.strip()[1:]
                    seqs[uname] = ""
                else:
                    seqs[uname] += line.strip()
    return seqs


def filter_by_selection(nt_fa: Path, outpath: Path, subsample: list[str], exclude: list[str]) -> tuple[Path, list[str]]:
    """Subsample or remove sequences before trimming."""
    # parse trimmed fasta file
    handle = outpath.with_suffix(outpath.suffix + ".tmp")
    filtered = []
    seqs = parse_fasta(nt_fa)

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

    # exclude if user-excluded
    if matched:
        with open(handle, 'w') as hout:
            for name in sorted(seqs):
                if name in matched:
                    logger.debug(f"{name} excluded by user args")
                    filtered.append(name)
                else:
                    hout.write(f">{name}\n{seqs[name]}\n")
        return handle, filtered
    return nt_fa, filtered


def parse_imap(seqs: dict[str, str], delim: str, delim_idxs: list[int], delim_join: str) -> dict[str, list[str]]:
    if delim is None:
        imap = {i: [i] for i in sorted(seqs)}
        # imap = {"ALL-SEQUENCES": sorted(seqs)}
    else:
        imap = {}
        for name in seqs:
            parts = name.split(delim)
            new_name = delim.join([parts[i] for i in delim_idxs])
            if new_name not in imap:
                imap[new_name] = [name]
            else:
                imap[new_name].append(name)
    # logger.debug(f"imap = {imap}")
    return imap


def call_macse_trim_non_homologous_fragments(
    cds_fasta: Path,
    min_homology_to_keep_seq: float,
    min_internal_homology_to_keep_seq: float,
    min_cov: int,
    min_trim_ext: int,
    min_trim_in: int,
    min_mem_length: int,
    outpath: str,
    verbose: bool,
    # kwargs: ...,
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
        "-out_NT", f"{outpath}.tmp.trimmed.nt.fa",
        "-out_AA", f"{outpath}.tmp.trimmed.aa.fa",
        "-out_trim_info", f"{outpath}.trim_info",
        "-out_mask_detail", f"{outpath}.tmp.trimmed.mask",
    ]
    logger.info(f"trimming {cds_fasta}")
    logger.debug(" ".join(cmd))
    if verbose:
        proc = subprocess.run(cmd, stderr=sys.stderr, check=True)
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode:
        raise subprocess.CalledProcessError(proc.stderr)
    return outpath.with_suffix(outpath.suffix + ".tmp.trimmed.nt.fa")


def filter_sequences(
    nt_fa: Path,
    outpath: Path,
    delim: str,
    delim_idxs: List[int],
    delim_join: str,
    split: tuple[str, int],
    split_file: Path,
    min_length: int,
    min_count: int,
    pre_filtered_count: int = 0,
    use_trim_info: bool = True,
    warn_on_existing_trim_info: bool = False,
) -> None:
    """Write fasta with only one isoform per gene. When multple are present
    the one with greatest homology to other sequences is retained, with
    ties broken by length, and then order.
    """
    # keep track of sequences filtered-by
    f = {"homology": 0, "min_length": 0, "isoform": 0, "user": pre_filtered_count}

    # -----------------------------------------------------------
    # parse input sequence fasta to dict
    # {taxonA-gene1: AAA..., taxonA-gene2: AAA..., ...}
    seqs = parse_fasta(nt_fa)

    # -----------------------------------------------------------
    # assign sequences to subgroups to either split/regex globally or to each taxon
    # parse taxon-gene names to dict {}
    # orig = [taxonA-gene1, taxonA-gene2, taxonB-gene1, taxonB-gene2]
    # imap = {taxonA: [taxonA-gene1, taxonA-gene2], taxonB: [...], ...}
    imap = parse_imap(seqs, delim, delim_idxs, delim_join)

    # -----------------------------------------------------------
    # parse sequence info (e.g., homology) if present, else use seq length
    # scores = {taxonA-gene1: 0.9, taxonA-gene2: 0.8, ...}
    scores = get_scores(
        nt_fa,
        use_trim_info=use_trim_info,
        warn_on_existing_trim_info=warn_on_existing_trim_info,
    )

    # -----------------------------------------------------------
    # build dict mapping {name: regex} for every name in imap
    groups = get_groups_from_split(split, split_file, imap)

    # select a single best in each isoform group
    collapsed = {}
    for group_key, members in groups.items():

        # get scores for group members
        sgroup = [scores[i] for i in members]

        # sort by these keys in this order where higher is better for all
        best_iso = max(
            sgroup,
            key=lambda x: (
                x.passed,
                x.effective_informative_homolog_bp,
                x.homology_internal,
                x.informative_bp,
                x.bp_kept,
                x.name,
            )
        )

        # report if iso's were collapsed ---------------------------
        kept = best_iso.name[len(group_key):]
        pruned = [i for i in members if i != best_iso.name]

        # get the leftover parts of names
        rests = [p[len(group_key):] for p in pruned]

        # shorten rests if very long
        if len(rests) > 5:
            rests = rests[:5] + ["..."]

        # report
        if len(members) > 1:
            logger.info(f"pruned {len(members) - 1} isoforms: group='{group_key}', kept='{kept}', pruned={rests}")
        else:
            logger.debug(f"pruned {len(members) - 1} isoforms: group='{group_key}', kept='{kept}', pruned={rests}")
        # -----------------------------------------------------------
        f["isoform"] += len(members) - 1

        name = best_iso.name
        seq = seqs[name]

        # min length filter
        if len(seq) < min_length:
            logger.info(f"{name} excluded by min_length ({len(seq)})")
            f["min_length"] += 1
            continue

        # min homology filter
        if not scores[name].passed:
            logger.info(f"{name} excluded by min_homology ({(scores[name].homology_total, scores[name].homology_internal)})")
            f["homology"] += 1
            continue

        # store passed result for this group
        collapsed[name] = seq

    # filter locus by min count
    if len(collapsed) < min_count:
        logger.warning(f"locus has fewer than min_count sequences ({len(collapsed)}). No result written.")
        return

    # if not filtered:
    with open(outpath, 'w') as hout:
        for uname in sorted(collapsed):
            hout.write(f">{uname}\n{seqs[uname]}\n")

    # report filtering stats
    mean_length = np.mean([scores[i].bp_kept for i in collapsed])
    mean_trimmed = np.mean([scores[i].bp_trim for i in collapsed])
    mean_homology = np.mean([scores[i].homology_internal for i in collapsed])
    logger.info(f"{len(seqs)} seqs -> {len(collapsed)} seqs, filtered by [min_homology={f['homology']}, min_length={f['min_length']}, user={f['user']}])")
    logger.info(f"stats of retained sequences: mean_nt_length={mean_length:.2f}; mean_nt_trimmed={mean_trimmed:.2f}; mean_internal_homology={mean_homology:.2f}")
    # if filtered:
    #     logger.info(f"[{outpath.name}] locus did not pass 'min-count' filter ({len(passed)} < {min_count})")
    #     outpath.unlink()  # keep zero size file?


def run_align_pre(args):
    """..."""
    set_log_level(args.log_level)#, args.log_file)

    # check that macse is in PATH
    assert Path(BIN_MACSE).exists(), f"macse binary not found. Checked: {BIN_MACSE}"

    # check for infile
    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")

    # check mutually exclusive args
    if args.exclude and args.subsample:
        raise ValueError("choose one of --exclude or --subample, but not both")
    if args.isoform_group_split and args.isoform_group_rules:
        raise TwigError("choose one of --isoform-group-split or --isoform-group-rules, but not both")

    # check for outpath outpath and pdir exists
    if args.outpath.is_dir():
        raise IOError("outpath should be a file path not dir")
    args.outpath.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    if args.outpath.exists() and not args.force:
        logger.warning(f"[skipping] {args.outpath} already exists. Use --force to overwrite")
        return

    # filter sequences by user selection
    nt_fa, f_user = filter_by_selection(args.input, args.outpath, args.subsample, args.exclude)

    # trim sequences
    if not args.skip_macse_trim:
        nt_fa = call_macse_trim_non_homologous_fragments(
            nt_fa,
            args.min_homology_full,
            args.min_homology_internal,
            args.min_homology_coverage,
            args.min_trim_length_homology_external,
            args.min_trim_length_homology_internal,
            args.mem_length,
            args.outpath,
            args.verbose,
        )

    # filter sequences by prune isoforms, min homology, min length, min cov
    filter_sequences(
        nt_fa,
        args.outpath,
        args.taxon_delim,
        args.taxon_delim_idxs,
        args.taxon_delim_join,
        args.isoform_group_split,
        args.isoform_group_rules,
        args.min_retained_nt_length,
        args.min_retained_seqs,
        len(f_user),
        use_trim_info=True,
        warn_on_existing_trim_info=args.skip_macse_trim,
    )

    # clean up tmp files
    suffices = [
        ".tmp",
        ".tmp.trimmed.nt.fa",
        ".tmp.trimmed.aa.fa",
        ".tmp.trimmed.mask"
    ]
    if not args.keep:
        for suffix in suffices:
            path = args.outpath.parent / (args.outpath.name + f"{suffix}")
            if path.exists():
                path.unlink()
    logger.info(f"trimmed/filtered sequences written to {args.outpath}")


def main():
    from twig.cli.subcommands import get_parser_align_pre
    parser = get_parser_align_pre()
    args = parser.parse_args()
    run_align_pre(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
