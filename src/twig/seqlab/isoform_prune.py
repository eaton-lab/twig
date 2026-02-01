#!/usr/bin/env python

"""Select ISOFORMS for downstream analysis.

For isoform selection we care first about how much usable homologous
signal a sequence will contribute to the downstream alignment/tree,
moreso than how homologous any remaining fragment is.

We pick the isoform that maximizes effective informative homologous bp, i.e.

score = (length − bp_trim_i) × homology_internal

where:
    - (length - bp_trim_i) = how many informative bases remain (non-N)
    - homology_internal (1.0 best) downweights sequences that only retain a small homologous fraction

"""

from typing import List
import re
from pathlib import Path
from loguru import logger
import pandas as pd
from dataclasses import dataclass
from twig.utils.logger_setup import set_log_level
from twig.utils.exceptions import TwigError


@dataclass
class TrimInfo:
    name: str
    homology_total: float
    homology_internal: float
    length: int
    bp_kept: int
    bp_trim: int
    bp_trim_i: int
    discarded: bool

    @property
    def informative_bp(self) -> int:
        return max(0, self.length - self.bp_trim_i)

    @property
    def effective_informative_homolog_bp(self) -> float:
        # primary score (bigger is better)
        return self.informative_bp * float(self.homology_internal)


def get_patterns_dict_from_regex(regex: str, regex_file: Path, imap: dict[str, list[str]]):
    """Return dict mapping {group-name: regex-pattern} for each imap group."""
    patterns = {}
    if regex_file:
        if regex_file.is_file() and regex_file.exists():
            with open(regex_file, 'r') as datain:
                for line in datain.readlines():
                    if line:

                        # [TODO] detect and handle header??
                        # it should be fine if present, just a key with no matches...

                        # parse regex file line
                        name, pattern, *x = line.strip().split("\t")
                        # check for optional group key
                        gkey = x[0] if x else 1
                        # try to cast to int if possible
                        try:
                            gkey = int(gkey)
                        except Exception:
                            gkey = str(gkey)

                        # store the {name: (regex, key)}
                        patterns[str(name)] = (re.compile(pattern), gkey)

            # if using regex_file then a pattern must be present for every sequence
            if len(imap) > 1:
                diff = set(imap) - set(patterns)
                if any(diff):
                    raise ValueError(f"regex-file tsv does not contain the following (-d) delimited names from the input: [{diff}]")
            else:
                diff = set(imap['ALL-SEQUENCES']) - set(patterns)
                if any(diff):
                    raise ValueError(f"regex-file tsv does not contain the following sequence names from the input: [{diff}]")
        else:
            raise ValueError(f"regex file is malformed: {regex_file}")
    else:
        patterns = {'ALL-SEQUENCES': (regex, 1)}
    return patterns


def get_groups_from_regex(regex: str, regex_file: Path, imap: dict[str, list[str]]):

    # get patterns dict {str: re.compile} of {group: regex}
    patterns = get_patterns_dict_from_regex(regex, regex_file, imap)

    # iterate over (taxon-group, regular expression) tuples
    observed = set()
    groups = {}
    for taxon_group, (pattern, gkey) in patterns.items():

        # select only the sequences in this taxon_group
        subgroup = imap.get(taxon_group)

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
            m = pattern.match(sname)
            if m:
                group_key = m.group(gkey)  # which capture group from regex pattern
            else:
                # no match, treat as its own group
                group_key = sname

            # store grouped data
            # logger.debug(f"{sname} -- {subgroup} -- {group_key}")
            if group_key in groups:
                groups[group_key].append(sname)
            else:
                groups[group_key] = [sname]
    return groups


def get_patterns_dict_from_split(split: str, split_file: Path, imap: dict[str, list[str]]):
    """Return dict mapping {group-name: regex-pattern} for each imap group."""
    patterns = {}
    if split_file:
        if split_file.is_file() and split_file.exists():
            with open(split_file, 'r') as datain:
                for line in datain.readlines():
                    if line:

                        # [TODO] detect and handle header??
                        # it should be fine if present, just a key with no matches...

                        # parse regex file line
                        name, delim, *x = line.strip().split("\t")
                        idx = int(x[0]) if x else 1

                        # store the {name: (regex, key)}
                        patterns[str(name)] = (delim, idx)

            # if using regex_file then a pattern must be present for every sequence
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
        patterns = {'ALL-SEQUENCES': (str(split[0]), int(split[1]))}
    logger.debug(patterns)
    # value check
    if any(i[1] == 0 for i in patterns.values()):
        raise TwigError("split index selector(s) cannot be zero. Use >1 or <1 to select the nth occurrence from start or end")
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
            # logger.debug(sname)
            # logger.debug((delim, idx))
            # logger.debug(sname.split(delim))
            parts = sname.split(delim)
            parts_left_of_nth = parts[:idx]
            group_key = delim.join(parts_left_of_nth)
            if group_key in groups:
                groups[group_key].append(sname)
            else:
                groups[group_key] = [sname]
    return groups


def get_scores(seqs: dict[str, str], info_file: Path):
    """Return dict mapping {sequence: dict} with trim info data"""
    scores = {}
    if info_file:
        if info_file.is_file() and info_file.exists():
            with open(info_file, 'r') as datain:
                for line in datain.readlines():
                    if line:
                        # use try/except to skip headers if present. Membership is checked below.
                        try:
                            score = line.strip().split("\t")
                            info = TrimInfo(
                                str(score[0]),
                                float(score[1]),
                                float(score[2]),
                                int(score[3]),
                                int(score[4]),
                                int(score[5]),
                                bool(score[6]),
                            )
                            scores[str(info.name)] = info
                        except KeyError:
                            pass

            # if using scores then a score must be present for every sequence
            diff = set(seqs) - set(scores)
            if any(diff):
                raise ValueError(f"scores tsv does not contain the following sequence names from the input: [{diff}]")
        else:
            raise ValueError(f"info file is malformed: {info_file}")
    else:
        for name, seq in seqs.items():
            length = len(seq.replace("-", ""))
            scores[name] = TrimInfo(
                name=name,
                homology_total=1.0,
                homology_internal=1.0,
                length=length,
                bp_kept=length,
                bp_trim=0,
                bp_trim_i=0,
                discarded=False,
            )
    # logger.debug(f"scores={scores}")
    return scores


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


def parse_imap(seqs: dict[str, str], delim: str, delim_idxs: list[int], delim_join: str) -> dict[str, list[str]]:
    if delim is None:
        imap = {"ALL-SEQUENCES": sorted(seqs)}
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


def prune_sequences(
    inpath: Path,
    outpath: Path,
    delim: str,
    delim_idxs: List[int],
    delim_join: str,
    info_file: Path,
    regex: re.compile,
    regex_file: Path,
    split: tuple[str, int],
    split_file: Path,
) -> None:
    """Write fasta with only one isoform per gene. When multple are present
    the one with greatest homology to other sequences is retained, with
    ties broken by length, and then order.
    """
    # -----------------------------------------------------------
    # parse input sequence fasta to dict
    # {taxonA-gene1: AAA..., taxonA-gene2: AAA..., ...}
    seqs = parse_fasta(inpath)

    # -----------------------------------------------------------
    # assign sequences to subgroups to either split/regex globally or to each taxon
    # parse taxon-gene names to dict {}
    # orig = [taxonA-gene1, taxonA-gene2, taxonB-gene1, taxonB-gene2]
    # imap = {taxonA: [taxonA-gene1, taxonA-gene2], taxonB: [...], ...}
    imap = parse_imap(seqs, delim, delim_idxs, delim_join)

    # -----------------------------------------------------------
    # parse sequence info (e.g., homology) if present, else use seq length
    # scores = {taxonA-gene1: 0.9, taxonA-gene2: 0.8, ...}
    scores = get_scores(seqs, info_file)

    # -----------------------------------------------------------
    # build dict mapping {name: regex} for every name in imap
    if regex or regex_file:
        groups = get_groups_from_regex(regex, regex_file, imap)
    else:
        groups = get_groups_from_split(split, split_file, imap)

    # select a single best in each isoform group
    collapsed = {}
    for group_key, members in groups.items():

        # get scores for group members
        sgroup = [scores[i] for i in members]

        # sort by these keys in this order
        best_iso = max(
            sgroup,
            key=lambda x: (
                x.effective_informative_homolog_bp,
                x.homology_internal,
                x.informative_bp,
                x.bp_kept,
                x.name,
            )
        )

        #
        collapsed[best_iso.name] = seqs[best_iso.name]

        # report if iso's were collapsed
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

    # filter the locus if too many were collapsed
    # filtered = len(collapsed) < min_count

    # write output
    # if not filtered:
    with open(outpath, 'w') as hout:
        for uname in sorted(collapsed):
            hout.write(f">{uname}\n{seqs[uname]}\n")


def run_isoform_prune(args):
    """..."""
    set_log_level(args.log_level)

    # check in files
    if not (args.input.exists() and args.input.is_file()):
        raise TwigError(f"{args.input} not found")

    # require only one
    group_args = [args.regex, args.regex_file, args.split, args.split_file]
    group_args = [i for i in group_args if i is not None]
    if not any(group_args):
        raise TwigError("must use one of: --regex, --regex-file, --split, --split-file")
    if len(group_args) > 1:
        raise TwigError("can select only one of: --regex, --regex-file, --split, or --split-file")

    # ensure outpath and pdir exists
    if args.outpath.is_dir():
        raise TwigError("outpath should be a file path not dir")
    args.outpath.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    if args.outpath.exists() and not args.force:
        logger.warning(f"[skipping] {args.outpath} already exists. Use --force to overwrite")
        return

    # filter by minimum length
    prune_sequences(
        args.input,
        args.outpath,
        #
        args.delim,
        args.delim_idxs,
        args.delim_join,
        #
        args.info_file,
        args.regex,
        args.regex_file,
        args.split,
        args.split_file,
    )


def main():
    from twig.cli.subcommands import get_parser_macse_prep
    parser = get_parser_macse_prep()
    args = parser.parse_args()
    run_isoform_prune(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
