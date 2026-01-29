#!/usr/bin/env python

"""Trim and filter sequences for macse alignment.

Command
-------
$ twig macse-prep -i CDS -o OUT/ID.nt.fa

Output
------
- OUT/ID.trim.nt.fa
- OUT/ID.trim.info.tsv
"""

from typing import List
import re
from pathlib import Path
from loguru import logger
# import pandas as pd
from twig.utils.logger_setup import set_log_level


def prune_sequences(
    inpath: Path,
    outpath: Path,
    regex: re.compile,
    regex_file: Path,
    delim: str,
    delim_idxs: List[int],
    delim_join: str,
    scores_file: Path,
    # exclude: List[str],
    # subsample: List[str],
) -> None:
    """Write fasta with only one isoform per gene. When multple are present
    the one with greatest homology to other sequences is retained, with
    ties broken by length, and then order.
    """
    # parse input sequence fasta to dict
    # {taxonA-gene1: AAA..., taxonA-gene2: AAA..., ...}
    seqs = {}
    with open(inpath, 'r') as datain:
        for line in datain.readlines():
            if line:
                if line.startswith(">"):
                    uname = line.strip()[1:]
                    seqs[uname] = ""
                else:
                    seqs[uname] += line.strip()
    # logger.debug(f"seqs keys = {list(seqs)}")

    # assign sequences to subgroups to each use their own regex pattern,
    # or apply one regex pattern to everyone.
    # parse taxon-gene names to dict {}
    # orig = [taxonA-gene1, taxonA-gene2, taxonB-gene1, taxonB-gene2]
    # imap = {taxonA: [taxonA-gene1, taxonA-gene2], taxonB: [...], ...}
    if delim is None:
        imap = {"SHARED-REGEX": sorted(seqs)}
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

    # parse sequence scores (e.g., homology) if present, else use seq length
    # scores = {taxonA-gene1: 0.9, taxonA-gene2: 0.8, ...}
    scores = {}
    if scores_file:
        if scores_file.is_file() and scores_file.exists():
            with open(scores_file, 'r') as datain:
                for line in datain.readlines():
                    if line:
                        # use try/except to skip headers if present. Membership is checked below.
                        try:
                            name, score, *_ = line.split("\t")
                            length = len(seqs[name].replace("-", ""))
                            scores[str(name)] = (float(score), int(length))
                        except KeyError:
                            pass

            # if using scores then a score must be present for every sequence
            diff = set(seqs) - set(scores)
            if any(diff):
                raise ValueError(f"scores tsv does not contain the following sequence names from the input: [{diff}]")
        else:
            raise ValueError(f"scores file is malformed: {scores_file}")
    else:
        for name, seq in seqs.items():
            length = len(seq.replace("-", ""))
            scores[name] = (length, length)
    # logger.debug(f"scores={scores}")

    # build dict mapping {name: regex} for every name in imap
    patterns = {}
    if regex_file:
        if regex_file.is_file() and regex_file.exists():
            with open(regex_file, 'r') as datain:
                for line in datain.readlines():
                    if line:

                        # [TODO]
                        # parse regex file line
                        name, pattern, *x = line.split("\t")
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
            if delim:
                diff = set(imap) - set(patterns)
                if any(diff):
                    raise ValueError(f"regex-file tsv does not contain the following (-d) delimited names from the input: [{diff}]")
            else:
                diff = set(seqs) - set(patterns)
                if any(diff):
                    raise ValueError(f"regex-file tsv does not contain the following sequence names from the input: [{diff}]")
        else:
            raise ValueError(f"regex file is malformed: {regex_file}")
    else:
        patterns = {'SHARED-REGEX': (regex, 1)}

    # report patterns to logger
    for key in patterns:
        pat, grp = patterns[key]
        logger.debug(f"group={key}, regex={pat}, select={grp}")

    # keep a global set of observed sequences and raise an exception
    # if a sequence lands in more than one regex group
    observed = set()

    # iterate over (taxon-group, regular expression) tuples
    groups = {}
    for taxon_group, (pattern, gkey) in patterns.items():

        # select only the sequences in this taxon_group
        subgroup = imap.get(taxon_group)

        # skip if this locus does not contain this subgroup
        if not subgroup:
            continue

        # apply this groups specific regular expression to find isoforms
        for sname in subgroup:

            # get score for this sequence
            score = scores[sname]

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
            logger.debug(f"{sname} -- {subgroup} -- {group_key}")
            if group_key in groups:
                groups[group_key].append((sname, score))
            else:
                groups[group_key] = [(sname, score)]

    # select a single best in each isoform group
    collapsed = {}
    for group_key, members in groups.items():
        # members: list of (name, idict)
        # keep only one per group by [highest homology_total, longest_length, or first]
        # logger.info((group_key, members))
        best_name, _ = max(
            members,
            key=lambda x: (x[1][0], x[1][1])  # (score, length)
        )
        collapsed[best_name] = seqs[best_name]

        # report if iso's were collapsed
        kept = best_name[len(group_key):]
        others = [i[0] for i in members if i[0] != best_name]

        # get the leftover parts of names
        rests = [s[len(group_key):] for s in others]

        # shorten rests if very long
        if len(rests) > 5:
            rests = rests[:5] + ["..."]

        # report
        logger.debug(f"pruned {len(members) - 1} isoforms: group='{group_key}', kept='{kept}', pruned={rests}")

    # filter the locus if too many were collapsed
    # filtered = len(collapsed) < min_count

    # write output
    # if not filtered:
    with open(outpath, 'w') as hout:
        for uname in sorted(collapsed):
            hout.write(f">{uname}\n{seqs[uname]}\n")

    # report filtering stats
    # keys = list(collapsed)
    # mean_length = data.loc[keys, "bp_kept"].mean()
    # mean_trimmed = data.loc[keys, "bp_trim"].mean()
    # mean_homology = data.loc[keys, "homology_total"].mean()
    logger.info(f"isoform pruning reduced {len(seqs)} seqs -> {len(collapsed)} seqs -> [{outpath}]")
    # if filtered:
    #     logger.info(f"[{outpath.name}] locus did not pass 'min-count' filter ({len(collapsed)} < {min_count})")


def run_isoform_prune(args):
    """..."""
    set_log_level(args.log_level)

    # check in files
    if not (args.input.exists() and args.input.is_file()):
        raise IOError(f"{args.input} not found")

    # only one or the other allowed
    # if args.exclude and args.subsample:
    #     raise ValueError("choose one of --exclude or --subample, but not both")

    # require one or the other
    if not (args.regex or args.regex_file):
        raise ValueError("must enter either --regex or --regex-file")

    # ensure outpath and pdir exists
    if args.outpath.is_dir():
        raise IOError("outpath should be a file path not dir")
    args.outpath.parent.mkdir(exist_ok=True, parents=True)

    # bail out if final file exists
    if args.outpath.exists() and not args.force:
        logger.warning(f"[skipping] {args.outpath} already exists. Use --force to overwrite")
        return

    # if skip isoform collapse then set isoform grouper to arbitrary str
    # if args.skip_isoform_collapse:
    #     args.isoform_regex = re.compile("@@@@@")

    # filter by minimum length
    prune_sequences(
        args.input,
        args.outpath,
        args.regex,
        args.regex_file,
        args.delim,
        args.delim_idxs,
        args.delim_join,
        args.scores,
        # args.exclude,
        # args.subsample,
        # args.min_length,
        # args.min_count,
        # args.force,
    )

    # clean up tmp files
    # suffices = [
        # ".trim_info",
        # ".tmp.aa.fa",
        # ".tmp.trim.mask",
    # ]
    # if not args.keep:
    #     for suffix in suffices:
    #         path = args.outpath.parent / (args.outpath.name + f"{suffix}")
    #         if path.exists():
    #             path.unlink()
    # logger.info(f"isoform pruned sequences written to {args.outpath}")


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
