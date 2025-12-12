#!/usr/bin/env python

"""Filter a set of gene trees based on various criteria.

This is useful to filter tree sets for species tree and network analysis,
among others.

>>> tree-filter -t TREES -o OUT -c 1
>>> tree-filter -t TREES -o OUT -I IMAP -M MINMAP
>>> tree-filter -t TREES -o OUT -m 20

TODO
----
Work over a generator of trees
"""

import sys
import textwrap
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
import numpy as np
import toytree
from twig.utils.logger_setup import set_log_level


KWARGS = dict(
    prog="tree-filter",
    usage="tree-filter [options]",
    help="filter and relabel a set of trees for downstream analyses",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | tree-filter: Filter/Relabel a set of trees for downstream analyses
        -------------------------------------------------------------------
        | The order of [optional] operations is (1) parse names from delim
        | split tip labels; (2) subsample and/or relabel names by imap;
        | (3) exclude outlier edges; (4) collapse outgroups; (5) exclude
        | trees w/ < min-tips; (6) exclude trees w/ > min-copies in any name.
        | Note that if names are parsed by delim, then the imap should contain
        | the parsed names. The labels on output trees will be the same as in
        | the inputs unless --relabel-delim or --relabel-imap is used. But
        | parsing shorter names can be useful to match names to imap even
        | if keeping full names. Filtering stats are written to stderr, which
        | can be redirected to a log file with `2> log.txt`.
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        # relabel by split-select-join on tip names
        $ twig tree-filter -i NWK -d '-' -di 0 2 -dj '-rd' > relabeled-trees.nwk

        # get only single-copy gene trees
        $ twig tree-filter -i NWK -c 1 > single-copy-trees.nwk

        # get only trees w/ >20 tips
        $ twig tree-filter -i NWK -m 20 > min20-trees.nwk

        # exclude outlier edges
        $ twig tree-filter -i NWK --exclude -ei 4 > cleaned-trees.nwk

        Example IMAP
        -------------
        sampleA     pop1
        sampleB     pop1
        sampleC     pop2
        sampleD     pop2
        sampleE     outgroup
        sampleF     outgroup

        Example MINMAP
        --------------
        pop1        1
        pop2        1
        outgroup    1

        Examples using IMAP/MINMAP
        --------------------------
        # subsample to names in imap
        $ twig tree-filter -i NWK -I IMAP > subsample-trees.nwk

        # subsample and relabel to pop names in imap
        $ twig tree-filter -i NWK -I IMAP -ri > relabeled-trees.nwk

        # subsample, relabel, and keep only one outgroup (best to root before)
        $ twig tree-filter -i NWK -I IMAP -ri --collapse > final-trees.nwk

        Example on a many trees in different paths
        -------------------------------------------
        $ parallel "twig tree-filter -i {} ... > {}.filtered" ::: TREES/*.nwk
    """)
)


def get_parser_tree_filter(parser: ArgumentParser | None = None) -> ArgumentParser:
    """Return a parser for relabel tool.
    """
    # create parser or connect as subparser to cli parser
    if parser:
        KWARGS['name'] = KWARGS.pop("prog")
        parser = parser.add_parser(**KWARGS)
    else:
        KWARGS.pop("help")
        parser = ArgumentParser(**KWARGS)

    # path i/o args
    parser.add_argument("-i", "--input", type=Path, metavar="path", required=True, help="newick or multi-newick trees file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # filter on names
    # parser.add_argument("-i", "--include", type=str, metavar="str", nargs="+", help="subset of names to keep")
    parser.add_argument("-I", "--imap", type=Path, metavar="path", help=r"filepath listing names to keep, or assigning (name population)")
    parser.add_argument("-M", "--minmap", type=Path, metavar="path", help=r"filepath listing (population mincov) for filters")

    # filter on tree data
    parser.add_argument("-m", "--min-tips", type=int, metavar="int", default=4, help="min tips after pruning [4]")
    parser.add_argument("-c", "--max-copies", type=int, metavar="int", help="filter trees with >c gene copies from a taxon [None]")
    parser.add_argument("-eo", "--edge-outlier-outgroup", type=float, metavar="float", default=10, help="exclude 'outgroup' population edges if >eo stdev from mean [10]")
    parser.add_argument("-ei", "--edge-outlier-ingroup", type=float, metavar="float", default=5, help="exclude non 'outgroup' population edges if >ei stdev from mean [5]")

    # actions
    parser.add_argument("-rd", "--relabel-delim", action="store_true", help="relabel tips by their delim parsed names")
    parser.add_argument("-ri", "--relabel-imap", action="store_true", help="relabel tips to their imap mapped names")
    # parser.add_argument("-x", "--nexus", action="store_true", help="export in NEXUS format. Retains path/tree names")
    parser.add_argument("--subsample", action="store_true", help="subsample to include only tips in imap")
    parser.add_argument("--exclude-outliers", action="store_true", help="exclude tips with outlier edge lengths (>ei or >eo)")
    parser.add_argument("--require-outgroups", action="store_true", help="require at least one 'outgroup' sample")
    parser.add_argument("--collapse-outgroups", action="store_true", help="keep only the most distant 'outgroup' (assumes rooted trees)")

    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="INFO", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def set_delim_and_imap_labels(tree, delim, idxs, join, imap):
    """Strip names to keep only accession IDs"""
    for node in tree[:tree.ntips]:
        items = node.name.split(delim)
        label = join.join([items[i] for i in idxs])
        node.delim = label
        node.imap = imap.get(node.delim)
    return tree


def subsample_by_imap(tree, imap, minmap, subsample):
    # subsample to keep only tips in imap
    if subsample:
        keep = []
        for node in tree[:tree.ntips]:
            if node.imap:
                keep.append(node.name)
        tree = tree.mod.prune(*keep)

    # filter by mincov
    if minmap:
        counts = {}
        for node in tree[:tree.ntips]:
            pop = imap.get(node.delim)
            if pop in counts:
                counts[pop] += 1
            else:
                counts[pop] = 1
        diff = set(counts) - set(minmap)
        if diff:
            raise KeyError(f"{diff} in imap but not in minmap")
        for key in counts:
            if counts[key] < minmap[key]:
                return tree, True        
    return tree, False


def exclude_long_tips(tree, ingroup_z, outgroup_z):
    # get dict mapping tips to their list of tips to consider dropping
    tips = tree[:tree.ntips]

    # get tip len distribution
    ingroup_dists = [i.dist for i in tips if i.imap != "outgroup"]
    mean = np.mean(ingroup_dists)
    stdev = np.std(ingroup_dists)

    droplist = []
    for node in tips:
        zi = abs(node.dist - mean) / stdev
        if node.imap != "outgroup":
            if zi > ingroup_z:
                droplist.append(node.name)
        else:
            if zi > outgroup_z:
                droplist.append(node.name)

    # return tree with outlier tips dropped
    if droplist:
        tree.mod.drop_tips(*droplist, inplace=True)
        return exclude_long_tips(tree, ingroup_z, outgroup_z)
    return tree


def collapse_and_require_outgroups(tree, imap, require_outgroups, collapse_outgroups):
    # get the delim labels from imap that are in the tree
    onodes = [i for i in tree[:tree.ntips] if i.imap == "outgroup"]

    # return options for no outgroups present
    if not onodes:
        if require_outgroups:
            return tree, True
        else:
            return tree, False

    # return option for >=1 outgroup present
    if (len(onodes) == 1) or (not collapse_outgroups):
        return tree, False

    # outgroups present and in need of collapse
    tree.ladderize(inplace=True)
    onodes = sorted(onodes, key=lambda x: x.idx)
    tree.mod.drop_tips(*onodes[1:], inplace=True)
    return tree, False


def filter_by_max_copies(tree, max_copies):
    delims = [i.delim for i in tree[:tree.ntips]]
    if any(delims.count(i) > max_copies for i in set(delims)):
        return True
    return False


def run_tree_filter(args):
    set_log_level(args.log_level)

    # parse the imap
    imap = {}
    if args.imap:
        with args.imap.open() as hin:
            for line in hin.readlines():
                data = line.strip().split()
                if data:
                    if len(data) == 1:
                        label = pop = data[0]
                    else:
                        label = data[0]
                        pop = data[1]
                    imap[label] = pop

    # parse the minmap
    minmap = {}
    if args.minmap:
        with args.minmap.open() as hin:
            for line in hin.readlines():
                pop, mincov, *other = line.strip().split()
                minmap[pop] = mincov

    # check required args
    if args.collapse_outgroups or args.require_outgroups:
        if not imap:
            raise ValueError("must provide an --imap when using --collapse-outgroups or --require-outgroups")
        if "outgroup" not in imap.values():
            raise ValueError("imap must contain a population named 'outgroup' when using --collapse-outgroups or --require-outgroups")

    if not args.max_copies:
        args.max_copies = float('inf')

    # store filtering stats
    ntrees_start = 0
    ntrees_end = 0
    outlier_tips_removed = 0
    filters = {
        "minmap": 0,
        "require-outgroup": 0,
        "min-tips": 0,
        "max-copies": 0,
    }

    if args.out:
        out_handle = open(args.out, "w")

    # loop over trees
    assert args.input.is_file, f"input file '{args.input}' not found"
    for newick in args.input.open().readlines():
        if not newick:
            continue
        tree = toytree.tree(newick)
        ntrees_start += 1
        tree = set_delim_and_imap_labels(tree, args.delim, args.delim_idxs, args.delim_join, imap)

        # imap methods
        if imap:
            tree, filt = subsample_by_imap(tree, imap, minmap, args.subsample)
            filters['minmap'] += int(filt)
        if args.exclude_outliers:
            npre = tree.ntips
            tree = exclude_long_tips(tree, args.edge_outlier_ingroup, args.edge_outlier_outgroup)
            npost = tree.ntips
            outlier_tips_removed += npre - npost
        if args.collapse_outgroups or args.require_outgroups:
            tree, filt = collapse_and_require_outgroups(tree, imap, args.require_outgroups, args.collapse_outgroups)
            if filt:
                filters['require-outgroup'] += 1
                continue

        # filters
        if tree.ntips < args.min_tips:
            filters["min-tips"] += 1
            continue
        if filter_by_max_copies(tree, args.max_copies):
            filters["max-copies"] += 1
            continue

        # relabel
        if args.relabel_imap or args.relabel_delim:
            for node in tree[:tree.ntips]:
                node.name = node.delim
        # write
        ntrees_end += 1
        if args.out:
            out_handle.write(f"{tree.write()}\n")
        else:
            print(tree.write(), file=sys.stdout)

    # print stats to stderr
    print(f"CMD: {' '.join(sys.argv[:])}", file=sys.stderr)
    print(f"ntrees start = {ntrees_start}", file=sys.stderr)
    for key in filters:
        print(f"trees filtered by {key} = {filters[key]}", file=sys.stderr)
    print(f"ntrees end = {ntrees_end}", file=sys.stderr)
    print(f"outlier tips removed = {outlier_tips_removed}", file=sys.stderr)    


def main():
    parser = get_parser_tree_filter()
    args = parser.parse_args()
    run_tree_filter(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
        raise
