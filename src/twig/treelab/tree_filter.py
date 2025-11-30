#!/usr/bin/env python

"""Filter a set of gene trees based on various criteria.

This is useful to filter tree sets for species tree and network analysis,
among others.

>>> tree-filter -t TREES -o OUT -c 1
>>> tree-filter -t TREES -o OUT -I IMAP -M MINMAP
>>> tree-filter -t TREES -o OUT -m 20
"""

from typing import List
import sys
import textwrap
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
import numpy as np
import toytree


KWARGS = dict(
    prog="tree-filter",
    usage="tree-filter [options]",
    help="filter and relabel a set of trees for downstream analyses",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | tree-filter: Filter a set of trees for downstream analyses
        -------------------------------------------------------------------
        | ...
        | The order of [optional] operations is (1) relabel by delim; (2)
        | subsample and/or relabel by imap; (3) exclude outlier; (4) collapse
        | outgroups; (5) exclude trees w/ < min-tips; (6) exclude trees w/
        | > min-copies in any name. Note that if names are relabed by delim,
        | then the imap should contain the relabeled names to subsample or
        | relabel further. Filtering stats are written to stderr, which can
        | be piped to a log file with `2> log.txt`.
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        # relabel by split-select-join on tip names
        $ twig tree-filter -t NWK -d '-' -di 0 2 -dj '-' > relabeled-trees.nwk

        # get only single-copy gene trees
        $ twig tree-filter -t NWK -c 1 > single-copy-trees.nwk

        # get only trees w/ >20 tips
        $ twig tree-filter -t NWK -m 20 > min20-trees.nwk

        # exclude outlier edges
        $ twig tree-filter -t NWK --exclude -ei 4 > cleaned-trees.nwk

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
        $ twig tree-filter -t NWK -i IMAP > subsample-trees.nwk

        # subsample and relabel to pop names in imap
        $ twig tree-filter -t NWK -i IMAP --relabel > relabeled-trees.nwk

        # subsample, relabel, and keep only one outgroup (best to root before)
        $ twig tree-filter -t NWK -i IMAP --relabel --collapse > final-trees.nwk
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
    parser.add_argument("-t", "--trees", type=Path, metavar="path", required=True, help="newick or multi-newick trees file")
    parser.add_argument("-o", "--out", type=Path, metavar="path", help="outfile name else printed to stdout")

    # parsing names
    parser.add_argument("-d", "--delim", type=str, metavar="str", help="delimiter to split tip labels")
    parser.add_argument("-di", "--delim-idxs", type=int, metavar="int", nargs="+", default=[0], help="index of delimited name items to keep")
    parser.add_argument("-dj", "--delim-join", type=str, metavar="str", default="-", help="join character on delimited name items")

    # filter on names
    # parser.add_argument("-i", "--include", type=str, metavar="str", nargs="+", help="subset of names to keep")
    parser.add_argument("-i", "--imap", type=Path, metavar="path", help=r"filepath listing names to keep, or assigning (name population)")
    parser.add_argument("-m", "--minmap", type=Path, metavar="path", help=r"filepath listing (population mincov) for filters")

    # filter on tree data
    parser.add_argument("-s", "--min-tips", type=int, metavar="int", default=4, help="min tips after pruning [4]")
    parser.add_argument("-c", "--max-copies", type=int, metavar="int", help="filter trees with >c gene copies from a taxon [None]")
    parser.add_argument("-eo", "--edge-outlier-outgroup", type=float, metavar="float", default=10, help="exclude 'outgroup' population edges if >eo stdev from mean [10]")
    parser.add_argument("-ei", "--edge-outlier-ingroup", type=float, metavar="float", default=5, help="exclude non 'outgroup' population edges if >ei stdev from mean [5]")

    # actions
    parser.add_argument("--subsample", action="store_true", help="subsample to include only tips in imap")
    parser.add_argument("--relabel", action="store_true", help="relabel tips to their imap population names")
    parser.add_argument("--exclude-outliers", action="store_true", help="exclude tips with outlier edge lengths (>ei or >eo)")
    parser.add_argument("--require-outgroups", action="store_true", help="require at least one 'outgroup' sample")
    parser.add_argument("--collapse-outgroups", action="store_true", help="keep only the most distant 'outgroup' (assumes rooted trees)")

    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="CRITICAL", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")
    return parser


def relabel_tips_by_delim(tree, delim, idxs, join):
    """Strip names to keep only accession IDs"""

    for node in tree[:tree.ntips]:
        items = node.name.split(delim)
        label = join.join([items[i] for i in idxs])
        node.name = label
    return tree


def relabel_and_subsample_by_imap(tree, imap, minmap, subsample, relabel):

    # subsample to keep only tips in imap
    if subsample:
        keep = []
        for node in tree[:tree.ntips]:
            if imap.get(node.name):
                keep.append(node.name)
        tree = tree.mod.prune(*keep)

    # subsample to keep only tips in imap
    if relabel:
        for node in tree[:tree.ntips]:
            node.name = imap[node.name]

    # filter by mincov
    if minmap:
        counts = {}
        for node in tree[:tree.ntips]:
            pop = imap[node.name]
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
    tip_dists = {i.name: i.dist for i in tree[:tree.ntips]}

    # get tip len distribution
    ingroup_dists = [j for (i, j) in tip_dists.items() if i != "outgroup"]
    mean = np.mean(ingroup_dists)
    stdev = np.std(ingroup_dists)

    droplist = []
    for name, dist in tip_dists.items():
        zi = abs(dist - mean) / stdev
        if name != "outgroup":
            if zi > ingroup_z:
                droplist.append(name)
        else:
            if zi > outgroup_z:
                droplist.append(name)

    # return tree with outlier tips dropped
    if droplist:
        tree.mod.drop_tips(*droplist, inplace=True)
        return exclude_long_tips(tree, ingroup_z, outgroup_z)
    return tree


def collapse_and_require_outgroups(tree, imap, relabel, require_outgroups, collapse_outgroups):
    if relabel:
        outgroups = [i for (i, j) in imap.items() if j == "outgroup"]
    else:
        outgroups = ["outgroup"]

    try:
        onodes = tree.get_nodes(*outgroups)
    except ValueError:
        onodes = []

    if not onodes:
        if require_outgroups:
            return tree, True
        else:
            return tree, False

    if (len(onodes) == 1) or (not collapse_outgroups):
        return tree, False

    # outgroups present and in need of collapse
    tree.ladderize(inplace=True)
    onodes = sorted(onodes, key=lambda x: x.idx)
    tree.mod.drop_tips(*onodes[1:], inplace=True)
    return tree, False


def filter_by_max_copies(tree, max_copies):
    tips = tree.get_tip_labels()
    if any(tips.count(i) > max_copies for i in set(tips)):
        return tree, True
    return tree, False


def run_tree_filter(args):
    toytree.set_log_level(args.log_level)

    # parse the imap
    imap = {}
    if args.imap:
        with args.imap.open() as hin:
            for line in hin.readlines():
                label, pop, *other = line.strip().split()
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
        assert imap, "must provide an --imap when using --collapse-outgroups or --require-outgroups"
        assert "outgroup" in imap.values(), "imap must contain a population named 'outgroup' when using --collapse-outgroups or --require-outgroups"

    # load the trees [todo: support, etc. here?]
    trees = toytree.mtree(args.trees).treelist

    # store filtering stats
    ntrees = len(trees)
    filters = {
        "minmap": 0,
        "require-outgroup": 0,
        "min-tips": 0,
        "max-copies": 0,
    }

    # [1] relabel by split-select-rejoin
    if args.delim:
        trees = [relabel_tips_by_delim(i, args.delim, args.delim_idxs, args.delim_join) for i in trees]

    # [2] subsample and/or relabel by imap
    if imap:
        results = [relabel_and_subsample_by_imap(i, imap, minmap, args.subsample, args.relabel) for i in trees]
        filters["minmap"] = sum(j for (i, j) in results)
        trees = [i for (i, j) in results if not j]

    # [3] exclude outlier edges
    if args.exclude_outliers:
        trees = [exclude_long_tips(i, args.edge_outlier_ingroup, args.edge_outlier_outgroup) for i in trees]

    # [4] collapse outgroups (warn if no samples named outgroup?)
    if args.collapse_outgroups or args.require_outgroups:
        results = [collapse_and_require_outgroups(i, imap, args.relabel, args.require_outgroups, args.collapse_outgroups) for i in trees]
        filters["require-outgroup"] = sum(j for (i, j) in results)
        trees = [i for (i, j) in results if not j]

    # [5] filter by min-tips
    if args.min_tips:
        results = [(i, i.ntips < args.min_tips) for i in trees]
        filters["min-tips"] = sum(j for (i, j) in results)
        trees = [i for (i, j) in results if not j]

    # [6] filter by max-copies
    if args.max_copies:
        results = [filter_by_max_copies(i, args.max_copies) for i in trees]
        filters["max-copies"] = sum(j for (i, j) in results)
        trees = [i for (i, j) in results if not j]

    # print stats to stderr
    print(f"ntrees start = {ntrees}", file=sys.stderr)
    for key in filters:
        print(f"trees filtered by {key} = {filters[key]}", file=sys.stderr)
    print(f"ntrees end = {len(trees)}", file=sys.stderr)

    # print to stdout or write to file
    if not trees:
        logger.warning("no trees passed filtering")
    newicks = "\n".join(i.write() for i in trees)
    if args.out:
        print(newicks, file=args.out)
    else:
        print(newicks, file=sys.stdout)


def main():
    parser = get_parser_tree_filter()
    args = parser.parse_args()
    run_tree_filter(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    # except TwigError as exc:
    #     logger.error(exc)
    except Exception as exc:
        logger.error(exc)
        raise
