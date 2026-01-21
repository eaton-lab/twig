#!/usr/bin/env python

"""Filter a set of gene trees based on various criteria.

This is useful to filter tree sets for species tree and network analysis,
among others.

>>> tree-filter -t TREES -o OUT -c 1
>>> tree-filter -t TREES -o OUT -I IMAP -M MINMAP
>>> tree-filter -t TREES -o OUT -m 20

"""

import sys
# from pathlib import Path
from loguru import logger
import numpy as np
import toytree
from twig.utils.logger_setup import set_log_level


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


def collapse_and_require_outgroups(tree, outgroups, require_outgroups, collapse_outgroups):
    """Return Tree and boolean where True means it failed to pass"""
    # get the delim labels from imap that are in the tree
    onodes = [i for i in tree[:tree.ntips] if i.delim in outgroups]

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


def collapse_by_min_support_and_require_min_splits(tree, min_support: int, min_splits: int):
    """Return tree with low support nodes collapsed and boolean of whether min_splits are still present"""
    tree = tree.mod.collapse_nodes(min_support=min_support)
    n_internal_edges = tree.nedges - tree.ntips
    if n_internal_edges > min_splits:
        return tree, False
    return tree, True


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

    # get list of outgroups from 'outgroups' or imap
    if args.outgroups:
        with args.outgroups.open() as hin:
            outgroups = [i.strip().split()[0] for i in hin]
    else:
        outgroups = [i for (i, j) in imap.items() if j == "outgroup"]
    logger.debug(f"outgroups = {outgroups}")

    # check required args
    if args.collapse_outgroups or args.require_outgroups:
        if not outgroups:
            raise ValueError("must designate outgroups in --imap or --outgroups to use --collapse-outgroups or --require-outgroups")
        # if "outgroup" not in imap.values():
        #     raise ValueError("imap must contain a population named 'outgroup' when using --collapse-outgroups or --require-outgroups")

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
            tree, filt = collapse_and_require_outgroups(tree, outgroups, args.require_outgroups, args.collapse_outgroups)
            if filt:
                filters['require-outgroup'] += 1
                continue

        # filters
        if tree.ntips < args.min_tips:
            # logger.debug(f"tree filtered ntips={tree.ntips}; nwk={tree.write()}")
            filters["min-tips"] += 1
            continue
        if filter_by_max_copies(tree, args.max_copies):
            filters["max-copies"] += 1
            continue

        # relabel
        if args.relabel_delim:
            for node in tree[:tree.ntips]:
                node.name = node.delim
        if args.relabel_imap:
            for node in tree[:tree.ntips]:
                node.name = node.imap
                if node.imap is None:
                    logger.warning(f"node {node.name} ({node.delim}) not in imap")

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
    from ..cli.subcommands import get_parser_tree_filter
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
