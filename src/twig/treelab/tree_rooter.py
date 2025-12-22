#!/usr/bin/env python

"""Root gene tree given a rooted species tree. 

...
"""

# from typing import List
import sys
# from pathlib import Path
from loguru import logger
import toytree
from twig.utils.logger_setup import set_log_level


def get_rooting_clades(sptree, outgroups):
    """Return list of rooting clades traversing down sptree"""
    tips = sptree.get_tip_labels()
    ingroup = set(tips) - set(outgroups)
    
    # traverse tree from root to tips selecting smaller child clade
    root_clades = []
    for node in sptree:
        if node.children:
            children = sorted(node.children, key=lambda x: len(x.get_leaves()))
            root_clade = children[0].get_leaf_names()
            if not any(i in root_clade for i in ingroup):
                root_clades.append(root_clade)
    return root_clades[::-1]


def set_delim_labels(tree, delim, idxs, join):
    """Strip names to keep only accession IDs"""
    for node in tree[:tree.ntips]:
        items = node.name.split(delim)
        label = join.join([items[i] for i in idxs])
        node.delim = label
    return tree


def run_tree_rooter(args):
    # require -s -r or -R
    set_log_level(args.log_level)
    assert args.input.exists(), f"trees file {args.input} not found."

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
        outgroups = [i for (i, j) in imap.items() if j == "outgroup"]
    else:
        outgroups = args.outgroups
    assert outgroups, "no outgroups specified with -r or -I"

    # get root_clades from sptree or outgroups
    if args.sptree:
        sptree = toytree.tree(args.sptree, internal_labels="name")
        root_clades = get_rooting_clades(sptree, outgroups)
    else:
        root_clades = [i.split(",") for i in outgroups] if outgroups else []
    logger.debug(f"ordered rooting clades: {root_clades}")
    assert root_clades, "no clades in sptree match to specified outgroups. Check delim args?"

    # track success
    count = {"rerooted": 0, "not-rerooted": 0}
    rtrees = []
    utrees = []

    # iterate over newicks in treefile
    with args.input.open() as datain:
        for tidx, nwk in enumerate(datain.readlines()):
            if not nwk:
                continue
            tree = toytree.tree(nwk, internal_labels="name")
            tree = set_delim_labels(tree, args.delim, args.delim_idxs, args.delim_join)
            rooted = False
            ntips = tree.ntips
            for rc in root_clades:
                try:
                    # get tipnodes in gtree that are in rootclade
                    onodes = [i for i in tree[:tree.ntips] if i.delim in rc]
                    if not onodes:
                        logger.debug(f"tree cannot be rooted on root clade {rc}")
                        continue
                    # try rooting the tree on a set of rootclade tips
                    _tree = tree.root(*onodes)

                    # validate that ingroup is monophyletic (i.e., there are not
                    # other outgroup samples nested in it).
                    inodes = [i for i in _tree[:_tree.ntips] if i.delim not in outgroups]
                    if not _tree.is_monophyletic(*inodes):
                        logger.debug("rooted tree does not contain monophlyetic ingroup")
                        continue

                    # success
                    rooted = True
                    break

                # errors that arise on bad .root() options
                except (AttributeError, ValueError):
                    pass
                except toytree.utils.ToytreeError:
                    pass

            # accepted rooted tree
            tree = _tree
            if tree.ntips != ntips:
                logger.warning(f"ntips changed, {ntips}->{tree.ntips}, {nwk}")

            # set labels on tree
            if args.relabel_delim:
                for node in tree[:tree.ntips]:
                    node.name = node.delim

            if rooted:
                count['rerooted'] += 1
                rtrees.append(tree)
            else:
                count['not-rerooted'] += 1
                utrees.append(tree)

    # use MAD to root not-rerooted trees
    if args.mad:
        # for t in utrees:
        #     try:
        #         t.mod.root_on_minimal_ancestor_deviation()
        #     except IndexError:
        #         print(t.write())
        #         raise
        rtrees += [i.mod.root_on_minimal_ancestor_deviation() for i in utrees]
        count['mad-rooted'] = count['not-rerooted']
        count['not-rerooted'] = 0

    # print stats
    for key in count:
        print(f"'{key}': {count[key]}", file=sys.stderr)

    # select return set
    if args.not_rooted:
        newicks = "\n".join(i.write() for i in utrees)
    else:
        newicks = "\n".join(i.write() for i in rtrees)        

    # print results
    if args.out:
        with args.out.open("w") as hout:
            hout.write(newicks)
    else:
        print(newicks, file=sys.stdout)



def main():
    from twig.cli.subcommands import get_parser_tree_rooter
    parser = get_parser_tree_rooter()
    args = parser.parse_args()
    run_tree_rooter(args)


if __name__ == "__main__":
    try:
        main()        
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
        raise

