#!/usr/bin/env python

"""Force a gene tree topology to match a species tree topology, with
gene duplicates ... grouped as sister.

"""

import sys
from loguru import logger
from toytree import tree
from twig.utils.logger_setup import set_log_level


def run_tree_skeleton(args):
    """The names in sptree should match those parsed from the gtree
    using delim/idx/join.
    """
    set_log_level(args.log_level)

    gtree = tree(args.input, internal_labels="name")
    sptree = tree(args.sptree, internal_labels="name")

    # map {spp-label: [g-label, ...], ...}
    copies = {}
    for tip in gtree.get_tip_labels():
        parts = tip.split(args.delim)
        name = args.delim_join.join([parts[i] for i in args.delim_idxs])
        if name in copies:
            copies[name].append(tip)
        else:
            copies[name] = [tip]

    # prune to only keep edges connecting the tips in 'copies'
    sptree = sptree.mod.prune(*list(copies))

    # iterate over tips in species tree adding extra copies to tips as needed
    for sname in copies:

        node = sptree.get_nodes(sname)[0]
        gnames = copies[sname]
        lidx = len(gnames)

        # set labels to existing or new nodes depending on copy number
        if lidx == 1:
            node.name = sname if args.relabel_delim else gnames[0]
        elif lidx == 2:
            logger.debug(f"two copies of {sname}")
            node.name = sname if args.relabel_delim else gnames[0]
            sptree.mod.add_internal_node_and_child(node, name=sname if args.relabel_delim else gnames[1], parent_name="", inplace=True)
        else:
            logger.debug(f"{len(gnames)} copies of {sname}")
            node.name = sname if args.relabel_delim else gnames[0]
            sptree.mod.add_internal_node_and_child(node, name=sname if args.relabel_delim else gnames[1], parent_name="", inplace=True)
            for sidx in range(2, len(gnames)):
                sptree.mod.add_sister_node(node, name=gnames[sidx], inplace=True)

    # write to file or stdout
    if args.out:
        sptree.write(args.out)
    else:
        print(sptree.write(None), file=sys.stdout)


def main():
    from twig.cli.subcommands import get_parser_tree_skeleton
    parser = get_parser_tree_skeleton()
    args = parser.parse_args()
    run_tree_skeleton(args)



if __name__ == "__main__":
    pass
