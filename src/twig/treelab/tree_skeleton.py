#!/usr/bin/env python

"""Force a gene tree topology to match a species tree topology, with
gene duplicates ... grouped as sister.

"""


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
        logger.debug(parts)
        name = args.delim_join.join([parts[i] for i in args.delim_idxs])
        logger.debug(name)
        if name in copies:
            copies[name].append(name)
        else:
            copies[name] = [name]
    logger.debug(copies)

    # iterate over species in the species tree
    for sname in copies:

        node = sptree.get_nodes(sname)[0]
        gnames = copies[sname]
        lidx = len(gnames)

        # set labels to existing or new nodes depending on copy number
        if lidx == 1:
            node.name = sname if args.relabel_delim else gnames[0]
        elif lidx == 2:
            node.name = sname if args.relabel_delim else gnames[0]
            sptree.mod.add_internal_node_and_child(node, name=sname if args.relabel_delim else gnames[1], parent_name="", inplace=True)
        else:
            node.name = sname if args.relabel_delim else gnames[0]
            sptree.mod.add_internal_node_and_child(node, name=sname if args.relabel_delim else gnames[1], parent_name="", inplace=True)
            sptree.mod.add_sister_node(node, name=gnames[2], inplace=True)
            # sptree.mod.add_internal_node_and_child(node, name=sname if args.relabel_delim else gnames[2], parent_name="", inplace=True)

    if args.out:
        sptree.write(args.out)
    else:
        sptree.write()


def main():
    from twig.cli.subcommands import get_parser_tree_skeleton
    parser = get_parser_tree_skeleton()
    args = parser.parse_args()
    run_tree_skeleton(args)



if __name__ == "__main__":
    pass
