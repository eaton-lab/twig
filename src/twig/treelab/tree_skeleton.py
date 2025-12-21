#!/usr/bin/env python

"""Force a gene tree topology to match a species tree topology, with
gene duplicates ... grouped as sister.

"""


from toytree import tree


def run_tree_skeleton(args):
    """The names in sptree should match those parsed from the gtree
    using delim/idx/join.
    """
    gtree = tree(args.input)
    sptree = tree(args.sptree)

    # map {spp-label: [g-label, ...], ...}
    copies = {}
    for tip in gtree.get_tip_labels():
        parts = tip.split(args.delim)
        name = args.delim_join.join([parts[i] for i in args.delim_idxs])
        if name in copies:
            copies[name].append(tip)
        else:
            copies[name] = [tip]

    for sname in copies:
        for idx, gname in enumerate(copies[sname]):
            # optionally use species labels
            if args.relabel_delim:
                gname = sname
            # set labels to existing or new nodes depending on copy number
            if idx == 0:
                sptree.get_nodes(sname)[0].name = gname
            elif idx == 1:
                sptree.mod.add_internal_node_and_child(sname, name=gname, parent_name="", inplace=True)
            elif idx > 1:
                sptree.mod.add_sister_node(sname, name=gname, inplace=True)

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
