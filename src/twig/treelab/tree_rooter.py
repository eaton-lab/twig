#!/usr/bin/env python

"""Root gene tree given a rooted species tree. 

...
"""

from typing import List
import sys
import textwrap
from pathlib import Path
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from loguru import logger
import toytree


KWARGS = dict(
    prog="tree-rooter",
    usage="tree-rooter [options]",
    help="Re-root gene trees by outgroup, species tree, and/or MAD",
    formatter_class=lambda prog: RawDescriptionHelpFormatter(prog, width=140, max_help_position=140),
    description=textwrap.dedent("""
        -------------------------------------------------------------------
        | tree-rooter: ...
        -------------------------------------------------------------------
        | Root trees by designating an outgroup, or a set of potential
        | outgroups in the case that a set of trees may or may not include
        | a single outgroup. The ...
        -------------------------------------------------------------------
    """),
    epilog=textwrap.dedent("""
        Examples
        --------
        # root trees given ordered options of outgroup clades
        $ twig tree-rooter -t NWK -r A > rooted-trees.nwk
        $ twig tree-rooter -t NWK -r A B A,B > rooted-trees.nwk

        # root trees using minimal ancestor deviation (MAD)
        $ twig tree-rooter -t NWK --mad > rooted-trees.nwk

        # root trees on outgroups if present and mad root the rest
        $ twig tree-rooter -t NWK -r A B --mad > rooted-trees.nwk

        # root trees on outgroups in order of a rooted species tree
        $ twig tree-rooter -t NWK -s SPTREE > rooted-trees.nwk

        # root trees on a restricted set of outgroups from a rooted species tree
        $ twig tree-rooter -t NWK -s SPTREE -r Z Y > rooted-trees.nwk

        # return trees that could not be rooted given the options
        $ twig tree-rooter -t NWK -s SPTREE -x > unrooted-trees.nwk
    """)
)


def get_parser_tree_rooter(parser: ArgumentParser | None = None) -> ArgumentParser:
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

    # rooting options
    parser.add_argument("-s", "--sptree", type=Path, metavar="path", help="rooted species tree")
    parser.add_argument("-r", "--outgroups", type=str, metavar="str", nargs="+", help="list outgroup tip labels")
    parser.add_argument("-R", "--outgroups-file", type=Path, metavar="path", help="file listing outgroup tip labels")

    # return option
    parser.add_argument("-m", "--mad", action="store_true", help="use minimal ancestor deviation to estimate root of remaining unrooted trees")
    parser.add_argument("-x", "--not-rooted", action="store_true", help="return the trees that could not be rooted given the options")

    # logging
    parser.add_argument("-l", "--log-level", type=str, metavar="level", default="CRITICAL", help="stderr logging level (DEBUG, [INFO], WARNING, ERROR)")
    # parser.add_argument("-L", "--log-file", type=Path, metavar="path", help="append stderr log to a file")    
    return parser


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


def run_tree_rooter(args):
    # require -s -r or -R
    toytree.set_log_level(args.log_level)#, args.log_file)

    # parse outfile args to list
    if args.outgroups_file:
        outgroups = [i.strip() for i in args.outgroups.open().readlines()]
    else:
        outgroups = args.outgroups

    if args.sptree:
        sptree = toytree.tree(args.sptree)
        root_clades = get_rooting_clades(sptree, args.outgroups)
    else:
        root_clades = [i.split(",") for i in outgroups] if outgroups else []
    logger.debug(f"ordered rooting clades: {root_clades}")

    # track success
    count = {"rerooted": 0, "not-rerooted": 0}
    rtrees = []
    utrees = []

    # iterate over newicks in treefile
    with args.trees.open() as datain:
        for tidx, nwk in enumerate(datain.readlines()):
            logger.debug(f"{tidx}...")            
            tree = toytree.tree(nwk)
            tips = tree.get_tip_labels()
            rooted = False
            for rc in root_clades:
                try:
                    rc = set(i for i in rc if i in tips)
                    if not rc:
                        continue
                    tree = tree.root(*rc)
                    rooted = True
                    break
                except (AttributeError, ValueError):
                    pass
                except toytree.utils.ToytreeError:
                    pass
            if rooted:
                count['rerooted'] += 1
                rtrees.append(tree)
            else:
                count['not-rerooted'] += 1
                utrees.append(tree)

    # use MAD to root not-rerooted trees
    if args.mad:
        for t in utrees:
            try:
                t.mod.root_on_minimal_ancestor_deviation()
            except IndexError:
                print(t.write())
                sys.exit(1)
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

