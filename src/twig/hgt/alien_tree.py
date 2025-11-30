#!/usr/bin/env python

"""Tree-based statistics to detect HGT.


"""

from toytree.core import ToyTree, Node
from toytree.infer.src.mul_trees import get_duplication_loss_coalescence



def ntaxa_in_hgt_containing_clade(species_tree: ToyTree, gene_tree: ToyTree, node: Node) -> int:
    """Return num taxa in subtree containing HGT clade.
    
    This is done by (1) pruning the HGT node out of the gene tree; 
    (2) identifying the subtree it is part of by assiging gene tree
    nodes as duplicated or coalesced in the absence of the HGT node,
    and then (3) putting the HGT back in, and tracing back to the last
    duplication event, or root species tree node. Finally, the number
    of descending species from that node is counted.

    Parameters
    ----------
    species_tree: ToyTree
        A rooted species tree...
    gene_tree: ToyTree
        A rooted gene tree...
    node: Node
        A Node of the gene tree putatively introduced by HGT.

    Examples
    --------
    ...
    """
    # label the sister of the input gtree
    tree = gene_tree.copy()
    tree[node.idx].get_sisters()[0].xxx = True

    # drop the HGT clade
    clade, clean_tree = tree.mod.bisect(node.idx)

    # label gene tree nodes with boolean feature 'duplication'
    clean_dup_tree = get_duplication_loss_coalescence(clean_tree, species_tree)['tree']

    # get sister node on the clean tree and add data back in
    sister = [i for i in clean_dup_tree if hasattr(i, 'xxx')][0]
    
    # add the HGT clade back to the tree
    ftree = clean_dup_tree.mod.add_internal_node_and_subtree(sister, subtree=clade)

    # trace back to last duplication or species tree node.
    xnode = [i for i in ftree if hasattr(i, 'xxx')][0].up
    for anode in xnode.iter_ancestors():
        if anode.duplication:
            break
        if hasattr(node, 'ns'):
            if anode.ns == species_tree.treenode.idx:
                break
        top_node = anode

    # count ntaxa
    a = set(top_node.get_leaf_names())
    b = set(node.get_leaf_names())
    c = set(species_tree.get_tip_labels())
    clade_size = len(c.intersection(a - b))
    return clade_size


def alien_tree_metric(species_tree: ToyTree, gene_tree: ToyTree, topology_only: bool) -> float:
    """

    """
    for node in gene_tree[:-1]:

        # get the sptree node that is mrca to this taxon set
        set_n = set(node.get_leaf_names())
        mrca_n = species_tree.get_mrca_node(*set_n)

        # get the sptree node that is mrca to the taxon set one gtree node up
        set_u = set(node.up.get_leaf_names())
        mrca_u = species_tree.get_mrca_node(*set_u)

        # get the sptree node that is mrca to the taxon set one gtree node up without this node
        set_o = set_u - set_n
        mrca_o = species_tree.get_mrca_node(*set_o)
        
        # 1: how much farther up the species tree does it cause to go up one gtree node?
        # higher values suggest this node is discordant with sptree (add this to score)
        dist_1 = species_tree.distance.get_node_distance(mrca_u, mrca_n, topology_only=topology_only)

        # 2: how much less far would it be if the orig gtree node was not here?
        # higher values suggest this node is not that far out of place (subtract from score)
        dist_2 = species_tree.distance.get_node_distance(mrca_u, mrca_o, topology_only=topology_only)
        
        # record the species tree nodes 
        node.mst = (mrca_n.name, mrca_u.name)

        # store index to node
        node.alien = max(0, dist_1 - dist_2)
    return gene_tree        



def alien_tree_metric2(species_tree, gene_tree, topology_only: bool = False) -> float:
    """...
    """
    for node in gene_tree[:-1]:

        # get the sptree node that is mrca to this taxon set
        mrca_st = species_tree.get_mrca_node(*node.get_leaf_names())

        # get the sptree node that is mrca to the taxon set one gtree node up
        mrca_up_st = species_tree.get_mrca_node(*node.up.get_leaf_names())
        node.mst = (mrca_st.name, mrca_up_st.name)

        # get the dist between the sptree mrca nodes
        dist_st = max([0.001, mrca_up_st.height - mrca_st.height])
        dist_st = species_tree.distance.get_node_distance(mrca_up_st, mrca_st, topology_only=topology_only)

        # get the dist between node and node.up in gene tree
        dist_gt = 1. if topology_only else node.dist

        # measure clade alien index
        node.alien = dist_st / dist_gt
    return gene_tree
        


def main():

    import toytree
    # create species tree
    stree = toytree.rtree.unittree(ntips=20, seed=1234).ladderize()

    # create gene tree with a duplication
    clade = stree.mod.extract_subtree("r0", "r8")
    gtree = stree.mod.add_internal_node_and_subtree("r0", "r8", subtree=clade)

    # insert a node as HGT into distant clade
    # hclade = gtree.mod.extract_subtree("r18", "r19")
    # gtree = gtree.mod.add_internal_node_and_subtree("r0", "r8", subtree=hclade)
    # stree.treenode.draw_ascii()
    # gtree.treenode.draw_ascii()



if __name__ == "__main__":

    main()
