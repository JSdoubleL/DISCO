import csv
import os
from typing import Any, Callable, Dict, Iterable, List, Set, Tuple, Union
import treeswift as ts
import argparse

def unroot(tree: ts.Tree) -> ts.Tree:
    """
    Unroots treeswift tree. Adapted from treeswift 'deroot' function.
    This one doesn't contract (A,B); to A;

    Parameters
    ----------
    tree: treeswift tree

    Returns unrooted treeswift tree
    """
    if tree.root == None:
        return tree
    if tree.root.num_children() == 2:
        [left, right] = tree.root.child_nodes()
        if not right.is_leaf():
            right.contract()
        elif not left.is_leaf():
            left.contract()
    tree.is_rooted = False
    return tree


def reroot_on_edge(tree: ts.Tree, node: ts.Node) -> None:
    """
    In order for DISCO to behave as expected, the root must bisect the edge
    above the node we reroot on.

    Parameters
    ----------
    tree: treeswift tree
    node: node from safe treeswift tree
    """
    tree.root.edge_length = None # prevent creation of leaf called ROOT
    if not node.is_root(): 
        if not hasattr(node, 'edge_length') or node.edge_length is None or node.edge_length == 0:
            node.edge_length = 1
        tree.reroot(node, length=node.edge_length / 2)


def remove_in_paralogs(tree: ts.Tree, gene_to_species: Callable[[str], str] = lambda x:x) -> int:
    """
    Removes in-paralogs from unrooted tree.

    Parameters
    ----------
    tree: treeswift tree
    gene_to_species: map from gene labels to species labels

    Returns number of in-paralogs removed
    """
    # root tree if not rooted
    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])

    num_paralogs = 0
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([gene_to_species(node.get_label())])
        else:
            node.s = set([])
            for child in node.child_nodes():
                node.s = node.s.union(child.s)
            # collapse if paralogs
            if len(node.s) == 1:
                for child in node.child_nodes()[1:]:
                    node.remove_child(child)
                    num_paralogs += 1

    # check over top of root
    if any(len(child.s) == 1 for child in tree.root.child_nodes()):
        for node in tree.traverse_preorder():
            if not node.is_root():
                parent = node.get_parent()
                node.up = set([]) if parent.is_root() else parent.up
                for sibl in parent.child_nodes():
                    if sibl != node:
                        node.up.union(sibl.s)
                if len(node.up) == 1:
                    for sibl in parent.child_nodes():
                        if sibl != node:
                            parent.remove_child(sibl)
                            num_paralogs += 1
    tree.suppress_unifurcations()
    return num_paralogs


def get_min_root(tree: ts.Tree, gene_to_species: Callable[[str], str] = lambda x:x, verbose: bool = False
                ) -> Tuple[ts.Node, int, List[Tuple[ts.Node, ts.Node]]]:
    """
    Calculates the root with the minimum score.

    Parameters
    ----------
    tree: treeswift tree
    gene_to_species: map from gene labels to species labels

    Returns vertex corresponding to best edge to root tree on
    """
    def score(total_set, set1, set2):
        if not len(set1.intersection(set2)) == 0:
                if total_set == set1 or total_set == set2:
                    if set1 == set2:
                        return 1
                    else:
                        return 2
                else:
                    return 3 
        return 0

    # check if tree is single leaf
    if tree.root.num_children() == 0:
        tree.root.s = set([tree.root.get_label()])
        return tree.root, 0, [] 

    # root tree if not rooted
    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])
    tree.resolve_polytomies()

    # Get down scores pass
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = set([gene_to_species(node.get_label())])
            node.d_score = 0
        else:
            if node.num_children() != 2:
                raise Exception("Vertex has more than 2 children")

            [left, right] = node.child_nodes()
            node.down = left.down.union(right.down)
            node.d_score = left.d_score + right.d_score + score(node.down, left.down, right.down)

    min_score, best_root, ties = float("inf"), None, []

    # Get scores above edge pass
    for node in tree.traverse_preorder():
        if node.is_root():
            root = node
            root.skip = True
        else:
            node.skip = False

    # Get 'up' set for children of root
    [left, right] = root.child_nodes()
    left.up = right.down
    left.u_score = right.d_score
    right.up = left.down
    right.u_score = left.d_score
    left.skip = True
    right.skip = True

    min_score = left.u_score + left.d_score + score(left.up.union(left.down), left.up, left.down)
    # we don't want to root at a leaf
    if not left.is_leaf():
        best_root = left 
    elif not right.is_leaf():
        best_root = right
    # if both are leaves (i.e. two leaf tree), we want to keep the same rooting
    else:
        best_root = root
    ties = [(best_root, left)]

    for node in tree.traverse_preorder(leaves=False):
        if not node.skip:
            
            parent = node.get_parent()
            if parent.child_nodes()[0] != node:
                other = parent.child_nodes()[0]
            else: 
                other = parent.child_nodes()[1]

            node.up = parent.up.union(other.down)
            node.u_score = parent.u_score + other.d_score + score(node.up, parent.up, other.down)

            total_score = node.u_score + node.d_score + score(node.up.union(node.down), node.up, node.down)

            if total_score == min_score:
                ties.append((node, parent))
                
            if total_score < min_score:
                min_score = total_score
                best_root = node
                ties = [(node, parent)]

    if verbose:            
        print('Best root had score', min_score, 'there were', len(ties), 'ties.')
        
    return best_root, min_score, ties


def tag(tree: ts.Tree, gene_to_species: Callable[[str], str] = lambda x:x) -> None:
    """
    Tags tree according to its current rooting.

    Parameters
    ----------
    tree: treeswift tree
    gene_to_species: map from gene labels to species labels
    """
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([gene_to_species(node.get_label())])
            node.n_dup = 0
        else:
            [left, right] = node.child_nodes()

            node.s = left.s.union(right.s)
            node.n_dup = left.n_dup + right.n_dup
            if len(left.s.intersection(right.s)) == 0:
                node.tag = 'S'
            else: 
                node.tag = 'D'
                node.n_dup += 1
    tree.n_dup = tree.root.n_dup


def decompose(tree: ts.Tree, single_tree: bool = False) -> List[ts.Tree]:
    """
    Decomposes a tagged tree, by separating clades at duplication vertices

    NOTE: must be run after 'tag()'

    Parameters
    ----------
    tree: tagged treeswift tree
    single_tree: return only the single large

    Returns result of the decomposition as a list of trees
    """
    out = []
    for node in tree.traverse_postorder(leaves=False):
        if node.tag == 'D':
            # trim off smallest subtree (i.e. subtree with least species)
            [left, right] = node.child_nodes()
            delete = left if len(left.s) < len(right.s) else right
            if not single_tree:
                out.append(tree.extract_subtree(delete))
                out[-1].suppress_unifurcations()
            node.remove_child(delete)
    tree.suppress_unifurcations() # all the duplication nodes will be unifurcations
    out.append(tree)
    return out


def relabel(tree: ts.Tree, gene_to_species: Callable[[str], str] = lambda x:x) -> ts.Tree:
    """
    Relabels leaf labels given map. Used to relabel leaf labels 
    corresponding to genes into leaf labels corresponding to species

    Parameters
    ----------
    tree: treeswift tree
    gene_to_species: mapping function used to relabel

    Returns relabeled treeswift tree
    """
    for l in tree.traverse_postorder(internal=False):
        l.set_label(gene_to_species(l.get_label()))
    return tree


def print_stats(i: int, stats: Dict[str, Any]) -> None:
    """Given number (i) and stats dictionary, prints stats for a tree"""
    print(f"""
Tree:                   {i}
Number of Species:      {stats["num_species"]}
Number of Duplications: {stats["num_duplications"]}
Best Rooting Score:     {stats["best_score"]}
Number of Rooting Ties: {stats["num_ties"]}
Outgroup:
{", ".join(stats["outgroup"])}
""")


def write_stats(stat_list: List[Dict[str, Any]], path: str) -> None:
    """Write list of tree stats to csv file"""
    with open(path, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["tree","num_species", "num_duplications", "best_score", "num_ties", "outgroup"])
        writer.writerows([[i] + [v if v != set() else "{}" for v in row.values() ] for i, row in enumerate(stat_list)])


def calcualte_outgroups(tree: Union[str, ts.Tree], leaf_to_species: Callable[[str], str]) -> Tuple[str, List[Set[str]]]:
    """Calculates outgroups (including ties) and outputs them to a text file"""
    if isinstance(tree, str):
        tree = ts.read_tree_newick(tree)
    if isinstance(tree, list): # check for unexpected treeswift outputs
        raise ValueError(f"{tree} could not be interpreted as a tree")

    root, _, ties = get_min_root(tree, leaf_to_species)
    reroot_on_edge(tree, root)
    tag(tree, leaf_to_species)

    outgroups = []
    if tree.n_dup == 0:
        tree_type = "Single-Copy"
        outgroups.append(set())
    elif len(tree.root.s) < 2:
        tree_type = "Uninformative"
        outgroups.append(set(tree.root.s.pop()))
    else:
        tree_type = "Multicopy-Tree"
        for n1, n2 in ties:
            # tied roots are recored as two nodes as we don't know the orintation of the tree
            # when they were recorded. We need the lowest of the two.
            new_root = n1 if n1 in n2.child_nodes() else n2
            reroot_on_edge(tree, new_root)
            tag(tree, leaf_to_species)
            outgroup = min((len(child.s), child.s) for child in tree.root.child_nodes())[1]
            outgroups.append(outgroup)
    return tree_type, outgroups


def _disco(input_tree: Union[str, ts.Tree], callback: Callable[[str], Any] = lambda x:x, minimum: int = 4, 
           leaf_to_species: Callable[[str], str] = lambda x:x, 
           keep_labels: bool = False, single_tree: bool = False, no_decomp: bool = False, verbose: bool = False
          ) -> Tuple[List[ts.Tree], Dict[str, Any]]:
    """Decompose single tree with DISCO algorithm. Returns list of trees and dictionary containing stats."""

    # stats we want to record for each tree
    stats = {"num_species":-1,
            "num_duplications":-1,
            "best_score":-1,
            "num_ties":-1,
            "outgroup":[]}

    # if tree is a newick string, make it a treesiwft object
    if isinstance(input_tree, str):
        tree = ts.read_tree_newick(input_tree)
    if isinstance(tree, list): # check for unexpected treeswift outputs
        raise ValueError(f"{input_tree} could not be interpreted as a tree")

    # Root and tag steps
    root, score, ties = get_min_root(tree, leaf_to_species)
    reroot_on_edge(tree, root)
    tag(tree, leaf_to_species)

    # Record stats
    stats["num_species"] = len(tree.root.s)
    stats["num_duplications"] = tree.n_dup
    stats["best_score"] = score
    stats["num_ties"] = len(ties)
    stats["outgroup"] = min([(len(child.s), child.s) for child in tree.root.child_nodes()], default=(None,set()))[1]

    if no_decomp:
        out = [tree]
    else:
        out = list(filter(lambda x:x.num_nodes(internal=False) >= minimum, decompose(tree, single_tree)))

    # Clean output trees
    for t in out:
        if not no_decomp: unroot(t)
        if not keep_labels: relabel(t, leaf_to_species)
        t.suppress_unifurcations()

    return callback([t.newick() for t in out]), stats


def disco(input_trees: Union[str, Iterable[str]], minimum: int = 4, leaf_to_species: Callable[[str], str] = lambda x:x,
          out_path: str = None, keep_labels: bool = False, single_tree: bool = False, no_decomp: bool = False, verbose: bool = False
         ) -> Tuple[Union[List[List[str]], None], List[Dict[str, Any]]]:
    """
    Runs DISCO algorithm on tree(s).

    Decomposition Into Single-COpy gene trees (DISCO) is a method for decomposing multi-copy gene-family 
    trees while attempting to preserve orthologs and discard paralogs. These single-copy gene trees can 
    be subsequently used by methods that can estimate species trees from single-copy gene trees such as 
    ASTRAL or ASTRID in order to obtain an accurate estimation of the species tree.

    Parameters
    ----------
    - input: Either file path, newick tree string, or treeswift tree to be used as algorithm input.
    - minimum: Trees (after algorithm is run) must have >= leaves to this value to be included in the output. 
        Default 4, as that is the minimum number of leaves that are topologically meaningful if unrooted.
    - leaf_to_species: Function or other callable mapping the leaf labels to the species labels. Default (x) -> x.
    - keep_labels: Keep original tree labels; otherwise these are replaced with species labels. Default False.
    - single_tree: Just output the single largest tree instead of the full output. Default False.
    - verbose: Enables verbose output

    Returns DISCO algorithm output as a tuple with two values:
    - 2d list of newick strings; the first index corresponding to the index of the original tree
    - A list of dictionaries containing stats - one for each original tree 

    Raises TypeError if input_trees type does not correct.
    """
    pass_args = {"minimum":minimum, "leaf_to_species":leaf_to_species, "keep_labels":keep_labels, "single_tree":single_tree, 
        "no_decomp":no_decomp, "verbose":verbose}
    
    # simple - single tree - case
    if isinstance(input_trees, str) and not os.path.isfile(input_trees):
        return _disco(input_trees, **pass_args)
    elif not os.path.isfile(input_trees) and not isinstance(input_trees, Iterable):
        raise TypeError(f"{type(input_trees).__name__} is not a valid input type for disco()")
    
    # set input iterator and callback
    reader = input_trees if not os.path.isfile(input_trees) else open(input_trees)
    output_trees = []
    output_file = open(out_path, "w")
    writer = output_trees.append if out_path is None else lambda x: output_file.writelines([l + '\n' for l in x])

    stats = []
    for i, tree in enumerate(reader):
        _, stat = _disco(tree, writer, **pass_args)
        if verbose:
            print_stats(i, stat)
        stats.append(stat)
    
    if os.path.isfile(input_trees):
        reader.close()
    if out_path is not None:
        output_file.close()
    
    return None if len(output_trees) == 0 else output_trees, stats


def run_console_disco(args: argparse.Namespace) -> None:
    """Run command line DISCO"""
    # set up gene -> species map if delimiter is used
    if args.delimiter is not None:
        gene_to_species = lambda x : args.delimiter.join(x.split(args.delimiter)[:args.nth_delimiter])
    else:
        gene_to_species = lambda x : x

    # set default name for output
    input_file = args.input
    prefix, postfix = input_file.rsplit('.', 1)
    if args.output is None: 
        output_file = prefix + '-decomp.' + postfix
    else:
        prefix, postfix = args.output.rsplit('.', 1)
        output_file = args.output

    # remove imparalogs 
    if args.remove_in_paralogs:
        # add "no-inparalogs" to all file names
        prefix += '_no-inparalogs'
        no_paralogs_file = prefix + '.' + postfix
        with open(input_file, "r") as fi, open(no_paralogs_file, "w") as fo:
            for line in fi:
                if line.strip() != '':
                    tree = ts.read_tree_newick(line)
                    remove_in_paralogs(tree, gene_to_species)
                    unroot(tree)
                    fo.write(tree.newick() + '\n')
        input_file = no_paralogs_file

    # run DISCO
    _, stats = disco(input_file, args.minimum, gene_to_species, output_file, args.keep_labels, args.single_tree, 
        args.no_decomp, args.verbose)

    # write stats
    if args.stats:
        write_stats(stats, prefix + '_stats.csv')

    # calculate and write outgroups
    if args.outgroups:
        with open(input_file, "r") as fi, open(prefix + '_outgroups.txt', "w") as fo:
            for i, line in enumerate(fi):
                tree_type, outgroups = calcualte_outgroups(line, gene_to_species)
                fo.write(f"Tree {i} is {tree_type}:\n" + "\n".join("    {" + ", ".join(list(map(str, outgroup))) + "}" for outgroup in outgroups) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='====================== DISCO v1.3.2 ======================')

    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree list file")
    parser.add_argument('-d', "--delimiter", type=str, 
                        help="Delimiter separating species name from rest of leaf label")
    parser.add_argument('-n', '--nth-delimiter', type=int, # Default is 1 -- set below
                        help="Split on nth delimiter (only works with -d)")
    parser.add_argument('-m', "--minimum", type=int, 
                        help="Minimum tree size outputted", default=4)
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enables verbose output")
    parser.add_argument('--stats', action='store_true',
                        help="Write stats about decomposed trees to csv file")
    parser.add_argument("--keep-labels", action='store_true', 
                        help="Keep original leaf labels instead of relabeling them with their species labels (only relevant with delimiter)")
    parser.add_argument('--single_tree', action='store_true',
                        help="Only output single large tree")
    parser.add_argument('--no-decomp', action='store_true', 
                        help="Outputs rooted trees without decomposition")
    parser.add_argument("--outgroups", action='store_true',
                        help="Output outgroups to file (including ties)")
    parser.add_argument("--remove_in_paralogs", action='store_true',
                        help="Remove in-paralogs before rooting/scoring tree.")

    args = parser.parse_args()
    if args.delimiter is None:
        if args.nth_delimiter is not None:
            parser.error("Cannot set -n without a delimiter")
        if args.keep_labels:
            parser.error("Cannot use --keep-labels without a delimiter")
    elif args.nth_delimiter is None:
        args.nth_delimiter = 1
    if args.single_tree and args.no_decomp:
        parser.error("Cannot combine --single_tree and --no-decomp")
    if not args.verbose and args.remove_in_paralogs:
        print("--remove_in_paralogs is meaningless without --verbose, as it does not change the optimal rooting. " + 
            "It may also slow the program.")

    run_console_disco(args)
