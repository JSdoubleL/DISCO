import treeswift
import argparse
import os

def unroot(tree):
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


def remove_in_paralogs(tree, delimiter=None):
    """
    Removes in-paralogs from unrooted tree.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

    Returns number of in-paralogs removed
    """
    # root tree if not rooted
    if tree.root.num_children() != 2:
        tree.reroot(tree.root)

    num_paralogs = 0
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([node.get_label().split(delimiter)[0]])
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


def get_min_root(tree, delimiter=None, verbose=False):
    """
    Calculates the root with the minimum score.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

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
        tree.reroot(tree.root)
    tree.resolve_polytomies()

    # Get down scores pass
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = set([node.get_label().split(delimiter)[0]])
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
    ties = [best_root]

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
                ties.append(node)
                
            if total_score < min_score:
                num_ties = 0
                min_score = total_score
                best_root = node
                ties = [node]

    if verbose:            
        print('Best root had score', min_score, 'there were', len(ties), 'ties.')
        
    return best_root, min_score, ties


def tag(tree, delimiter=None):
    """
    Tags tree according to its current rooting.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label
    """
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = set([node.get_label().split(delimiter)[0]])
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


def decompose(tree, single_tree=False):
    """
    Decomposes a tagged tree, by separating clades at duplication vetices

    NOTE: must be run after 'tag()'

    Parameters
    ----------
    tree: tagged treeswift tree
    single_tree: return only the single large

    Returns result of the decomposition as a list of trees
    """
    out = []
    root = tree.root
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


def parse_notung_gtree(gtree_file):
    tree = treeswift.read_tree_newick(gtree_file)
    tree.n_dup = 0
    for u in tree.traverse_postorder():
        if u.is_leaf():
            u.s = set([u.get_label()])
        else: 
            [left, right] = u.child_nodes()
            u.s =  left.s.union(right.s)
            if hasattr(u, 'edge_params') and 'D=Y' in u.edge_params:
                u.tag = 'D' 
                tree.n_dup += 1
            else: 
                u.tag = 'S'
                #if hasattr(u, 'edge_params') and 'D=Y' in u.edge_params else 'S'
    return tree


def run_notung(gtree, stree_path, notung_path, dup_cost=1.5,  loss_cost=1):
    outdir = stree_path.rsplit(os.sep, 1)[0] + os.sep if os.sep in stree_path else ''
    tmp_file = outdir + 'tmp.tree'
    notung_output = tmp_file + '.rooting.0'

    # confirm species tree is rooted.
    stree = treeswift.read_tree_newick(stree_path)
    if stree.root.num_children() != 2:
        raise Exception("Species tree must be rooted.")

    # clean species tree
    for u in stree.traverse_postorder(leaves=False):
        u.label = None
        u.edge_length = None

    prefix, postfix = stree_path.rsplit('.', 1)
    clean_stree_path = "{}-cleaned.{}".format(prefix, postfix)
    with open(clean_stree_path, 'w') as f:
        f.write(stree.newick())

    # Write gene tree to temp file so Notung can read it 
    gtree.resolve_polytomies()
    with open(tmp_file, 'w') as f:
        f.write(gtree.newick())

    # run Notung
    command = 'java -jar {} {} -s {} --root --infertransfers false --log --treeoutput ' \
        'nhx --speciestag prefix --costdup {} --costloss {} --nolosses' \
        .format(notung_path, tmp_file, clean_stree_path, dup_cost, loss_cost)
    if outdir != '': command += ' --outputdir {}'.format(outdir)
    os.system(command)

    return parse_notung_gtree(notung_output)


def relabel(tree, delimiter=None):
    if delimiter is None:
        return tree
    for l in tree.traverse_postorder(internal=False):
        l.set_label(l.get_label().split(delimiter)[0])
    return tree


def main(args):

    if args.output is None:
        split = args.input.rsplit('.', 1)
        output = split[0] + '-decomp.' + split[1]
    else:
        output = args.output

    # delete existing outgroup file (so you don't append to it)
    outgroup_file_name = args.input.rsplit('.', 1)[0] + '_outgroups.txt'
    if args.outgroups:
        open(outgroup_file_name, 'w').close()

    dup_stats = 'num-dups,num-outputted,mean-size,max-size\n'

    with open(args.input, 'r') as fi, open(output, 'w') as fo:
        for i, line in enumerate(fi, 1):
            tree = treeswift.read_tree_newick(line)

            if args.remove_in_paralogs:
                num_paralogs = remove_in_paralogs(tree, args.delimiter)

            if args.species_tree is None:
                root, score, ties = get_min_root(tree, args.delimiter)
                tree.reroot(root)
                tag(tree, args.delimiter)

                n_dups = tree.n_dup # record number of duplication events

                if args.verbose:
                    print('Tree ', i, ': Tree has ', len(tree.root.s), ' species.', sep='')
                    if args.remove_in_paralogs:
                        print(num_paralogs, 'in-paralogs removed prior to rooting/scoring.')  
                    if len(tree.root.s) < 2:
                        print('Uninformative')
                    elif tree.n_dup == 0:
                        print('Single-Copy')                 
                    else:
                        outgroup = min((len(child.s), child.s) for child in tree.root.child_nodes())                    
                        print('Best root had score ', score, ' with ', tree.n_dup, ' non-terminal' if args.remove_in_paralogs else '',
                            ' duplications; there were ', len(ties), ' ties.\nOutgroup: {',','.join(outgroup[1]),'}', sep='')
            else: # Notung rooting
                tree = run_notung(tree, args.species_tree, args.notung_path, args.dup_cost, args.loss_cost)
                if args.dup_stats:
                    n_dups = tree.n_dup

            # Choose modes
            if args.no_decomp:
                out = [tree]
            else:
                out = list(filter(lambda x:x.num_nodes(internal=False) >= args.minimum, decompose(tree, args.single_tree)))

            # Output trees
            for t in out:
                if not args.no_decomp: unroot(t)
                if args.relabel: relabel(t, args.delimiter)
                t.suppress_unifurcations()
                fo.write(t.newick() + '\n')
            
            if args.verbose:
                print('Decomposition strategy outputted', len(out), 'tree(s) with minimum size', args.minimum, '.\n')

            if args.dup_stats:
                sizes = [t.num_nodes(internal=False) for t in out]
                dup_stats += '{},{},{},{}\n'.format(n_dups, len(out), sum(sizes) / len(out) if len(out) != 0 else 0, 
                                                    max(sizes) if len(out) != 0 else 0)

            # output outgroups
            if args.outgroups:
                og_tree = treeswift.read_tree_newick(line)
                root, score, ties = get_min_root(og_tree, args.delimiter)
                og_tree.reroot(root)
                tag(og_tree, args.delimiter)
                if len(og_tree.root.s) >= 2 and og_tree.n_dup >= 1:
                    with open(outgroup_file_name, 'a') as outgfile:
                        outgfile.write('Tree ' + str(i) + ':\n')
                        for t in ties:
                            og_tree.reroot(t)
                            tag(og_tree, args.delimiter)
                            outgroup = min((len(child.s), child.s) for child in og_tree.root.child_nodes())
                            outgfile.write('{' + ','.join(outgroup[1]) + '}\n')

        open(args.input.rsplit('.', 1)[0] + '_stats_{}.csv'.format('notung' if args.species_tree is not None else 'apro'), 'w').write(dup_stats)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree list file")
    parser.add_argument('-d', "--delimiter", type=str, 
                        help="Delimiter separating species name from rest of leaf label")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enables verbose output")
    parser.add_argument('-m', "--minimum", type=int, 
                        help="Minimum tree size outputted", default=4)
    parser.add_argument('-s', '--species-tree', type=str, 
                        help="Species tree for reconciliation rooting using Notung")
    parser.add_argument('--relabel', action='store_true', 
                        help="Overwrite current labels with species label in output (used with delimiter).")
    parser.add_argument('--dup-cost', type=float, default=1.5,
                        help="Duplication cost (only works with Notung rooting for now).")
    parser.add_argument('--loss-cost', type=float, default=1,
                        help="Loss cost (only works with Notung rooting for now).")
    parser.add_argument('--notung-path', type=str, default='./Notung-2.9.1.5.jar',
                        help="Path to Notung jar file")
    parser.add_argument("--outgroups", action='store_true',
                        help="Output outgroups to file (including ties)")
    parser.add_argument("--dup-stats", action='store_true',
                        help="Output file listing number of duplication events detected per gene tree and number of trees outputted.")
    parser.add_argument("--remove_in_paralogs", action='store_true',
                        help="Remove in-paralogs before rooting/scoring tree.")
    parser.add_argument('--no-decomp', action='store_true', 
                        help="Outputs rooted trees without decomposition")
    parser.add_argument('--single-tree', action='store_true',
                        help="Only output single large tree")

    main(parser.parse_args())
