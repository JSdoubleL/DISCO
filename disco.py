from collections import deque
import treeswift
import argparse
from itertools import combinations

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


def reroot_on_edge(tree, node):
    if not hasattr(node, 'edge_length') or node.edge_length is None or node.edge_length == 0:
        node.edge_length = 1
    tree.reroot(node, length=node.edge_length / 2)
    

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


def preprocess_tree(tree, delimiter=None, rand_resolve=False):
    max_degree = 0
    if rand_resolve:
        tree.resolve_polytomies()
    tree.suppress_unifurcations()

    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])

    assert tree.root.num_children() == 2

    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = set([node.get_label().split(delimiter)[0]])
        else:
            node.down = set().union(*[u.down for u in node.child_nodes()])
            if node.num_children() + 1 > max_degree:
                max_degree = node.num_children() + 1

    assert tree.root.num_children() == 2
    [left, right] = tree.root.child_nodes()
    left.up, right.up = right.down, left.down

    for node in tree.traverse_preorder():
        parent = node.get_parent()
        if parent is not None and parent is not tree.root:
            node.up = parent.up.union(*[c.down for c in parent.child_nodes() if c is not node])

    return max_degree


def get_min_root_old(tree, loss_cost=1, delimiter=None, verbose=False):
    """
    Calculates the root with the minimum score.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

    Returns vertex corresponding to best edge to root tree on
    """

    def score(total_set, set1, set2, loss_cost):
        """DL-score function"""
        if not len(set1.intersection(set2)) == 0:
                if total_set == set1 or total_set == set2:
                    if set1 == set2:
                        return 1
                    else:
                        return 1 + loss_cost
                else:
                    return 1 + 2 * loss_cost
        return 0

    # check if tree is single leaf
    if tree.root.num_children() == 0:
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
            node.d_score = left.d_score + right.d_score + score(node.down, left.down, right.down, loss_cost)

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

    min_score = left.u_score + left.d_score + score(left.up.union(left.down), left.up, left.down, loss_cost)
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
            node.u_score = parent.u_score + other.d_score + score(node.up, parent.up, other.down, loss_cost)

            total_score = node.u_score + node.d_score + score(node.up.union(node.down), node.up, node.down, loss_cost)

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


def get_min_root(tree, loss_cost=1, constraint_clades=None, verbose=False):
    """
    Calculates the root with the minimum score.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label

    Returns vertex corresponding to best edge to root tree on
    """

    def score(total_set, set1, set2, loss_cost):
        """DL-score function"""
        if not len(set1.intersection(set2)) == 0:
                if total_set == set1 or total_set == set2:
                    if set1 == set2:
                        return 1
                    else:
                        return 1 + loss_cost
                else:
                    return 1 + 2 * loss_cost
        return 0

    def score_polytomy(v, constraint_clades, loss_cost):
        """Dynamic Programing polytomy scoring algorithm"""
        if v.parent.is_root(): # no real value for 'up' at the root; the correct value in this case should be the other side of tree
            v.parent.up = v.child_nodes()[0].down if v.child_nodes()[1] is v else v.child_nodes()[1].down

        neighbors = v.child_nodes() + [v.parent]
        m = len(neighbors)
        #print(m)
        assert m > 3, '{} is not a polytomy'.format(v.get_label())

        C = {}
        A_v = []
        #solutions = [frozenset(frozenset(neighbors) - {u}) for u in neighbors] # force all possible solutions
        for i in range(1, m):
            for A in combinations(neighbors, i):
                temp = frozenset().union(*[u.down if u is not v.parent else u.up for u in A])         
                # TODO: implement constraints correctly       
                if True: #constraint_clades is None or temp in constraint_clades or frozenset(A) in solutions:
                    A_v.append(frozenset(A))
                    C[A_v[-1]] = temp

        #A_v = [frozenset(A) for i in range(1, m) for A in combinations(neighbors, i) if frozenset().union(u.s for u in A) in constraint_clades]
        v.M = {frozenset({u}):0 for u in neighbors}
        v.backtrace = {}

        for i, A in enumerate(A_v):
            if i >= len(neighbors): 
                min_score = float('inf')
                for A_1 in A_v[:i]:
                    A_2 = A - A_1
                    if A_2 in A_v[:i] and A_1.union(A_2) == A:                                      
                        cur_score = score(C[A], C[A_1], C[A_2], loss_cost) + v.M[A_1] + v.M[A_2] 
                        if cur_score < min_score:
                            min_score = cur_score
                            v.backtrace[A] = A_1
                v.M[A] = min_score
        return v.M[frozenset(neighbors) - {v.parent}] + sum(u.d_score for u in v.child_nodes())


    def polytomy_child_score(v):
        """Calculate the "up score" for the child of a polytomy"""
        poly = v.parent
        neighbors = poly.child_nodes() + [poly.parent]
        assert poly.M[frozenset(neighbors) - {v}] != float('inf')
        return poly.M[frozenset(neighbors) - {v}] + sum(u.d_score for u in poly.child_nodes() if u is not v) + poly.u_score


    """Main "get_min_root" code"""
    # check if tree is single leaf
    if tree.root.num_children() == 0:
        return tree.root, 0, [] 

    assert tree.root.num_children() == 2
    # root tree if not rooted
    #if tree.root.num_children() != 2:
    #    reroot_on_edge(tree, tree.root.get_children[0])

    # Get down scores pass
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.d_score = 0
        elif node.num_children() == 2:
            [left, right] = node.child_nodes()
            node.d_score = left.d_score + right.d_score + score(node.down, left.down, right.down, loss_cost)
        else:
            node.d_score = score_polytomy(node, constraint_clades, loss_cost)
    
    [left, right] = tree.root.child_nodes()
    left.u_score, right.u_score = right.d_score, left.d_score
    min_score = left.u_score + left.d_score + score(left.up.union(left.down), left.up, left.down, loss_cost)
    # we don't want to root at a leaf
    if not left.is_leaf():
        best_root = left 
    elif not right.is_leaf():
        best_root = right
    # if both are leaves (i.e. two leaf tree), we want to keep the same rooting
    else:
        best_root = tree.root
    ties = [best_root]

    for node in tree.traverse_preorder():
        parent = node.get_parent()
        if parent is not None and parent is not tree.root:           
            if parent.num_children() == 2: # no polytomies
                # find the other child of parent
                other = parent.child_nodes()[0] if parent.child_nodes()[0] != node  else parent.child_nodes()[1]
                node.u_score = parent.u_score + other.d_score + score(node.up, parent.up, other.down, loss_cost)
                total_score = node.u_score + node.d_score + score(node.up.union(node.down), node.up, node.down, loss_cost)                
            else: # polytomies
                min_total = min_score + 1
                for v in parent.child_nodes():
                    v.u_score = polytomy_child_score(v)
                    total_score = v.u_score + v.d_score + score(v.up.union(v.down), v.up, v.down, loss_cost)
                    if total_score < min_total:
                        min_total = total_score
                total_score = min_total

            if total_score == min_score:
                    ties.append(node)
            if total_score < min_score:
                min_score = total_score
                best_root = node
                ties = [node]

    if verbose:            
        print('Best root had score', min_score, 'there were', len(ties), 'ties.')
        
    return best_root, min_score, ties


def backtrace_polytomies(tree):
    def backtrace(polytomy):
        q = deque(); q.append(polytomy)
        while len(q) != 0:
            v = q.popleft()
            A_1 = polytomy.backtrace[frozenset(v.child_nodes())]
            A_2 = frozenset(v.child_nodes()) - A_1
            n1, n2 = treeswift.Node(), treeswift.Node()
            if len(A_1) == 1: 
                n1, = A_1 # comma unpacks the first (and only) element
            else:
                for a in A_1:
                    n1.add_child(a)
            if len(A_2) == 1: 
                n2, = A_2
            else:
                for a in A_2:
                    n2.add_child(a)
            if len(A_1) > 2:
                q.append(n1)
            if len(A_2) > 2:
                q.append(n2)
            v.children = []
            v.add_child(n1)
            v.add_child(n2)
        #assert all(u == u.parent.child_nodes()[0] for u in tree.traverse_postorder() if u.num_children() == 1)

    for u in tree.traverse_preorder():
        if u.num_children() > 2:
            if not hasattr(u, 'backtrace'):
                t = tree.extract_subtree(u.parent) 
                print(t.newick())
            backtrace(u)


def tag(tree, delimiter=None):
    """
    Tags tree according to its current rooting.

    Parameters
    ----------
    tree: treeswift tree
    delimiter: delimiter separating species name from rest of leaf label
    """
    assert all([u.num_children() == 2 for u in tree.traverse_postorder(leaves=False)])
    #tree.suppress_unifurcations()
    #tree.resolve_polytomies()
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


def trivial(newick_str):
    """
    Determines if a newick string represents a trivial tree (tree containing no quartets).

    Parameters
    ----------
    newick_str: newick string

    Returns True if tree contains less than two '('
    """
    count = 0
    for c in newick_str:
        if c == '(':
            count += 1
        if count > 1:
            return False
    return True


def get_tree_clades(tree, delimiter=None):
    """
    Lists all clades in a tree. Works with multi-copy trees.

    Parameters
    ----------
    tree: treeswift tree

    Returns set of clades (as frozenset objects) found in the input gene tree
    """
    assert max(u.num_children() for u in tree.traverse_postorder()) <= 2, \
        "Input trees contain polytomies. Cannot extract clades" 

    clades = set()
    for u in tree.traverse_postorder():
        if u.is_leaf():
            u.down = frozenset({u.get_label().split(delimiter)[0]})            
        else:
            [left, right] = u.child_nodes()
            u.down = left.down.union(right.down)
        clades.add(frozenset(u.down))

    [left, right] = tree.root.child_nodes()
    left.up, right.up = right.down, left.down

    for u in tree.traverse_preorder():
        parent = u.get_parent()
        if parent is not None and not parent.is_root():            
            [left, right] = parent.child_nodes()
            u.up = parent.up.union(left.down) if u is right else parent.up.union(right.down)
            clades.add(frozenset(u.up))
        
    return clades


def relabel(tree, delimiter=None):
    if delimiter is None:
        return tree
    for l in tree.traverse_postorder(internal=False):
        l.set_label(l.get_label().split(delimiter)[0])
    return tree


def contract_low_support_with_max(tree, threshold, max_degree):
    low_support = []
    for u in tree.traverse_postorder(leaves=False):
        try:
            if float(str(u)) < threshold:
                low_support.append(u)
        except:
            pass
    low_support.sort(key=lambda x:float(str(x)))
    for u in low_support:
        if u.parent is not None and u.parent.num_children() + u.num_children() <= max_degree:
            u.contract()


def max_degree(tree):
    return max(u.num_children() + 1 if not u.is_root() else u.num_children() for u in tree.traverse_postorder())


def main(args):
    # set output file name
    if args.output is None:
        split = args.input.rsplit('.', 1)
        output = split[0] + '-decomp.' + split[1]
    else:
        output = args.output

    # delete existing outgroup file (so you don't append to it)
    outgroup_file_name = args.input.rsplit('.', 1)[0] + '_outgroups.txt'
    if args.outgroups:
        open(outgroup_file_name, 'w').close()

    # create constraint set of clades from input gene trees
    clades = set()
    if not args.classic:
        # if the input trees are not fully resolved, randomly resolve them to get constraints (to ensure a solution will exist)
        for gtree in treeswift.read_tree_newick(args.input):
            gtree.resolve_polytomies()
            clades.update(get_tree_clades(gtree, args.delimiter))

    if args.verbose and not args.classic:
        print('constraint set', len(clades))

    with open(args.input, 'r') as fi, open(output, 'w') as fo:
        for i, line in enumerate(fi, 1):
            tree = treeswift.read_tree_newick(line)

            if args.remove_in_paralogs:
                num_paralogs = remove_in_paralogs(tree, args.delimiter)

            if max_degree(tree) > args.max_degree:
                tree.resolve_polytomies() # will set support of new branches to 0

            if args.threshold is not None:
                contract_low_support_with_max(tree, args.threshold, args.max_degree)

            if args.classic:
                tree.resolve_polytomies()
                root, score, ties = get_min_root_old(tree, args.loss_cost, args.delimiter)
                reroot_on_edge(tree, root)
            else:
                max_poly = preprocess_tree(tree, args.delimiter, args.random)
                root, score, ties = get_min_root(tree, args.loss_cost, clades if len(clades) != 0 else None)                

                reroot_on_edge(tree, root)
                backtrace_polytomies(tree)
                tree.suppress_unifurcations() # remove previous root node

            tag(tree, args.delimiter)

            if args.verbose:
                print('Tree ', i, ': Tree has ', len(tree.root.s), ' species.', sep='')
                if not args.classic:
                    print('Max polytomy size is', max_poly)
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

            # Choose modes
            if args.no_decomp:
                out = [tree]
            else:
                out = list(filter(lambda x:x.num_nodes(internal=False) >= args.minimum, decompose(tree, args.single_tree)))

            # Output trees
            for t in out:
                if not args.no_decomp: unroot(t)
                if not args.keep_original_labels: relabel(t)
                t.suppress_unifurcations()
                fo.write(t.newick() + '\n')
            
            if args.verbose:
                print('Decomposition strategy outputted ', len(out), ' tree(s) with minimum size ', args.minimum, '.\n', sep='')

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
            

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output tree list file")
    parser.add_argument('-d', "--delimiter", type=str, 
                        help="Delimiter separating species name from rest of leaf label")
    parser.add_argument('-l', '--loss-cost', type=float, default=1,
                        help="Lost cost relative to duplication cost")
    parser.add_argument('-n', '--no-decomp', action='store_true', 
                        help="Outputs rooted trees without decomposition")
    parser.add_argument('-s', '--single_tree', action='store_true',
                        help="Only output single large tree")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enables verbose output")
    parser.add_argument('-m', '--max-degree', type=int, default=10,
                        help="Max allowed polytomy degree")
    parser.add_argument('-t', '--threshold', type=float, 
                        help="Support threshold for collapsing gene tree branches")
    parser.add_argument('-r', "--random", action='store_true',
                        help="Resolve polytomies randomly (faster)")
    parser.add_argument("--minimum", type=int, 
                        help="Minimum tree size outputted", default=4)
    parser.add_argument('-k', "--keep-original-labels", action='store_true', 
                        help="Keep original leaf labels instead of relabling them with their species labels (only relevent with delimiter)")
    parser.add_argument("--outgroups", action='store_true',
                        help="Output outgroups to file (including ties)")
    parser.add_argument('-rp', "--remove_in_paralogs", action='store_true',
                        help="Remove in-paralogs before rooting/scoring tree.")
    parser.add_argument('--classic', action='store_true', 
                        help="Use the old version of DISCO")

    main(parser.parse_args())
