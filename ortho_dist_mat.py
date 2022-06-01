import treeswift
from itertools import product

class OrthoDistMat:

    def _init_matrix(self):
        return {t1:{t2:0 for t2 in self.taxaset} for t1 in self.taxaset}


    def _mats_from_tree(self, tree, delimiter):
        assert all(hasattr(u, 'tag') for u in tree.traverse_postorder(leaves=False)), "Tree must be rooted and tagged"
        # Memoize depth of each vertex
        # (here we say that duplication vertices don't contribute to internode distance)
        tree.root.depth = 0
        for u in tree.traverse_preorder():
            if not u.is_root():
                u.depth = u.parent.depth + 1 if u.is_leaf() or u.tag == 'S' else u.parent.depth

        M, N = self._init_matrix(), self._init_matrix() # Dist Matrix and Number of Paths Matrix

        for u in tree.traverse_postorder():
            if u.is_leaf():
                u.below = [u]
            else:
                [left, right] = u.child_nodes()
                if u.tag == "S":
                    for l1, l2 in product(left.below, right.below):
                        t1, t2 = l1.get_label().split(delimiter)[0], l2.get_label().split(delimiter)[0]
                        M[t1][t2] += l1.depth + l2.depth - 2 * u.depth
                        N[t1][t2] += 1
                        M[t2][t1] = M[t1][t2]
                        N[t2][t1] = N[t1][t2]
                u.below = left.below + right.below
                delattr(left, 'below'); delattr(right, 'below')

        return M, N


    def __init__(self, tree=None, delimiter=None):
        self.taxaset = set()
        self.dist_mat = {}
        self.n_paths_mat = {}

        if not tree is None:
            self.taxaset = {l.get_label().split(delimiter)[0] for l in tree.traverse_leaves()}
            self.dist_mat, self.n_paths_mat = self._mats_from_tree(tree, delimiter)    


    def dist(self, taxa1, taxa2):
        return self.dist_mat[taxa1][taxa2], self.n_paths_mat[taxa1][taxa2]

    
    def average_mat(self):
        div = lambda x:x[0]/x[1] if x[1] != 0 else 0
        return {t1:{t2:div(self.dist(t1,t2)) for t2 in self.taxaset} for t1 in self.taxaset}


    def phylip(self):
        avg_mat = self.average_mat()
        result = "{}\n".format(len(self.taxaset))
        for t in self.taxaset:
            #print(t, avg_mat[t].values())
            result += "{} {}\n".format(t, " ".join(map(str, avg_mat[t].values())))
        return result


    def __add__(self, other):
        result = OrthoDistMat()
        result.taxaset = self.taxaset.union(other.taxaset)
        result.dist_mat = {t1:{t2:self.dist_mat.get(t1, {}).get(t2, 0) + other.dist_mat.get(t1, {}).get(t2, 0) for t2 in result.taxaset} for t1 in result.taxaset}
        result.n_paths_mat = {t1:{t2:self.n_paths_mat.get(t1, {}).get(t2, 0) + other.n_paths_mat.get(t1, {}).get(t2, 0) for t2 in result.taxaset} for t1 in result.taxaset}
        return result
