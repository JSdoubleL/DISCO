import disco
from ortho_dist_mat import OrthoDistMat
import argparse
import treeswift

def main(args):
    output_name = args.output
    result = OrthoDistMat()
    if args.output is None:
        output_name = "{}.phy".format(args.input.rsplit('.', 1)[0])
    with open(args.input, "r") as f:        
        for line in f:
            tree = treeswift.read_tree_newick(line)
            best_root, *_ = disco.get_min_root(tree, args.delimiter)
            tree.reroot(best_root)
            disco.tag(tree, args.delimiter)
            result += OrthoDistMat(tree, args.delimiter)
    with open(output_name, "w") as f:
        f.write(result.phylip())

            
if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', "--input", type=str, 
                        help="Input newick gene trees", required=True)
    parser.add_argument('-o', "--output", type=str, 
                        help="Output matrix file name")
    parser.add_argument('-d', "--delimiter", type=str, 
                        help="Delimiter separating species name from rest of leaf label")

    main(parser.parse_args())
