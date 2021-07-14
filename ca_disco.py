from disco import *
import treeswift as ts
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import argparse

from argparse import ArgumentParser, ArgumentTypeError
import re

def fileList(path):
    return [line.strip() for line in open(path).readlines()]

def parseTaxonList(string):
    if os.path.isfile(string): 
        return [line.strip() for line in open(string).readlines()]
    # https://stackoverflow.com/a/6512463
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise ArgumentTypeError("'" + string + "' is not a range of number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start,10), int(end,10)+1))


def retrieve_alignment(tre, alnpath, taxonset=range(0,101), delimiter='_'):
    """
    Parameters
    ----------------
    tre : single-copy treeswift tree generated from James's code
    alnpath : path to the phylip formatted alignment of the genes. The row labels should be a superset of the leafset of 'tre'
    seqlen : sequence length parameter, only the first seqlen columns are taken from the MSA
    taxonset: set, the taxon set of the entire dataset

    Returns the MSA that corresponds to the input tree.
    """
    aln = AlignIO.read(open(alnpath), "phylip") 
    seqlen = len(aln[0].seq)
    blank = "-" * seqlen
    whitelist = set(tre.labels(True, False))
    rest = set(taxonset)
    #print(rest)
    res = MultipleSeqAlignment([])
    for r in aln[:,:seqlen]:
        if r.id in whitelist:
            rid = r.id.split(delimiter)[0]
            rid_i = rid
            res.append(SeqRecord(r.seq, id=rid))
            rest.remove(rid_i)
    for rst in rest:
        res.append(SeqRecord(Seq(blank), id=str(rst)))
    res.sort()
    return res

def format_phy(i, n):
    return str(i).zfill(n) + ".phy"

if __name__=="__main__":

    parser = argparse.ArgumentParser(description = "generate concatenation files from gene-family trees using decomposition strategies")

    # some of these arguments are directly copied from James's code
    parser.add_argument("-i", "--input", type=str,
                        help="Input tree list file", required=True)
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output tree list file")
    parser.add_argument("-a", "--alignment", type=fileList, required=True,
                        help="File containing paths to all alignment files in order the genes are found in the input newick file")
    parser.add_argument('-d', '--delimiter', type=str, default='_',
                        help="Delimiter separating taxon label from the rest of the leaf label.")
    parser.add_argument('-m', '--filter', type=int, default=4,
                        help="Exclude decomposed trees with less then X taxa")
    #parser.add_argument("-k", "--ngenes", type=int, default = math.inf, help="maximum number of input gene trees to use")
    #parser.add_argument("-d", "--decomp", type=str, default = "s", help = "decomposition method, either s or d for sampling (linear) or decomposition")
    #parser.add_argument("-l", "--seqln", type=int, default = math.inf, help="maximum kept sequence length of alignments")
    parser.add_argument("-t", "--taxonset", type=parseTaxonList, required=True, 
                        help="taxon set (a range of numbers in the format of a-b) or file listing taxa (each label alone on one line)")
    args = parser.parse_args()


    # arguments: n (seqlength), alnpath, inpath, k (genetree limit), decompm: either "s" or "d"
    #N = args.seqln # seqlength
    ALNROOT = args.alignment # alignment directory, containing zero-padded files like 0001.phy ... 0030.phy that corresponds to the alignment of the genes
    INPATH = args.input # input multi-copy gene family trees path
    #K = args.ngenes # maximum number of multi-copy trees to use from INPATH
    #DECOMPM = args.decomp # decomposition method, either "s" or "d" for sampling or decomposition
    LEAFSET = args.taxonset
    OUTPUT = args.output #or INPATH + f".{K}.{DECOMPM}.aln"
    delim = args.delimiter
    F = args.filter

    num_genes = sum(1 for _ in open(INPATH))

    with open(INPATH) as fh:
        aln = None
        for i, l in enumerate(fh, start=1):
            t = ts.read_tree_newick(l)
            t.reroot(get_min_root(t, delim)[0])
            tag(t, delim)
            out = list(filter(lambda x:x.num_nodes(internal=False) >= F, decompose(t)))
            #phy = os.path.join(ALNROOT, format_phy(i, len(str(num_genes))))
            phy = ALNROOT[i - 1]
            for ot in out:
                if aln == None:
                    aln = retrieve_alignment(ot, phy, LEAFSET, delimiter=delim)
                else:
                    aln += retrieve_alignment(ot, phy, LEAFSET, delimiter=delim)
        AlignIO.write(aln, OUTPUT, "phylip")
