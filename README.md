# DISCO

Decomposition Into Single-COpy gene trees (DISCO) is a method for decomposing multi-copy gene-family trees while attempting to preserve orthologs and discard paralogs. These single-copy gene trees can be subsequently used by methods that can estimate species trees from single-copy gene trees such as [ASTRAL](https://github.com/smirarab/ASTRAL) or [ASTRID](https://github.com/pranjalv123/ASTRID) in order to obtain an accurate estimation of the species tree. 

## Algorithm

Given a list of multi-copy gene trees, DISCO does the following for each tree:

1. Root the tree and tag each internal vertex as either a duplication event or a speciation event in such a way that minimizes the total number of duplications and losses. We do this with the ASTRAL-Pro rooting and tagging algorithm ([Zhang et. al. 2020](https://doi.org/10.1093/molbev/msaa139)).
2. Decompose gene tree by splitting off the smallest subtree under every vertex tagged as a duplication from the bottom up until all duplication events are resolved; it returns the set of single-copy trees produced.

## Dependencies

- Python 3
- [TreeSwift](https://github.com/niemasd/TreeSwift)

Treeswift can be installed with: `pip install treeswift`

## Usage

### disco.py

**Input**: File containing list of multi-copy trees in newick format

**Output**: File containing resulting list of single-copy trees after decomposition in newick format

```
python3 disco.py -i <input_file> -o <ouput_file> -d <delimiter>
```

#### Arguments

- **Required**
  - `-i`: Input newick tree file
- **Optional**
  - `-o`: Output newick tree file
  - `-d`: Delimiter separating species name from rest of leaf label. Default None.
  - `-s`: Output only single tree (discarding smallest duplicate clades).
  - `-m`: Minimum number of taxa required for tree to be outputted. Default 4.
  - `-n`: No decomposition (outputs rooted gene trees).
  - `-v`: Enable verbose output
  - `-k`: Keep original leaf labels (otherwise leaves are relabeled with species label; only relevant with delimiter option)
  - `-rp`: Remove in-paralogs before rooting/scoring (does not affect output, only reported score)
  - `--outgroups`: Write outgroups (including ties) to txt file. (Might make program slower).

#### Example

```cmd
python3 tag_decomp.py -i example/gtrees-mult.trees
```

### ca_disco.py

**Input**: File containing list of multi-copy trees in newick format and set of alignment files in phylip format corresponding to the gene families.

**Output**: Concatenated alignment file in the phylip format

```
python3 ca_disco.py -i <input_trees> -a <alignments_list> -t <taxa_list> -o <output> -d <delimiter> -m <n> 
```

`disco.py` must be present in the same directory as `ca_disco.py` in order for it to run. Also, unlike `disco.py`, it is necessary for the input newick trees given to `ca_disco.py` to have unique leaf labels where the taxon name comes first and is separated from the rest of the name by some delimiter. 

#### Arguments

- **Required**
  - `-i`: Input newick tree file
  - `-a`: Text file containing paths to alignment files (one path for line, each path corresponding to gene-family tree on the same line in teh input tree file)
  - `-t`: Text file containing taxa list (one taxon per line)
  - `-o`: Output concatenated alignment file
- **Optional**
  - `-m`: Minimum number of taxa required for tree to be outputted. Default 4.
  - `-d`: Delimiter separating species name from rest of leaf label. Default _.

#### Example

```cmd
python3 ca_disco.py -i example/g_100.trees -o example.phy -a example/seq_list.txt -t example/taxa_list.txt
```
