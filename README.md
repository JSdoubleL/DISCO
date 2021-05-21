# DISCO

Decomposition Into Single-COpy gene trees (DISCO) is a method for decomposing multi-copy gene-family trees while attempting to preserve orthologs and discard paralogs.

## Dependencies

- Python 3
- [TreeSwift](https://github.com/niemasd/TreeSwift)

## Usage

**Input**: File containing list of multi-copy trees in newick format

**Output**: File containing list resulting list of trees after decomposition in newick format

```
python3 tag_decomp.py -i <input_file> -o <ouput_file> -d <delimiter>
```

### Arguments

#### Required

- `-i`: Input newick tree file

#### Optional

- `-o`: Output newick tree file
- `-d`: Delimiter separating species name from rest of leaf label. Default None.
- `-m`: Output only single tree (discarding smallest duplicate clades).
- `-n`: No decomposition (outputs rooted gene trees).
- `-v`: Enable verbose output
- `-rp`: Remove in-paralogs before rooting/scoring (does not affect output, only reported score)
- `--trivial`: Includes trivial trees in decomposition output (by default trees not containing a quartet or rooted triple are discarded).
- `--outgroups`: Write outgroups (including ties) to txt file. (Might make program slower).

### Example

```cmd
python tag_decomp.py -i example/gtrees-mult.trees
```
