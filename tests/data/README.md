# Data files

The `tests/data` directory contains the following files:

1. `README.md`: A markdown file with project documentation
2. `example_alignment.fa`: A FASTA file with sequence data
3. `example_tree.tre`: A Newick tree file with polytomies, some tips occur in 2, is a clade from 8
4. `output_in_tree.fasta`: A FASTA file with sequence data, matches tips in 3
5. `output_not_in_tree.fasta`: A FASTA file with sequence data, does not match tips in 3
6. `example_tree_resolved.tre`: A resolved Newick tree file, based on 3 but sparser
7. `example_tree_iqtree.tre`: A Newick tree file with branch lengths optimized by IQTree, based on 4 and 6
8. `outfile-labels-internals.tre`: A Newick tree file, based on taxonomy from backbone alignment
