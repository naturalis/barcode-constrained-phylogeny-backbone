# Polytomy Resolution Tool

A Python tool for resolving polytomies in large phylogenetic trees using the Open Tree of Life (OpenToL) API.

## Overview

This tool processes phylogenetic trees with polytomies (multifurcating nodes) and resolves them into fully bifurcating trees. It uses a combination of:

1. OpenToL-based resolution to apply evolutionary knowledge
2. Minimal-loss pruning for nodes that cannot be resolved through OpenToL
3. Branch length optimization using IQTree or RAxML-NG
4. Sequence placement for pruned tips and additional sequences

Designed specifically for handling large phylogenies (~50,000 tips) efficiently.

## Key Features

- **Memory-efficient tree processing** for large phylogenies
- **OpenToL integration** to resolve polytomies based on established phylogenetic knowledge
- **Minimal-loss pruning** strategy for remaining polytomies
- **Branch length optimization** using IQTree (preferred for large trees) or RAxML-NG
- **Sequence placement** for pruned tips and additional sequences
- **Subtree grafting** compatible with the Bactria pipeline's `graft_clades` step

## Installation

```bash
# Clone the repository
git clone https://github.com/naturalis/barcode-constrained-phylogeny-backbone.git
cd barcode-constrained-phylogeny-backbone

# Install dependencies
conda env create -f environment.yml
```

## Usage

Basic usage:

```bash
python resolve_polytomies.py --input tree.newick --output resolved_tree.newick
```

Advanced usage:

```bash
python resolve_polytomies.py \
  --input tree.newick \
  --output resolved_tree.newick \
  --optimization-tool iqtree \
  --sequences additional_seqs.fasta \
  --alignment alignment.fasta \
  --family-subtrees family_dir/ \
  --log-level INFO
```

## Dependencies

- Python 3.8+
- DendroPy
- Requests
- IQTree and/or RAxML-NG (for branch length computation)
- RAxML (for phylogenetic placement)

## Project Structure

```
barcode-constrained-phylogeny-backbone/
├── resolve_polytomies.py     # Main script
├── polytomy/                 # Core package
│   ├── __init__.py
│   ├── tree_parser.py        # Tree parsing functionality
│   ├── polytomy_finder.py    # Polytomy identification
│   ├── opentol_client.py     # OpenToL API client
│   ├── polytomy_resolver.py  # Polytomy resolution logic
│   ├── branch_optimizer.py   # Branch length optimization
│   ├── sequence_placer.py    # Sequence placement
│   └── pipeline.py           # Pipeline orchestration
├── tests/                    # Test files
│   ├── data/                 # Test data
│   └── ...                   # Test modules
├── examples/                 # Example files and tutorials
├── requirements.txt          # Python dependencies
└── README.md                 # This file
```

## Testing

Run the test suite:

```bash
pytest
```

## License

[Apache License](LICENSE)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
