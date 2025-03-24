# Polytomy Resolution Tool

A Python tool for resolving polytomies in large phylogenetic trees using the Open Tree of Life (OpenToL) API.

## Overview

This tool processes phylogenetic trees with polytomies (multifurcating nodes) and resolves them into fully bifurcating trees. It offers the exemplar pair handling capabilities, which allow placing both first and second exemplars in a controlled manner. The pipeline combines:

1. OpenToL-based resolution to apply evolutionary knowledge
2. Minimal-loss pruning for nodes that cannot be resolved through OpenToL
3. Branch length optimization using IQTree or RAxML-NG
4. Sequence placement with EPA-RAxML
5. Exemplar pair placement strategy

Designed specifically for handling large phylogenies (~50,000 tips) efficiently.

## Key Features

- **Memory-efficient tree processing** for large phylogenies
- **OpenToL integration** to resolve polytomies based on established phylogenetic knowledge
- **Minimal-loss pruning** strategy for remaining polytomies
- **Branch length optimization** using IQTree or RAxML-NG
- **Sequence placement** for pruned tips and additional sequences
- **Exemplar pair handling** with two-phase placement:
  - First place one exemplar from each taxon
  - Then graft the second exemplar next to the first
- **Alignment compression** to speed up placement with most informative sites
- **Parallelization support** for multi-threaded execution
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

Complete Exemplar Pair usage:

```bash
python resolve_polytomies.py \
  --input tree.newick \
  --output final_tree_optimized.tre \
  --alignment alignment.fa \
  --sequences alignment.fa \
  --exemplar-pairs exemplars_table.txt \
  --place-first-exemplars \
  --graft-second-exemplars \
  --optimization-tool raxml-ng \
  --model "GTR+G" \
  --threads 52 \
  --compress-alignment \
  --compress-columns 700 \
  --log-file full_pipeline.log \
  --keep-files
```


# Key Command Line Arguments

- Input/Output  
--input: Input tree file in Newick format  
--output: Output resolved tree file  
--alignment: Alignment file for branch length optimization  
--sequences: FASTA file with sequences to place onto the backbone  
--exemplar-pairs: Path to file with exemplar pair information (format: taxon<tab>exemplar1<tab>exemplar2)  
  
- Workflow Options  
--skip-opentol: Skip OpenToL resolution for trees that are already resolved  
--filter-exemplars: Filter tree to keep at most one exemplar per pair before optimization  
--place-first-exemplars: Place only the first exemplar for each missing taxon   
--graft-second-exemplars: Graft the second exemplar of each pair next to the first one  
  
- Optimization Options  
--optimization-tool: Tool to use (iqtree, raxml-ng)  
--model: Sequence evolution model ("GTR+G" for IQTree, "GTRCAT" for RAxML)  
--threads: Number of threads for parallel processing  
--max-memory: Maximum memory usage in MB  
  
- Performance Options  
--compress-alignment: Compress alignment to most informative columns  
--compress-columns: Number of columns to keep when compressing (default: 700)  
  
- Utility Options  
--log-level: Set logging level (debug, info, warning, error, critical)  
--log-file: Path to output log file  
--keep-files: Keep temporary files generated during sequence placement  
  
  
## Dependencies

- Python 3.8+
- DendroPy
- Requests
- IQTree and/or RAxML-NG (for branch length computation)
- RAxML (for phylogenetic placement)

## Project Structure

```
barcode-constrained-phylogeny-backbone/
├── resolve_polytomies.py          # Main script
├── polytomy/                      # Core package
│   ├── __init__.py
│   ├── tree_parser.py             # Tree parsing functionality
│   ├── opentol_client.py          # OpenToL API client
│   ├── polytomy_resolver.py       # Polytomy resolution logic
│   ├── branch_optimizer.py        # Branch length optimization
│   ├── sequence_placer.py         # Sequence placement with exemplar support
│   ├── tree_alignment_matcher.py  # Tree-alignment compatibility
│   └── pipeline.py                # Pipeline orchestration
├── tests/                         # Test files
│   ├── data/                      # Test data
│   └── ...                        # Test modules
├── environment.yml                # Conda dependencies
├── LICENSE                        # Apache 2.0 license
├── pyproject.toml                 # Python package metadata
├── requirements.txt               # Python dependencies
└── README.md                      # This file
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
