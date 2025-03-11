# Bactria backbone resolution

## Project Overview
- This project aims to tooling for resolving polytomies (multifurcating nodes) in BOLD taxonomic trees using the Open Tree of Life (OpenToL) API.
- The pipeline processes a Newick format tree, identifies polytomies, and attempts to resolve them through a combination of OpenToL-based resolution and weighted minimal-loss pruning.
- Once resolved, the tree undergoes branch length computation and can incorporate additional sequences through phylogenetic placement.
- The tool is designed to handle large-scale phylogenies with approximately 50,000 tips, requiring highly optimized processing strategies.
- This solution specifically replaces the problematic `prep_raxml_backbone`, `run_raxml_backbone`, and `reroot_backbone` steps in the Bactria pipeline when dealing with the full BOLD database (approximately 10M records).

## Core Components

### Tree Processing
- Parse Newick format trees using the DendroPy library
- Traverse tree structures efficiently in post-order (children before parents)
- Identify interior nodes with more than two children (polytomies)
- Record node information including OTT IDs using DendroPy annotations
- Implement memory-efficient algorithms to handle trees with ~50,000 tips

### OpenToL API Integration
- Resolve taxon names to OpenToL Taxonomy IDs (OTT IDs) using the TNRS API
- Request induced subtrees from the OpenToL API based on OTT IDs
- Extract branching order information from OpenToL subtrees
- Apply this branching information to resolve polytomies in the original tree

### Polytomy Resolution Algorithm

This is discussed in more detail in the [Polytomy Resolution Algorithm](resolution-pruning.md) document.

### Branch Length Optimization
- Compute branch lengths on the resolved topology using IQTree (preferred for large trees) or RAxML-NG
- Ensure all branch lengths are proportional and biologically plausible
- Investigate outliers, maybe contaminants (like Wolbachia); consider pruning

### Sequence Integration
- Place pruned tips and unplaced sequences onto the backbone tree
- Support for grafting family-level subtrees using the Bactria pipeline's `graft_clades` step
- Consider using just one of an exemplar pair in initial resolved backbone and placement
- Graft other member of pair later, to guarantee exemplar monophyly, then re-optimize
- Get exemplair pair list from MaaS C11G

## Implementation Details

### Dependencies
- DendroPy: for phylogenetic tree manipulation
- Requests: for API communication with OpenToL
- IQTree: primary tool for branch length computation and possible branch swapping on large trees
- RAxML-NG: alternative for branch length computation on smaller trees or subtrees
- RAxML (classic): for phylogenetic placement of sequences

### Architecture Design
We've defined a modular class structure with 6 core components:
1. `TreeParser`: For efficiently parsing large Newick trees with DendroPy
2. `PolytomyResolver`: For detecting and resolving polytomies using OpenToL and weighted pruning strategies
3. `OpenToLClient`: For interacting with OpenToL APIs with caching
4. `BranchLengthOptimizer`: For computing branch lengths using IQTree or RAxML-NG
5. `SequencePlacer`: For placing sequences onto backbone trees
6. `PolytomyResolutionPipeline`: For orchestrating the entire workflow

### Design Patterns
- Facade Pattern in the pipeline class for overall orchestration
- Strategy Pattern for branch length optimization selection
- Adapter Pattern for OpenToL API integration
- Cache Pattern for API response storage
- Factory Method for tree creation
- Composite Pattern for tree structure representation

### Performance Optimization
- Strategic caching of API responses to reduce redundant network calls
- Use of IQTree for faster branch length calculation on large phylogenies
- Efficient propagation of annotations through tree transformations

### Robustness Features
- Comprehensive logging throughout the codebase
- Detailed error handling and fallback mechanisms
- Performance optimizations for large trees
- Data caching to minimize API calls
- Configuration flexibility for all components
- Special handling for MRCA nodes and polyphyletic subtrees

### Testing Strategy
- Comprehensive test coverage using pytest
- Unit tests for each component:
  - `tree_parser.py`: Tests for parsing Newick trees from files and strings
  - `polytomy_resolver.py`: Tests for polytomy detection and resolution logic
  - `opentol_client.py`: Tests for taxon name resolution and API interactions- 
  - `branch_optimizer.py`: Tests for branch length optimization with IQTree/RAxML-NG  - 
  - `sequence_placer.py`: Tests for sequence placement algorithms
  - `pipeline.py`: Integration tests for the complete workflow
- Real-world test data from the example tree included in the project
- Tests designed to work even when external tools or APIs are unavailable

### Code Structure
- Following nbitk coding style with its configuration and logging system
- Function documentation using pydoc
- PEP8 compliance with snake_case naming conventions
- Command-line interface using argparse
- Modular design to minimize git diff impacts when modifying existing code

### Project Structure
We've established a clear project structure with:
- A main command-line entry point script (`resolve_polytomies.py`)
- Modular Python files for each component in the `polytomy` package
- Directory structure for tests and examples
- Documentation for usage and API

```
polytomy-resolution/
├── resolve_polytomies.py     # Main script
├── polytomy/                 # Core package
│   ├── __init__.py
│   ├── tree_parser.py        # Tree parsing functionality
│   ├── opentol_client.py     # OpenToL API client
│   ├── polytomy_resolver.py  # Polytomy detection and resolution
│   ├── branch_optimizer.py   # Branch length optimization
│   ├── sequence_placer.py    # Sequence placement
│   └── pipeline.py           # Pipeline orchestration
├── tests/                    # Test files
│   ├── data/                 # Test data
│   │   ├── example_tree.tre  # Example tree with polytomies
│   │   └── example_alignment.fa # Alignment for testing
│   ├── test_tree_parser.py   # Tests for tree parser
│   ├── test_opentol_client.py # Tests for OpenToL client
│   └── ...                   # Additional test modules
├── examples/                 # Example files and tutorials
├── requirements.txt          # Python dependencies
└── README.md                 # Project documentation
```

## Workflow

1. **Input**: Newick format tree (potentially with polytomies)
2. **Processing**:
   - Parse tree using DendroPy with optimized memory settings
   - Traverse tree structure post-order
   - For each polytomy:
     - Extract child node names
     - Resolve names to OTT IDs using OpenToL TNRS API (batched where possible)
     - Request induced subtree from OpenToL
     - Apply subtree topology to resolve polytomy
     - If OpenToL resolution fails, apply weighted pruning strategy
   - Annotations (tip counts, pruned lists, OTT IDs) are propagated throughout
3. **Optimization**:
   - Compute branch lengths using IQTree (preferred for 50k tip trees)
   - Consider limited branch swapping if needed using IQTree
4. **Integration**:
   - Place pruned tips and additional sequences onto the backbone
   - Graft family-level subtrees using Bactria pipeline's `graft_clades` function

## Challenges and Considerations
- API rate limiting and error handling for OpenToL requests
- Handling cases where OpenToL lacks information for certain taxa
- Managing computational efficiency and memory usage for trees with ~50,000 tips
- Minimizing tip loss when pruning is necessary
- Special handling for MRCA nodes and polyphyletic subtrees
- Ensuring compatibility with the existing Bactria pipeline
- Optimizing the transition between different computational tools (DendroPy, IQTree, RAxML)
- Scaling to handle the full BOLD database with approximately 10M records

## Future Extensions
- Parallelization of API requests and tree operations
- Caching of OpenToL responses to improve efficiency
- Improved handling of polyphyletic groups in OpenToL data
- Refinement of the weighted pruning strategy based on empirical results