I've reviewed the additional document on the polytomy resolution algorithm and will update the project documentation to include these details, particularly around the weighted pruning approach and the full algorithm steps.

Here's the updated project documentation:

# Polytomy Resolution Project Documentation

## Project Overview
- This project aims to develop a Python tool for resolving polytomies (multifurcating nodes) in phylogenetic trees using the Open Tree of Life (OpenToL) API.
- The pipeline processes a Newick format tree, identifies polytomies, and attempts to resolve them through a combination of OpenToL-based resolution and weighted minimal-loss pruning.
- Once resolved, the tree undergoes branch length computation and can incorporate additional sequences through phylogenetic placement.
- The tool is designed to handle large-scale phylogenies with approximately 50,000 tips, requiring highly optimized processing strategies.
- This solution specifically replaces the problematic `prep_raxml_backbone`, `run_raxml_backbone`, and `reroot_backbone` steps in the Bactria pipeline when dealing with the full BOLD database (approximately 10M records).

## Core Components

### Tree Processing
- Parse Newick format trees using the DendroPy library
- Traverse tree structures efficiently in post-order (children before parents)
- Identify interior nodes with more than two children (polytomies)
- Record node information including taxon names and relationships
- Implement memory-efficient algorithms to handle trees with ~50,000 tips

### OpenToL API Integration
- Resolve taxon names to OpenToL Taxonomy IDs (OTT IDs) using the TNRS API
- Request induced subtrees from the OpenToL API based on OTT IDs
- Extract branching order information from OpenToL subtrees
- Apply this branching information to resolve polytomies in the original tree
- Batch API requests where possible to minimize network overhead

### Polytomy Resolution Algorithm
1. **Tree Annotation**:
   - Annotate tree with tip counts, pruned tips lists, and OTT IDs
   - Calculate tip counts in post-order traversal (children before parents)

2. **OpenToL-based Resolution**:
   - Extract taxon names from polytomy child nodes
   - Resolve names to OTT IDs using OpenToL TNRS API
   - Fetch induced subtree from OpenToL based on these OTT IDs
   - Apply OpenToL topology to resolve the polytomy
   - Special handling for MRCA nodes and polyphyletic subtrees

3. **Weighted Pruning Strategy**:
   - When OpenToL resolution fails or is partial:
     - Sort child nodes by tip count (ascending)
     - Keep the two largest children (those with most tips)
     - Prune the smaller children
     - Track pruned tips by their process IDs (leaf labels)
   - This ensures minimal information loss while achieving bifurcation

4. **Post-Resolution Processing**:
   - Compute branch lengths using BranchOptimizer
   - Place pruned sequences back onto the backbone
   - Graft family-level subtrees

### Branch Length Optimization
- Compute branch lengths on the resolved topology using IQTree (preferred for large trees) or RAxML-NG
- Ensure all branch lengths are proportional and biologically plausible

### Sequence Integration
- Place pruned tips and unplaced sequences onto the backbone tree
- Support for grafting family-level subtrees using the Bactria pipeline's `graft_clades` step

## Implementation Details

### Dependencies
- DendroPy: for phylogenetic tree manipulation
- Requests: for API communication with OpenToL
- IQTree: primary tool for branch length computation and possible branch swapping on large trees
- RAxML-NG: alternative for branch length computation on smaller trees or subtrees
- RAxML (classic): for phylogenetic placement of sequences
- Bactria pipeline: for the final grafting of family-level subtrees

### Architecture Design
We've defined a modular class structure with 7 core components:
- `TreeParser`: For efficiently parsing large Newick trees with DendroPy
- `PolytomyFinder`: For identifying polytomies using tree traversal
- `OpenToLClient`: For interacting with OpenToL APIs with caching
- `PolytomyResolver`: For resolving polytomies using OpenToL and weighted pruning strategies
- `BranchLengthOptimizer`: For computing branch lengths using IQTree or RAxML-NG
- `SequencePlacer`: For placing sequences onto backbone trees
- `PolytomyResolutionPipeline`: For orchestrating the entire workflow

### Design Patterns
We've incorporated several design patterns:
- Facade Pattern in the pipeline class for overall orchestration
- Strategy Pattern for branch length optimization selection
- Adapter Pattern for OpenToL API integration
- Cache Pattern for API response storage
- Factory Method for tree creation
- Composite Pattern for tree structure representation

### Performance Optimization
- Chunked tree processing to manage memory consumption
- Strategic caching of API responses to reduce redundant network calls
- Parallel processing of independent subtrees where applicable
- Memory profiling and optimization for large tree structures
- Use of IQTree for faster branch length calculation on large phylogenies
- Efficient propagation of annotations through tree transformations

### Robustness Features
We've added numerous features for production robustness:
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
  - `polytomy_finder.py`: Tests for identifying polytomies in tree structures
  - `opentol_client.py`: Tests for taxon name resolution and API interactions
  - `branch_optimizer.py`: Tests for branch length optimization with IQTree/RAxML-NG
  - `polytomy_resolver.py`: Tests for polytomy resolution logic
  - `sequence_placer.py`: Tests for sequence placement algorithms
  - `pipeline.py`: Integration tests for the complete workflow
- Real-world test data from the example tree included in the project
- Tests designed to work even when external tools or APIs are unavailable

### Key Algorithm Steps
1. **Tree Annotation**:
   - Use DendroPy annotations to store tip counts, pruned lists, and OTT IDs
   - Calculate tip counts in post-order (children before parents)

2. **Polytomy Identification**:
   - Find all nodes with more than two children during tree traversal
   - Process in post-order sequence (children before parents)

3. **OpenToL Resolution**:
   - Extract taxon names from child nodes
   - Resolve names to OTT IDs via OpenToL TNRS API
   - Fetch induced subtree from OpenToL API
   - Apply OpenToL topology to resolve polytomy
   - Special handling for MRCA nodes (polyphyletic subtrees)

4. **Weighted Pruning**:
   - Sort child nodes by tip count (ascending)
   - Keep the two largest children (most tips)
   - Prune smaller children to achieve bifurcation
   - Track pruned tips for later placement

5. **Post-Resolution Processing**:
   - Compute branch lengths
   - Place pruned sequences back
   - Graft family-level subtrees

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
│   ├── polytomy_finder.py    # Polytomy identification
│   ├── opentol_client.py     # OpenToL API client
│   ├── polytomy_resolver.py  # Polytomy resolution logic
│   ├── branch_optimizer.py   # Branch length optimization
│   ├── sequence_placer.py    # Sequence placement
│   └── pipeline.py           # Pipeline orchestration
├── tests/                    # Test files
│   ├── data/                 # Test data
│   │   ├── example_tree.tre  # Example tree with polytomies
│   │   └── example_alignment.fa # Alignment for testing
│   ├── test_tree_parser.py   # Tests for tree parser
│   ├── test_polytomy_finder.py # Tests for polytomy finder
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