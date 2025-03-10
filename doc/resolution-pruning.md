# Polytomy Resolution Algorithm

## Overview
This document describes a comprehensive algorithm for resolving polytomies (nodes with more than two children) in 
phylogenetic trees. The algorithm uses a combination of the Open Tree of Life (OpenToL) API for taxonomically informed 
resolution and a weighted pruning strategy for cases where OpenToL data is insufficient. The approach is designed to 
minimize tip loss and preserve as much taxonomic structure as possible.

## Core Concepts

### Tree Representation
- Trees are represented using DendroPy's tree structures
- Annotations are stored directly on nodes using DendroPy's annotation system
- Post-order traversal ensures child nodes are processed before their parents

### Types of Resolution
1. **OpenToL-based Resolution**: Using the OpenToL API to obtain taxonomically informed topologies
2. **Weighted Pruning**: When OpenToL resolution fails or is partial, using tip counts to determine which nodes to retain

## Algorithm Details

### Step 1: Tree Annotation
The first phase involves annotating the tree with essential information:

1. Traverse the tree in post-order (children before parents)

2. For each node:
   - Use DendroPy annotations to store additional information, whose keys are:
     - `size`: Number of tips in the subtree rooted at this node
     - `pruned`: List of pruned tip labels (BOLD process IDs) for this node
     - `ott_id`: OTT ID of the taxon represented by this node (if known)
   - If leaf node: set `size = 1`, `pruned = []`
   - If internal node: set `size = sum(child.annotations['size'].value for all children)`
   - If internal node: make pruned the union of all children's pruned lists
   - If internal node is polytomous, apply resolution, i.e. Step 3 during this traversal

This annotation process efficiently computes tip counts for all nodes and prepares the tree for resolution.

### Step 2: Polytomy Identification
1. Find all nodes with more than two children (polytomies)
2. Process these polytomies in the natural post-order sequence (already guaranteed by the traversal)

### Step 3: Polytomy Resolution
For each polytomy encountered during traversal:

1. **Extract Taxonomic Information**
   - Get taxon names from all child nodes
   - Resolve names to OTT IDs using OpenToL TNRS API

2. **OpenToL Resolution Attempt**
   - Fetch induced subtree from OpenToL API using OTT IDs
   - If successful, process the returned subtree:
     - Special handling for MRCA nodes (see Step 4)
     - Graft OpenToL structure onto the original polytomy
     - Remove extraneous tips added by OpenToL
     - Propagate annotations through the new structure using post-order traversal
     - If the resulting structure still contains polytomies, apply weighted pruning

3. **Weighted Pruning (when OpenToL resolution fails or is partial)**
   - Sort child nodes by tip count (ascending)
   - Keep the two largest children (those with most tips)
   - Prune the smaller children
   - Track pruned tips by their process IDs (leaf labels)

4. **Resolution Tracking**
   - Update node annotations to track resolution progress
   - Counts of nodes resolved by OpenToL is in `size` and by pruning in `len(pruned)`
   - We also keep list of pruned tips for later placement

### Step 4: MRCA Node Handling
MRCA nodes from OpenToL (e.g., "mrcaott170987ott201497") require special handling:

1. These represent polyphyletic subtrees and are "broken" in the returned JSON
2. First, graft the existing subtrees at the appropriate point
3. Remove extraneous tips that OpenToL added to describe the polyphyly
4. If still polytomous, apply weighted pruning (see 3.3)
5. Propagate annotations through the new structure

**Note**: This is so tricky due to the different polyphyly patterns that we skip this for now.

### Step 5: Post-Resolution Processing
After all polytomies are resolved:

1. **Branch Length Computation**
   - Use BranchOptimizer to compute branch lengths on the resolved tree
   - This is done after topology is fully resolved

2. **Sequence Placement**
   - Place pruned sequences back onto the backbone
   - This can be done using evolutionary placement algorithms

3. **Family-Level Subtree Grafting**
   - Graft family-level subtrees onto the backbone
   - Use existing Bactria pipeline functionality

This algorithm provides a comprehensive approach to resolving polytomies in phylogenetic trees, leveraging both the 
taxonomic knowledge from OpenToL and information-preserving weighted pruning strategies.