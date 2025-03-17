#!/usr/bin/env python
"""
Polytomy Resolution Pipeline - Main orchestration module

This module provides a pipeline that coordinates the complete polytomy resolution workflow,
from tree parsing to polytomy resolution, branch length optimization, and sequence placement.
"""

import os
import logging
import time
from pathlib import Path

# Import pipeline components
from polytomy.tree_parser import TreeParser
from polytomy.opentol_client import OpenToLClient
from polytomy.polytomy_resolver import PolytomyResolver
from polytomy.branch_optimizer import BranchLengthOptimizer
from polytomy.sequence_placer import SequencePlacer
from polytomy.tree_alignment_matcher import TreeAlignmentMatcher


class PolytomyResolutionPipeline:
    """Orchestrates the complete polytomy resolution workflow."""

    def __init__(self, config=None):
        """
        Initialize with optional configuration.

        Args:
            config (dict, optional): Configuration options for the pipeline.
        """
        self.config = config or {}
        self.logger = logging.getLogger(__name__)

        # Initialize pipeline state
        self.tree = None
        self.original_tree = None
        self.resolved_tree = None
        self.optimized_tree = None
        self.final_tree = None

        # Initialize pipeline components
        self.parser = TreeParser(config=self.config.get('parser', {}))

        # Initialize OpenToL client if not disabled
        if not self.config.get('opentol', {}).get('skip_opentol', False):
            self.opentol_client = OpenToLClient(config=self.config.get('opentol', {}))
        else:
            self.opentol_client = None
            self.logger.warning("OpenToL integration disabled")

        # Track pipeline execution stats
        self.stats = {
            'start_time': None,
            'end_time': None,
            'elapsed_time': None,
            'tree_size': None,
            'polytomy_count': None,
            'resolved_count': None,
            'pruned_tips': None,
        }

        self.logger.info("Polytomy resolution pipeline initialized")

    def load_tree(self, tree_source):
        """
        Load a tree from file or string.

        Args:
            tree_source (str): Path to a tree file or a Newick string.

        Returns:
            bool: True if tree was loaded successfully, False otherwise.
        """
        self.logger.info(f"Loading tree from {tree_source}")

        try:
            # Check if tree_source is a file path
            if os.path.exists(tree_source):
                self.tree = self.parser.parse_from_file(tree_source)
            else:
                # Assume it's a Newick string
                self.tree = self.parser.parse_from_string(tree_source)

            # Make a deep copy of the original tree for reference
            import dendropy
            self.original_tree = dendropy.Tree(self.tree)

            # Record tree size
            self.stats['tree_size'] = len(self.tree.leaf_nodes())
            self.logger.info(f"Tree loaded with {self.stats['tree_size']} tips")

            return True

        except Exception as e:
            self.logger.error(f"Failed to load tree: {str(e)}")
            return False

    def resolve_polytomies(self):
        """
        Execute the complete polytomy resolution process.

        Returns:
            bool: True if polytomies were resolved successfully, False otherwise.
        """
        if not self.tree:
            self.logger.error("No tree loaded. Call load_tree() first.")
            return False

        self.logger.info("Starting polytomy resolution")
        self.stats['start_time'] = time.time()

        # Resolve polytomies
        resolver = PolytomyResolver(self.tree, self.opentol_client)
        resolver.resolve_all_polytomies()

        # Store the resolved tree
        self.resolved_tree = self.tree

        # Calculate elapsed time
        self.stats['end_time'] = time.time()
        self.stats['elapsed_time'] = self.stats['end_time'] - self.stats['start_time']

        self.logger.info(f"Polytomy resolution completed in {self.stats['elapsed_time']:.2f} seconds")

        return True

    def match_tree_and_alignment(self, tree_path, alignment_path, output_dir=None):
        """
        Match tree and alignment files to ensure they have the same taxa.
        
        Args:
            tree_path: Path to the tree file
            alignment_path: Path to the alignment file
            output_dir: Directory for output files (default: self.output_dir)
        
        Returns:
            tuple: (filtered_tree_path, filtered_alignment_path) or (None, None) if failed
        """
        self.logger.info(f"Matching tree and alignment taxa")
        
        # Initialize matcher
        matcher = TreeAlignmentMatcher(config={'output_dir': output_dir or '.'})
        
        # Match tree and alignment
        filtered_tree, filtered_alignment = matcher.match_tree_and_alignment(
            tree_path, 
            alignment_path
        )
        
        if filtered_tree and filtered_alignment:
            self.logger.info(f"Successfully matched tree and alignment")
            
            # Load the filtered tree
            self.tree = self.parser.parse_from_file(filtered_tree)
            
            return filtered_tree, filtered_alignment
        else:
            self.logger.error(f"Failed to match tree and alignment")
            return None, None

    def filter_exemplar_pairs(self, pairs_file):
        """Filter tree to keep at most one exemplar per pair"""
        self.logger.info(f"Filtering tree to keep at most one exemplar per pair")
        
        # Read pairs file
        pairs = {}
        with open(pairs_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    taxon = parts[0]
                    exemplars = [parts[1], parts[2] if len(parts) > 2 else None]
                    pairs[taxon] = [ex for ex in exemplars if ex]
        
        # Get current tips
        current_tips = {leaf.taxon.label for leaf in self.tree.leaf_node_iter()}
        
        # Find tips to prune
        to_prune = []
        for taxon, exemplars in pairs.items():
            present = [ex for ex in exemplars if ex in current_tips]
            if len(present) > 1:  # If both exemplars present, keep only first
                to_prune.extend(present[1:])
        
        # Prune identified tips
        if to_prune:
            self.logger.info(f"Removing {len(to_prune)} tips to ensure at most one exemplar per pair")
            self.tree.prune_taxa_with_labels(to_prune)
        
        return len(to_prune)

    def optimize_tree(self, alignment=None):
        """
        Optimize the tree structure and branch lengths.

        Args:
            alignment (str, optional): Path to alignment file for branch length optimization.
                                     If not specified, uses alignment from config.

        Returns:
            bool: True if tree was optimized successfully, False otherwise.
        """
        if not self.resolved_tree:
            self.logger.error("No resolved tree available. Call resolve_polytomies() first.")
            return False

        self.logger.info("Starting tree optimization")

        # Configure the branch length optimizer
        optimizer_config = self.config.get('optimizer', {})
        if alignment:
            optimizer_config['alignment'] = alignment

        # Check if alignment is available
        if 'alignment' not in optimizer_config:
            self.logger.warning("No alignment provided for branch length optimization")
            self.optimized_tree = self.resolved_tree
            return True

        # Run branch length optimization
        optimizer = BranchLengthOptimizer(
            self.resolved_tree,
            tool=optimizer_config.get('tool', 'iqtree'),
            config=optimizer_config
        )

        self.optimized_tree = optimizer.optimize_branch_lengths()

        if not self.optimized_tree:
            self.logger.error("Branch length optimization failed")
            # Fall back to resolved tree without optimized branch lengths
            self.optimized_tree = self.resolved_tree
            return False

        self.logger.info("Tree optimization completed successfully")
        return True

    def place_additional_sequences(self, sequences, pairs_file=None):
        """
        Place additional sequences onto the optimized tree.
        
        If a pairs file is provided, only the first exemplar from each pair
        that isn't already represented in the tree will be placed.
        
        Args:
            sequences (str): Path to FASTA file with sequences to place.
            pairs_file (str, optional): Path to file with exemplar pair information.
                Format: taxon_name<tab>exemplar1<tab>exemplar2
                
        Returns:
            bool: True if sequences were placed successfully, False otherwise.
        """
        if not self.optimized_tree:
            self.logger.error("No optimized tree available. Call optimize_tree() first.")
            return False

        # Skip if sequences is None or empty or placement is disabled
        if not sequences or self.config.get('placer', {}).get('skip_placement', False):
            self.logger.info("Skipping sequence placement (no sequences or disabled)")
            self.final_tree = self.optimized_tree
            return True

        # Configure the sequence placer
        placer_config = self.config.get('placer', {})
        placer = SequencePlacer(self.optimized_tree, config=placer_config)
        
        # Handle exemplar pairs if provided
        if pairs_file:
            self.logger.info(f"Placing first exemplars from pairs in {pairs_file}")
            self.final_tree = placer.place_first_exemplars(sequences, pairs_file)
        else:
            # Standard placement behavior (no pairs file)
            self.logger.info(f"Placing all additional sequences from {sequences}")
            
            # If sequences is a list, run batch placement
            if isinstance(sequences, list):
                self.final_tree = placer.batch_place_sequences(sequences)
            else:
                self.final_tree = placer.place_sequences(sequences)

        if not self.final_tree:
            self.logger.error("Sequence placement failed")
            # Fall back to optimized tree without additional sequences
            self.final_tree = self.optimized_tree
            return False

        # Count how many sequences were placed
        placed_count = len(self.final_tree.leaf_nodes()) - len(self.optimized_tree.leaf_nodes())
        self.logger.info(f"Placed {placed_count} additional sequences")

        return True

    def graft_second_exemplars(self, alignment, pairs_file):
        """
        Graft second exemplars from each pair next to the first one.
        
        Args:
            alignment (str): Path to alignment with all exemplar sequences
            pairs_file (str): Path to exemplar pair information
        
        Returns:
            bool: True if grafting succeeded, False otherwise
        """
        if not self.optimized_tree and not self.final_tree:
            self.logger.error("No tree available for grafting. Run resolve_polytomies() first.")
            return False
        
        # Use the most advanced tree available
        tree_to_use = self.final_tree if self.final_tree else self.optimized_tree
        
        self.logger.info(f"Grafting second exemplars from {pairs_file}")
        
        # Initialize placer with the current tree
        placer = SequencePlacer(tree_to_use, config=self.config.get('placer', {}))
        
        # Graft second exemplars
        grafted_tree = placer.graft_second_exemplars(alignment, pairs_file)
        
        if not grafted_tree:
            self.logger.error("Failed to graft second exemplars")
            return False
        
        # Update the final tree
        self.final_tree = grafted_tree
        
        self.logger.info("Second exemplars grafted successfully")
        return True

    def write_tree(self, output_path, format='newick'):
        """
        Write the final tree to file.

        Args:
            output_path (str): Path to output file.
            format (str, optional): Output format ('newick' or 'nexus').

        Returns:
            bool: True if tree was written successfully, False otherwise.
        """
        if not self.final_tree and not self.optimized_tree and not self.resolved_tree:
            self.logger.error("No tree available to write")
            return False

        # Use the most advanced tree available
        if self.final_tree:
            tree_to_write = self.final_tree
        elif self.optimized_tree:
            tree_to_write = self.optimized_tree
        else:
            tree_to_write = self.resolved_tree

        # Ensure output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Write tree
        try:
            tree_to_write.write(
                path=output_path,
                schema=format.lower()
            )
            self.logger.info(f"Tree written to {output_path}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to write tree: {str(e)}")
            return False
