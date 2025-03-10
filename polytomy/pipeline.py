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

    def place_additional_sequences(self, sequences):
        """
        Place additional sequences onto the optimized tree.

        Args:
            sequences (str): Path to FASTA file with sequences to place.

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

        self.logger.info(f"Placing additional sequences from {sequences}")

        # Configure the sequence placer
        placer_config = self.config.get('placer', {})

        # Run sequence placement
        placer = SequencePlacer(self.optimized_tree, config=placer_config)

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
