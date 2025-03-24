#!/usr/bin/env python
"""
Polytomy Resolution Pipeline - Main orchestration module

This module provides a pipeline that coordinates the complete polytomy resolution workflow,
from tree parsing to polytomy resolution, branch length optimization, and sequence placement.
"""

import os
import logging
import time
import dendropy
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
            
            # Ensure newly placed exemplars are retained in the tree
            self.tree = self.parser.parse_from_file(filtered_tree)
            self.tree.retain_taxa_with_labels([leaf.taxon.label for leaf in self.final_tree.leaf_node_iter()])
            
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
        # Find the most advanced tree to optimize
        tree_to_optimize = None
        if hasattr(self, 'final_tree') and self.final_tree:
            tree_to_optimize = self.final_tree
            self.logger.info(f"Optimizing final tree with {len(tree_to_optimize.leaf_nodes())} tips")
        elif hasattr(self, 'optimized_tree') and self.optimized_tree:
            tree_to_optimize = self.optimized_tree
            self.logger.info(f"Optimizing optimized tree with {len(tree_to_optimize.leaf_nodes())} tips") 
        elif hasattr(self, 'resolved_tree'):
            tree_to_optimize = self.resolved_tree
            self.logger.info(f"Optimizing resolved tree with {len(tree_to_optimize.leaf_nodes())} tips")
        else:
            self.logger.error("No tree available to optimize")
            return False

        # Configure the branch length optimizer
        optimizer_config = self.config.get('optimizer', {})
        if alignment:
            optimizer_config['alignment'] = alignment

        # Check if alignment is available
        if 'alignment' not in optimizer_config:
            self.logger.warning("No alignment provided for branch length optimization")
            self.optimized_tree = tree_to_optimize
            return True
            
        # Run branch length optimization
        optimizer = BranchLengthOptimizer(
            tree_to_optimize,
            tool=optimizer_config.get('tool', 'iqtree'),
            config=optimizer_config
        )

        self.optimized_tree = optimizer.optimize_branch_lengths()

        if not self.optimized_tree:
            self.logger.error("Branch length optimization failed")
            # Fall back to original tree without optimization
            self.optimized_tree = tree_to_optimize
            return False
            
        # If we were optimizing the final tree, update it with the optimized version
        if hasattr(self, 'final_tree') and self.final_tree is tree_to_optimize:
            self.final_tree = self.optimized_tree

        self.logger.info("Tree optimization completed successfully")
        return True

    def place_additional_sequences(self, sequences, pairs_file=None):
        """
        Place additional sequences onto the optimized tree.
        
        Args:
            sequences (str): Path to FASTA file with sequences to place.
            pairs_file (str, optional): Path to exemplar pairs file in format:
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
            
            # Store the alignment created for placement - this includes the first exemplars
            self.placement_alignment = placer.last_placement_alignment
        else:
            # Standard placement behavior (no pairs file)
            self.logger.info(f"Placing all additional sequences from {sequences}")
            
            # If sequences is a list, run batch placement
            if isinstance(sequences, list):
                self.final_tree = placer.batch_place_sequences(sequences)
            else:
                self.final_tree = placer.place_sequences(sequences)
                
            # Store the alignment used for placement
            self.placement_alignment = placer.last_placement_alignment

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
        # First debug tree references
        self.logger.info(f"DEBUG tree references before selecting grafting tree:")
        for attr in ['final_tree', 'optimized_tree', 'resolved_tree']:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                tree_obj = getattr(self, attr)
                self.logger.info(f"  {attr}: {len(tree_obj.leaf_nodes())} tips")
        
        # Find the tree with the most tips
        tree_sizes = []
        if hasattr(self, 'final_tree') and self.final_tree:
            final_tips = len(self.final_tree.leaf_nodes())
            tree_sizes.append((self.final_tree, final_tips, "final"))
            self.logger.info(f"Final tree has {final_tips} tips")

        if hasattr(self, 'optimized_tree') and self.optimized_tree:
            optimized_tips = len(self.optimized_tree.leaf_nodes())
            tree_sizes.append((self.optimized_tree, optimized_tips, "optimized"))
            self.logger.info(f"Optimized tree has {optimized_tips} tips")
            
        if hasattr(self, 'resolved_tree') and self.resolved_tree:
            resolved_tips = len(self.resolved_tree.leaf_nodes())
            tree_sizes.append((self.resolved_tree, resolved_tips, "resolved"))
            self.logger.info(f"Resolved tree has {resolved_tips} tips")
            
        if not tree_sizes:
            self.logger.error("No trees available for grafting second exemplars")
            return False
            
        # Sort by number of tips (descending)
        tree_sizes.sort(key=lambda x: x[1], reverse=True)
        
        # Use the tree with the most tips
        tree_to_use, tip_count, tree_name = tree_sizes[0]
        self.logger.info(f"Using {tree_name} tree with {tip_count} tips for grafting (has most tips)")
        
        # Configure the sequence placer
        placer_config = self.config.get('placer', {})
        placer = SequencePlacer(tree_to_use, config=placer_config)
        
        # Graft second exemplars
        self.logger.info(f"Grafting second exemplars from {pairs_file}")
        self.final_tree = placer.graft_second_exemplars(alignment, pairs_file)
        
        # Store the alignment used for grafting - important for subsequent optimization 
        self.placement_alignment = placer.last_placement_alignment
        self.logger.info(f"Updated placement alignment path to: {self.placement_alignment}")
        
        if not self.final_tree:
            self.logger.error("Second exemplar grafting failed")
            # Fall back to previous tree
            self.final_tree = tree_to_use
            return False
            
        self.logger.info("Second exemplars grafted successfully")
        return True
        
    def optimize_tree_after_placement(self, alignment=None):
        """
        Optimize branch lengths on the tree after sequence placement.
        
        Args:
            alignment (str, optional): Path to alignment file. If None, uses the 
                                    placement alignment stored during placement/grafting.
        
        Returns:
            bool: True if optimization succeeded, False otherwise.
        """
        if not self.final_tree:
            self.logger.error("No tree available for optimization after placement")
            return False
            
        # Use placement alignment if available and no specific alignment provided
        if alignment is None and hasattr(self, 'placement_alignment') and self.placement_alignment:
            alignment = self.placement_alignment
            self.logger.info(f"Using placement alignment for optimization: {alignment}")
        
        if not alignment:
            self.logger.warning("No alignment available for branch length optimization after placement")
            return False
        
        # Configure the optimizer
        optimizer_config = self.config.get('optimizer', {})
        optimizer = BranchLengthOptimizer(self.final_tree, 
                                        tool=self.config.get('optimization_tool', 'iqtree'),
                                        config=optimizer_config)
        
        # Optimize the tree
        self.logger.info(f"Optimizing final tree with {len(self.final_tree.leaf_nodes())} tips")
        optimized_tree = optimizer.optimize_branch_lengths(alignment)
        
        if optimized_tree:
            self.final_tree = optimized_tree
            self.logger.info("Tree optimization completed successfully")
            return True
        else:
            self.logger.error("Tree optimization failed")
            return False

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
