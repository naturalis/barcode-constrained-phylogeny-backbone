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
from polytomy.polytomy_finder import PolytomyFinder
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

        # Find polytomies
        finder = PolytomyFinder(self.tree)
        polytomies = finder.find_all_polytomies()
        self.stats['polytomy_count'] = len(polytomies)

        if self.stats['polytomy_count'] == 0:
            self.logger.info("No polytomies found in tree")
            self.resolved_tree = self.tree
            return True

        self.logger.info(f"Found {self.stats['polytomy_count']} polytomies")

        # Resolve polytomies
        resolver = PolytomyResolver(self.tree, self.opentol_client)
        resolved_count, failed_count, pruned_tips = resolver.resolve_all_polytomies()

        self.stats['resolved_count'] = resolved_count
        self.stats['pruned_tips'] = pruned_tips

        # Check if resolution was successful
        if failed_count > 0:
            self.logger.warning(f"Failed to resolve {failed_count} polytomies")

        if resolved_count == 0:
            self.logger.error("No polytomies were resolved")
            return False

        # Store the resolved tree
        self.resolved_tree = self.tree

        # Calculate elapsed time
        self.stats['end_time'] = time.time()
        self.stats['elapsed_time'] = self.stats['end_time'] - self.stats['start_time']

        self.logger.info(f"Polytomy resolution completed in {self.stats['elapsed_time']:.2f} seconds")
        self.logger.info(f"Resolved {resolved_count}/{self.stats['polytomy_count']} polytomies")
        self.logger.info(f"Pruned {len(pruned_tips)} tips")

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

    def graft_family_subtrees(self, subtree_files):
        """
        Graft family-level subtrees using the Bactria pipeline.

        Args:
            subtree_files (str or list): Path to subtree file or directory, or list of files.

        Returns:
            bool: True if subtrees were grafted successfully, False otherwise.
        """
        if not self.final_tree and not self.optimized_tree:
            self.logger.error("No backbone tree available. Resolve polytomies first.")
            return False

        # Use final_tree if available, otherwise use optimized_tree
        backbone_tree = self.final_tree if self.final_tree else self.optimized_tree

        # Skip if subtree_files is None or empty or grafting is disabled
        if not subtree_files or self.config.get('subtrees', {}).get('skip_grafting', False):
            self.logger.info("Skipping subtree grafting (no subtrees or disabled)")
            self.final_tree = backbone_tree
            return True

        self.logger.info("Grafting family-level subtrees")

        # Configure the sequence placer for grafting
        placer_config = self.config.get('subtrees', {})
        placer = SequencePlacer(backbone_tree, config=placer_config)

        # Collect subtree files
        if isinstance(subtree_files, str):
            # If it's a directory, find all tree files
            if os.path.isdir(subtree_files):
                files = [
                    os.path.join(subtree_files, f)
                    for f in os.listdir(subtree_files)
                    if f.endswith(('.tree', '.newick', '.nwk'))
                ]
            else:
                # Single file
                files = [subtree_files]
        else:
            # Assume it's a list of files
            files = subtree_files

        if not files:
            self.logger.warning("No subtree files found")
            self.final_tree = backbone_tree
            return True

        self.logger.info(f"Found {len(files)} subtree files")

        # Process each subtree file
        import dendropy
        grafted_tree = backbone_tree

        for subtree_file in files:
            try:
                # Extract family name from filename
                family_name = os.path.splitext(os.path.basename(subtree_file))[0]

                # Load subtree
                subtree = dendropy.Tree.get(path=subtree_file, schema="newick")

                # Find attachment point based on family name
                attachment_point = self._find_attachment_point(grafted_tree, family_name)

                if not attachment_point:
                    self.logger.warning(f"Could not find attachment point for {family_name}")
                    continue

                # Graft subtree
                grafted_tree = placer.graft_subtree(subtree, attachment_point)

                self.logger.info(f"Grafted subtree for {family_name}")

            except Exception as e:
                self.logger.error(f"Error grafting subtree {subtree_file}: {str(e)}")

        self.final_tree = grafted_tree

        # Count how many tips were added
        tip_count = len(self.final_tree.leaf_nodes())
        orig_count = len(backbone_tree.leaf_nodes())
        added_count = tip_count - orig_count

        self.logger.info(f"Grafted {len(files)} subtrees adding {added_count} tips")

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

    def _find_attachment_point(self, tree, family_name):
        """
        Find attachment point for a family subtree.

        Args:
            tree (dendropy.Tree): Tree to search.
            family_name (str): Family name to match.

        Returns:
            str: Taxon label of attachment point, or None if not found.
        """
        # Look for nodes with matching family name
        for node in tree.preorder_node_iter():
            # Check taxon label
            if hasattr(node, 'taxon') and node.taxon:
                label = node.taxon.label
                if label and family_name.lower() in label.lower():
                    return label

            # Check node label
            if hasattr(node, 'label') and node.label:
                if family_name.lower() in node.label.lower():
                    return node.label

        return None