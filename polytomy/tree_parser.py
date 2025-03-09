#!/usr/bin/env python
"""
Tree Parser Module - Parses Newick format trees efficiently

This module provides functionality for parsing phylogenetic trees
in Newick format, with optimizations for handling large trees.
"""

import os
import logging
import dendropy
from dendropy.datamodel import treemodel


class TreeParser:
    """Parses Newick format trees into DendroPy tree objects with memory optimization."""

    def __init__(self, config=None):
        """
        Initialize the tree parser with memory optimization settings.

        Args:
            config (dict, optional): Configuration dictionary with memory settings.
                                    Can include 'max_memory' in MB.
        """
        self.config = config or {}
        self.tree = None
        self.logger = logging.getLogger(__name__)

        # Configure memory settings
        self.max_memory = self.config.get('max_memory', 4000)  # Default 4GB
        self.logger.debug(f"Tree parser initialized with max_memory={self.max_memory}MB")

    def parse_from_file(self, filepath):
        """
        Parse a Newick tree from a file path.

        Args:
            filepath (str): Path to the Newick tree file.

        Returns:
            dendropy.Tree: The parsed tree object.

        Raises:
            FileNotFoundError: If the file doesn't exist.
            ValueError: If the file cannot be parsed as a Newick tree.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Tree file not found: {filepath}")

        self.logger.info(f"Parsing tree from file: {filepath}")

        try:
            # Configure memory-optimized schema settings
            schema_kwargs = self._get_schema_kwargs()

            # Parse the tree
            self.tree = dendropy.Tree.get(
                path=filepath,
                schema="newick",
                **schema_kwargs
            )

            # Log tree statistics
            self._log_tree_stats()

            return self.tree

        except Exception as e:
            self.logger.error(f"Failed to parse tree file: {str(e)}")
            raise ValueError(f"Could not parse tree file: {str(e)}")

    def parse_from_string(self, newick_string):
        """
        Parse a Newick tree from a string.

        Args:
            newick_string (str): Newick tree string.

        Returns:
            dendropy.Tree: The parsed tree object.

        Raises:
            ValueError: If the string cannot be parsed as a Newick tree.
        """
        self.logger.info("Parsing tree from string")

        try:
            # Configure memory-optimized schema settings
            schema_kwargs = self._get_schema_kwargs()

            # Parse the tree
            self.tree = dendropy.Tree.get(
                data=newick_string,
                schema="newick",
                **schema_kwargs
            )

            # Log tree statistics
            self._log_tree_stats()

            return self.tree

        except Exception as e:
            self.logger.error(f"Failed to parse tree string: {str(e)}")
            raise ValueError(f"Could not parse tree string: {str(e)}")

    def get_optimized_tree(self):
        """
        Return the parsed tree with optimized memory settings.

        Returns:
            dendropy.Tree: The optimized tree object.

        Raises:
            ValueError: If no tree has been parsed yet.
        """
        if self.tree is None:
            raise ValueError("No tree has been parsed yet")

        return self.tree

    def _get_schema_kwargs(self):
        """
        Get schema-specific keyword arguments for memory optimization.

        Returns:
            dict: Schema-specific keyword arguments.
        """
        # Settings for optimizing memory usage with large trees
        # These settings are based on DendroPy recommendations
        schema_kwargs = {
            'preserve_underscores': True,
            'suppress_internal_node_taxa': True,
            'suppress_leaf_node_taxa': False,
            'case_sensitive_taxon_labels': False,
            #'suppress_annotations': True,
        }

        # Add any schema-specific settings from config
        if 'schema' in self.config:
            schema_kwargs.update(self.config['schema'])

        return schema_kwargs

    def _log_tree_stats(self):
        """Log statistics about the parsed tree."""
        if self.tree is None:
            return

        try:
            num_tips = len(self.tree.leaf_nodes())
            num_internal = len(self.tree.internal_nodes())
            num_edges = len(self.tree.edges())

            self.logger.info(f"Tree parsed successfully with {num_tips} tips, "
                             f"{num_internal} internal nodes, and {num_edges} edges")
        except Exception as e:
            self.logger.warning(f"Could not compute tree statistics: {str(e)}")