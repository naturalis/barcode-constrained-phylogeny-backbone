#!/usr/bin/env python
"""
Polytomy Finder Module - Identifies polytomies in phylogenetic trees

This module traverses tree structures and identifies nodes with more than
two child nodes (polytomies).
"""

import logging
from collections import namedtuple

# Define a named tuple to represent polytomy information
Polytomy = namedtuple('Polytomy', ['node', 'degree', 'depth', 'children', 'is_internal'])


class PolytomyFinder:
    """Traverses tree structure and identifies polytomies."""

    def __init__(self, tree):
        """
        Initialize with a DendroPy tree object.

        Args:
            tree (dendropy.Tree): The tree to analyze for polytomies.
        """
        self.tree = tree
        self.polytomies = []
        self.logger = logging.getLogger(__name__)

        # Cache for node depths to avoid recalculating
        self._node_depths = {}

    def find_all_polytomies(self):
        """
        Traverse the tree and identify all polytomies.

        Returns:
            list: List of Polytomy objects representing identified polytomy nodes.
        """
        self.logger.info("Searching for polytomies in tree")
        self.polytomies = []

        # Calculate node depths first to avoid repeated traversals
        self._calculate_node_depths()

        # Use post-order traversal to visit children before parents
        for node in self.tree.postorder_node_iter():
            # Skip leaf nodes
            if node.is_leaf():
                continue

            # Check if node has more than 2 children (polytomy)
            child_nodes = list(node.child_node_iter())
            if len(child_nodes) > 2:
                # It's a polytomy
                depth = self._node_depths.get(node, 0)
                is_internal = not any(child.is_leaf() for child in child_nodes)

                # Create a Polytomy object
                polytomy = Polytomy(
                    node=node,
                    degree=len(child_nodes),
                    depth=depth,
                    children=child_nodes,
                    is_internal=is_internal
                )

                self.polytomies.append(polytomy)
                self.logger.debug(f"Found polytomy with {len(child_nodes)} children at depth {depth}")

        # Sort polytomies by depth (deepest first) for efficient resolution
        self.polytomies.sort(key=lambda p: (-p.depth, -p.degree))

        self.logger.info(f"Found {len(self.polytomies)} polytomies in tree")
        return self.polytomies

    def get_polytomies(self):
        """
        Return list of identified polytomy nodes.

        Returns:
            list: List of Polytomy objects.
        """
        if not self.polytomies:
            self.logger.warning("get_polytomies() called but no polytomies have been found yet")
            self.find_all_polytomies()

        return self.polytomies

    def get_polytomy_stats(self):
        """
        Return statistics about identified polytomies.

        Returns:
            dict: Statistics about identified polytomies.
        """
        if not self.polytomies:
            self.find_all_polytomies()

        # Count polytomies by degree
        degree_counts = {}
        for polytomy in self.polytomies:
            degree = polytomy.degree
            degree_counts[degree] = degree_counts.get(degree, 0) + 1

        # Count internal vs. terminal polytomies
        internal_count = sum(1 for p in self.polytomies if p.is_internal)
        terminal_count = len(self.polytomies) - internal_count

        # Calculate max depth of polytomies
        max_depth = max((p.depth for p in self.polytomies), default=0)

        # Build statistics dictionary
        stats = {
            'total_polytomies': len(self.polytomies),
            'by_degree': degree_counts,
            'internal_polytomies': internal_count,
            'terminal_polytomies': terminal_count,
            'max_depth': max_depth,
            'max_degree': max((p.degree for p in self.polytomies), default=0)
        }

        return stats

    def _calculate_node_depths(self):
        """Calculate the depth of each node in the tree."""
        self._node_depths = {}

        # Use level-order traversal to calculate depths efficiently
        # Root node has depth 0
        if not self.tree.seed_node:
            return

        self._node_depths[self.tree.seed_node] = 0

        # Traverse tree in level-order
        for node in self.tree.preorder_node_iter():
            parent_depth = self._node_depths.get(node, 0)

            # Set depth of each child
            for child in node.child_node_iter():
                self._node_depths[child] = parent_depth + 1