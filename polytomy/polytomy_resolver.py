#!/usr/bin/env python
"""
Polytomy Resolver Module - Resolves polytomies in phylogenetic trees

This module resolves polytomies using a combination of OpenToL topology
information and minimal-loss pruning strategies.
"""

import logging
import random
from collections import defaultdict, deque

from polytomy.tree_parser import TreeParser
from polytomy.polytomy_finder import PolytomyFinder


class PolytomyResolver:
    """Handles resolution of polytomies using various strategies."""

    def __init__(self, tree, opentol_client=None):
        """
        Initialize with a tree and optional OpenToL client.

        Args:
            tree (dendropy.Tree): The tree containing polytomies to resolve.
            opentol_client (OpenToLClient, optional): Client for OpenToL API interactions.
        """
        self.tree = tree
        self.opentol_client = opentol_client
        self.logger = logging.getLogger(__name__)

        # Track resolution statistics
        self.resolved_count = 0
        self.opentol_resolved_count = 0
        self.pruning_resolved_count = 0
        self.pruned_tips = []

        # Store node mappings
        self.node_to_taxon_map = {}
        self.node_to_ott_id_map = {}

        # Initialize mapping from nodes to taxon names and OTT IDs
        self._initialize_taxon_maps()

    def resolve_with_opentol(self, polytomy_node):
        """
        Resolve a polytomy using OpenToL topology.

        Args:
            polytomy_node (Polytomy): The polytomy to resolve.

        Returns:
            bool: True if polytomy was resolved, False otherwise.
        """
        if not self.opentol_client:
            self.logger.warning("No OpenToL client provided for resolve_with_opentol")
            return False

        # Get child nodes from the polytomy
        nodes = polytomy_node.children
        node_count = len(nodes)

        self.logger.info(f"Attempting to resolve polytomy with {node_count} children using OpenToL")

        # Extract taxon names and OTT IDs for the child nodes
        ott_ids = []
        node_to_ott_map = {}

        for node in nodes:
            node_id = id(node)  # Use node object ID as a key

            # Check if node has an OTT ID directly
            if node_id in self.node_to_ott_id_map:
                ott_id = self.node_to_ott_id_map[node_id]
                if ott_id:
                    ott_ids.append(ott_id)
                    node_to_ott_map[ott_id] = node
                    continue

            # If not, try to resolve it
            taxon_name = self._get_taxon_name(node)
            if not taxon_name:
                self.logger.debug(f"No taxon name found for node {node_id}")
                continue

            # Resolve taxon name to OTT ID if needed
            if node_id not in self.node_to_ott_id_map or not self.node_to_ott_id_map[node_id]:
                # If we have a name but need to get the OTT ID
                if self.opentol_client:
                    result = self.opentol_client.resolve_names([taxon_name])
                    if result and taxon_name in result and result[taxon_name]:
                        ott_id = result[taxon_name]['ott_id']
                        self.node_to_ott_id_map[node_id] = ott_id
                        ott_ids.append(ott_id)
                        node_to_ott_map[ott_id] = node
            else:
                # Already have the OTT ID
                ott_id = self.node_to_ott_id_map[node_id]
                ott_ids.append(ott_id)
                node_to_ott_map[ott_id] = node

        # If we don't have at least 3 OTT IDs, we can't resolve the polytomy using OpenToL
        if len(ott_ids) < 3:
            self.logger.debug(f"Insufficient OTT IDs ({len(ott_ids)}) to resolve polytomy")
            return False

        # Get induced subtree from OpenToL
        subtree_data = self.opentol_client.get_induced_subtree(ott_ids)
        if not subtree_data or 'newick' not in subtree_data:
            self.logger.warning("Failed to get valid induced subtree data from OpenToL")
            return False

        # Parse the subtree topology
        newick = subtree_data['newick']

        # Extract the topology information
        try:
            # Use TreeParser to parse the induced subtree
            parser = TreeParser(config={
                'schema': {'preserve_underscores': True, 'case_sensitive_taxon_labels': True }
            })
            opentol_tree = parser.parse_from_string(newick)

            # Build a resolution plan based on the OpenToL tree
            resolution_plan = self._extract_topology_from_opentol_tree(opentol_tree, node_to_ott_map)

            if not resolution_plan:
                self.logger.warning("Could not extract valid topology from OpenToL tree")
                return False

            # Apply the resolution plan
            success = self._apply_resolution_plan(polytomy_node.node, resolution_plan)

            if success:
                self.opentol_resolved_count += 1
                self.resolved_count += 1
                self.logger.info(f"Successfully resolved polytomy with {node_count} children using OpenToL")
                return True
            else:
                self.logger.warning("Failed to apply resolution plan from OpenToL")
                return False

        except Exception as e:
            self.logger.error(f"Error parsing OpenToL subtree: {str(e)}")
            return False

    def resolve_by_pruning(self, polytomy_node, max_loss=None):
        """
        Resolve a polytomy by pruning to minimize tip loss.

        Args:
            polytomy_node (Polytomy): The polytomy to resolve.
            max_loss (int, optional): Maximum number of tips to prune.

        Returns:
            bool: True if polytomy was resolved, False otherwise.
        """
        original_node = polytomy_node.node
        children = list(original_node.child_node_iter())

        # Need at least 3 children for a polytomy
        if len(children) < 3:
            return False

        self.logger.info(f"Resolving polytomy with {len(children)} children by pruning")

        # Count tips under each child
        tip_counts = {}
        for child in children:
            tip_counts[child] = len([n for n in child.leaf_nodes()])

        # Sort children by tip count (ascending)
        sorted_children = sorted(children, key=lambda c: tip_counts[c])

        # If max_loss is specified, check if we can resolve without exceeding it
        total_tips = sum(tip_counts.values())
        if max_loss is not None:
            # Calculate minimum tips we'd need to keep
            min_keep = total_tips - max_loss

            # Check if we can keep enough tips
            cumulative_keep = 0
            for i, child in enumerate(sorted_children):
                if i >= 2:  # We need at least 2 children to form a bifurcating tree
                    cumulative_keep += tip_counts[child]
                    if cumulative_keep >= min_keep:
                        break

            if cumulative_keep < min_keep:
                self.logger.warning(
                    f"Cannot resolve polytomy without exceeding max_loss={max_loss} (would lose {total_tips - cumulative_keep} tips)"
                )
                return False

        # Resolve the polytomy by keeping the two largest children
        # and creating a ladder with the remaining children
        self._create_ladder_topology(original_node, sorted_children)

        # Keep track of pruned tips
        pruned_tip_count = 0
        for i, child in enumerate(sorted_children):
            if i < 2:  # The two smallest children
                for tip in child.leaf_nodes():
                    if hasattr(tip, 'taxon') and tip.taxon is not None:
                        self.pruned_tips.append(tip.taxon.label)
                    pruned_tip_count += 1

        self.pruning_resolved_count += 1
        self.resolved_count += 1
        self.logger.info(f"Resolved polytomy by pruning {pruned_tip_count} tips")

        return True

    def resolve_all_polytomies(self):
        """
        Attempt to resolve all polytomies in the tree.

        Returns:
            tuple: (resolved_count, failed_count, pruned_tips)
        """
        # Find all polytomies
        finder = PolytomyFinder(self.tree)
        polytomies = finder.find_all_polytomies()

        # Reset counters
        self.resolved_count = 0
        self.opentol_resolved_count = 0
        self.pruning_resolved_count = 0
        self.pruned_tips = []

        polytomy_count = len(polytomies)
        if polytomy_count == 0:
            self.logger.info("No polytomies found in tree")
            return (0, 0, [])

        self.logger.info(f"Attempting to resolve {polytomy_count} polytomies")

        # Resolve each polytomy
        failed_count = 0
        for polytomy in polytomies:
            # First try to resolve using OpenToL
            if self.opentol_client:
                resolved = self.resolve_with_opentol(polytomy)
                if resolved:
                    continue

            # If OpenToL resolution fails, use pruning
            resolved = self.resolve_by_pruning(polytomy)
            if not resolved:
                failed_count += 1
                self.logger.warning(f"Failed to resolve polytomy with {polytomy.degree} children")

        # Log resolution statistics
        self.logger.info(f"Resolved {self.resolved_count}/{polytomy_count} polytomies")
        self.logger.info(f"  - {self.opentol_resolved_count} resolved using OpenToL")
        self.logger.info(f"  - {self.pruning_resolved_count} resolved by pruning")
        self.logger.info(f"  - {failed_count} failed to resolve")
        self.logger.info(f"  - {len(self.pruned_tips)} tips pruned")

        return (self.resolved_count, failed_count, self.pruned_tips)

    def get_resolution_stats(self):
        """
        Return statistics about polytomy resolution.

        Returns:
            dict: Statistics about the resolution process.
        """
        return {
            'total_resolved': self.resolved_count,
            'resolved_by_opentol': self.opentol_resolved_count,
            'resolved_by_pruning': self.pruning_resolved_count,
            'pruned_tip_count': len(self.pruned_tips),
            'pruned_tips': self.pruned_tips
        }

    def _initialize_taxon_maps(self):
        """Initialize mappings from nodes to taxon names and OTT IDs."""
        for node in self.tree.preorder_node_iter():
            node_id = id(node)

            # Get taxon name
            taxon_name = self._get_taxon_name(node)
            if taxon_name:
                self.node_to_taxon_map[node_id] = taxon_name

    def _get_taxon_name(self, node):
        """
        Get taxon name for a node.

        Args:
            node (dendropy.Node): The node to get taxon name for.

        Returns:
            str or None: Taxon name if available, None otherwise.
        """
        # Check if node has taxon attribute
        if hasattr(node, 'taxon') and node.taxon is not None:
            return node.taxon.label

        # Check if node has label
        if hasattr(node, 'label') and node.label:
            return node.label

        # If node has no taxon or label, check children if it's an internal node
        if not node.is_leaf():
            # For internal nodes, get the taxon name from a potential genus-level label
            child_taxa = []
            for child in node.child_node_iter():
                if hasattr(child, 'taxon') and child.taxon is not None:
                    child_taxa.append(child.taxon.label)

            # Check if all children share a common genus
            genera = set()
            for taxon in child_taxa:
                if taxon and ' ' in taxon:
                    genus = taxon.split(' ')[0]
                    genera.add(genus)

            if len(genera) == 1:
                return list(genera)[0]

        return None

    def _extract_topology_from_opentol_tree(self, opentol_tree, node_to_ott_map):
        """
        Extract topology information from OpenToL tree.

        Args:
            opentol_tree (dendropy.Tree): The tree from OpenToL.
            node_to_ott_map (dict): Mapping from OTT IDs to nodes.

        Returns:
            list: Resolution plan as a list of node groups.
        """
        # Create a mapping from OTT IDs in the OpenToL tree to leaf nodes
        ott_to_leaf_map = {}
        for node in opentol_tree.leaf_node_iter():
            if hasattr(node, 'taxon') and node.taxon is not None:
                label = node.taxon.label

                # Extract OTT ID from label (format: name_ottXXXXXX)
                if '_ott' in label:
                    parts = label.split('_ott')
                    if len(parts) > 1 and parts[1].isdigit():
                        ott_id = int(parts[1])
                        ott_to_leaf_map[ott_id] = node

        # Build a list of internal nodes and their descendant OTT IDs
        internal_nodes = []
        for node in opentol_tree.preorder_internal_node_iter():
            descendants = []
            for leaf in node.leaf_nodes():
                if hasattr(leaf, 'taxon') and leaf.taxon is not None:
                    label = leaf.taxon.label

                    # Extract OTT ID from label
                    if '_ott' in label:
                        parts = label.split('_ott')
                        if len(parts) > 1 and parts[1].isdigit():
                            ott_id = int(parts[1])
                            descendants.append(ott_id)

            if descendants:
                internal_nodes.append((node, descendants))

        # Sort internal nodes by depth (deepest first)
        internal_nodes.sort(key=lambda x: len(x[1]))

        # Build a resolution plan
        resolution_plan = []
        for _, descendants in internal_nodes:
            if len(descendants) < 2:
                continue

            # Map OTT IDs to original tree nodes
            nodes = []
            for ott_id in descendants:
                if ott_id in node_to_ott_map:
                    nodes.append(node_to_ott_map[ott_id])

            if len(nodes) >= 2:
                resolution_plan.append(nodes)

        return resolution_plan

    def _apply_resolution_plan(self, polytomy_node, resolution_plan):
        """
        Apply a resolution plan to a polytomy.

        Args:
            polytomy_node (dendropy.Node): The polytomy node to resolve.
            resolution_plan (list): List of node groups to combine.

        Returns:
            bool: True if resolution was successful, False otherwise.
        """
        # Get the set of all children
        children = set(polytomy_node.child_node_iter())

        # Make a copy of the children set
        remaining_children = children.copy()

        # Start building the tree from the smallest groups
        for group in resolution_plan:
            # Check if all nodes in the group are still available
            if not all(node in remaining_children for node in group):
                continue

            # Create a new internal node
            new_node = polytomy_node.__class__()

            # Remove the group nodes from the polytomy
            for node in group:
                polytomy_node.remove_child(node)
                remaining_children.remove(node)

                # Add them as children of the new internal node
                new_node.add_child(node)

            # Add the new internal node as a child of the polytomy
            polytomy_node.add_child(new_node)

            # Add the new node to the remaining set
            remaining_children.add(new_node)

        # Check if we've resolved the polytomy
        return len(list(polytomy_node.child_node_iter())) < 3

    def _create_ladder_topology(self, polytomy_node, children):
        """
        Create a ladder topology from a list of children.

        Args:
            polytomy_node (dendropy.Node): The polytomy node to resolve.
            children (list): Sorted list of child nodes (smallest first).

        Returns:
            bool: True if successful, False otherwise.
        """
        # Need at least 3 children for a polytomy
        if len(children) < 3:
            return False

        # Remove all children from the polytomy
        for child in children:
            polytomy_node.remove_child(child)

        # Keep the two smallest children for later
        small_children = children[:2]

        # Start building a ladder with the remaining children
        current_parent = polytomy_node

        # Add children in reverse order (largest first)
        for child in reversed(children[2:]):
            current_parent.add_child(child)

            # Create a new internal node for the next level
            new_node = polytomy_node.__class__()
            current_parent.add_child(new_node)
            current_parent = new_node

        # Add the two smallest children to the final internal node
        for child in small_children:
            current_parent.add_child(child)

        return True