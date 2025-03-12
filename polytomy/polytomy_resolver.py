#!/usr/bin/env python
"""
Polytomy Resolver Module - Resolves polytomies in phylogenetic trees

This module resolves polytomies using a combination of OpenToL topology
information and minimal-loss pruning strategies.
"""

import re
import logging
from polytomy.tree_parser import TreeParser
from dendropy.datamodel.treemodel import Node

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

    def resolve_all_polytomies(self):
        """
        Attempt to resolve all polytomies in the tree.
        """

        # Traverse the tree from tips to root
        for node in self.tree.postorder_node_iter():

            # Propagate annotations from children to parent node
            self.propagate_annotations(node)

            # Node is polytomy, resolve it
            if len(node.child_nodes()) > 2:
                self.resolve_polytomy(node)

        # Return value?
        # resolved_count, failed_count, pruned_tips

    def resolve_polytomy(self, polytomy: Node) -> bool:
        """
        Resolve a polytomy node using opentol and pruning strategies.

        :param polytomy: The polytomy node to resolve.
        :return: True if the polytomy was resolved, False otherwise.
        """

        # Nothing to be done if focal node is None, is not internal, or is bifurcating
        if polytomy is None or not polytomy.is_internal() or len(polytomy.child_nodes()) < 3:
            return False

        self.logger.info(f"Resolving polytomy for {polytomy.label}")

        # Step 1: Get OTT IDs for all immediate children
        ott_ids = self.map_opentol_children(polytomy)

        # Step 2: Integrate induced subtree from OpenToL
        # Need to be defensive about the TNRS results:
        # if we have fewer than 3 OTT IDs, there's no subtree
        if len(ott_ids) > 2:
            self.graft_induced_subtree(polytomy, ott_ids)

        # Step 3: Remove MRCA leaves
        # If a clade is polyphyletic, OpenToL marks this as 'broken' and
        # returns a leaf with a name like 'mrcaott1234ott5678'. We need to
        # remove these leaves from the tree because fetching the subtree
        # (different endpoint) gives a mess that is very hard to resolve in
        # any general way.
        count = self.remove_mrca_tips(polytomy)
        self.logger.info(f"Removed {count} MRCA leaves from {polytomy}")

        # Step 4: Propagate annotations from tips to focal node
        self.propagate_all_annotations(polytomy)

        # Step 5: Weighted prune
        for n in polytomy.postorder_iter():
            if len(n.child_nodes()) > 2:
                self.weighted_prune(n)

        return True

    def remove_mrca_tips(self, polytomy: Node) -> int:
        """
        Remove MRCA leaves from a polytomy node.

        :param polytomy: The polytomy node to prune.
        :return: The number of MRCA leaves removed.
        """
        self.logger.debug(f"Removing MRCA leaves from {polytomy}")
        pattern = r"mrcaott(\d+)ott(\d+)"
        count = 0
        for leaf in polytomy.leaf_nodes():

            # There is probably no taxon if this is an mrca node
            if leaf.taxon:
                taxon_name = leaf.taxon.label
            else:
                taxon_name = leaf.label

            # If the taxon name matches the pattern, remove the leaf recursively
            if re.match(pattern, taxon_name):
                node = leaf
                count += 1

                # We traverse from the leaf up to the polytomy node
                # We break when the parent still has children after the prune
                while node != polytomy:
                    parent = node.parent_node
                    parent.remove_child(node)
                    if parent.is_internal():
                        break
                    node = parent
        return count

    def propagate_annotations(self, node: Node) -> None:
        """
        Propagate subtended tip number from direct children to focal node.

        :param node: The focal node to propagate annotations to.
        """
        self.logger.debug(f"Propagating annotations for {node}")

        # Node is a tip, set size to 1
        if node.is_leaf():
            node.annotations['size'] = 1

        # Node is internal, sum sizes of children
        else:
            size = 0
            for c in node.child_nodes():
                size += c.annotations['size'].value
            node.annotations['size'] = size

    def propagate_all_annotations(self, root: Node = None) -> None:
        """
        Propagate subtended tip number from all tips to focal node or root.

        :param root: The root node to propagate annotations to.
        """

        # Propagate to the root if no node is specified, otherwise propagate to the focal node
        if root is None:
            root = self.tree.seed_node

        # Note that this is potentially a costly operation
        self.logger.info(f"Propagating all annotations from tips to {root}")
        for node in root.postorder_iter():
            self.propagate_annotations(node)

    def graft_induced_subtree(self, polytomy: Node, ott_ids: dict) -> None:
        """
        Get an induced subtree from OpenToL for a set of OTT IDs.

        :param polytomy: The polytomy node to graft the subtree onto.
        :param ott_ids: A dictionary of taxon names and their corresponding OTT IDs.
        """
        subtree_data = self.opentol_client.get_induced_subtree(list(ott_ids.values()))

        # Parse the subtree topology
        newick = subtree_data['newick']
        parser = TreeParser(config={
            'schema': {'preserve_underscores': True, 'case_sensitive_taxon_labels': True}
        })
        opentol_tree = parser.parse_from_string(newick)

        # Iterate over the tips of the opentol tree. Clean up labels. Match with the polytomy's children.
        matches = []
        for opentol_leaf in opentol_tree.leaf_node_iter():

            # If leaf label matches '_ott' pattern, clean it up for matching
            if '_ott' in opentol_leaf.taxon.label:
                parts = opentol_leaf.taxon.label.split('_ott')
                if len(parts) > 1 and parts[1].isdigit():
                    opentol_leaf.taxon.label = parts[0]

                    # Name match the node's children and graft their children on the leaf
                    for polytomy_child in polytomy.child_nodes():
                        if polytomy_child.label == opentol_leaf.taxon.label:
                            matches.append([opentol_leaf, polytomy_child])
                            break

        # Graft the matched nodes
        for match in matches:
            opentol_leaf, polytomy_child = match

            # Remove the current child from the polytomy
            polytomy.remove_child(polytomy_child)

            # Replace the matched opentol leaf by the polytomy child
            parent = opentol_leaf.parent_node
            parent.remove_child(opentol_leaf)
            parent.add_child(polytomy_child)

        # Here we graft the opentol root's children onto the polytomy.
        for opentol_root_child in opentol_tree.seed_node.child_nodes():
            polytomy.add_child(opentol_root_child)

    def map_opentol_children(self, node: Node) -> dict:
        """
        Resolve taxon names to OTT IDs for all children of a node.

        :param node: The node to resolve children for.
        :return: A dictionary of taxon names and their corresponding OTT IDs.
        """
        ott_ids = {}
        taxon_names = []

        # Aggregate the names of the direct children
        for c in node.child_nodes():

            # Sometimes we arrive at a leaf, after pruning. It probably
            # still doesn't have a taxon, but we can check.
            if c.taxon:
                taxon_name = c.taxon.label
            else:
                taxon_name = c.label
            taxon_names.append(taxon_name)

        # Resolve the taxon name to an OTT ID
        # We need at least 2 names. If we have 2, we can get a split that combines them,
        # taking them out of the pool of polytomy children.
        if len(taxon_names) > 1:
            result = self.opentol_client.resolve_names(taxon_names)
            for taxon_name in taxon_names:
                if result and taxon_name in result and result[taxon_name]:
                    ott_id = result[taxon_name]['ott_id']
                    ott_ids[taxon_name] = ott_id

        # Return the list of OTT IDs
        return ott_ids

    def weighted_prune(self, polytomy: Node):
        """
        Prune a polytomy node to minimize tip loss.

        :param polytomy: The polytomy node to prune.
        """

        # We'll do a Schwartzian transform to sort the children by weight
        child_list = []
        for c in polytomy.child_nodes():

            # The weight is the number of tips in the subtree
            size = int(c.annotations['size'].value)
            child_list.append((c, size))

        # Sort children by weight in descending order
        child_list.sort(key=lambda x: x[1], reverse=True)

        # Keep the two largest children, so start iteration at 3rd element
        total_tips_pruned = 0
        for i, (child, size) in enumerate(child_list[2:], start=2):

            # Remove any other children and add up their sizes
            total_tips_pruned += size
            polytomy.remove_child(child)

        polytomy.annotations['size'] = int(polytomy.annotations['size'].value) - total_tips_pruned
        current_size = int(polytomy.annotations['size'].value)
        self.logger.info(f"Pruned {total_tips_pruned} tips from {polytomy.label} (remaining: {current_size})")