#!/usr/bin/env python
"""
Polytomy Resolver Module - Resolves polytomies in phylogenetic trees

This module resolves polytomies using a combination of OpenToL topology
information and minimal-loss pruning strategies.
"""

import re
import logging
from polytomy.tree_parser import TreeParser

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
            self._propagate_annotations(node)

            # Node is polytomy, resolve it
            if len(node.child_nodes()) > 2:
                self.resolve_polytomy(node)

        # Return value?
        # resolved_count, failed_count, pruned_tips

    def resolve_polytomy(self, polytomy):
        """
        Resolve a polytomy node using opentol and pruning strategies.

        """
        if polytomy is None or not polytomy.is_internal() or len(polytomy.child_nodes()) < 3:
            return False
        leaves = [c for c in polytomy.leaf_nodes()]

        # Step 1: Get OTT IDs for all immediate children
        ott_ids = self._tnrs_children(polytomy)

        # Step 2: Integrate induced subtree from OpenToL
        self._opentol_subtree(polytomy, ott_ids)

        # Step 3: Handle MRCA leaves
        pattern = r"mrcaott(\d+)ott(\d+)"
        for leaf in polytomy.leaf_nodes():
            if re.match(pattern, leaf.taxon.label):
                # self._handle_mrca_leaf(leaves, leaf)
                leaf.parent_node.remove_child(leaf)

        # Step 4: Propagate annotations
        self.propagate_all_annotations(polytomy)

        # Step 5: Weighted prune
        for n in polytomy.postorder_iter():
            if len(n.child_nodes()) > 2:
                self.weighted_prune(n)

    def _propagate_annotations(self, node):
        """
        Propagate annotations from children to parent node.
        """
        if node.is_leaf():
            node.annotations['size'] = 1
            node.annotations['pruned'] = set()
        else:
            size = 0
            pruned = set()
            for c in node.child_nodes():
                size += c.annotations['size'].value
                pruned = pruned.union(c.annotations['pruned'].value)
            node.annotations['size'] = size
            node.annotations['pruned'] = pruned

    def propagate_all_annotations(self, root = None):
        """
        Propagate annotations from tips to root.
        """
        if root is None:
            root = self.tree.seed_node
        for node in root.postorder_iter():
            self._propagate_annotations(node)

    def _opentol_subtree(self, polytomy, ott_ids):
        """
        Get an induced subtree from OpenToL for a set of OTT IDs.
        """
        subtree_data = self.opentol_client.get_induced_subtree(ott_ids)

        # Parse the subtree topology
        newick = subtree_data['newick']
        parser = TreeParser(config={
            'schema': {'preserve_underscores': True, 'case_sensitive_taxon_labels': True}
        })
        opentol_tree = parser.parse_from_string(newick)

        # Iterate over the tips of the opentol tree. Clean up labels. Match with the polytomy's children.
        matches = []
        for leaf in opentol_tree.leaf_node_iter():

            # If leaf label matches '_ott' pattern, clean it up for matching
            if '_ott' in leaf.taxon.label:
                parts = leaf.taxon.label.split('_ott')
                if len(parts) > 1 and parts[1].isdigit():
                    ott_id = int(parts[1])
                    leaf.annotations['ott_id'] = ott_id
                    leaf.taxon.label = parts[0]

                    # Name match the node's children and graft their children on the leaf
                    for c in polytomy.child_nodes():
                        if c.label == leaf.taxon.label:
                            matches.append([leaf, c])
                            break

        # Graft the matched nodes
        for match in matches:
            leaf, node = match
            polytomy.remove_child(node)
            leaf.add_child(node)

        # Graft the opentol tree to the input node
        for c in opentol_tree.seed_node.child_nodes():
            polytomy.add_child(c)

    def _tnrs_children(self, node):
        """
        Resolve taxon names to OTT IDs for all children of a node.
        """
        ott_ids = []
        for c in node.child_nodes():
            if c.is_leaf():

                # Should be impossible to get here because we should only
                # have terminal cherries due to exemplar selection
                raise RuntimeError

            else:
                ott_id = c.annotations['ott_id'].value
                if ott_id:
                    ott_ids.append(ott_id)
                    continue

                # We need to resolve the OTT ID for this node
                taxon_name = c.label
                if taxon_name:
                    result = self.opentol_client.resolve_names([taxon_name])
                    if result and taxon_name in result and result[taxon_name]:
                        ott_id = result[taxon_name]['ott_id']
                        c.annotations['ott_id'] = ott_id
                        ott_ids.append(ott_id)
        return ott_ids

    def weighted_prune(self, polytomy):
        """
        Prune a polytomy node to minimize tip loss.
        """

        # We'll do a Schwartzian transform to sort the children by weight
        child_list = []
        for c in polytomy.child_nodes():

            # The weight is the number of tips in the subtree minus the number of pruned tips
            weight = int(c.annotations['size'].value) - len(c.annotations['pruned'].value)
            child_list.append((c, weight))

        # Sort children by weight in descending order
        child_list.sort(key=lambda x: x[1], reverse=True)

        # Keep the two largest children, so start iteration at 3rd element
        pruned = set()
        for i, (child, _) in enumerate(child_list[2:], start=2):

            # Child may have already been propagated previous pruning results
            pruned = pruned.union(child.annotations['pruned'].value)

            # Child may subtend leaves whose labels we want to record
            for leaf in child.leaf_nodes():
                pruned.add(leaf.taxon.label)

            # Now we can remove the child
            polytomy.remove_child(child)

        polytomy.annotations['pruned'] = pruned

    def _handle_mrca_leaf(self, leaves, opentol_leaf):
        """
        Fetch the polyphyletic ("broken") subtree rooted at the OpenToL leaf
        """

        subtree_data = self.opentol_client.get_subtree(opentol_leaf.taxon.label)

        # Parse the subtree topology
        newick = subtree_data['newick']
        parser = TreeParser(config={
            'schema': {'preserve_underscores': False, 'case_sensitive_taxon_labels': True}
        })
        opentol_tree = parser.parse_from_string(newick)

        # Clean leaf labels of opentol tree to extract OTT ID, copy it to annotation, and map nodes to dict
        leaf_map = {}
        for leaf in opentol_tree.leaf_node_iter():
            if ' ott' in leaf.taxon.label:
                parts = leaf.taxon.label.split(' ott')
                if len(parts) > 1 and parts[1].isdigit():
                    ott_id = int(parts[1])
                    leaf.annotations['ott_id'] = ott_id
                    leaf.taxon.label = parts[0]
            leaf_map[leaf.taxon.label] = leaf

        # Iterate over polytomy leaves, traverse to nearest internal node match, and graft there
        is_mapped = {}
        for leaf in leaves:

            # Tips are process IDs so need to traverse up
            for parent in leaf.ancestor_iter():
                if parent.label in leaf_map:
                    otol_node = leaf_map[parent.label]

                    # Clear existing children in the first pass
                    if parent.label not in is_mapped:
                        otol_node.clear_child_nodes()
                        is_mapped[parent.label] = True

                    # Graft the leaf to the OpenToL node
                    otol_node.add_child(leaf)

        # Do a postorder traversal on the opentol subtree to identify unmapped tips
        tips_to_prune = []
        for onode in opentol_tree.postorder_node_iter():
            if onode.is_leaf():
                label = onode.taxon.label
                if label not in is_mapped:
                    tips_to_prune.append(onode)

        # Prune the tips
        for tip in tips_to_prune:
            tip.parent_node.remove_child(tip)

        # Graft the opentol tree to the input node
        for c in opentol_tree.seed_node.child_nodes():
            opentol_leaf.add_child(c)