#!/usr/bin/env python
"""
Unit tests for the opentol_client module focusing on real polytomy resolution.

These tests use the example tree to test real-world polytomy resolution scenarios
by resolving the names of polytomy descendants and retrieving induced subtrees.
"""

import os
import pytest
import tempfile
import shutil
import dendropy
from pathlib import Path

# Import modules to test
from polytomy.opentol_client import OpenToLClient
from polytomy.tree_parser import TreeParser
from polytomy.polytomy_finder import PolytomyFinder


# Fixtures
@pytest.fixture
def data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def example_tree_path(data_dir):
    """Path to example tree file."""
    return data_dir / "example_tree.tre"


@pytest.fixture
def example_tree(example_tree_path):
    """Parse the example tree."""
    parser = TreeParser()
    return parser.parse_from_file(str(example_tree_path))


@pytest.fixture
def temp_cache_dir():
    """Create a temporary directory for cache files."""
    temp_dir = tempfile.mkdtemp(prefix="opentol_cache_")
    yield temp_dir
    # Clean up
    shutil.rmtree(temp_dir)


@pytest.fixture
def opentol_client(temp_cache_dir):
    """Create an OpenToLClient instance with a temporary cache directory."""
    config = {
        'cache_dir': temp_cache_dir,
        'rate_limit': 0.5  # Shorter rate limit for tests
    }
    return OpenToLClient(config=config)


# Helper function to extract taxonomy name from a node
def get_taxon_name(node):
    """Extract a clean taxonomic name from a node."""
    # First try direct taxon label
    if hasattr(node, 'taxon') and node.taxon is not None and node.taxon.label:
        label = node.taxon.label
        # Remove any process IDs (usually in parentheses)
        if "'" in label:
            # Extract the part between quotes
            parts = label.split("'")
            if len(parts) >= 3:
                return parts[1]
        return label

    # Next try node label
    if hasattr(node, 'label') and node.label:
        return node.label

    # For internal nodes, try to get a genus name
    if not node.is_leaf():
        # Check if this node represents a named taxon (often labeled as genus)
        for child in node.child_node_iter():
            if child.is_leaf() and hasattr(child, 'taxon') and child.taxon is not None:
                child_label = child.taxon.label
                if child_label and "'" in child_label:
                    parts = child_label.split("'")
                    if len(parts) >= 3:
                        # Get genus from species name
                        species_name = parts[1]
                        if " " in species_name:
                            return species_name.split(" ")[0]

    return None


# Tests
def test_find_polytomies(example_tree):
    """Verify that the example tree contains polytomies and print their taxonomic context."""
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    # Verify we have polytomies
    assert len(polytomies) > 0

    # Print info about the polytomies for debugging
    print(f"\nFound {len(polytomies)} polytomies in the example tree:")
    for i, p in enumerate(polytomies[:5]):  # Show the first 5
        print(f"Polytomy {i + 1}: {p.degree} children at depth {p.depth}")

        # Get taxonomic context
        parent_taxa = []
        current = p.node
        while current.parent_node:
            current = current.parent_node
            taxon_name = get_taxon_name(current)
            if taxon_name:
                parent_taxa.append(taxon_name)
                break

        # Get direct child taxa
        child_taxa = []
        for child in p.node.child_node_iter():
            taxon_name = get_taxon_name(child)
            if taxon_name:
                child_taxa.append(taxon_name)

        # Print context
        parent_context = parent_taxa[0] if parent_taxa else "Unknown"
        print(f"  Parent taxon: {parent_context}")
        print(f"  Direct child taxa: {', '.join(child_taxa[:10])}...")


def test_resolve_polytomy_taxa(opentol_client, example_tree):
    """Test resolving taxon names from polytomies in the example tree."""
    # Skip test if OpenToL is disabled
    if opentol_client.skip_opentol:
        pytest.skip("OpenToL API queries are disabled")

    # Find polytomies
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    # Process the first 3 polytomies (or fewer if there aren't 3)
    test_count = min(3, len(polytomies))
    for i in range(test_count):
        polytomy = polytomies[i]
        print(f"\nProcessing polytomy {i + 1} with {polytomy.degree} children")

        # Collect taxonomic names from the direct children
        taxon_names = []
        for child in polytomy.node.child_node_iter():
            name = get_taxon_name(child)
            if name:
                taxon_names.append(name)

        # If we have at least 3 names, try to resolve them
        if len(taxon_names) >= 3:
            print(f"  Found {len(taxon_names)} taxonomic names to resolve")
            print(f"  Names: {', '.join(taxon_names)}")

            # Resolve names
            results = opentol_client.resolve_names(taxon_names)

            # Check results
            resolved_count = sum(1 for r in results.values() if r is not None)
            print(f"  Successfully resolved {resolved_count}/{len(taxon_names)} names")

            # At least some names should resolve if we're not in skip mode
            if not opentol_client.skip_opentol:
                assert resolved_count > 0, f"Failed to resolve any names for polytomy {i + 1}"

            # Print sample of resolved names
            resolved_names = [(name, result['ott_id']) for name, result in results.items()
                              if result is not None]
            print(f"  Resolved: {resolved_names}")
        else:
            print(f"  Insufficient taxonomic names ({len(taxon_names)}) to test resolution")


def test_get_mrca_for_polytomies(opentol_client, example_tree):
    """Test getting MRCA for polytomy descendants in the example tree."""
    # Skip test if OpenToL is disabled
    if opentol_client.skip_opentol:
        pytest.skip("OpenToL API queries are disabled")

    # Find polytomies
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    # Process the first 3 polytomies (or fewer if there aren't 3)
    test_count = min(3, len(polytomies))
    for i in range(test_count):
        polytomy = polytomies[i]
        print(f"\nProcessing polytomy {i + 1} with {polytomy.degree} children for MRCA")

        # Collect taxonomic names from the direct children
        taxon_names = []
        for child in polytomy.node.child_node_iter():
            name = get_taxon_name(child)
            if name:
                taxon_names.append(name)

        # If we have at least 3 names, try to resolve them and find MRCA
        if len(taxon_names) >= 3:
            print(f"  Finding MRCA for taxa: {', '.join(taxon_names)}")

            # Resolve names first
            name_results = opentol_client.resolve_names(taxon_names)

            # Extract OTT IDs from resolved names
            ott_ids = [result['ott_id'] for result in name_results.values()
                       if result is not None]

            if len(ott_ids) >= 2:
                print(f"  Got {len(ott_ids)} OTT IDs for MRCA calculation")
                print(f"  OTT IDs: {ott_ids}")

                # Get MRCA
                mrca_result = opentol_client.get_mrca(ott_ids)

                # Check result
                if mrca_result:
                    assert 'mrca' in mrca_result, "MRCA result missing 'mrca' field"

                    # Print MRCA info
                    mrca_info = mrca_result['mrca']

                    # The OTT ID might be directly in the mrca_info dict or in a nested 'taxon' field
                    ott_id = None
                    name = None

                    if 'taxon' in mrca_info and isinstance(mrca_info['taxon'], dict):
                        ott_id = mrca_info['taxon'].get('ott_id')
                        name = mrca_info['taxon'].get('name')
                    else:
                        ott_id = mrca_info.get('ott_id')
                        name = mrca_info.get('name')

                    print(f"  MRCA: {name} (OTT ID: {ott_id})")

                    # Either directly or in the taxon field, we should have MRCA info
                    assert ott_id is not None or name is not None, "No valid MRCA information found"
                else:
                    print("  No MRCA result returned")
            else:
                print(f"  Insufficient OTT IDs ({len(ott_ids)}) to find MRCA")
        else:
            print(f"  Insufficient taxonomic names ({len(taxon_names)}) to test MRCA")


def test_get_induced_subtree_for_polytomies(opentol_client, example_tree):
    """Test getting induced subtree for polytomy descendants in the example tree."""
    # Skip test if OpenToL is disabled
    if opentol_client.skip_opentol:
        pytest.skip("OpenToL API queries are disabled")

    # Find polytomies
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    # Process the first 3 polytomies (or fewer if there aren't 3)
    test_count = min(3, len(polytomies))
    for i in range(test_count):
        polytomy = polytomies[i]
        print(f"\nProcessing polytomy {i + 1} with {polytomy.degree} children for induced subtree")

        # Collect taxonomic names from the direct children
        taxon_names = []
        for child in polytomy.node.child_node_iter():
            name = get_taxon_name(child)
            if name:
                taxon_names.append(name)

        # If we have at least 3 names, try to resolve them and get induced subtree
        if len(taxon_names) >= 3:
            print(f"  Getting induced subtree for taxa: {', '.join(taxon_names)}")

            # Resolve names first
            name_results = opentol_client.resolve_names(taxon_names)

            # Extract OTT IDs from resolved names
            ott_ids = [result['ott_id'] for result in name_results.values()
                       if result is not None]

            if len(ott_ids) >= 3:
                print(f"  Got {len(ott_ids)} OTT IDs for induced subtree")
                print(f"  OTT IDs: {ott_ids}")

                # Get induced subtree
                subtree_result = opentol_client.get_induced_subtree(ott_ids)

                # Check result
                if subtree_result:
                    assert 'newick' in subtree_result, "Subtree result missing 'newick' field"

                    # Newick string should be non-empty
                    newick = subtree_result['newick']
                    assert len(newick) > 0

                    # Print subtree info
                    print(f"  Induced subtree newick length: {len(newick)}")
                    print(f"  Newick preview: {newick[:100]}...")

                    # Try to parse the newick string with dendropy
                    try:
                        tree = dendropy.Tree.get(data=newick, schema="newick")
                        print(f"  Successfully parsed newick into tree with {len(tree.leaf_nodes())} leaves")
                        assert len(tree.leaf_nodes()) > 0
                    except Exception as e:
                        print(f"  Failed to parse newick: {str(e)}")
                else:
                    print("  No subtree result returned")
            else:
                print(f"  Insufficient OTT IDs ({len(ott_ids)}) to get induced subtree")
        else:
            print(f"  Insufficient taxonomic names ({len(taxon_names)}) to test induced subtree")