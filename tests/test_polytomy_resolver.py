#!/usr/bin/env python
"""
Unit tests for the polytomy_resolver module.

These tests verify that the PolytomyResolver can correctly resolve polytomies
in phylogenetic trees using both OpenToL and pruning strategies. Tests use
the real example_tree.tre to ensure realistic test conditions.
"""

import os
import pytest
import logging
import dendropy
from pathlib import Path

# Import modules to test
from polytomy.tree_parser import TreeParser
from polytomy.opentol_client import OpenToLClient
from polytomy.polytomy_resolver import PolytomyResolver

# Set up logging
logging.basicConfig(level=logging.INFO)


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
    parser = TreeParser(config={'schema': {'preserve_underscores': True, 'case_sensitive_taxon_labels': True}})
    return parser.parse_from_file(str(example_tree_path))


@pytest.fixture
def temp_cache_dir(tmp_path):
    """Create a temporary directory for cache files."""
    cache_dir = tmp_path / "opentol_cache"
    cache_dir.mkdir()
    return str(cache_dir)


@pytest.fixture
def opentol_client(temp_cache_dir):
    """Create an OpenToLClient instance with a temporary cache directory."""
    config = {
        'cache_dir': temp_cache_dir,
        'rate_limit': 0.5  # Shorter rate limit for tests
    }
    return OpenToLClient(config=config)


@pytest.fixture
def polytomy_resolver(example_tree, opentol_client):
    """Create a PolytomyResolver instance with the example tree."""
    return PolytomyResolver(example_tree, opentol_client)


def test_example_tree(example_tree):
    """Test that the example tree was parsed correctly."""
    assert example_tree is not None

    # Explore the basal polytomy
    assert example_tree.seed_node.label == 'Psychodidae', "Root node is incorrect"
    assert len(example_tree.seed_node.child_nodes()) == 5, "Incorrect number of children"
    known_labels = {'Phlebotominae', 'Psychodinae', 'Sycoracinae', 'Trichomyiinae', 'Bruchomyiinae'}
    seen_labels = {node.label for node in example_tree.seed_node.child_nodes()}
    assert known_labels == seen_labels, "Incorrect child nodes"

def test_tnrs(polytomy_resolver, example_tree):
    known_otts = {
        369395, # Phlebotominae
        199307, # Psychodinae
        886451, # Sycoracinae
        5036144, # Trichomyiinae
        767849, # Bruchomyiinae
    }
    root = example_tree.seed_node

    # Returns a dictionary of taxon names to OTT IDs
    ott_ids = polytomy_resolver.map_opentol_children(root)
    for name in ott_ids:
        ott_id = ott_ids[name]
        assert ott_id in known_otts

def test_weighted_prune(polytomy_resolver, example_tree):

    # Propagate annotations
    polytomy_resolver.propagate_all_annotations()

    # Prune the root
    polytomy_resolver.weighted_prune(example_tree.seed_node)

    # Root should be bifurcating
    assert len(example_tree.seed_node.child_nodes()) == 2

def test_resolve(polytomy_resolver):
    """Test that the PolytomyResolver can resolve the example tree."""
    polytomy_resolver.resolve_all_polytomies()

    # Check that the tree is now binary
    root = polytomy_resolver.tree.seed_node
    for node in root.postorder_iter():
        assert len(node.child_nodes()) <= 2

    newick = polytomy_resolver.tree.as_string(schema='newick')
    assert newick is not None

