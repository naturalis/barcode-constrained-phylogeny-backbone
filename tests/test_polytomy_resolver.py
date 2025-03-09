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
from polytomy.polytomy_finder import PolytomyFinder
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
    parser = TreeParser()
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


@pytest.fixture
def offline_polytomy_resolver(example_tree):
    """Create a PolytomyResolver without OpenToL client for offline tests."""
    return PolytomyResolver(example_tree, None)


# Test helper function
def count_polytomies(tree):
    """Count the number of polytomies in a tree."""
    finder = PolytomyFinder(tree)
    polytomies = finder.find_all_polytomies()
    return len(polytomies)


# Tests
def test_initialize_taxon_maps(polytomy_resolver, example_tree):
    """Test that taxon maps are initialized correctly."""
    # Run the initialization
    polytomy_resolver._initialize_taxon_maps()

    # Check that the taxon map has entries
    assert len(polytomy_resolver.node_to_taxon_map) > 0

    # Verify some key taxa from the example tree are in the map
    found_taxon = False
    for node in example_tree.preorder_node_iter():
        if hasattr(node, 'taxon') and node.taxon and 'Bichromomyia' in str(node.taxon.label):
            node_id = id(node)
            if node_id in polytomy_resolver.node_to_taxon_map:
                found_taxon = True
                break

    assert found_taxon, "Failed to map known taxon from example tree"


def test_get_taxon_name(polytomy_resolver, example_tree):
    """Test extraction of taxon names from nodes."""
    # Test with a leaf node that has a taxon
    leaf_nodes = example_tree.leaf_nodes()
    if leaf_nodes:
        leaf_node = leaf_nodes[0]
        if hasattr(leaf_node, 'taxon') and leaf_node.taxon:
            name = polytomy_resolver._get_taxon_name(leaf_node)
            assert name is not None
            assert isinstance(name, str)
            assert len(name) > 0

    # Test with internal nodes
    for node in example_tree.internal_nodes():
        name = polytomy_resolver._get_taxon_name(node)
        # Not all internal nodes might have a name, but if they do, it should be a string
        if name is not None:
            assert isinstance(name, str)
            assert len(name) > 0


def test_resolve_by_pruning(polytomy_resolver, example_tree):
    """Test resolving polytomies by pruning."""
    # Find a polytomy to resolve
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    if not polytomies:
        pytest.skip("No polytomies found in example tree")

    # Get the first polytomy
    polytomy = polytomies[0]

    # Get the number of tips before resolution
    total_tips_before = len(example_tree.leaf_nodes())

    # Resolve the polytomy by pruning
    result = polytomy_resolver.resolve_by_pruning(polytomy)

    # Check that resolution was successful
    assert result is True

    # Check that the polytomy node now has at most 2 children
    resolved_children = list(polytomy.node.child_node_iter())
    assert len(resolved_children) <= 2

    # Check that some tips were pruned
    total_tips_after = len(example_tree.leaf_nodes())
    assert total_tips_after <= total_tips_before

    # Check that pruned tips were recorded
    assert len(polytomy_resolver.pruned_tips) > 0


def test_resolve_by_pruning_with_max_loss(polytomy_resolver, example_tree):
    """Test resolving polytomies by pruning with a maximum loss constraint."""
    # Find a polytomy to resolve
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    if not polytomies:
        pytest.skip("No polytomies found in example tree")

    # Get the first polytomy
    polytomy = polytomies[0]

    # Get the number of tips before resolution
    total_tips_before = len(example_tree.leaf_nodes())

    # Set a very low max_loss to force failure
    result = polytomy_resolver.resolve_by_pruning(polytomy, max_loss=1)

    # If the polytomy is large, this should fail because max_loss is too small
    if polytomy.degree > 3:
        assert result is False

    # Try with a more reasonable max_loss
    max_loss = polytomy.degree - 1
    result = polytomy_resolver.resolve_by_pruning(polytomy, max_loss=max_loss)

    # This should succeed
    assert result is True

    # Check that the number of pruned tips is within max_loss
    total_tips_after = len(example_tree.leaf_nodes())
    tips_lost = total_tips_before - total_tips_after
    assert tips_lost <= max_loss


def test_create_ladder_topology(offline_polytomy_resolver):
    """Test creating a ladder topology from child nodes."""
    # Create a test tree with a polytomy
    newick = "((A,B,C,D,E)F)G;"
    tree = dendropy.Tree.get(data=newick, schema="newick")

    # Get the polytomy node (F)
    polytomy_node = tree.seed_node.child_nodes()[0]
    children = list(polytomy_node.child_nodes())

    # Create ladder topology
    result = offline_polytomy_resolver._create_ladder_topology(polytomy_node, children)

    # Check success
    assert result is True

    # Check that polytomy is resolved
    new_children = list(polytomy_node.child_nodes())
    assert len(new_children) == 2

    # The structure should now be binary throughout
    for node in tree.preorder_node_iter():
        assert len(list(node.child_node_iter())) <= 2


def test_resolve_all_polytomies_offline(offline_polytomy_resolver, example_tree):
    """Test resolving all polytomies in offline mode (pruning only)."""
    # Count polytomies before resolution
    polytomies_before = count_polytomies(example_tree)

    if polytomies_before == 0:
        pytest.skip("No polytomies found in example tree")

    # Resolve all polytomies
    resolved_count, failed_count, pruned_tips = offline_polytomy_resolver.resolve_all_polytomies()

    # Check that some polytomies were resolved
    assert resolved_count > 0

    # Check that some tips were pruned
    assert len(pruned_tips) > 0

    # Count polytomies after resolution
    polytomies_after = count_polytomies(example_tree)

    # There should be fewer polytomies after resolution
    assert polytomies_after <= polytomies_before - resolved_count


@pytest.mark.skipif(os.environ.get('SKIP_OPENTOL', 'false').lower() == 'true',
                    reason="Skip tests requiring OpenToL API")
def test_resolve_with_opentol(polytomy_resolver, example_tree):
    """Test resolving a polytomy using OpenToL topology."""
    # Skip if OpenToL client is set to skip
    if not polytomy_resolver.opentol_client or polytomy_resolver.opentol_client.skip_opentol:
        pytest.skip("OpenToL client not available or skipping OpenToL")

    # Find a polytomy to resolve
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    if not polytomies:
        pytest.skip("No polytomies found in example tree")

    # Try to resolve each polytomy until one succeeds
    resolved = False
    for polytomy in polytomies[:5]:  # Try first 5 polytomies
        result = polytomy_resolver.resolve_with_opentol(polytomy)
        if result:
            resolved = True
            break

    # Check if any resolution succeeded
    # This is a soft assertion because it depends on OpenToL having data for the taxa
    if not resolved:
        pytest.skip("Could not resolve any polytomies with OpenToL, possibly due to taxa not in OpenToL")

    # If resolution succeeded, check OpenToL resolution count
    assert polytomy_resolver.opentol_resolved_count > 0


@pytest.mark.skipif(os.environ.get('SKIP_OPENTOL', 'false').lower() == 'true',
                    reason="Skip tests requiring OpenToL API")
def test_resolve_all_polytomies_with_opentol(polytomy_resolver, example_tree):
    """Test resolving all polytomies using OpenToL and pruning."""
    # Skip if OpenToL client is set to skip
    if not polytomy_resolver.opentol_client or polytomy_resolver.opentol_client.skip_opentol:
        pytest.skip("OpenToL client not available or skipping OpenToL")

    # Count polytomies before resolution
    polytomies_before = count_polytomies(example_tree)

    if polytomies_before == 0:
        pytest.skip("No polytomies found in example tree")

    # Resolve all polytomies
    resolved_count, failed_count, pruned_tips = polytomy_resolver.resolve_all_polytomies()

    # Check that some polytomies were resolved
    assert resolved_count > 0

    # Get resolution stats
    stats = polytomy_resolver.get_resolution_stats()

    # Check that the stats match the results
    assert stats['total_resolved'] == resolved_count
    assert stats['pruned_tip_count'] == len(pruned_tips)

    # We should have at least some OpenToL resolutions if the API is working
    # This is a soft expectation since it depends on the OpenToL data
    if stats['resolved_by_opentol'] == 0:
        print("Warning: No polytomies resolved by OpenToL. Check API access or taxa presence in OpenToL.")

    # Count polytomies after resolution
    polytomies_after = count_polytomies(example_tree)

    # There should be fewer polytomies after resolution
    assert polytomies_after <= polytomies_before - resolved_count


def test_extract_topology_from_opentol_tree(polytomy_resolver):
    """Test extracting topology information from an OpenToL tree."""
    # Create a simple OpenToL-like tree
    newick = "((A_ott123,B_ott456)node1,(C_ott789,D_ott321)node2)root;"
    opentol_tree = dendropy.Tree.get(data=newick, schema="newick")

    # Create a node to OTT map
    node_to_ott_map = {
        123: dendropy.Node(label="A"),
        456: dendropy.Node(label="B"),
        789: dendropy.Node(label="C"),
        321: dendropy.Node(label="D")
    }

    # Extract topology
    resolution_plan = polytomy_resolver._extract_topology_from_opentol_tree(opentol_tree, node_to_ott_map)

    # Check that a plan was created
    assert resolution_plan is not None
    assert isinstance(resolution_plan, list)

    # The plan should include groupings of nodes
    assert len(resolution_plan) > 0


def test_apply_resolution_plan(polytomy_resolver):
    """Test applying a resolution plan to a polytomy."""
    # Create a test tree with a polytomy
    newick = "((A,B,C,D)E)F;"
    tree = dendropy.Tree.get(data=newick, schema="newick")

    # Get the polytomy node (E)
    polytomy_node = tree.seed_node.child_nodes()[0]
    children = list(polytomy_node.child_nodes())

    # Create a resolution plan
    resolution_plan = [
        [children[0], children[1]],  # Group A and B together
        [children[2], children[3]]  # Group C and D together
    ]

    # Apply the resolution plan
    result = polytomy_resolver._apply_resolution_plan(polytomy_node, resolution_plan)

    # Check success
    assert result is True

    # Check that polytomy is resolved
    new_children = list(polytomy_node.child_nodes())
    assert len(new_children) == 2  # Should now have 2 internal nodes as children


def test_get_resolution_stats(polytomy_resolver, example_tree):
    """Test getting resolution statistics."""
    # Resolve some polytomies first
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    if not polytomies:
        pytest.skip("No polytomies found in example tree")

    # Resolve a few polytomies by pruning
    for polytomy in polytomies[:3]:
        polytomy_resolver.resolve_by_pruning(polytomy)

    # Get stats
    stats = polytomy_resolver.get_resolution_stats()

    # Check that stats include the expected keys
    expected_keys = [
        'total_resolved',
        'resolved_by_opentol',
        'resolved_by_pruning',
        'pruned_tip_count',
        'pruned_tips'
    ]

    for key in expected_keys:
        assert key in stats

    # Check that the stats match the actual resolutions
    assert stats['total_resolved'] == stats['resolved_by_opentol'] + stats['resolved_by_pruning']
    assert stats['pruned_tip_count'] == len(stats['pruned_tips'])