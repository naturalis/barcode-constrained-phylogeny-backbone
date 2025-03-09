#!/usr/bin/env python
"""
Unit tests for the polytomy_finder module.

These tests verify that the PolytomyFinder class correctly identifies
polytomies in phylogenetic trees and provides accurate information about them.
"""

import pytest
import logging
import dendropy
from pathlib import Path

# Import the module to test
from polytomy.polytomy_finder import PolytomyFinder
from polytomy.tree_parser import TreeParser

# Set up logging
logging.basicConfig(level=logging.ERROR)


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
def basic_polytomy_tree():
    """Create a simple tree with known polytomies."""
    newick = "((A,B,C,D)E,(F,G)H,(I,J,K)L)M;"
    return dendropy.Tree.get(data=newick, schema="newick")


@pytest.fixture
def nested_polytomy_tree():
    """Create a tree with nested polytomies."""
    newick = "((A,B,C)D,(E,F,G,(H,I,J,K)L)M)N;"
    return dendropy.Tree.get(data=newick, schema="newick")


@pytest.fixture
def binary_tree():
    """Create a fully binary tree without polytomies."""
    newick = "(((A,B)C,(D,E)F)G,((H,I)J,(K,L)M)N)O;"
    return dendropy.Tree.get(data=newick, schema="newick")


# Tests
def test_find_all_polytomies_basic(basic_polytomy_tree):
    """Test finding polytomies in a simple tree."""
    finder = PolytomyFinder(basic_polytomy_tree)
    polytomies = finder.find_all_polytomies()

    # Should find 2 polytomies (E and L nodes, with 4 and 3 children respectively)
    assert len(polytomies) == 3  # Root node M is also a polytomy with 3 children

    # Check for the specific polytomies
    degrees = sorted([p.degree for p in polytomies])
    assert degrees == [3, 3, 4]  # 3 children for L and M, 4 children for E


def test_find_all_polytomies_nested(nested_polytomy_tree):
    """Test finding polytomies in a tree with nested polytomies."""
    finder = PolytomyFinder(nested_polytomy_tree)
    polytomies = finder.find_all_polytomies()

    # Should find 3 polytomies
    assert len(polytomies) == 3

    # Verify polytomy degrees
    degrees = sorted([p.degree for p in polytomies])
    assert degrees == [3, 4, 4]  # 3 children for D, 4 children for L and M


def test_find_all_polytomies_binary(binary_tree):
    """Test finding polytomies in a fully binary tree (should find none)."""
    finder = PolytomyFinder(binary_tree)
    polytomies = finder.find_all_polytomies()

    # Should find no polytomies
    assert len(polytomies) == 0


def test_find_all_polytomies_real_tree(example_tree):
    """Test finding polytomies in a real tree."""
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    # Should find at least one polytomy in the example tree
    assert len(polytomies) > 0

    # Check if polytomies have the expected properties
    for polytomy in polytomies:
        # Each polytomy should have at least 3 children
        assert polytomy.degree >= 3
        # Each polytomy should have a node reference
        assert polytomy.node is not None
        # Each polytomy should have a depth value
        assert isinstance(polytomy.depth, int)
        # Each polytomy should have a list of children
        assert len(polytomy.children) >= 3


def test_get_polytomies_without_finding(basic_polytomy_tree):
    """Test that get_polytomies() runs find_all_polytomies() if needed."""
    finder = PolytomyFinder(basic_polytomy_tree)
    # Should automatically call find_all_polytomies()
    polytomies = finder.get_polytomies()

    assert len(polytomies) == 3


def test_get_polytomy_stats_basic(basic_polytomy_tree):
    """Test getting statistics about polytomies in a basic tree."""
    finder = PolytomyFinder(basic_polytomy_tree)
    stats = finder.get_polytomy_stats()

    # Check statistics structure
    assert stats['total_polytomies'] == 3
    assert 3 in stats['by_degree']
    assert 4 in stats['by_degree']
    assert stats['by_degree'][3] == 2  # Two polytomies with degree 3
    assert stats['by_degree'][4] == 1  # One polytomy with degree 4
    assert stats['max_degree'] == 4


def test_get_polytomy_stats_empty(binary_tree):
    """Test getting statistics when there are no polytomies."""
    finder = PolytomyFinder(binary_tree)
    stats = finder.get_polytomy_stats()

    assert stats['total_polytomies'] == 0
    assert stats['by_degree'] == {}
    assert stats['max_degree'] == 0


def test_calculate_node_depths(basic_polytomy_tree):
    """Test that node depths are calculated correctly."""
    finder = PolytomyFinder(basic_polytomy_tree)
    finder._calculate_node_depths()

    # Root node should have depth 0
    assert finder._node_depths[basic_polytomy_tree.seed_node] == 0

    # Check some leaf nodes depths
    for node in basic_polytomy_tree.leaf_nodes():
        # All leaf nodes should have depth > 0
        if node in finder._node_depths:
            assert finder._node_depths[node] > 0


def test_polytomy_sorting(nested_polytomy_tree):
    """Test that polytomies are sorted by depth (deepest first)."""
    finder = PolytomyFinder(nested_polytomy_tree)
    polytomies = finder.find_all_polytomies()

    # Check that polytomies are sorted by depth (descending)
    for i in range(1, len(polytomies)):
        assert polytomies[i - 1].depth >= polytomies[i].depth


def test_internal_vs_terminal_polytomies(basic_polytomy_tree):
    """Test identification of internal vs. terminal polytomies."""
    finder = PolytomyFinder(basic_polytomy_tree)
    stats = finder.get_polytomy_stats()

    # Check counts of internal and terminal polytomies
    assert stats['internal_polytomies'] + stats['terminal_polytomies'] == stats['total_polytomies']


def test_example_tree_specific_polytomies(example_tree):
    """Test for specific known polytomies in the example tree."""
    finder = PolytomyFinder(example_tree)
    polytomies = finder.find_all_polytomies()

    # Find a large polytomy in the example tree (if one exists)
    large_polytomies = [p for p in polytomies if p.degree >= 5]

    if large_polytomies:
        # If we have large polytomies, verify they're structured correctly
        for polytomy in large_polytomies:
            assert len(polytomy.children) == polytomy.degree