#!/usr/bin/env python
"""
Unit tests for the tree_parser module.

These tests verify that the TreeParser class correctly parses Newick trees
from files and strings, with appropriate error handling and configuration.
"""

import os
import tempfile
import logging
import pytest
from pathlib import Path

# Import the module to test
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
def default_parser():
    """Create a TreeParser with default settings."""
    return TreeParser()


@pytest.fixture
def custom_parser():
    """Create a TreeParser with custom memory settings."""
    return TreeParser(config={'max_memory': 2000})


# Tests
def test_parse_from_file_success(default_parser, example_tree_path):
    """Test parsing a tree from a file."""
    # Parse the example tree
    tree = default_parser.parse_from_file(str(example_tree_path))

    # Basic checks for correct parsing
    assert tree is not None

    # Check tree structure properties
    assert len(tree.leaf_nodes()) > 0

    # Check if taxa are preserved properly
    assert len(tree.taxon_namespace) > 0

    # Find a specific taxon that should be in the example tree
    found = False
    for taxon in tree.taxon_namespace:
        if "AFBR541-14" in taxon.label:
            found = True
            break
    assert found, "Expected taxon not found in parsed tree"


def test_parse_from_file_nonexistent(default_parser):
    """Test parsing a non-existent file raises appropriate error."""
    with pytest.raises(FileNotFoundError):
        default_parser.parse_from_file("nonexistent_file.tre")


def test_parse_from_file_invalid(default_parser):
    """Test parsing an invalid tree file raises appropriate error."""
    # Create a temporary file with invalid Newick content
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
        temp_file.write("This is not a valid Newick string")
        temp_path = temp_file.name

    try:
        with pytest.raises(ValueError):
            default_parser.parse_from_file(temp_path)
    finally:
        # Clean up the temporary file
        os.unlink(temp_path)


def test_parse_from_string_success(default_parser):
    """Test parsing a tree from a Newick string."""
    # Simple Newick string for testing
    newick_string = "(A,B,(C,D));"

    # Parse the string
    tree = default_parser.parse_from_string(newick_string)

    # Basic checks for correct parsing
    assert tree is not None
    assert len(tree.leaf_nodes()) == 4

    # Check taxon labels
    leaf_labels = [leaf.taxon.label for leaf in tree.leaf_nodes()]
    assert "A" in leaf_labels
    assert "B" in leaf_labels
    assert "C" in leaf_labels
    assert "D" in leaf_labels


def test_parse_from_string_invalid(default_parser):
    """Test parsing an invalid Newick string raises appropriate error."""
    with pytest.raises(ValueError):
        default_parser.parse_from_string("This is not a valid Newick string")


def test_memory_settings(custom_parser, example_tree_path):
    """Test memory optimization settings are correctly applied."""
    # Verify the memory setting was applied
    assert custom_parser.max_memory == 2000

    # Parse the example tree with custom memory settings
    tree = custom_parser.parse_from_file(str(example_tree_path))

    # Basic checks that the tree was parsed correctly
    assert tree is not None
    assert len(tree.leaf_nodes()) > 0


def test_schema_kwargs(example_tree_path):
    """Test schema kwargs are correctly applied."""
    # Create parser with custom schema settings
    parser = TreeParser(config={
        'schema': {
            'preserve_underscores': False,
            'case_sensitive_taxon_labels': True
        }
    })

    # Get schema kwargs and check settings
    schema_kwargs = parser._get_schema_kwargs()
    assert not schema_kwargs['preserve_underscores']
    assert schema_kwargs['case_sensitive_taxon_labels']

    # Parse example tree with these settings
    tree = parser.parse_from_file(str(example_tree_path))

    # Basic checks that the tree was parsed correctly
    assert tree is not None
    assert len(tree.leaf_nodes()) > 0


def test_get_optimized_tree_without_parsing():
    """Test getting optimized tree before parsing raises error."""
    fresh_parser = TreeParser()
    with pytest.raises(ValueError):
        fresh_parser.get_optimized_tree()


def test_real_tree_statistics(default_parser, example_tree_path):
    """Test that realistic tree statistics match expectations."""
    tree = default_parser.parse_from_file(str(example_tree_path))

    # Check specific properties of our known example tree
    tip_count = len(tree.leaf_nodes())
    internal_count = len(tree.internal_nodes())

    # These values should match the known structure of example_tree.tre
    # Adjust these numbers to match your actual example tree
    expected_tip_count = 78  # Based on example_tree.txt

    assert tip_count == expected_tip_count, f"Expected {expected_tip_count} tips, got {tip_count}"
    assert internal_count > 0, "Expected at least one internal node"