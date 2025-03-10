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


