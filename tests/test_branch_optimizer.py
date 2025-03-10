#!/usr/bin/env python
"""
Unit tests for the branch_optimizer module.

These tests verify that the BranchLengthOptimizer class correctly optimizes
branch lengths using external tools like IQTree and RAxML-NG.
"""

import os
import pytest
import tempfile
import dendropy
import subprocess
from pathlib import Path

# Import the module to test
from polytomy.branch_optimizer import BranchLengthOptimizer
from polytomy.tree_parser import TreeParser


# Helper function to check if tools are available (not a fixture)
def is_tool_available(tool_name):
    """Check if a command-line tool is available."""
    try:
        subprocess.run(
            [tool_name, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False
        )
        return True
    except FileNotFoundError:
        return False


# Skip markers - using the helper function directly
skip_if_no_iqtree = pytest.mark.skipif(
    not is_tool_available("iqtree"),
    reason="IQTree not available"
)

skip_if_no_raxml = pytest.mark.skipif(
    not is_tool_available("raxml-ng"),
    reason="RAxML-NG not available"
)

skip_if_no_tools = pytest.mark.skipif(
    not is_tool_available("iqtree") or not is_tool_available("raxml-ng"),
    reason="IQTree or RAxML-NG not available"
)


# Fixtures
@pytest.fixture
def data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def example_tree_path(data_dir):
    """Path to example tree file."""
    return data_dir / "example_tree_resolved.tre"


@pytest.fixture
def example_alignment_path(data_dir):
    """Path to example alignment file."""
    return data_dir / "output_in_tree.fasta"


@pytest.fixture
def example_tree(example_tree_path):
    """Parse the example tree."""
    parser = TreeParser()
    return parser.parse_from_file(str(example_tree_path))


# Tests
def test_initialization(example_tree):
    """Test initializing the optimizer with different settings."""
    # Default settings
    optimizer = BranchLengthOptimizer(example_tree)
    assert optimizer.tool == "iqtree"
    assert optimizer.model == "GTR+G"

    # Custom settings
    config = {
        'threads': 4,
        'memory': "8G",
        'model': "GTR+I+G",
        'prefix': "test_opt"
    }
    optimizer = BranchLengthOptimizer(example_tree, tool="raxml-ng", config=config)
    assert optimizer.tool == "raxml-ng"
    assert optimizer.threads == 4
    assert optimizer.memory == "8G"
    assert optimizer.model == "GTR+I+G"
    assert optimizer.prefix == "test_opt"


def test_invalid_tool(example_tree):
    """Test initializing with an invalid tool."""
    optimizer = BranchLengthOptimizer(example_tree, tool="invalid_tool")
    # Should fall back to iqtree
    assert optimizer.tool == "iqtree"


def test_is_program_available(example_tree):
    """Test checking if a program is available."""
    optimizer = BranchLengthOptimizer(example_tree)

    # Test with IQTree
    iqtree_available = is_tool_available("iqtree")
    result = optimizer._is_program_available("iqtree")
    assert result == iqtree_available

    # Test with RAxML-NG
    raxml_available = is_tool_available("raxml-ng")
    result = optimizer._is_program_available("raxml-ng")
    assert result == raxml_available

    # Test with nonexistent program
    result = optimizer._is_program_available("nonexistent_program")
    assert result is False


@skip_if_no_iqtree
def test_run_iqtree(example_tree, example_alignment_path):
    """Test running IQTree on the example tree."""
    optimizer = BranchLengthOptimizer(example_tree)

    with tempfile.TemporaryDirectory() as temp_dir:
        # Save tree to file
        tree_path = os.path.join(temp_dir, "input.tree")
        example_tree.write(path=tree_path, schema="newick")

        # Run IQTree
        result = optimizer.run_iqtree(str(example_alignment_path), tree_path, temp_dir)

        # Check if output file exists
        expected_output = os.path.join(temp_dir, f"{optimizer.prefix}.treefile")
        assert os.path.exists(expected_output)

        # Verify result path
        assert result == expected_output

        # Verify the tree can be parsed
        result_tree = dendropy.Tree.get(path=result, schema="newick")
        assert result_tree is not None

        # Check that the resulting tree has some taxa
        assert len(result_tree.taxon_namespace) > 0


@skip_if_no_raxml
def test_run_raxml(example_tree, example_alignment_path):
    """Test running RAxML-NG on the example tree."""
    optimizer = BranchLengthOptimizer(example_tree, tool="raxml-ng")

    with tempfile.TemporaryDirectory() as temp_dir:
        # Save tree to file
        tree_path = os.path.join(temp_dir, "input.tree")
        example_tree.write(path=tree_path, schema="newick")

        # Run RAxML-NG
        result = optimizer.run_raxml(str(example_alignment_path), tree_path, temp_dir)

        # Check if output file exists
        expected_output = os.path.join(temp_dir, f"{optimizer.prefix}.raxml.bestTree")
        assert os.path.exists(expected_output)

        # Verify result path
        assert result == expected_output

        # Verify the tree can be parsed
        result_tree = dendropy.Tree.get(path=result, schema="newick")
        assert result_tree is not None

        # Check that the resulting tree has some taxa
        assert len(result_tree.taxon_namespace) > 0


@skip_if_no_iqtree
def test_optimize_branch_lengths_iqtree(example_tree, example_alignment_path):
    """Test full branch length optimization using IQTree with example data."""
    optimizer = BranchLengthOptimizer(example_tree)

    # Run optimization
    optimized_tree = optimizer.optimize_branch_lengths(str(example_alignment_path))

    # Check that a tree was returned
    assert optimized_tree is not None

    # Check that branch lengths are set
    branch_lengths = [edge.length for edge in optimized_tree.edges() if edge.length is not None]
    assert len(branch_lengths) > 0
    assert all(length >= 0 for length in branch_lengths)


@skip_if_no_raxml
def test_optimize_branch_lengths_raxml(example_tree, example_alignment_path):
    """Test full branch length optimization using RAxML-NG with example data."""
    optimizer = BranchLengthOptimizer(example_tree, tool="raxml-ng")

    # Run optimization
    optimized_tree = optimizer.optimize_branch_lengths(str(example_alignment_path))

    # Check that a tree was returned
    assert optimized_tree is not None

    # Check that branch lengths are set
    branch_lengths = [edge.length for edge in optimized_tree.edges() if edge.length is not None]
    assert len(branch_lengths) > 0
    assert all(length >= 0 for length in branch_lengths)


def test_optimize_branch_lengths_no_alignment(example_tree):
    """Test handling missing alignment file."""
    optimizer = BranchLengthOptimizer(example_tree)
    result = optimizer.optimize_branch_lengths()

    # Should return None when no alignment is provided
    assert result is None


def test_optimize_branch_lengths_nonexistent_alignment(example_tree):
    """Test handling nonexistent alignment file."""
    optimizer = BranchLengthOptimizer(example_tree)
    result = optimizer.optimize_branch_lengths("nonexistent_alignment.fasta")

    # Should return None when alignment file doesn't exist
    assert result is None


@skip_if_no_iqtree
def test_keep_files(example_tree, example_alignment_path):
    """Test keeping output files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir)

        # Configure optimizer to keep files
        config = {
            'keep_files': True,
            'output_dir': output_dir
        }

        optimizer = BranchLengthOptimizer(example_tree, config=config)
        optimizer.optimize_branch_lengths(str(example_alignment_path))

        # Check if files were copied to output directory
        output_files = os.listdir(output_dir)
        assert len(output_files) > 0
        assert any(file.startswith(optimizer.prefix) for file in output_files)


@skip_if_no_iqtree
def test_model_parameter(example_tree, example_alignment_path):
    """Test optimization with different substitution models."""
    # Test with JC model (simplest model)
    config = {'model': 'JC'}
    optimizer = BranchLengthOptimizer(example_tree, config=config)

    # Run optimization
    optimized_tree = optimizer.optimize_branch_lengths(str(example_alignment_path))

    # Check that a tree was returned
    assert optimized_tree is not None

    # Check that branch lengths are set
    branch_lengths = [edge.length for edge in optimized_tree.edges() if edge.length is not None]
    assert len(branch_lengths) > 0
    assert all(length >= 0 for length in branch_lengths)


@skip_if_no_tools
def test_compare_tools(example_tree, example_alignment_path):
    """Compare branch length optimization results between IQTree and RAxML-NG."""
    # Run optimization with IQTree
    iqtree_optimizer = BranchLengthOptimizer(example_tree, tool="iqtree")
    iqtree_tree = iqtree_optimizer.optimize_branch_lengths(str(example_alignment_path))

    # Run optimization with RAxML-NG
    raxml_optimizer = BranchLengthOptimizer(example_tree, tool="raxml-ng")
    raxml_tree = raxml_optimizer.optimize_branch_lengths(str(example_alignment_path))

    # Both should return a valid tree
    assert iqtree_tree is not None
    assert raxml_tree is not None

    # Compare branch lengths - they won't be identical but should be reasonable
    iqtree_lengths = [edge.length for edge in iqtree_tree.edges() if edge.length is not None]
    raxml_lengths = [edge.length for edge in raxml_tree.edges() if edge.length is not None]

    # Just check that both have branch lengths
    assert len(iqtree_lengths) > 0
    assert len(raxml_lengths) > 0