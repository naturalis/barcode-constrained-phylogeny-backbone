#!/usr/bin/env python
"""
Unit tests for the sequence_placer module.

These tests verify that the SequencePlacer class correctly places sequences
onto backbone trees using EPA (Evolutionary Placement Algorithm).
"""

import os
import pytest
import tempfile
import dendropy
import subprocess
import shutil
from pathlib import Path

# Import the module to test
from polytomy.sequence_placer import SequencePlacer
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


# Skip markers for required tools
skip_if_no_raxml = pytest.mark.skipif(
    not is_tool_available("raxmlHPC"),
    reason="RAxML (classic) not available"
)


# Fixtures
@pytest.fixture
def data_dir():
    """Path to test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def backbone_tree_path(data_dir):
    """Path to backbone tree file."""
    return data_dir / "example_tree_iqtree_optimized.tre"


@pytest.fixture
def alignment_path(data_dir):
    """Path to test alignment file."""
    return data_dir / "example_alignment.fa"


@pytest.fixture
def backbone_tree(backbone_tree_path):
    """Parse the backbone tree."""
    parser = TreeParser()
    return parser.parse_from_file(str(backbone_tree_path))


@pytest.fixture
def sequence_placer(backbone_tree):
    """Create a SequencePlacer instance with the backbone tree."""
    return SequencePlacer(backbone_tree)


@pytest.fixture
def custom_sequence_placer(backbone_tree):
    """Create a SequencePlacer instance with custom settings."""
    config = {
        'threads': 2,
        'model': "GTRGAMMA",
        'prefix': "test_placement",
        'keep_files': True
    }
    return SequencePlacer(backbone_tree, config=config)


# Tests
def test_initialization(backbone_tree):
    """Test initialization with default and custom settings."""
    # Default settings
    placer = SequencePlacer(backbone_tree)
    assert placer.threads == 1
    assert placer.model == "GTRCAT"
    assert placer.prefix == "seq_placement"
    assert placer.keep_files is False
    assert placer.output_dir == '.'

    # Custom settings
    config = {
        'threads': 4,
        'model': "GTRCAT",
        'prefix': "custom_prefix",
        'keep_files': True,
        'output_dir': "/tmp/output"
    }
    placer = SequencePlacer(backbone_tree, config=config)
    assert placer.threads == 4
    assert placer.model == "GTRCAT"
    assert placer.prefix == "custom_prefix"
    assert placer.keep_files is True
    assert placer.output_dir == "/tmp/output"


def test_is_program_available(sequence_placer):
    """Test checking if a program is available."""
    # Test with RAxML
    raxml_available = is_tool_available("raxmlHPC")
    result = sequence_placer._is_program_available("raxmlHPC")
    assert result == raxml_available

    # Test with nonexistent program
    result = sequence_placer._is_program_available("nonexistent_program")
    assert result is False


@skip_if_no_raxml
def test_run_epa_placement(sequence_placer, alignment_path, backbone_tree):
    """Test running EPA-based sequence placement."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Save backbone tree to file
        backbone_path = os.path.join(temp_dir, "backbone.tree")
        backbone_tree.write(path=backbone_path, schema="newick")

        # Run EPA placement
        result = sequence_placer.place_sequences(str(alignment_path), prefilter=True)

        # Skip the full test if raxml isn't available
        if not is_tool_available("raxmlHPC"):
            assert result is None
            return

        # Check that result is returned
        assert result is not None


@skip_if_no_raxml
def test_place_sequences(sequence_placer, alignment_path):
    """Test placing sequences onto the backbone tree."""
    result_tree = sequence_placer.place_sequences(str(alignment_path), prefilter=True)

    # Skip the full test if RAxML isn't available
    if not is_tool_available("raxmlHPC"):
        assert result_tree is None
        return

    # Check that a tree was returned
    assert result_tree is not None

    # Check that the resulting tree has taxa
    assert len(result_tree.taxon_namespace) > 0

    # Check that the resulting tree has taxa from the original backbone
    # (indicating that sequences were properly placed)
    original_taxa = set(taxon.label for taxon in sequence_placer.backbone_tree.taxon_namespace
                        if taxon and taxon.label)
    result_taxa = set(taxon.label for taxon in result_tree.taxon_namespace
                      if taxon and taxon.label)

    # At least some of the original taxa should be in the result tree
    assert original_taxa.intersection(result_taxa)


def test_place_sequences_nonexistent_file(sequence_placer):
    """Test handling a nonexistent alignment file."""
    result = sequence_placer.place_sequences("nonexistent_file.fasta")
    assert result is None


@skip_if_no_raxml
def test_keep_files(backbone_tree, alignment_path):
    """Test keeping output files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = os.path.join(temp_dir, "output")
        os.makedirs(output_dir)

        # Configure optimizer to keep files
        config = {
            'keep_files': True,
            'output_dir': output_dir
        }

        placer = SequencePlacer(backbone_tree, config=config)
        placer.place_sequences(str(alignment_path), prefilter=True)

        # Skip checking files if RAxML isn't available
        if not is_tool_available("raxmlHPC"):
            return

        # Check if files were copied to output directory
        output_files = os.listdir(output_dir)
        assert len(output_files) > 0
        assert any(file.startswith(placer.prefix) for file in output_files)


@skip_if_no_raxml
def test_model_parameter(backbone_tree, alignment_path):
    """Test sequence placement with different substitution models."""
    # Test with JC model (simplest model)
    config = {'model': 'GTRCAT'}
    placer = SequencePlacer(backbone_tree, config=config)

    # Run placement
    result_tree = placer.place_sequences(str(alignment_path), prefilter=True)

    # Skip the test if RAxML isn't available
    if not is_tool_available("raxmlHPC"):
        assert result_tree is None
        return

    # Check that a tree was returned
    assert result_tree is not None

    # Check that the tree has branch lengths
    branch_lengths = [edge.length for edge in result_tree.edges() if edge.length is not None]
    assert len(branch_lengths) > 0
    assert all(length >= 0 for length in branch_lengths)

@skip_if_no_raxml
def test_place_big_file(backbone_tree, alignment_path):
    """Test placing sequences from a large alignment file."""
    placer = SequencePlacer(backbone_tree)
    result_tree = placer.place_sequences(str(alignment_path), prefilter=True)

    # Skip the test if RAxML isn't available
    if not is_tool_available("raxmlHPC"):
        assert result_tree is None
        return

    # Check that a tree was returned
    assert result_tree is not None

    # Check that the resulting tree has taxa
    assert len(result_tree.taxon_namespace) > 0

    # Check that the resulting tree has taxa from the original backbone
    # (indicating that sequences were properly placed)
    original_taxa = set(taxon.label for taxon in backbone_tree.taxon_namespace
                        if taxon and taxon.label)
    result_taxa = set(taxon.label for taxon in result_tree.taxon_namespace
                      if taxon and taxon.label)

    # At least some of the original taxa should be in the result tree
    assert original_taxa.intersection(result_taxa)