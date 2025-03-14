#!/usr/bin/env python
"""
Sequence Placer Module - Places sequences onto backbone trees

This module handles placement of pruned tips and additional sequences
onto a backbone tree using EPA (Evolutionary Placement Algorithm).
"""

import os
import sys
import time
import logging
import tempfile
import subprocess
import shutil
from pathlib import Path
from dendropy import Tree
from Bio import SeqIO
from polytomy.tree_parser import TreeParser


class SequencePlacer:
    """Places sequences onto a backbone tree."""

    def __init__(self, backbone_tree, config=None):
        """
        Initialize with a backbone tree.

        Args:
            backbone_tree (dendropy.Tree): The backbone tree to place sequences on.
            config (dict, optional): Configuration options for sequence placement.
        """
        self.backbone_tree = backbone_tree
        self.config = config or {}
        self.logger = logging.getLogger(__name__)

        # Configure placement parameters
        self.threads = self.config.get('threads', 6)
        self.model = self.config.get('model', "GTRCAT")
        self.prefix = self.config.get('prefix', "seq_placement")
        self.keep_files = self.config.get('keep_files', False)
        self.output_dir = self.config.get('output_dir', '.')

        # Track placed sequences
        self.placed_sequences = []

        self.logger.info("Sequence placer initialized")

    def _filter_tree(self, alignment_path) -> Tree:
        """
        Filter the backbone tree to only include tips present in the alignment.

        Args:
            alignment_path (str): Path to reference alignment.

        Returns:
            dendropy.Tree: The filtered tree.
        """
        self.logger.info("Filtering tree to include only sequences in alignment")

        # Load alignment
        alignment = SeqIO.to_dict(SeqIO.parse(alignment_path, "fasta"))

        # Filter tree tips
        filtered_tree = self.backbone_tree.extract_tree_with_taxa_labels(
            [leaf.taxon.label for leaf in self.backbone_tree.leaf_node_iter() if leaf.taxon.label in alignment]
        )

        self.logger.info("Tree filtering completed")
        return filtered_tree

    def place_sequences(self, alignment: str, prefilter: bool = False) -> Tree:
        """
        Place sequences onto the backbone.

        :param alignment: Path to reference alignment.
        :param prefilter: Whether to prefilter the tree by removing all tips not in the alignment.
        :return: dendropy.Tree: The tree with placed sequences.
        """

        self.logger.info(f"Placing sequences from {alignment} onto backbone tree")

        # Create temporary directory
        with (tempfile.TemporaryDirectory() as temp_dir):

            if not alignment or not os.path.exists(alignment):
                self.logger.error("Failed to create or find reference alignment")
                return None

            tree = self.backbone_tree
            if prefilter:
                # Filter the backbone tree to only include tips present in the alignment.
                tree = self._filter_tree(alignment)

            # Write backbone tree to temporary file
            backbone_path = os.path.join(temp_dir, "backbone.tree")
            tree.write(path=backbone_path, schema="newick")

            # Run sequence placement
            placement_results = self._run_epa_placement(
                backbone_path,
                alignment,
                temp_dir
            )

            if not placement_results:
                self.logger.error("Sequence placement failed")
                return None

            # Parse placement results
            result_tree_path = os.path.join(temp_dir, f"{self.prefix}.raxml.bestTree")
            if not os.path.exists(result_tree_path):
                self.logger.error(f"Placement result tree not found: {result_tree_path}")
                return None

            # Read the tree with placed sequences
            parser = TreeParser(config={
                'schema': {'preserve_underscores': True, 'case_sensitive_taxon_labels': True}
            })
            result_tree = parser.parse_from_file(result_tree_path)

            # Copy output files if keep_files is True
            if self.keep_files:
                if not os.path.exists(self.output_dir):
                    os.makedirs(self.output_dir)

                # Copy all relevant output files
                for file in os.listdir(temp_dir):
                    if file.startswith(self.prefix):
                        src = os.path.join(temp_dir, file)
                        dst = os.path.join(self.output_dir, file)
                        shutil.copy2(src, dst)

            self.logger.info("Sequence placement completed")
            return result_tree

    def _run_epa_placement(self, backbone_path, alignment_path, temp_dir):
        """
        Run EPA-based sequence placement using RAxML.

        Args:
            backbone_path (str): Path to backbone tree.
            alignment_path (str): Path to reference alignment.
            temp_dir (str): Temporary directory for files.

        Returns:
            dict: Placement results, or None if failed.
        """
        self.logger.info("Running EPA sequence placement with RAxML")

        # Build command
        cmd = [
            "raxmlHPC-PTHREADS",
            "-f", "v",  # EPA placement algorithm
            "-t", backbone_path,  # Reference tree
            "-s", alignment_path,  # Alignment with query sequences
            "-m", self.model,  # Model
            "-n", self.prefix,  # Output prefix
            "-T", str(self.threads),  # Threads
            "-w", temp_dir  # Working directory
        ]

        # Add any additional options from config
        if 'raxml_options' in self.config:
            cmd.extend(self.config['raxml_options'])

        # Execute command
        try:
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                cwd=temp_dir
            )

            # Check for errors
            if result.returncode != 0:
                self.logger.error(f"EPA placement failed with exit code {result.returncode}")
                self.logger.error(f"EPA placement stderr: {result.stderr}")
                return None

            # Check for output files
            result_tree = os.path.join(temp_dir, f"RAxML_labelledTree.{self.prefix}")
            if not os.path.exists(result_tree):
                self.logger.error(f"EPA placement result tree not found: {result_tree}")
                return None

            # Copy the result tree to the expected location
            shutil.copy2(result_tree, os.path.join(temp_dir, f"{self.prefix}.raxml.bestTree"))

            # Return a simple result object
            return {
                'tree_path': os.path.join(temp_dir, f"{self.prefix}.raxml.bestTree"),
                'epa_result': result_tree
            }

        except Exception as e:
            self.logger.error(f"Error running EPA placement: {str(e)}")
            return None

    def _is_program_available(self, program):
        """
        Check if a program is available in PATH.

        Args:
            program (str): Program name to check.

        Returns:
            bool: True if program is available, False otherwise.
        """
        try:
            subprocess.run(
                [program, "--version"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False
            )
            return True
        except FileNotFoundError:
            return False
