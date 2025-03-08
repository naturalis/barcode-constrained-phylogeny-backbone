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
        self.threads = self.config.get('threads', 1)
        self.model = self.config.get('model', "GTR+G")
        self.prefix = self.config.get('prefix', "seq_placement")
        self.keep_files = self.config.get('keep_files', False)
        self.output_dir = self.config.get('output_dir', '.')

        # Track placed sequences
        self.placed_sequences = []

        self.logger.info("Sequence placer initialized")

    def place_sequences(self, sequences, alignment=None):
        """
        Place sequences onto the backbone.

        Args:
            sequences (str): Path to FASTA file with sequences to place.
            alignment (str, optional): Path to reference alignment.
                                      If not provided, sequences will be aligned first.

        Returns:
            dendropy.Tree: The tree with placed sequences.
        """
        if not sequences or not os.path.exists(sequences):
            self.logger.error(f"Sequences file not found: {sequences}")
            return None

        self.logger.info(f"Placing sequences from {sequences} onto backbone tree")

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Write backbone tree to temporary file
            backbone_path = os.path.join(temp_dir, "backbone.tree")
            self.backbone_tree.write(path=backbone_path, schema="newick")

            # If alignment is not provided, create one
            ref_alignment = alignment
            if not ref_alignment:
                ref_alignment = self._create_alignment(sequences, temp_dir)

            if not ref_alignment or not os.path.exists(ref_alignment):
                self.logger.error("Failed to create or find reference alignment")
                return None

            # Run sequence placement
            placement_results = self._run_epa_placement(
                backbone_path,
                ref_alignment,
                sequences,
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
            import dendropy
            result_tree = dendropy.Tree.get(path=result_tree_path, schema="newick", preserve_underscores=True)

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

    def graft_subtree(self, subtree, attachment_point):
        """
        Graft a subtree at the specified attachment point.

        Args:
            subtree (dendropy.Tree): The subtree to graft.
            attachment_point (str): Taxon label or node ID to attach subtree to.

        Returns:
            dendropy.Tree: The tree with grafted subtree.
        """
        if not subtree:
            self.logger.error("No subtree provided for grafting")
            return self.backbone_tree

        self.logger.info(f"Grafting subtree at attachment point: {attachment_point}")

        # Find attachment node
        attachment_node = None
        for node in self.backbone_tree.preorder_node_iter():
            # Check if node has matching taxon label
            if hasattr(node, 'taxon') and node.taxon and node.taxon.label == attachment_point:
                attachment_node = node
                break
            # Check if node has matching label
            if hasattr(node, 'label') and node.label == attachment_point:
                attachment_node = node
                break

        if not attachment_node:
            self.logger.error(f"Attachment point not found: {attachment_point}")
            return self.backbone_tree

        # Clone the backbone tree to avoid modifying the original
        import dendropy
        result_tree = dendropy.Tree(self.backbone_tree)

        # Find the corresponding node in the cloned tree
        for node in result_tree.preorder_node_iter():
            if (hasattr(node, 'taxon') and node.taxon and
                    hasattr(attachment_node, 'taxon') and attachment_node.taxon and
                    node.taxon.label == attachment_node.taxon.label):

                # Get the parent of the attachment node
                parent = node.parent_node
                if not parent:
                    self.logger.error("Cannot graft at root node")
                    return result_tree

                # Remove the attachment node
                parent.remove_child(node)

                # Create a new internal node
                new_node = result_tree.node_factory()
                parent.add_child(new_node)

                # Add the attachment node to the new internal node
                new_node.add_child(node)

                # Add the subtree to the new internal node
                for child in subtree.seed_node.child_node_iter():
                    # Clone the child node and its descendants
                    child_clone = dendropy.Node(child)
                    new_node.add_child(child_clone)

                self.logger.info("Subtree grafted successfully")
                return result_tree

        self.logger.error("Failed to find attachment node in cloned tree")
        return result_tree

    def batch_place_sequences(self, sequence_batches):
        """
        Place multiple batches of sequences.

        Args:
            sequence_batches (list): List of paths to FASTA files with sequences to place.

        Returns:
            dendropy.Tree: The tree with all placed sequences.
        """
        if not sequence_batches:
            self.logger.warning("No sequence batches provided")
            return self.backbone_tree

        current_tree = self.backbone_tree

        for i, batch in enumerate(sequence_batches):
            self.logger.info(f"Placing batch {i + 1}/{len(sequence_batches)}: {batch}")

            result_tree = self.place_sequences(batch)
            if result_tree:
                current_tree = result_tree
                self.backbone_tree = current_tree
            else:
                self.logger.warning(f"Failed to place batch {i + 1}: {batch}")

        return current_tree

    def _create_alignment(self, sequences, temp_dir):
        """
        Create alignment for sequences.

        Args:
            sequences (str): Path to sequences file.
            temp_dir (str): Temporary directory for files.

        Returns:
            str: Path to aligned sequences, or None if failed.
        """
        self.logger.info("Creating sequence alignment")

        # Extract reference sequences from backbone
        ref_seqs_path = os.path.join(temp_dir, "ref_seqs.fasta")
        self._extract_reference_sequences(ref_seqs_path)

        if not os.path.exists(ref_seqs_path):
            self.logger.error("Failed to extract reference sequences")
            return None

        # Combine reference and query sequences
        combined_seqs_path = os.path.join(temp_dir, "combined_seqs.fasta")
        with open(combined_seqs_path, 'w') as outfile:
            with open(ref_seqs_path, 'r') as infile:
                outfile.write(infile.read())
            with open(sequences, 'r') as infile:
                outfile.write(infile.read())

        # Build command for alignment (using MAFFT)
        alignment_output = os.path.join(temp_dir, "aligned.fasta")
        cmd = [
            "mafft",
            "--auto",
            "--thread", str(self.threads),
            "--add", sequences,
            ref_seqs_path
        ]

        # Execute command
        try:
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            with open(alignment_output, 'w') as outfile:
                result = subprocess.run(
                    cmd,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False
                )

            # Check for errors
            if result.returncode != 0:
                self.logger.error(f"Alignment failed with exit code {result.returncode}")
                self.logger.error(f"Alignment stderr: {result.stderr}")
                return None

            if not os.path.exists(alignment_output) or os.path.getsize(alignment_output) == 0:
                self.logger.error("Alignment output file is empty or doesn't exist")
                return None

            return alignment_output

        except Exception as e:
            self.logger.error(f"Error creating alignment: {str(e)}")
            return None

    def _extract_reference_sequences(self, output_path):
        """
        Extract reference sequences from backbone tree.

        Args:
            output_path (str): Path to write sequences to.

        Returns:
            bool: True if successful, False otherwise.
        """
        # For this function to work, we need a way to get sequences from the tree
        # This depends on how sequences are stored with the tree
        # Here we'll use a simple implementation that assumes sequences are available
        # in the taxon.annotations or need to be fetched from a database

        # Placeholder implementation
        self.logger.warning("Extracting reference sequences not fully implemented")

        # Create a dummy reference sequence for testing
        with open(output_path, 'w') as f:
            f.write(">reference_seq1\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")

        return True

    def _run_epa_placement(self, backbone_path, alignment_path, query_seqs_path, temp_dir):
        """
        Run EPA-based sequence placement using RAxML.

        Args:
            backbone_path (str): Path to backbone tree.
            alignment_path (str): Path to reference alignment.
            query_seqs_path (str): Path to query sequences.
            temp_dir (str): Temporary directory for files.

        Returns:
            dict: Placement results, or None if failed.
        """
        self.logger.info("Running EPA sequence placement with RAxML")

        # RAxML classic is preferred for EPA, as RAxML-NG's EPA implementation
        # might work differently or have different parameter names

        # Build command
        cmd = [
            "raxmlHPC",
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