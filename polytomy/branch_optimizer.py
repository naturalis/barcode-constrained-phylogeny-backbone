#!/usr/bin/env python
"""
Branch Length Optimizer Module - Computes and optimizes branch lengths

This module handles the computation of branch lengths on resolved trees
using external tools such as IQTree and RAxML-NG.
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


class BranchLengthOptimizer:
    """Computes and optimizes branch lengths on resolved trees."""

    def __init__(self, tree: Tree, tool: str="iqtree", config: dict=None):
        """
        Initialize with a tree and tool selection.

        :param tree: The tree to optimize branch lengths for.
        :param tool: Tool to use for optimization ("iqtree" or "raxml-ng").
        :param config: Configuration options for the optimization.
        """
        self.tree = tree
        self.tool = tool.lower()
        self.config = config or {}
        self.logger = logging.getLogger(__name__)

        # Validate tool selection
        if self.tool not in ["iqtree", "raxml-ng"]:
            self.logger.warning(f"Invalid tool '{self.tool}', falling back to iqtree")
            self.tool = "iqtree"

        # Configure optimization parameters
        self.threads = self.config.get('threads', 1)
        self.memory = self.config.get('memory', "4G")
        self.model = self.config.get('model', "GTR+G")
        self.prefix = self.config.get('prefix', "branch_opt")
        self.keep_files = self.config.get('keep_files', False)

        self.logger.info(f"Branch length optimizer initialized with tool={self.tool}, model={self.model}")

    def optimize_branch_lengths(self, alignment_path: str=None) -> Tree:
        """
        Compute optimal branch lengths for the tree.

        :param alignment_path: Path to alignment file. Required if not already set in config.
        :return: The tree with optimized branch lengths.
        """
        # Check if alignment is provided
        alignment = alignment_path or self.config.get('alignment')
        if not alignment:
            self.logger.error("No alignment provided for branch length optimization")
            return None

        # Check if alignment file exists
        if not os.path.exists(alignment):
            self.logger.error(f"Alignment file not found: {alignment}")
            return None

        self.logger.info(f"Optimizing branch lengths using {self.tool}")

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Write tree to temporary file
            tree_path = os.path.join(temp_dir, "input.tree")
            self.tree.write(path=tree_path, schema="newick")

            # Run the selected tool
            if self.tool == "iqtree":
                optimized_tree_path = self.run_iqtree(alignment, tree_path, temp_dir)
            else:  # raxml-ng
                optimized_tree_path = self.run_raxml(alignment, tree_path, temp_dir)

            if not optimized_tree_path or not os.path.exists(optimized_tree_path):
                self.logger.error("Branch length optimization failed")
                return None

            # Read the optimized tree
            import dendropy
            optimized_tree = dendropy.Tree.get(path=optimized_tree_path, schema="newick", preserve_underscores=True)

            # Copy output files if keep_files is True
            if self.keep_files:
                output_dir = self.config.get('output_dir', '.')
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                # Copy all relevant output files
                for file in os.listdir(temp_dir):
                    if file.startswith(self.prefix):
                        src = os.path.join(temp_dir, file)
                        dst = os.path.join(output_dir, file)
                        shutil.copy2(src, dst)

            self.logger.info("Branch length optimization completed")
            return optimized_tree

    # Fix for BranchLengthOptimizer.run_iqtree method
    def run_iqtree(self, alignment_path: str, tree_path: str, working_dir: str) -> str:
        """
        Run IQTree for branch length optimization.

        :param alignment_path: Path to alignment file.
        :param tree_path: Path to input tree file.
        :param working_dir: Directory for temporary files.
        :return: Path to the optimized tree file, or None if failed.
        """
        self.logger.info("Running IQTree for branch length optimization")

        # Build command
        cmd = [
            "iqtree",
            "-s", alignment_path,
            "-te", tree_path,
            "-m", self.model,
            "-nt", str(self.threads),
            "-mem", self.memory,
            "-fixbr",  # Fix tree topology, optimize branch lengths only
            "-quiet",
            "-pre", os.path.join(working_dir, self.prefix)  # Fixed: -pre instead of -prefix
        ]

        # Add any additional options from config
        if 'iqtree_options' in self.config:
            cmd.extend(self.config['iqtree_options'])

        # Execute command
        try:
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )

            # Check for errors
            if result.returncode != 0:
                self.logger.error(f"IQTree failed with exit code {result.returncode}")
                self.logger.error(f"IQTree stderr: {result.stderr}")
                return None

            # Check for output tree
            output_tree = os.path.join(working_dir, f"{self.prefix}.treefile")
            if not os.path.exists(output_tree):
                self.logger.error(f"IQTree output tree not found: {output_tree}")
                return None

            return output_tree

        except Exception as e:
            self.logger.error(f"Error running IQTree: {str(e)}")
            return None

    def run_raxml(self, alignment_path: str, tree_path: str, working_dir: str) -> str:
        """
        Run RAxML-NG for branch length optimization.

        :param alignment_path: Path to alignment file.
        :param tree_path: Path to input tree file.
        :param working_dir: Directory for temporary files.
        :return: Path to the optimized tree file, or None if failed.
        """
        self.logger.info("Running RAxML-NG for branch length optimization")

        # Create an unrooted version of the tree for RAxML-NG
        import dendropy
        unrooted_tree_path = os.path.join(working_dir, "unrooted_input.tree")

        try:
            # Read the original tree
            tree = dendropy.Tree.get(path=tree_path, schema="newick")

            # Unroot the tree
            tree.is_rooted = False

            # Suppress unifurcations (collapsing unbranched interior nodes)
            tree.suppress_unifurcations()

            # Write the unrooted tree
            tree.write(
                path=unrooted_tree_path,
                schema="newick",
                suppress_rooting=True,
                suppress_internal_node_labels=True,
                suppress_leaf_taxon_labels=False,
                suppress_leaf_node_labels=False,
                suppress_edge_lengths=False,
                unquoted_underscores=True
            )
        except Exception as e:
            self.logger.error(f"Error creating unrooted tree: {str(e)}")
            return None

        # Build command using the unrooted tree
        cmd = [
            "raxml-ng",
            "--evaluate",
            "--msa", alignment_path,
            "--tree", unrooted_tree_path,
            "--model", self.model,
            "--threads", str(self.threads),
            "--prefix", os.path.join(working_dir, self.prefix)
        ]

        # Execute command
        try:
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False
            )

            # Log the full output for debugging
            self.logger.debug(f"RAxML-NG stdout: {result.stdout}")
            self.logger.debug(f"RAxML-NG stderr: {result.stderr}")

            # Check for errors
            if result.returncode != 0:
                self.logger.error(f"RAxML-NG failed with exit code {result.returncode}")
                self.logger.error(f"RAxML-NG stderr: {result.stderr}")
                return None

            # Check for output tree
            output_tree = os.path.join(working_dir, f"{self.prefix}.raxml.bestTree")
            if not os.path.exists(output_tree):
                self.logger.error(f"RAxML-NG output tree not found: {output_tree}")
                return None

            return output_tree

        except Exception as e:
            self.logger.error(f"Error running RAxML-NG: {str(e)}")
            return None

    def _is_program_available(self, program: str) -> bool:
        """
        Check if a program is available in PATH.

        :param program: Program name to check.
        :return: True if program is available, False otherwise.
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