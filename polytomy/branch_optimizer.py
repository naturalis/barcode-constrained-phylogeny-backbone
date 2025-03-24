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
import dendropy
from pathlib import Path
from dendropy import Tree
from Bio import SeqIO


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
        
        # Validate and convert model format based on selected tool
        self.model = self._validate_model_format(self.model)
        
        self.logger.info(f"Branch length optimizer initialized with tool={self.tool}, model={self.model}")

    def _create_sequence_id_map(self, tree_tips, alignment_seqs):
        """Create a mapping between truncated tree IDs and full alignment IDs."""
        id_map = {}
        
        for tip_id in tree_tips:
            # First try exact match
            if tip_id in alignment_seqs:
                id_map[tip_id] = tip_id
                continue
                
            # Try to find if tip_id is a truncated version of any alignment ID
            for seq_id in alignment_seqs:
                if seq_id.endswith(tip_id):
                    self.logger.info(f"Mapped truncated ID {tip_id} to full ID {seq_id}")
                    id_map[tip_id] = seq_id
                    break
        
        return id_map

    def _validate_model_format(self, model):
        """
        Validate the model format based on the selected optimization tool.
        
        Different tools require different model string formats:
        - IQTree: Uses formats like "GTR+G", "GTR+I+G"
        - RAxML-NG: Uses similar format to IQTree
        - Standard RAxML: Uses formats like "GTRGAMMA", "GTRCATI"
        
        Args:
            model (str): The input model string
            
        Returns:
            str: The validated model string
            
        Raises:
            ValueError: If model format is incompatible with the selected tool
        """
        # For RAxML-NG
        if self.tool == "raxml-ng":
            if model in ["GTRGAMMA", "GTRCAT", "GTRCATI"]:
                raise ValueError(
                    f"Invalid model format '{model}' for RAxML-NG. "
                    f"Use 'GTR+G' instead of 'GTRGAMMA', or 'GTR+I+G' instead of 'GTRCATI'."
                )
        
        # For IQTree
        elif self.tool == "iqtree":
            if model in ["GTRGAMMA", "GTRCAT", "GTRCATI"]:
                raise ValueError(
                    f"Invalid model format '{model}' for IQTree. "
                    f"Use 'GTR+G' instead of 'GTRGAMMA', or 'GTR+I+G' instead of 'GTRCATI'."
                )
        
        # Model seems appropriate for the selected tool
        return model

    def prepare_tree_for_optimization(self, tree):
        """Remove QUERY___ prefixes to prevent IQTree from filtering them out"""
        self.logger.info("Preparing tree for optimization by handling QUERY___ prefixes")
        
        # Count prefixed taxa before modification
        prefixed_count = sum(1 for taxon in tree.taxon_namespace if taxon.label.startswith("QUERY___"))
        self.logger.info(f"Found {prefixed_count} taxa with QUERY___ prefix")
        
        # Create a mapping between old and new names
        rename_map = {}
        for taxon in tree.taxon_namespace:
            if taxon.label.startswith("QUERY___"):
                rename_map[taxon.label] = taxon.label[9:]  # Remove prefix
        
        # Apply renaming
        for old_name, new_name in rename_map.items():
            for taxon in tree.taxon_namespace:
                if taxon.label == old_name:
                    taxon.label = new_name
                    break
        
        self.logger.info(f"Renamed {len(rename_map)} taxa by removing QUERY___ prefix")
        return tree

    def _check_tree_alignment_congruence(self, alignment_path: str, workdir: str) -> [str,str]:
        # Get tip labels from the tree
        tree_tips = set([tip.taxon.label for tip in self.tree.leaf_node_iter()])
        
        # Load alignment
        alignment = SeqIO.to_dict(SeqIO.parse(alignment_path, "fasta"))
        alignment_ids = set(alignment.keys())
        
        # Create ID mapping between tree tips and alignment sequences
        id_map = self._create_sequence_id_map(tree_tips, alignment)
        mapped_tips = set(id_map.keys())
        
        # Find missing sequences (tips without a mapping)
        missing_tips = tree_tips - mapped_tips
        
        self.logger.info(f"Tree has {len(tree_tips)} tips, alignment has {len(alignment_ids)} sequences")
        self.logger.info(f"{len(missing_tips)} tips have no matching sequences in alignment")
        
        # Write the alignment with only mapped sequences
        congruent_alignment_path = os.path.join(workdir, "congruent_alignment.fasta")
        with open(congruent_alignment_path, "w") as out:
            for tip_id, aln_id in id_map.items():
                if aln_id in alignment:
                    SeqIO.write(alignment[aln_id], out, "fasta")
        
        # Write the tree with only the mapped tips
        congruent_tree_path = os.path.join(workdir, "congruent_tree.nwk")
        congruent_tree = self.tree.extract_tree_with_taxa_labels(mapped_tips)
        congruent_tree.write(path=congruent_tree_path, schema="newick")
        
        return congruent_tree_path, congruent_alignment_path

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
        
        # Process tree to remove QUERY___ prefixes
        self.tree = self.prepare_tree_for_optimization(self.tree)

        self.logger.info(f"Optimizing branch lengths using {self.tool}")
        
        # Create a copy of the full tree to preserve structure including tips not in alignment
        full_tree = dendropy.Tree(self.tree)
        
        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            # Get tip labels that exist in the alignment
            tree_tips = [tip.taxon.label for tip in self.tree.leaf_node_iter()]
            alignment_seqs = SeqIO.to_dict(SeqIO.parse(alignment, "fasta"))
            
            # Create ID mapping between tree tips and alignment sequences
            id_map = self._create_sequence_id_map(tree_tips, alignment_seqs)
            mapped_tips = set(id_map.keys())
            
            # Create a modified alignment with tree tip IDs
            mapped_alignment_path = os.path.join(temp_dir, "mapped_alignment.fa")
            with open(mapped_alignment_path, "w") as f:
                # Write records for mapped IDs
                for tree_id, aln_id in id_map.items():
                    if aln_id in alignment_seqs:
                        seq = alignment_seqs[aln_id]
                        f.write(f">{tree_id}\n{str(seq.seq)}\n")
            
            alignment_path = mapped_alignment_path
            self.logger.info(f"Created alignment with renamed sequences to match tree IDs")
            
            # Check for missing tips after mapping
            missing_tips = set(tree_tips) - mapped_tips
            missing_count = len(missing_tips)
            
            self.logger.info(f"Tree has {len(tree_tips)} tips, alignment has {len(mapped_tips)} sequences after mapping")
            self.logger.info(f"{missing_count} tips have no matching sequences in alignment after mapping")
            
            # If more than 50% of tips are missing from alignment, we should warn and return original tree
            if missing_count > 0 and missing_count / len(tree_tips) > 0.5:
                self.logger.warning(f"Over 50% of tree tips ({missing_count}/{len(tree_tips)}) not in alignment")
                self.logger.warning("Skipping branch length optimization to preserve tree structure")
                return full_tree
                
            # Create a filtered tree and alignment for optimization
            if missing_count > 0:
                # Extract tree with only tips that have matching sequences
                optimization_tree = self.tree.extract_tree_with_taxa_labels(list(mapped_tips))
                tree_path = self._preserve_taxon_labels(optimization_tree)
            else:
                # No filtering needed
                tree_path = self._preserve_taxon_labels(self.tree)

            # Run the selected tool
            if self.tool == "iqtree":
                optimized_tree_path = self.run_iqtree(alignment_path, tree_path, temp_dir)
            else:  # raxml-ng
                optimized_tree_path = self.run_raxml(alignment_path, tree_path, temp_dir)

            if not optimized_tree_path or not os.path.exists(optimized_tree_path):
                self.logger.error("Branch length optimization failed")
                return full_tree

            # Read the optimized tree
            optimized_tree = dendropy.Tree.get(path=optimized_tree_path, schema="newick", preserve_underscores=True)

            # If we filtered tips, we need to copy optimized branch lengths to the full tree
            if missing_count > 0:
                # Create a mapping of nodes in optimized tree to nodes in full tree by taxon label
                optimized_nodes = {}
                for node in optimized_tree.leaf_node_iter():
                    if node.taxon:
                        optimized_nodes[node.taxon.label] = node
                        
                # Copy optimized branch lengths to the full tree
                for node in full_tree.postorder_node_iter():
                    if node.is_leaf() and node.taxon and node.taxon.label in optimized_nodes:
                        # Find matching node and path back to root in optimized tree
                        opt_node = optimized_nodes[node.taxon.label]
                        
                        # Copy edge length for this tip
                        if opt_node.edge_length is not None:
                            node.edge_length = opt_node.edge_length
                    elif not node.is_leaf():
                        # For internal nodes, we can't easily map between trees
                        # We'll keep the original branch lengths for these
                        pass
                            
                self.logger.info(f"Merged optimized branch lengths into full tree with {len(tree_tips)} tips")
                return full_tree
            else:
                # No filtering was done, just return the optimized tree
                # Copy original taxon namespace to preserve any annotations
                optimized_tree.taxon_namespace = self.tree.taxon_namespace
                
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
        unrooted_tree_path = os.path.join(working_dir, "unrooted_input.tree")

        try:
            # Read the original tree
            tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)

            # Unroot the tree
            tree.is_rooted = False

            # Suppress unifurcations (collapsing unbranched interior nodes)
            tree.suppress_unifurcations()

            # Use _preserve_taxon_labels method instead of directly writing
            unrooted_tree_path = self._preserve_taxon_labels(tree)
                    
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

        # Add a check for None values in the command
        if None in cmd:
            self.logger.error(f"Command contains None value: {cmd}")
            return None

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
        
    def _preserve_taxon_labels(self, tree, schema="newick"):
        """Ensure taxon labels are preserved exactly when writing tree files"""
        output = tempfile.NamedTemporaryFile(delete=False, suffix='.nwk')
        output.close()
        tree.write(
            path=output.name,
            schema=schema,
            suppress_rooting=True,
            suppress_internal_node_labels=False, 
            suppress_edge_lengths=False,
            unquoted_underscores=True
        )
        return output.name

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