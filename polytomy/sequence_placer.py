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
import re
import numpy as np
import dendropy
from pathlib import Path
from dendropy import Tree
from Bio import SeqIO
from polytomy.tree_parser import TreeParser
from polytomy.branch_optimizer import BranchLengthOptimizer
from polytomy.tree_alignment_matcher import TreeAlignmentMatcher


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
        self.threads = self.config.get('threads', 52)
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
        Uses TreeAlignmentMatcher to ensure compatibility.

        Args:
            alignment_path (str): Path to reference alignment.

        Returns:
            dendropy.Tree: The filtered tree.
        """
        self.logger.info("Filtering tree to include only sequences in alignment")

        with tempfile.TemporaryDirectory() as temp_dir:
            # Write the backbone tree to a temporary file
            tree_path = os.path.join(temp_dir, "backbone.tre")
            self.backbone_tree.write(
                path=tree_path, 
                schema="newick",
                suppress_rooting=True,
                suppress_edge_lengths=False
            )
            
            # Initialize the matcher
            matcher = TreeAlignmentMatcher(config={'output_dir': temp_dir})
            
            # Filter the tree to match alignment
            filtered_tree_path = matcher.match_tree_to_alignment(
                tree_path=tree_path,
                alignment_path=alignment_path
            )
            
            if not filtered_tree_path:
                self.logger.error("Failed to filter tree based on alignment")
                return None
            
            # Parse the filtered tree
            parser = TreeParser()
            filtered_tree = parser.parse_from_file(filtered_tree_path)
            
            self.logger.info(f"Tree filtering completed: {len(filtered_tree.leaf_nodes())} tips remain")
            return filtered_tree

    def _compress_alignment(self, alignment_path: str, suffix: str = ".compressed", retain: int = 700) -> str:
        """
        Compress alignment by selecting the most informative columns.
        
        Args:
            alignment_path: Path to the alignment file
            suffix: Suffix to add to the compressed alignment file name
            retain: Number of columns to retain (default: 700)
            
        Returns:
            str: Path to the compressed alignment
        """
        self.logger.info(f"Compressing alignment to {retain} informative columns")
        
        # Create output path
        output_path = alignment_path + suffix
        
        try:
            # Read sequences
            sequences = list(SeqIO.parse(alignment_path, "fasta"))
            if not sequences:
                self.logger.error(f"No sequences found in {alignment_path}")
                return None
                
            # Convert to numpy array for efficient column operations
            seq_length = len(sequences[0].seq)
            alignment_array = np.array([list(str(seq.seq).upper()) for seq in sequences])
            
            # Calculate column coverage (non-gap characters)
            coverage = np.sum(alignment_array != '-', axis=0)
            
            # Get indices of best columns
            sorted_indices = np.argsort(-coverage)
            top_indices = sorted_indices[:retain]
            
            # Sort to maintain original column order
            top_indices.sort()
            
            # Extract selected columns
            compressed_seqs = alignment_array[:, top_indices]
            
            # Convert back to SeqIO format
            compressed_records = []
            for i, seq in enumerate(sequences):
                seq.seq = ''.join(compressed_seqs[i])
                compressed_records.append(seq)
                
            # Write compressed alignment
            SeqIO.write(compressed_records, output_path, "fasta")
            
            self.logger.info(f"Alignment compressed from {seq_length} to {retain} columns")
            return output_path
            
        except Exception as e:
            self.logger.error(f"Error compressing alignment: {e}")
            return None

    def place_sequences(self, alignment: str, prefilter: bool = False, compress: bool = True) -> Tree:
        """
        Place sequences onto the backbone.

        Args:
            alignment: Path to reference alignment.
            prefilter: Whether to prefilter the tree by removing all tips not in the alignment.
            compress: Whether to compress the alignment to most informative columns.
            
        Returns:
            dendropy.Tree: The tree with placed sequences.
        """

        self.logger.info(f"Placing sequences from {alignment} onto backbone tree")

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            try:
                if not alignment or not os.path.exists(alignment):
                    self.logger.error("Failed to create or find reference alignment")
                    return None

                # Compress alignment if requested
                alignment_to_use = alignment
                if compress:
                    compressed = self._compress_alignment(alignment)
                    if compressed:
                        alignment_to_use = compressed

                tree = self.backbone_tree
                if prefilter and not self.config.get('skip_prefiltering', False):
                    # Filter the backbone tree to only include tips present in the alignment.
                    tree = self._filter_tree(alignment_to_use)
                    if tree is None:
                        self.logger.error("Failed to filter tree based on alignment")
                        return None

                # Write backbone tree to temporary file
                backbone_path = os.path.join(temp_dir, "backbone.tree")
                try:
                    tree.write(path=backbone_path, schema="newick")
                except Exception as e:
                    self.logger.error(f"Failed to write backbone tree: {e}")
                    return None

                # Run sequence placement
                placement_results = self._run_epa_placement(
                    backbone_path,
                    alignment_to_use,
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
                try:
                    parser = TreeParser(config={
                        'schema': {'preserve_underscores': True, 'case_sensitive_taxon_labels': True}
                    })
                    result_tree = parser.parse_from_file(result_tree_path)
                except Exception as e:
                    self.logger.error(f"Failed to parse result tree: {e}")
                    return None

                # Copy output files if keep_files is True
                if self.keep_files:
                    try:
                        if not os.path.exists(self.output_dir):
                            os.makedirs(self.output_dir)

                        # Copy all relevant output files
                        for file in os.listdir(temp_dir):
                            if file.startswith(self.prefix):
                                src = os.path.join(temp_dir, file)
                                dst = os.path.join(self.output_dir, file)
                                shutil.copy2(src, dst)
                    except Exception as e:
                        self.logger.warning(f"Failed to copy output files: {e}")

                self.logger.info("Sequence placement completed successfully")
                return result_tree
                
            except Exception as e:
                self.logger.error(f"Unexpected error during sequence placement: {e}")
                return None

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

    def place_first_exemplars(self, alignment_path: str, pairs_file: str) -> Tree:
        """
        Place only the first exemplar from each pair that isn't represented in the tree.
        
        Args:
            alignment_path: Path to alignment with all sequences
            pairs_file: Path to exemplar pairs file
            
        Returns:
            dendropy.Tree: Tree with first exemplars placed
        """
        self.logger.info(f"Placing first exemplars from {pairs_file}")
        
        try:
            # Read the alignment
            self.logger.info(f"Reading alignment from {alignment_path}")
            try:
                all_sequences = SeqIO.to_dict(SeqIO.parse(alignment_path, "fasta"))
                self.logger.info(f"Read {len(all_sequences)} sequences from alignment")
            except Exception as e:
                self.logger.error(f"Failed to read alignment: {e}")
                return None
            
            # Read pairs file
            self.logger.info(f"Reading exemplar pairs from {pairs_file}")
            pairs = {}  # taxon -> [exemplar1, exemplar2]
            try:
                with open(pairs_file, 'r') as f:
                    for line in f:
                        if line.strip() and not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 3:
                                taxon = parts[0]
                                exemplar1 = parts[1]
                                exemplar2 = parts[2]
                                pairs[taxon] = [exemplar1, exemplar2]
                self.logger.info(f"Read {len(pairs)} exemplar pairs")
            except Exception as e:
                self.logger.error(f"Failed to read pairs file: {e}")
                return None
            
            # Get current tips in tree
            current_tips = {leaf.taxon.label for leaf in self.backbone_tree.leaf_node_iter()}
            self.logger.info(f"Found {len(current_tips)} tips in backbone tree")
            
            # For each pair, check if any representative is in the tree
            to_place = []
            for taxon, exemplars in pairs.items():
                exemplar1, exemplar2 = exemplars
                
                # If neither exemplar is in the tree, add the first one to place
                if exemplar1 not in current_tips and exemplar2 not in current_tips:
                    if exemplar1 in all_sequences:
                        to_place.append(exemplar1)
                    elif exemplar2 in all_sequences:
                        # If first exemplar is not in alignment, try second
                        to_place.append(exemplar2)
                        self.logger.warning(f"First exemplar {exemplar1} not in alignment, using {exemplar2}")
            
            self.logger.info(f"Found {len(to_place)} first exemplars to place")
            
            if not to_place:
                self.logger.info("No exemplars to place, all taxa already represented")
                return self.backbone_tree
            
            # Create temporary directory for placement
            with tempfile.TemporaryDirectory() as temp_dir:
                # Create filtered alignment with backbone sequences + first exemplars
                filtered_path = os.path.join(temp_dir, "first_exemplars.fa")
                
                # Write filtered alignment
                records = []
                for seq_id in list(current_tips) + to_place:
                    if seq_id in all_sequences:
                        records.append(all_sequences[seq_id])
                
                try:
                    SeqIO.write(records, filtered_path, "fasta")
                    self.logger.info(f"Created filtered alignment with {len(records)} sequences")
                except Exception as e:
                    self.logger.error(f"Failed to write filtered alignment: {e}")
                    return None
                
                # Use the existing place_sequences method
                return self.place_sequences(filtered_path, prefilter=True)
        
        except Exception as e:
            self.logger.error(f"Unexpected error during first exemplar placement: {e}")
            return None

    def graft_second_exemplars(self, alignment_path: str, pairs_file: str) -> Tree:
        """
        Graft the second exemplar of each pair next to the first one.
        
        This ensures that exemplars from the same taxon form monophyletic groups.
        The second exemplar is added with a branch length of 0 initially
        (to be optimized in a later step).
        
        Args:
            alignment_path: Path to alignment with all exemplar sequences
            pairs_file: Path to file with exemplar pair information
                        Format: taxon_name<tab>exemplar1<tab>exemplar2
        
        Returns:
            dendropy.Tree: The tree with grafted exemplars
        """
        self.logger.info(f"Grafting second exemplars from {pairs_file}")
        
        # Read the alignment
        alignment = SeqIO.to_dict(SeqIO.parse(alignment_path, "fasta"))
        
        # Read pairs file
        pairs = {}  # taxon -> [exemplar1, exemplar2]
        with open(pairs_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:  # Ensure we have taxon + at least 2 exemplars
                        taxon = parts[0]
                        exemplar1 = parts[1]
                        exemplar2 = parts[2]
                        pairs[taxon] = [exemplar1, exemplar2]
        
        # Get current tips in tree
        current_tips = {leaf.taxon.label for leaf in self.backbone_tree.leaf_node_iter()}
        
        # Clone the tree to avoid modifying the original
        tree = dendropy.Tree(self.backbone_tree)
        
        # Count how many exemplars we graft
        grafted_count = 0
        
        # For each pair, check if the first exemplar is in the tree
        for taxon, exemplars in pairs.items():
            exemplar1, exemplar2 = exemplars
            
            # Skip if second exemplar is already in the tree
            if exemplar2 in current_tips:
                continue
            
            # Skip if first exemplar is not in the tree
            if exemplar1 not in current_tips:
                self.logger.warning(f"First exemplar {exemplar1} for taxon {taxon} not found in tree")
                continue
                
            # Skip if second exemplar is not in the alignment
            if exemplar2 not in alignment:
                self.logger.warning(f"Second exemplar {exemplar2} for taxon {taxon} not found in alignment")
                continue
            
            # Find the node for the first exemplar
            exemplar1_node = None
            for leaf in tree.leaf_node_iter():
                if leaf.taxon and leaf.taxon.label == exemplar1:
                    exemplar1_node = leaf
                    break
            
            if not exemplar1_node:
                self.logger.warning(f"Node for exemplar {exemplar1} not found in tree")
                continue
            
            # Create a node for the second exemplar
            exemplar2_taxon = tree.taxon_namespace.new_taxon(label=exemplar2)
            exemplar2_node = tree.node_factory()
            exemplar2_node.taxon = exemplar2_taxon
            
            # Create a new internal node that will be the parent of both exemplars
            new_parent = tree.node_factory()
            
            # Get parent of exemplar1
            old_parent = exemplar1_node.parent_node
            
            # Remember the original branch length
            original_edge_length = exemplar1_node.edge_length
            
            # Remove exemplar1 from its parent
            old_parent.remove_child(exemplar1_node)
            
            # Add both exemplars to the new parent
            new_parent.add_child(exemplar1_node)
            new_parent.add_child(exemplar2_node)
            
            # Add the new parent to the old parent of exemplar1
            old_parent.add_child(new_parent)
            
            # Set branch length of the new parent to the original branch length of exemplar1
            new_parent.edge_length = original_edge_length
            
            # Set branch lengths of both exemplars to 0
            exemplar1_node.edge_length = 0.0
            exemplar2_node.edge_length = 0.0
            
            grafted_count += 1
        
        self.logger.info(f"Grafted {grafted_count} second exemplars")
        return tree
            
    def reoptimize_tree_after_grafting(self, alignment):
        """Re-optimize branch lengths after grafting second exemplars"""
        self.logger.info("Re-optimizing branch lengths after grafting second exemplars")
        
        # Configure the branch length optimizer
        optimizer_config = self.config.get('optimizer', {})
        if alignment:
            optimizer_config['alignment'] = alignment
        
        # Run branch length optimization
        optimizer = BranchLengthOptimizer(
            self.final_tree,
            tool=optimizer_config.get('tool', 'iqtree'),
            config=optimizer_config
        )
        
        self.final_tree = optimizer.optimize_branch_lengths()
        
        if not self.final_tree:
            self.logger.error("Branch length re-optimization failed")
            return False
        
        self.logger.info("Branch length re-optimization completed successfully")
        return True

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
