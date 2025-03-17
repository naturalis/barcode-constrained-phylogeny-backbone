#!/usr/bin/env python
"""
Tree-Alignment Matcher - Ensures tree and alignment files have matching taxa

This module provides functionality to match tree and alignment files by filtering
or mapping taxa between them.
"""

import os
import re
import csv
import logging
import dendropy
from pathlib import Path
from Bio import AlignIO, SeqIO


class TreeAlignmentMatcher:
    """Matches tree and alignment taxa to ensure compatibility."""
    
    def __init__(self, config=None):
        """
        Initialize the tree-alignment matcher with optional configuration.
        
        Args:
            config (dict, optional): Configuration options.
        """
        self.config = config or {}
        self.logger = logging.getLogger(__name__)
        
        # Configure parameters
        self.output_dir = self.config.get('output_dir', '.')

    def get_tree_taxa(self, tree_file, format="newick"):
        """Extract leaf node names from a tree file"""
        self.logger.info(f"Reading tree from: {tree_file}")
        tree = dendropy.Tree.get(path=tree_file, schema=format, preserve_underscores=True)
        taxa = set([leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon is not None])
        self.logger.info(f"Found {len(taxa)} taxa in tree")
        return taxa, tree

    def get_alignment_taxa(self, alignment_file, format="fasta"):
        """Extract sequence IDs from an alignment file"""
        self.logger.info(f"Reading alignment from: {alignment_file}")
        try:
            alignment = AlignIO.read(alignment_file, format)
            taxa = set([seq.id for seq in alignment])
            self.logger.info(f"Found {len(taxa)} sequences in alignment")
            return taxa, alignment
        except Exception as e:
            self.logger.error(f"Error reading alignment: {e}")
            return set(), None

    def filter_alignment(self, input_alignment, output_alignment, taxa_to_keep, format="fasta"):
        """Create a filtered alignment with only the specified taxa"""
        self.logger.info(f"Filtering alignment to include only {len(taxa_to_keep)} specified taxa")
        
        try:
            # Read all sequences
            records = list(SeqIO.parse(input_alignment, format))
            
            # Filter records
            filtered_records = [rec for rec in records if rec.id in taxa_to_keep]
            found_taxa = {rec.id for rec in filtered_records}
            
            # Check if all taxa were found
            missing_taxa = set(taxa_to_keep) - found_taxa
            if missing_taxa:
                self.logger.warning(f"Could not find {len(missing_taxa)} taxa in the alignment")
            
            # Write filtered alignment
            SeqIO.write(filtered_records, output_alignment, format)
            
            self.logger.info(f"Created filtered alignment with {len(filtered_records)} sequences: {output_alignment}")
            return True
        except Exception as e:
            self.logger.error(f"Error filtering alignment: {e}")
            return False

    def filter_tree(self, input_tree, output_tree, taxa_to_keep, format="newick"):
        """Create a filtered tree with only the specified taxa"""
        self.logger.info(f"Filtering tree to keep {len(taxa_to_keep)} taxa")
        
        try:
            # Read tree with conservative settings
            tree = dendropy.Tree.get(
                path=input_tree, 
                schema=format, 
                preserve_underscores=True,
                rooting="default-rooted"
            )
            
            # Count before filtering
            before_count = len(tree.leaf_nodes())
            
            # Convert set to list for retain_taxa_with_labels
            taxa_list = list(taxa_to_keep)
            
            # Filter taxa
            tree.retain_taxa_with_labels(taxa_list)
            
            # Count after filtering
            after_count = len(tree.leaf_nodes())
            
            if after_count == 0:
                self.logger.error("Filtering removed all taxa! Cannot write empty tree.")
                return False
            
            self.logger.info(f"Tree filtered: {before_count} → {after_count} tips")
            
            # Write the tree
            tree.write(
                path=output_tree,
                schema=format,
                suppress_rooting=True,
                suppress_edge_lengths=False,
                unquoted_underscores=True,
                suppress_leaf_node_labels=False,
                suppress_internal_node_labels=True
            )
            
            self.logger.info(f"Written filtered tree to: {output_tree}")
            
            return True
        except Exception as e:
            self.logger.error(f"Error filtering tree: {e}")
            return False

    def analyze_mismatches(self, tree_taxa, aln_taxa):
        """Analyze and report on mismatches between tree and alignment"""
        in_tree_only = tree_taxa - aln_taxa
        in_aln_only = aln_taxa - tree_taxa
        common = tree_taxa.intersection(aln_taxa)
        
        self.logger.info(f"Analysis of taxa mismatches:")
        self.logger.info(f"- Taxa in common: {len(common)}")
        self.logger.info(f"- Taxa only in tree: {len(in_tree_only)}")
        self.logger.info(f"- Taxa only in alignment: {len(in_aln_only)}")
        
        return {
            "common": common,
            "tree_only": in_tree_only,
            "aln_only": in_aln_only
        }

    def match_tree_to_alignment(self, tree_path, alignment_path, output_tree=None, 
                                tree_format="newick", aln_format="fasta"):
        """
        Filter tree to include only tips that are present in alignment.
        
        Args:
            tree_path: Path to input tree
            alignment_path: Path to alignment file
            output_tree: Path for filtered tree (if None, creates one based on input paths)
            tree_format: Format of the tree file
            aln_format: Format of the alignment file
            
        Returns:
            str: Path to filtered tree file, or None if failed
        """
        # Extract taxa
        tree_taxa, _ = self.get_tree_taxa(tree_path, tree_format)
        aln_taxa, _ = self.get_alignment_taxa(alignment_path, aln_format)
        
        # Analyze mismatches
        mismatches = self.analyze_mismatches(tree_taxa, aln_taxa)
        common_taxa = mismatches["common"]
        
        # Set output path if not specified
        if output_tree is None:
            tree_name = Path(tree_path).stem
            output_tree = os.path.join(self.output_dir, f"{tree_name}_filtered.tre")
        
        # Filter tree
        self.logger.info(f"Filtering tree to keep {len(common_taxa)} common taxa")
        success = self.filter_tree(tree_path, output_tree, common_taxa, tree_format)
        
        if success:
            return output_tree
        else:
            return None

    def match_alignment_to_tree(self, alignment_path, tree_path, output_alignment=None,
                              aln_format="fasta", tree_format="newick"):
        """
        Filter alignment to include only sequences matching tree tips.
        
        Args:
            alignment_path: Path to input alignment
            tree_path: Path to tree file
            output_alignment: Path for filtered alignment (if None, creates one based on input paths)
            aln_format: Format of the alignment file
            tree_format: Format of the tree file
            
        Returns:
            str: Path to filtered alignment file, or None if failed
        """
        # Extract taxa
        tree_taxa, _ = self.get_tree_taxa(tree_path, tree_format)
        aln_taxa, _ = self.get_alignment_taxa(alignment_path, aln_format)
        
        # Analyze mismatches
        mismatches = self.analyze_mismatches(tree_taxa, aln_taxa)
        
        # Set output path if not specified
        if output_alignment is None:
            aln_name = Path(alignment_path).stem
            output_alignment = os.path.join(self.output_dir, f"{aln_name}_filtered.fa")
        
        # Filter alignment
        self.logger.info(f"Filtering alignment to keep {len(tree_taxa)} taxa from tree")
        success = self.filter_alignment(alignment_path, output_alignment, tree_taxa, aln_format)
        
        if success:
            return output_alignment
        else:
            return None

    def match_tree_and_alignment(self, tree_path, alignment_path, output_dir=None,
                               tree_format="newick", aln_format="fasta"):
        """
        Filter both tree and alignment to include only common taxa.
        
        Args:
            tree_path: Path to input tree
            alignment_path: Path to alignment file
            output_dir: Directory for output files (if None, uses instance output_dir)
            tree_format: Format of the tree file
            aln_format: Format of the alignment file
            
        Returns:
            tuple: (filtered_tree_path, filtered_alignment_path) or (None, None) if failed
        """
        if output_dir is None:
            output_dir = self.output_dir
            
        # Make sure output dir exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Extract taxa
        tree_taxa, _ = self.get_tree_taxa(tree_path, tree_format)
        aln_taxa, _ = self.get_alignment_taxa(alignment_path, aln_format)
        
        # Analyze mismatches
        mismatches = self.analyze_mismatches(tree_taxa, aln_taxa)
        common_taxa = mismatches["common"]
        
        # Set output paths
        tree_name = Path(tree_path).stem
        aln_name = Path(alignment_path).stem
        output_tree = os.path.join(output_dir, f"{tree_name}_common.tre")
        output_alignment = os.path.join(output_dir, f"{aln_name}_common.fa")
        
        # Filter both files to common taxa
        self.logger.info(f"Filtering tree and alignment to keep {len(common_taxa)} common taxa")
        tree_success = self.filter_tree(tree_path, output_tree, common_taxa, tree_format)
        aln_success = self.filter_alignment(alignment_path, output_alignment, common_taxa, aln_format)
        
        if tree_success and aln_success:
            return (output_tree, output_alignment)
        else:
            return (None, None)