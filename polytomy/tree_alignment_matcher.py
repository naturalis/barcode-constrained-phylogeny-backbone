#!/usr/bin/env python
"""
Tree-Alignment Matcher - Ensures tree and alignment files have matching taxa

This script compares tree leaf names and alignment sequence IDs and creates
filtered versions where the taxa match perfectly for downstream analyses.
It can also create taxonomic mapping between sequence IDs and tree labels.
"""
import os
import sys
import argparse
import logging
import re
from pathlib import Path
import csv
import numpy as np
import dendropy
from Bio import AlignIO, SeqIO

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("tree_alignment_matcher")

def get_tree_taxa(tree_file, format="newick"):
    """Extract leaf node names from a tree file"""
    logger.info(f"Reading tree from: {tree_file}")
    tree = dendropy.Tree.get(path=tree_file, schema=format, preserve_underscores=True)
    taxa = set([leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon is not None])
    logger.info(f"Found {len(taxa)} taxa in tree")
    return taxa, tree

def get_alignment_taxa(alignment_file, format="fasta"):
    """Extract sequence IDs from an alignment file"""
    logger.info(f"Reading alignment from: {alignment_file}")
    try:
        alignment = AlignIO.read(alignment_file, format)
        taxa = set([seq.id for seq in alignment])
        logger.info(f"Found {len(taxa)} sequences in alignment")
        return taxa, alignment
    except Exception as e:
        logger.error(f"Error reading alignment: {e}")
        sys.exit(1)

def filter_alignment(input_alignment, output_alignment, taxa_to_keep, format="fasta"):
    """Create a filtered alignment with only the specified taxa"""
    logger.info(f"Filtering alignment to include only {len(taxa_to_keep)} specified taxa")
    
    try:
        # Read all sequences
        records = list(SeqIO.parse(input_alignment, format))
        
        # Filter records
        filtered_records = [rec for rec in records if rec.id in taxa_to_keep]
        found_taxa = {rec.id for rec in filtered_records}
        
        # Check if all taxa were found
        missing_taxa = set(taxa_to_keep) - found_taxa
        if missing_taxa:
            logger.warning(f"Could not find {len(missing_taxa)} taxa in the alignment")
            if len(missing_taxa) < 10:
                logger.warning(f"Missing taxa: {missing_taxa}")
            else:
                logger.warning(f"First 10 missing taxa: {list(missing_taxa)[:10]}")
        
        # Write filtered alignment
        SeqIO.write(filtered_records, output_alignment, format)
        
        logger.info(f"Created filtered alignment with {len(filtered_records)} sequences: {output_alignment}")
        return True
    except Exception as e:
        logger.error(f"Error filtering alignment: {e}")
        return False

def filter_tree(input_tree, output_tree, taxa_to_keep, format="newick"):
    """Create a filtered tree with only the specified taxa"""
    logger.info(f"Filtering tree to keep {len(taxa_to_keep)} taxa")
    
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
            logger.error("Filtering removed all taxa! Cannot write empty tree.")
            return False
        
        logger.info(f"Tree filtered: {before_count} tips → {after_count} tips (removed {before_count - after_count})")
        
        # Write to a string first to verify correctness
        newick_string = tree.as_string(
            schema="newick",
            suppress_rooting=True,  # To ensure compatibility with IQTree
            suppress_edge_lengths=False,
            unquoted_underscores=True,
            suppress_leaf_node_labels=False,
            suppress_internal_node_labels=True  # Remove internal labels that can cause parsing issues
        )
        
        # Ensure the string ends with a semicolon
        if not newick_string.strip().endswith(';'):
            newick_string = newick_string.strip() + ';'
        
        # Verify the tree is valid by parsing it again
        try:
            test_tree = dendropy.Tree.get_from_string(
                newick_string,
                schema="newick"
            )
            
            # Check that all taxa were retained
            test_count = len(test_tree.leaf_nodes())
            if test_count != after_count:
                logger.warning(f"Tree validation found inconsistency: {after_count} vs {test_count} tips")
        except Exception as e:
            logger.error(f"Generated tree fails validation: {e}")
            
            # If validation fails, try a simpler approach
            logger.info("Attempting to create simpler tree...")
            simple_tree = dendropy.Tree()
            
            # Create a star tree with just the tips
            for taxon_label in taxa_to_keep:
                if taxon_label in tree.taxon_namespace.labels():
                    taxon = tree.taxon_namespace.get_taxon(label=taxon_label)
                    if taxon:
                        simple_tree.create_node(taxon=taxon)
            
            newick_string = simple_tree.as_string(schema="newick")
        
        # Write the validated tree string
        with open(output_tree, 'w') as f:
            f.write(newick_string)
        
        logger.info(f"Created filtered tree with {after_count} tips (removed {before_count - after_count} tips)")
        logger.info(f"Output file: {output_tree}")
        
        return True
    except Exception as e:
        logger.error(f"Error filtering tree: {e}")
        return False

def validate_newick_file(newick_file):
    """Check if a Newick tree file is valid"""
    logger.info(f"Validating Newick file: {newick_file}")
    
    try:
        with open(newick_file, 'r') as f:
            newick_string = f.read()
        
        # Basic validation
        open_count = newick_string.count('(')
        close_count = newick_string.count(')')
        
        if open_count != close_count:
            logger.error(f"Tree structure is invalid! Unbalanced parentheses: {open_count} opening vs {close_count} closing")
            return False
        
        # Try parsing with DendroPy
        tree = dendropy.Tree.get_from_string(
            newick_string,
            schema="newick"
        )
        
        leaf_count = len(tree.leaf_nodes())
        logger.info(f"Tree validated successfully: {leaf_count} leaf nodes")
        return True
    except Exception as e:
        logger.error(f"Tree validation failed: {e}")
        return False

def analyze_mismatches(tree_taxa, aln_taxa):
    """Analyze and report on mismatches between tree and alignment"""
    in_tree_only = tree_taxa - aln_taxa
    in_aln_only = aln_taxa - tree_taxa
    common = tree_taxa.intersection(aln_taxa)
    
    logger.info(f"Analysis of taxa mismatches:")
    logger.info(f"- Taxa in common: {len(common)}")
    logger.info(f"- Taxa only in tree: {len(in_tree_only)}")
    logger.info(f"- Taxa only in alignment: {len(in_aln_only)}")
    
    if in_tree_only:
        logger.info(f"First 10 taxa only in tree: {list(in_tree_only)[:10]}")
    
    if in_aln_only:
        logger.info(f"First 10 taxa only in alignment: {list(in_aln_only)[:10]}")
    
    return {
        "common": common,
        "tree_only": in_tree_only,
        "aln_only": in_aln_only
    }

def create_taxonomy_mapping(tree_taxa, aln_taxa, output_file, tree_obj=None, aln_obj=None):
    """
    Create a mapping file between taxonomic names and sequence IDs
    
    This function uses several approaches to guess which sequence IDs might be 
    associated with taxonomic names in the tree.
    """
    logger.info("Creating taxonomy-to-sequence mapping file")
    
    # Initialize mapping dictionary
    mapping = {}
    mapped_count = 0
    
    # 1. Direct matches (exact sequence IDs that appear in both)
    direct_matches = tree_taxa.intersection(aln_taxa)
    for match in direct_matches:
        mapping[match] = match
    mapped_count += len(direct_matches)
    logger.info(f"Found {len(direct_matches)} direct matches")
    
    # 2. Look for BOLD IDs in both tree and alignment
    bold_pattern = re.compile(r'[A-Z]+\d+-\d+')
    tree_bold_ids = {taxon for taxon in tree_taxa if bold_pattern.match(taxon)}
    aln_bold_ids = {seq_id for seq_id in aln_taxa if bold_pattern.match(seq_id)}
    logger.info(f"Found {len(tree_bold_ids)} BOLD IDs in tree and {len(aln_bold_ids)} in alignment")
    
    # 3. Extract genus names from tree taxa (assuming binomial nomenclature)
    genera = {}
    for taxon in tree_taxa:
        if '_' in taxon:  # Likely a scientific name
            genus = taxon.split('_')[0]
            if genus not in genera:
                genera[genus] = []
            genera[genus].append(taxon)
    logger.info(f"Extracted {len(genera)} potential genera from tree")
    
    # 4. For remaining unmapped tree taxa, try to find potential matching sequences
    remaining_tree_taxa = tree_taxa - set(mapping.keys())
    remaining_aln_taxa = aln_taxa - set(mapping.values())
    
    # Write mapping file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['tree_label', 'sequence_id', 'match_type'])
        
        # First write direct matches
        for taxon in direct_matches:
            writer.writerow([taxon, taxon, 'exact_match'])
        
        # Write remaining tree taxa with no matches yet
        for taxon in sorted(remaining_tree_taxa):
            writer.writerow([taxon, '', 'no_match'])
        
        # Write alignment sequence IDs not matched to any tree taxon
        for seq_id in sorted(remaining_aln_taxa):
            writer.writerow(['', seq_id, 'unmatched_sequence'])
    
    logger.info(f"Created mapping file with {mapped_count} mappings: {output_file}")
    logger.info(f"Remaining unmapped tree taxa: {len(remaining_tree_taxa)}")
    logger.info(f"Remaining unmapped sequence IDs: {len(remaining_aln_taxa)}")
    
    return mapping

def create_mapped_alignment(input_alignment, output_alignment, mapping, format="fasta"):
    """Create a new alignment with sequences renamed according to mapping"""
    logger.info(f"Creating mapped alignment: {output_alignment}")
    
    # Read alignment
    records = list(SeqIO.parse(input_alignment, format))
    
    # Create reverse mapping (sequence_id -> tree_label)
    reverse_map = {v: k for k, v in mapping.items() if v}
    
    # Count matches
    matched = 0
    
    # Create new alignment
    new_records = []
    for record in records:
        if record.id in reverse_map:
            # This sequence has a mapping to a tree label
            record.id = reverse_map[record.id]
            record.description = f"mapped_from={record.id}"
            new_records.append(record)
            matched += 1
        else:
            # Keep original sequence (could be filtered later)
            new_records.append(record)
    
    # Write new alignment
    SeqIO.write(new_records, output_alignment, format)
    logger.info(f"Created mapped alignment with {matched} renamed sequences")
    
    return matched

def relabel_tree_with_mapping(input_tree, output_tree, mapping, format="newick"):
    """Create a new tree with taxon labels updated according to mapping"""
    logger.info(f"Creating relabeled tree: {output_tree}")
    
    # Read tree
    tree = dendropy.Tree.get(path=input_tree, schema=format, preserve_underscores=True)
    
    # Count relabeled nodes
    relabeled = 0
    
    # Update taxa labels
    for node in tree.leaf_node_iter():
        if node.taxon and node.taxon.label in mapping and mapping[node.taxon.label]:
            # This taxon has a mapping
            old_label = node.taxon.label
            node.taxon.label = mapping[old_label]
            relabeled += 1
    
    # Write new tree
    tree.write(path=output_tree, schema=format)
    logger.info(f"Created relabeled tree with {relabeled} renamed taxa")
    
    return relabeled

def main():
    parser = argparse.ArgumentParser(description='Match tree and alignment files by filtering or mapping taxa')
    parser.add_argument('--tree', '-t', required=True, help='Input tree file (Newick format)')
    parser.add_argument('--alignment', '-a', required=True, help='Input alignment file (FASTA format)')
    parser.add_argument('--output-dir', '-o', default=None, help='Output directory for filtered files')
    parser.add_argument('--tree-format', default='newick', help='Tree file format')
    parser.add_argument('--aln-format', default='fasta', help='Alignment file format')
    parser.add_argument('--filter-tree', action='store_true', help='Filter tree to include only taxa in alignment')
    parser.add_argument('--filter-alignment', action='store_true', help='Filter alignment to include only taxa in tree')
    parser.add_argument('--filter-both', action='store_true', help='Filter both to include only common taxa')
    parser.add_argument('--create-mapping', action='store_true', help='Create a mapping file between tree and alignment')
    parser.add_argument('--map-alignment', action='store_true', help='Create a new alignment with sequence IDs mapped to tree labels')
    parser.add_argument('--map-tree', action='store_true', help='Create a new tree with taxa labels mapped to sequence IDs')
    parser.add_argument('--mapping-file', help='Path to custom mapping file (CSV with tree_label,sequence_id columns)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Set up output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path('.')
    
    # Set filter flags if filter-both is selected
    if args.filter_both:
        args.filter_tree = True
        args.filter_alignment = True
    
    # Extract taxa from tree and alignment
    tree_taxa, tree_obj = get_tree_taxa(args.tree, args.tree_format)
    aln_taxa, aln_obj = get_alignment_taxa(args.alignment, args.aln_format)
    
    # Analyze mismatches
    mismatches = analyze_mismatches(tree_taxa, aln_taxa)
    
    # Create file base names
    tree_name = Path(args.tree).stem
    aln_name = Path(args.alignment).stem
    
    # Create mapping if requested
    mapping = {}
    if args.create_mapping:
        mapping_file = output_dir / f"{tree_name}_to_{aln_name}_mapping.csv"
        mapping = create_taxonomy_mapping(tree_taxa, aln_taxa, mapping_file, tree_obj, aln_obj)
    
    # Use provided mapping file if specified
    if args.mapping_file:
        logger.info(f"Using custom mapping file: {args.mapping_file}")
        try:
            with open(args.mapping_file, 'r') as f:
                reader = csv.DictReader(f)
                mapping = {row['tree_label']: row['sequence_id'] for row in reader if row['tree_label'] and row['sequence_id']}
            logger.info(f"Loaded {len(mapping)} mappings from file")
        except Exception as e:
            logger.error(f"Error reading mapping file: {e}")
            sys.exit(1)
    
    # Map alignment if requested
    if args.map_alignment and mapping:
        output_aln = output_dir / f"{aln_name}_mapped.{Path(args.alignment).suffix}"
        create_mapped_alignment(args.alignment, output_aln, mapping, args.aln_format)
    
    # Map tree if requested
    if args.map_tree and mapping:
        output_tree = output_dir / f"{tree_name}_mapped.{Path(args.tree).suffix}"
        relabel_tree_with_mapping(args.tree, output_tree, mapping, args.tree_format)
    
    # Filter alignment if requested (but not if filter-both is also requested)
    if args.filter_alignment and not args.filter_both:
        output_aln = output_dir / f"{aln_name}_filtered{Path(args.alignment).suffix}"
        filter_alignment(args.alignment, output_aln, tree_taxa, args.aln_format)
    
    # Filter tree if requested (but not if filter-both is also requested)
    if args.filter_tree and not args.filter_both:
        output_tree = output_dir / f"{tree_name}_filtered{Path(args.tree).suffix}"
        logger.info(f"Filtering tree to keep only the {len(mismatches['common'])} common taxa (removing {len(tree_taxa) - len(mismatches['common'])} taxa)")
        filter_tree(args.tree, output_tree, mismatches['common'], args.tree_format)
    
    # If both were filtered, create common filtered files
    if args.filter_both:
        output_aln = output_dir / f"{aln_name}_common{Path(args.alignment).suffix}"
        output_tree = output_dir / f"{tree_name}_common{Path(args.tree).suffix}"
        filter_alignment(args.alignment, output_aln, mismatches['common'], args.aln_format)
        filter_tree(args.tree, output_tree, mismatches['common'], args.tree_format)
    
    # Validate output files
    if args.filter_tree or args.filter_both:
        if args.filter_both:
            tree_output = output_dir / f"{tree_name}_common{Path(args.tree).suffix}"
        else:
            tree_output = output_dir / f"{tree_name}_filtered{Path(args.tree).suffix}"
            
        validate_newick_file(tree_output)
    
    # If no actions were requested
    if not (args.filter_tree or args.filter_alignment or args.filter_both or 
            args.create_mapping or args.map_alignment or args.map_tree):
        logger.info("No actions requested. Use --filter-tree, --filter-alignment, --filter-both, "
                   "--create-mapping, --map-alignment, or --map-tree")
        
    

if __name__ == '__main__':
    main()