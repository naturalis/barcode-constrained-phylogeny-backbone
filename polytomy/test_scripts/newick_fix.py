#!/usr/bin/env python
"""
Newick Tree Format Fix - Cleans up Newick tree files with format issues

This script fixes common formatting problems in Newick trees:
1. Redundant nested parentheses - ((node)) -> (node)
2. Empty nodes - () -> removed
3. Malformed branch length expressions
"""
import os
import re
import argparse
import logging
import dendropy

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("newick_fix")

def fix_newick_format(input_path, output_path=None):
    """
    Fix formatting issues in Newick tree files.
    
    Args:
        input_path (str): Path to input Newick tree file
        output_path (str, optional): Path to write fixed tree
        
    Returns:
        str: Path to the fixed Newick tree file
    """
    logger.info(f"Reading Newick tree from {input_path}")
    
    # First approach: Fix using dendropy
    try:
        tree = dendropy.Tree.get(path=input_path, schema="newick", preserve_underscores=True)
        
        # Clean up the tree
        logger.info("Cleaning tree structure")
        tree.suppress_unifurcations()
        
        # Generate output path if not provided
        if not output_path:
            base, ext = os.path.splitext(input_path)
            output_path = f"{base}.fixed{ext}"
        
        # Write the fixed tree
        logger.info(f"Writing fixed tree to {output_path}")
        tree.write(path=output_path, schema="newick", unquoted_underscores=True)
        
        return output_path
    
    except Exception as e:
        logger.warning(f"DendroPy approach failed: {e}. Trying direct string manipulation.")
    
    # Second approach: Direct string manipulation if dendropy fails
    with open(input_path, 'r') as f:
        newick = f.read().strip()
    
    # Keep count of initial issues
    original_length = len(newick)
    
    # Fix 1: Remove redundant parentheses
    while True:
        # Find instances of ((...)) with nothing in between the double parentheses
        modified = re.sub(r'\(\(([^()]*)\)\)', r'(\1)', newick)
        if modified == newick:
            break
        newick = modified

    # Fix 2: Remove empty nodes
    newick = re.sub(r'\(\)', '', newick)
    
    # Fix 3: Clean up any malformed branch lengths
    newick = re.sub(r':\.([0-9])', r':0.\1', newick)  # Fix things like :.5 -> :0.5
    newick = re.sub(r':,', r':0,', newick)  # Fix empty branch lengths
    
    # Generate output path if not provided
    if not output_path:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}.fixed{ext}"
    
    # Write the fixed tree
    with open(output_path, 'w') as f:
        f.write(newick)
    
    logger.info(f"Fixed tree written to {output_path}")
    logger.info(f"Reduced tree size from {original_length} to {len(newick)} characters")
    
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Fix formatting issues in Newick tree files')
    parser.add_argument('--input', '-i', required=True, help='Input Newick tree file')
    parser.add_argument('--output', '-o', help='Output fixed Newick tree file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    fix_newick_format(args.input, args.output)

if __name__ == '__main__':
    main()