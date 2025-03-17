#!/usr/bin/env python
"""
Empty Node Fix - Repairs Newick trees with empty node issues
"""
import os
import re
import argparse
import logging
import subprocess

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("empty_node_fix")

def fix_empty_nodes(newick_string):
    """Fix empty nodes in Newick format that appear as ,(, pattern"""
    
    # Track replacements
    replacements = 0
    
    # Fix empty nodes: ,(, → ,( 
    # This pattern shows up when a tip is pruned but its structure remains
    modified = newick_string
    while True:
        previous = modified
        modified = re.sub(r',\(,', r',(', modified)
        if modified == previous:
            break
        replacements += 1
    
    # Also fix pattern: ),) → ))
    while True:
        previous = modified
        modified = re.sub(r'\),\)', r'))', modified)
        if modified == previous:
            break
        replacements += 1
    
    # Fix other empty node patterns
    while True:
        previous = modified
        # Fix empty node: (,( → (
        modified = re.sub(r'\(,\(', r'(', modified)
        if modified == previous:
            break
        replacements += 1
    
    # Report how many replacements were made
    logger.info(f"Fixed {replacements} empty nodes")
    
    return modified

def fix_redundant_brackets(newick_string):
    """Fix redundant brackets in Newick string"""
    
    # Track replacements
    replacements = 0
    
    # Apply fixes repeatedly
    while True:
        previous = newick_string
        
        # Fix common pattern: ((content)) → (content)
        newick_string = re.sub(r'\(\(([^()]+)\)\)', r'(\1)', newick_string)
        
        # Fix redundant brackets with more complex content
        newick_string = re.sub(r'\(\(([^()]*(?:\([^()]*\)[^()]*)*)\)\)', r'(\1)', newick_string)
        
        # Check if we made changes
        if newick_string == previous:
            break
            
        replacements += 1
    
    # Report how many replacements were made
    logger.info(f"Fixed {replacements} redundant bracket patterns")
    
    return newick_string

def repair_newick_file(input_path, output_path=None, validate=False, alignment_path=None):
    """Repair Newick file with various fixes"""
    
    logger.info(f"Repairing Newick file: {input_path}")
    
    # Read the file
    with open(input_path, 'r') as f:
        content = f.read()
    
    # Remove any comment lines
    if content.startswith('//'):
        logger.info("Removing comment lines")
        content = re.sub(r'^//.*\n', '', content)
    
    # Save original length
    original_length = len(content)
    logger.info(f"Original content length: {original_length} characters")
    
    # Fix empty nodes
    content = fix_empty_nodes(content)
    
    # Fix redundant brackets
    content = fix_redundant_brackets(content)
    
    # Ensure the tree is valid
    # Must start with ( and end with ;
    if not content.startswith('('):
        logger.info("Adding opening bracket")
        content = '(' + content
    
    if not content.endswith(';'):
        logger.info("Adding closing semicolon")
        content = content + ';'
    
    # Set output path
    if not output_path:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}.fixed{ext}"
    
    # Write fixed content
    with open(output_path, 'w') as f:
        f.write(content)
    
    logger.info(f"Fixed tree written to {output_path} ({len(content)} characters)")
    logger.info(f"Reduced length by {original_length - len(content)} characters")
    
    # Validate with IQTree if requested
    if validate and alignment_path:
        logger.info("Validating with IQTree...")
        
        try:
            cmd = ['iqtree2', '-s', alignment_path, '-te', output_path, '-m', 'GTR+G', '-n', '0', '-quiet']
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info("✓ Tree validated successfully with IQTree")
                return output_path
            else:
                logger.warning(f"✗ IQTree reports issues: {result.stderr}")
                
                # Check for specific error patterns
                error_match = re.search(r'ending at \(line \d+ column (\d+)\)', result.stderr)
                if error_match:
                    problem_col = int(error_match.group(1))
                    logger.info(f"Problem detected at column {problem_col}")
                    
                    # Extract context around problem
                    start = max(0, problem_col - 20)
                    end = min(len(content), problem_col + 20)
                    context = content[start:end]
                    logger.info(f"Context: ...{context}...")
        except Exception as e:
            logger.error(f"Error during validation: {e}")
    
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Fix empty nodes in Newick tree files')
    parser.add_argument('--input', '-i', required=True, help='Input Newick tree file')
    parser.add_argument('--output', '-o', help='Output fixed tree file')
    parser.add_argument('--validate', '-v', action='store_true', help='Validate with IQTree')
    parser.add_argument('--alignment', '-a', help='Alignment file for validation')
    
    args = parser.parse_args()
    repair_newick_file(args.input, args.output, args.validate, args.alignment)

if __name__ == '__main__':
    main()