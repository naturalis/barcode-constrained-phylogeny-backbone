#!/usr/bin/env python
"""
Advanced Newick Bracket Fix - Resolves nested parenthesis issues in Newick trees

This script focuses on fixing redundant brackets in Newick tree files,
which is a common error after pruning operations.
"""
import os
import re
import argparse
import logging
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger("newick_bracket_fix")

def find_problematic_brackets(newick_str, column_hint=None):
    """
    Finds problematic double brackets in the Newick string.
    
    Args:
        newick_str: The Newick string to check
        column_hint: A hint about where the problem might be
        
    Returns:
        List of positions where redundant brackets were found
    """
    problematic_positions = []
    
    # Stack to track brackets
    stack = []
    
    # Search area around hint if provided
    search_start = max(0, column_hint - 100) if column_hint else 0
    search_end = min(len(newick_str), column_hint + 100) if column_hint else len(newick_str)
    
    # Process the string
    for i in range(search_start, search_end):
        char = newick_str[i]
        
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                opening_pos = stack.pop()
                
                # Check if this closing bracket forms a redundant pair
                # by looking at the content between brackets
                content = newick_str[opening_pos+1:i]
                
                # Check if content contains only another bracket pair and no other data
                if re.match(r'^\([^()]*\)$', content):
                    problematic_positions.append((opening_pos, i))
    
    return problematic_positions

def fix_newick_brackets(input_path, output_path=None, column_hint=None):
    """
    Fix redundant bracket issues in Newick files.
    
    Args:
        input_path: Path to input Newick file
        output_path: Path to write the fixed file
        column_hint: Approximate column where error was reported
        
    Returns:
        Path to fixed file
    """
    logger.info(f"Reading Newick tree from {input_path}")
    
    with open(input_path, 'r') as f:
        newick_str = f.read().strip()
    
    original_length = len(newick_str)
    logger.info(f"Original Newick string length: {original_length}")
    
    # Check for issues around the hint if provided
    if column_hint:
        logger.info(f"Looking for issues around column {column_hint}")
        
        # Extract context around the hint
        start = max(0, column_hint - 50)
        end = min(len(newick_str), column_hint + 50)
        context = newick_str[start:end]
        
        logger.info(f"Context around problematic area: {context}")
    
    # First approach: Direct string replacements
    iterations = 0
    fixed = False
    
    while iterations < 10 and not fixed:  # Limit iterations to avoid infinite loop
        iterations += 1
        logger.info(f"Iteration {iterations} - fixing nested brackets")
        
        # Replace redundant nested brackets: ((x)) -> (x)
        previous = newick_str
        newick_str = re.sub(r'\(\(([^()]*)\)\)', r'(\1)', newick_str)
        
        # Replace more complex nested brackets with content: ((x,y)) -> (x,y)
        newick_str = re.sub(r'\(\((.*?)\)\)', lambda m: '(' + m.group(1) + ')', newick_str)
        
        # Handle empty brackets
        newick_str = re.sub(r'\(\)', '', newick_str)
        
        # Check if we made any changes
        if newick_str == previous:
            fixed = True
            logger.info("No more redundant brackets found")
    
    # Second approach: Character-by-character analysis if first approach didn't work
    if not fixed or (column_hint and iterations >= 10):
        logger.warning("First approach didn't fully resolve issues. Using character-by-character analysis.")
        
        problematic_positions = find_problematic_brackets(newick_str, column_hint)
        
        if problematic_positions:
            logger.info(f"Found {len(problematic_positions)} problematic bracket pairs")
            
            # Process from end to start to avoid index shifting
            for opening_pos, closing_pos in sorted(problematic_positions, reverse=True):
                content = newick_str[opening_pos+1:closing_pos]
                logger.info(f"Fixing redundant brackets at positions {opening_pos}-{closing_pos}: {newick_str[opening_pos-10:closing_pos+10]}")
                
                # Replace ((content)) with (content)
                newick_str = newick_str[:opening_pos] + content + newick_str[closing_pos+1:]
    
    # Generate output path if not provided
    if not output_path:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}.bracket_fixed{ext}"
    
    # Write the fixed tree
    with open(output_path, 'w') as f:
        f.write(newick_str)
    
    new_length = len(newick_str)
    logger.info(f"Fixed tree written to {output_path}")
    logger.info(f"Newick string length reduced from {original_length} to {new_length} characters")
    
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Fix redundant bracket issues in Newick tree files')
    parser.add_argument('--input', '-i', required=True, help='Input Newick tree file')
    parser.add_argument('--output', '-o', help='Output fixed Newick tree file')
    parser.add_argument('--column', '-c', type=int, help='Approximate column where error was reported')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    fix_newick_brackets(args.input, args.output, args.column)

if __name__ == '__main__':
    main()