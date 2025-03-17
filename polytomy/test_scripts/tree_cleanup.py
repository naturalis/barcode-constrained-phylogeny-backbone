#!/usr/bin/env python
"""
Redundant-Bracket Fix - Specialized tool for fixing double-bracket issues in Newick trees
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
logger = logging.getLogger("bracket_fix")

def fix_redundant_brackets(newick_string):
    """
    Fix redundant brackets in Newick string
    
    Args:
        newick_string: The Newick tree string
        
    Returns:
        Fixed Newick string
    """
    logger.info("Fixing redundant brackets...")
    modified = True
    iteration = 0
    
    while modified and iteration < 10:  # Limit iterations to prevent infinite loop
        iteration += 1
        original = newick_string
        
        # Fix most common pattern: ((content)) -> (content)
        newick_string = re.sub(r'\(\(([^()]+)\)\)', r'(\1)', newick_string)
        
        # Handle more complex patterns with recursive replacement
        previous_len = len(newick_string) + 1  # Ensure we enter the loop
        
        # Keep applying fixes until no more changes
        while len(newick_string) < previous_len:
            previous_len = len(newick_string)
            
            # Fix complex nested patterns with content that might include commas
            newick_string = re.sub(r'\(\(([^()]*(?:\([^()]*\)[^()]*)*)\)\)', r'(\1)', newick_string)
            
            # Fix empty subgroups: () -> ''
            newick_string = re.sub(r'\(\)', '', newick_string)
            
            # Fix unnecessary brackets around single taxa: (taxon) -> taxon
            newick_string = re.sub(r'\(([^,():;]+)\)', r'\1', newick_string)
        
        # Check if we made any changes this iteration
        modified = (newick_string != original)
        
        logger.info(f"Iteration {iteration}: {'Modified' if modified else 'No change'} " + 
                   f"({len(original)} → {len(newick_string)} chars)")
    
    return newick_string

def deep_fix_specific_problem(newick_string, problem_col=9255):
    """
    Fix a specific problem at a given column by analyzing the local structure
    
    Args:
        newick_string: Newick string
        problem_col: Approximate column position of the problem
        
    Returns:
        Fixed Newick string
    """
    logger.info(f"Targeting specific problem around column {problem_col}")
    
    # Safety check - ensure the column is within bounds
    if problem_col >= len(newick_string):
        logger.warning(f"Problem column {problem_col} is beyond string length {len(newick_string)}")
        return newick_string
    
    # Extract context around the problem
    start = max(0, problem_col - 100)
    end = min(len(newick_string), problem_col + 100)
    context = newick_string[start:end]
    
    logger.info(f"Context: ...{context}...")
    
    # Find the problematic area - look for double brackets
    # Focus on the closing bracket area mentioned in error
    problem_area_search = re.search(r'\)\)', context)
    
    if not problem_area_search:
        logger.warning("Couldn't identify double closing brackets in context")
        return newick_string
    
    # Position within context
    problem_pos = problem_area_search.start() + 1  # +1 to point to second ')'
    
    # Find matching opening brackets - trace backwards counting brackets
    open_count = 2  # We need to find two opening brackets
    nested_level = 0
    
    # Start from the problem position and go backwards
    pos = problem_pos
    while pos >= 0 and open_count > 0:
        if context[pos] == ')':
            nested_level += 1
        elif context[pos] == '(':
            if nested_level == 0:
                open_count -= 1
            else:
                nested_level -= 1
        pos -= 1
    
    if open_count > 0:
        logger.warning("Couldn't find matching opening brackets")
        return newick_string
    
    # Extract the area with redundant brackets
    opening_pos = pos + 1
    closing_pos = problem_pos + 1
    
    # The area to fix
    area_to_fix = context[opening_pos:closing_pos+1]
    logger.info(f"Area to fix: {area_to_fix}")
    
    # Fix the specific problem - remove outer brackets
    if area_to_fix.startswith('((') and area_to_fix.endswith('))'):
        fixed_area = '(' + area_to_fix[2:-1] + ')'
        
        # Replace in the full string
        abs_opening = start + opening_pos
        abs_closing = start + closing_pos
        fixed_newick = newick_string[:abs_opening] + fixed_area + newick_string[abs_closing+1:]
        
        logger.info(f"Fixed specific problem: {area_to_fix} → {fixed_area}")
        return fixed_newick
    
    logger.warning(f"Area doesn't match the pattern we can fix: {area_to_fix}")
    return newick_string

def repair_tree(input_path, output_path=None, problem_col=None, test_iqtree=False, alignment_path=None):
    """
    Repair a Newick tree file
    
    Args:
        input_path: Path to input tree file
        output_path: Path to output fixed tree file
        problem_col: Column position where problem was reported
        test_iqtree: Whether to test with IQTree
        alignment_path: Path to alignment file for testing with IQTree
    
    Returns:
        Path to fixed tree file
    """
    logger.info(f"Loading tree from {input_path}")
    
    # Read the tree file
    with open(input_path, 'r') as f:
        newick = f.read().strip()
    
    original_length = len(newick)
    logger.info(f"Original tree length: {original_length} characters")
    
    # Apply general bracket fixes
    fixed_newick = fix_redundant_brackets(newick)
    
    # If a problem column is provided, apply targeted fix
    if problem_col is not None:
        fixed_newick = deep_fix_specific_problem(fixed_newick, problem_col)
    
    # Generate output path if not provided
    if output_path is None:
        base, ext = os.path.splitext(input_path)
        output_path = f"{base}.bracket_fixed{ext}"
    
    # Write the fixed tree
    with open(output_path, 'w') as f:
        f.write(fixed_newick)
    
    logger.info(f"Fixed tree written to {output_path} ({len(fixed_newick)} characters)")
    
    # Test with IQTree if requested
    if test_iqtree and alignment_path:
        logger.info(f"Testing with IQTree...")
        
        cmd = ['iqtree2', '-s', alignment_path, '-te', output_path, '-mredo', '-m', 'GTR+G', '-n', '0']
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                logger.info("IQTree validation passed!")
            else:
                logger.error(f"IQTree reported issues: {result.stderr}")
                
                # Extract position from error message
                error_match = re.search(r'ending at \(line \d+ column (\d+)\)', result.stderr)
                if error_match:
                    new_problem_col = int(error_match.group(1))
                    logger.info(f"New problem detected at column {new_problem_col}")
                    
                    # Recursively attempt one more fix
                    if new_problem_col != problem_col:
                        logger.info("Attempting to fix new problem...")
                        return repair_tree(output_path, None, new_problem_col, test_iqtree, alignment_path)
        except Exception as e:
            logger.error(f"Error during IQTree validation: {e}")
    
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Fix redundant brackets in Newick trees')
    parser.add_argument('--input', '-i', required=True, help='Input Newick tree file')
    parser.add_argument('--output', '-o', help='Output fixed tree file')
    parser.add_argument('--column', '-c', type=int, help='Problem column position from error message')
    parser.add_argument('--test', '-t', action='store_true', help='Test with IQTree')
    parser.add_argument('--alignment', '-a', help='Alignment file for IQTree testing')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    repair_tree(args.input, args.output, args.column, args.test, args.alignment)

if __name__ == '__main__':
    main()