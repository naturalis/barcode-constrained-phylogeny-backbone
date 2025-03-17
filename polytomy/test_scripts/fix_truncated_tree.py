#!/usr/bin/env python3
import dendropy
import re
import sys

def fix_tree(input_file, output_file):
    # First fix the truncated mrcaott labels
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Find and fix truncated mrcaott labels
    pattern = r'mrcaott(\d+)ott?(\d*)'
    
    # This replaces truncated labels like "mrcaott14848ot" with "mrcaott14848"
    fixed_content = re.sub(pattern, r'mrcaott\1', content)
    
    # Now process with DendroPy to ensure valid Newick format
    try:
        tree = dendropy.Tree.get_from_string(
            fixed_content,
            schema="newick", 
            preserve_underscores=True
        )
        
        # Write out clean tree
        tree.write(
            path=output_file,
            schema="newick",
            unquoted_underscores=True
        )
        
        print(f"Fixed tree written to {output_file}")
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fix_truncated_tree.py input.tre output.tre")
        sys.exit(1)
    
    success = fix_tree(sys.argv[1], sys.argv[2])
    sys.exit(0 if success else 1)