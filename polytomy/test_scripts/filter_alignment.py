#!/usr/bin/env python3
import csv
import os
import sys
from Bio import SeqIO

# Read exemplars table
exemplars = set()  # Use a set to store all exemplar IDs
with open('./tests/my_data/exemplars_table.txt', 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)  # Skip header
    
    # Assuming format: taxon_name, exemplar1_id, exemplar2_id
    for row in reader:
        if len(row) > 1:
            exemplar1 = row[1]
            if exemplar1:  # Add first exemplar if not empty
                exemplars.add(exemplar1)
        
        if len(row) > 2:
            exemplar2 = row[2]
            if exemplar2:  # Add second exemplar if not empty
                exemplars.add(exemplar2)

print(f"Found {len(exemplars)} total exemplar sequences in pairs file")

# Create filtered alignment with both exemplars from each taxon
records = []
for record in SeqIO.parse('./tests/my_data/raxml-ready.fa', 'fasta'):
    if record.id in exemplars:
        records.append(record)

# Write filtered alignment
output_file = './tests/my_data/exemplars_all.fa'
SeqIO.write(records, output_file, 'fasta')
print(f"Created alignment with {len(records)} exemplars in {output_file}")
print(f"Found {len(exemplars) - len(records)} exemplars from pairs file that weren't in the alignment")