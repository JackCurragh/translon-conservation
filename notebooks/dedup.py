#!/usr/bin/env python3

import sys
from collections import defaultdict

def make_names_unique(bed_file, output_file):
    # Dictionary to keep track of name occurrences
    name_counts = defaultdict(int)
    
    # Read the input file and store lines
    lines = []
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:  # Ensure there's a name field
                name = fields[3]
                name_counts[name] += 1
                if name_counts[name] > 1:
                    # Append dup number to duplicate names
                    fields[3] = f"{name}_dup{name_counts[name]-1}"
            lines.append('\t'.join(fields))
    
    # Write the modified lines to output file
    with open(output_file, 'w') as f:
        for line in lines:
            f.write(line + '\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: python make_unique_names.py input.bed output.bed")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        make_names_unique(input_file, output_file)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()