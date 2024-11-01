#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.bed>"
    echo "Script will create a new file named 'sorted_coordinates.bed'"
    exit 1
fi

input_file=$1
output_file="sorted_coordinates.bed"

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found"
    exit 1
fi

# Create temporary file
temp_file=$(mktemp)

# Process the BED file
while IFS=$'\t' read -r line; do
    # Skip header lines if they exist (lines starting with #)
    if [[ $line == \#* ]]; then
        echo "$line" >> "$temp_file"
        continue
    fi
    
    # Skip empty lines
    if [ -z "$line" ]; then
        continue
    fi
    
    # Split the line into an array
    IFS=$'\t' read -ra fields <<< "$line"
    
    # Check if we have at least the required chr, start, end fields
    if [ ${#fields[@]} -lt 3 ]; then
        continue
    fi
    
    chr="${fields[0]}"
    start="${fields[1]}"
    end="${fields[2]}"
    
    # Extract numeric values for comparison
    start_num=${start//[!0-9]/}
    end_num=${end//[!0-9]/}
    
    # Compare coordinates and rebuild the line
    if [ "$start_num" -gt "$end_num" ]; then
        # Swap coordinates
        fields[1]="$end"
        fields[2]="$start"
    fi
    
    # Print all fields that exist, tab-separated
    (IFS=$'\t'; echo "${fields[*]}") >> "$temp_file"
    
done < "$input_file"

# Move temporary file to output file
mv "$temp_file" "$output_file"

echo "Processed file saved as '$output_file'"
