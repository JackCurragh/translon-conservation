#!/usr/bin/env python3
"""
Add CDS overlap analysis to PhyloP results.
This will help distinguish truly novel conserved elements from those
that overlap existing protein-coding genes.
"""

import pandas as pd
import pyranges as pr
import subprocess
from pathlib import Path

# Load PhyloP results
phylop_results = pd.read_csv('notebooks/results/phylop/translon_phylop_v2_comprehensive.tsv', sep='\t')

print("="*80)
print("CDS OVERLAP ANALYSIS")
print("="*80)

# First, check if we have GENCODE CDS annotations
gencode_gtf = Path('data/gencode.v46.annotation.gtf.gz')

if not gencode_gtf.exists():
    print("\nGENCODE annotation not found. Downloading...")
    print("This will take a few minutes...")
    subprocess.run([
        'wget',
        '-O', str(gencode_gtf),
        'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz'
    ], check=True)

print(f"\nReading GENCODE annotations from {gencode_gtf}...")

# Load GENCODE and extract CDS features
import gzip

print("Extracting CDS features...")
cds_records = []

with gzip.open(gencode_gtf, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        if fields[2] == 'CDS':
            chrom = fields[0]
            start = int(fields[3]) - 1  # GTF is 1-based, convert to 0-based
            end = int(fields[4])
            strand = fields[6]

            # Parse attributes
            attrs = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr:
                    key_val = attr.split(' ', 1)
                    if len(key_val) == 2:
                        key, val = key_val
                        attrs[key] = val.strip('"')

            cds_records.append({
                'Chromosome': chrom,
                'Start': start,
                'End': end,
                'Strand': strand,
                'gene_id': attrs.get('gene_id', ''),
                'gene_name': attrs.get('gene_name', ''),
                'gene_type': attrs.get('gene_type', '')
            })

print(f"Found {len(cds_records):,} CDS features")

# Create PyRanges objects
cds_pr = pr.PyRanges(pd.DataFrame(cds_records))

# Create PyRanges for translons (using 470way results for specificity)
translon_470 = phylop_results[phylop_results['phylop_dataset'] == '470way'].copy()

translon_records = []
for _, row in translon_470.iterrows():
    translon_records.append({
        'Chromosome': row['chrom'],
        'Start': row['start'],
        'End': row['end'],
        'Strand': row['strand'],
        'translon_id': row['name'] if 'name' in row else row['translon_id'],
        'feature_mean': row['feature_mean'],
        'conservation_specificity': row['conservation_specificity'],
        'blockCount': row['blockCount']
    })

translon_pr = pr.PyRanges(pd.DataFrame(translon_records))

print("\nFinding overlaps between translons and CDS...")

# Find overlaps (any overlap, regardless of strand)
overlaps = translon_pr.join(cds_pr, how='left', suffix='_cds')

# Process results
overlap_results = []
for idx, row in translon_470.iterrows():
    translon_id = row['name'] if 'name' in row else row['translon_id']

    # Find if this translon overlaps any CDS
    translon_overlaps = overlaps.df[
        (overlaps.df['Chromosome'] == row['chrom']) &
        (overlaps.df['Start'] == row['start']) &
        (overlaps.df['End'] == row['end'])
    ]

    if len(translon_overlaps) > 0 and 'Start_cds' in translon_overlaps.columns:
        # Has overlap
        has_cds_overlap = translon_overlaps['Start_cds'].notna().any()

        if has_cds_overlap:
            # Calculate overlap details
            overlapping_genes = translon_overlaps['gene_name'].dropna().unique()
            overlapping_gene_types = translon_overlaps['gene_type'].dropna().unique()

            # Check if same strand
            same_strand = (translon_overlaps['Strand'] == translon_overlaps['Strand_cds']).any()

            overlap_results.append({
                'translon_id': translon_id,
                'has_cds_overlap': True,
                'same_strand': same_strand,
                'n_overlapping_genes': len(overlapping_genes),
                'overlapping_genes': ','.join(overlapping_genes[:5]),  # First 5
                'overlapping_gene_types': ','.join(overlapping_gene_types[:3])  # First 3
            })
        else:
            overlap_results.append({
                'translon_id': translon_id,
                'has_cds_overlap': False,
                'same_strand': False,
                'n_overlapping_genes': 0,
                'overlapping_genes': '',
                'overlapping_gene_types': ''
            })
    else:
        overlap_results.append({
            'translon_id': translon_id,
            'has_cds_overlap': False,
            'same_strand': False,
            'n_overlapping_genes': 0,
            'overlapping_genes': '',
            'overlapping_gene_types': ''
        })

overlap_df = pd.DataFrame(overlap_results)

# Merge with PhyloP results
translon_470_with_overlap = translon_470.merge(overlap_df, on='translon_id', how='left')

print("\n" + "="*80)
print("OVERALL CDS OVERLAP STATISTICS")
print("="*80)

total = len(translon_470_with_overlap)
with_overlap = translon_470_with_overlap['has_cds_overlap'].sum()
same_strand_overlap = (translon_470_with_overlap['has_cds_overlap'] & translon_470_with_overlap['same_strand']).sum()

print(f"\nTotal translons (470way): {total:,}")
print(f"Overlap any CDS: {with_overlap:,} ({100*with_overlap/total:.1f}%)")
print(f"  Same strand:   {same_strand_overlap:,} ({100*same_strand_overlap/total:.1f}%)")
print(f"  Opposite strand: {with_overlap - same_strand_overlap:,} ({100*(with_overlap - same_strand_overlap)/total:.1f}%)")
print(f"No CDS overlap:  {total - with_overlap:,} ({100*(total - with_overlap)/total:.1f}%)")

print("\n" + "="*80)
print("CDS OVERLAP BY CONSERVATION SPECIFICITY")
print("="*80)

for threshold, label in [(2.0, '>2.0 (highest)'), (1.0, '>1.0 (high)'), (0.5, '>0.5 (moderate)'), (0.0, '>0.0 (any)')]:
    subset = translon_470_with_overlap[translon_470_with_overlap['conservation_specificity'] > threshold]

    if len(subset) > 0:
        with_cds = subset['has_cds_overlap'].sum()
        same_strand = (subset['has_cds_overlap'] & subset['same_strand']).sum()

        print(f"\nSpecificity {label}: {len(subset):,} translons")
        print(f"  With CDS overlap:    {with_cds:,} ({100*with_cds/len(subset):5.1f}%)")
        print(f"    Same strand:       {same_strand:,} ({100*same_strand/len(subset):5.1f}%)")
        print(f"    Opposite strand:   {with_cds - same_strand:,} ({100*(with_cds - same_strand)/len(subset):5.1f}%)")
        print(f"  WITHOUT CDS overlap: {len(subset) - with_cds:,} ({100*(len(subset) - with_cds)/len(subset):5.1f}%)")

print("\n" + "="*80)
print("TRULY NOVEL CONSERVED ELEMENTS (No CDS overlap)")
print("="*80)

# High confidence candidates without CDS overlap
novel_high_conf = translon_470_with_overlap[
    (translon_470_with_overlap['conservation_specificity'] > 1.0) &
    (~translon_470_with_overlap['has_cds_overlap'])
].sort_values('conservation_specificity', ascending=False)

print(f"\nTranslons with specificity >1.0 AND no CDS overlap: {len(novel_high_conf):,}")

if len(novel_high_conf) > 0:
    print("\nTop 20 truly novel candidates:")
    cols = ['translon_id', 'chrom', 'start', 'end', 'blockCount', 'exonic_length',
            'feature_mean', 'conservation_specificity']
    print(novel_high_conf[cols].head(20).to_string(index=False))

    # Save to file
    output_file = 'notebooks/results/phylop/novel_conserved_no_cds_overlap.tsv'
    novel_high_conf.to_csv(output_file, sep='\t', index=False)
    print(f"\nFull list saved to: {output_file}")

print("\n" + "="*80)
print("CONSENSUS CANDIDATES (n=8) - CDS OVERLAP CHECK")
print("="*80)

consensus_ids = ['c4norep39', 'c16norep108', 'c3norep228', 'c1riboseqorf33',
                 'c16riboseqorf85', 'c10norep84', 'c2riboseqorf156', 'c1norep63']

consensus_overlap = translon_470_with_overlap[translon_470_with_overlap['translon_id'].isin(consensus_ids)]

print(f"\nConsensus candidates:")
cols = ['translon_id', 'feature_mean', 'conservation_specificity', 'has_cds_overlap',
        'same_strand', 'overlapping_genes']
print(consensus_overlap[cols].to_string(index=False))

novel_consensus = consensus_overlap[~consensus_overlap['has_cds_overlap']]
print(f"\nConsensus candidates WITHOUT CDS overlap: {len(novel_consensus)}")
if len(novel_consensus) > 0:
    print(novel_consensus['translon_id'].tolist())

# Save full results with CDS overlap annotations
output_full = 'notebooks/results/phylop/translon_phylop_v2_with_cds_overlap.tsv'
translon_470_with_overlap.to_csv(output_full, sep='\t', index=False)
print(f"\n\nFull results with CDS overlap saved to: {output_full}")
