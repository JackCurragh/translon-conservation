#!/usr/bin/env python3
"""
Identify translons with RESTRICTED conservation:
- High feature conservation (PhyloP >1.5)
- No CDS overlap
- LOW flanking conservation (feature conservation is specific, not regional)
"""

import pandas as pd
import numpy as np

# Load corrected results
df = pd.read_csv('notebooks/results/phylop/translon_phylop_v2_with_cds_overlap_CORRECTED.tsv', sep='\t')

print("="*80)
print("RESTRICTED CONSERVATION ANALYSIS")
print("="*80)

# Start with the 1,565 conserved translons without CDS overlap
novel_conserved = df[
    (df['feature_mean'] > 1.5) &
    (~df['has_cds_overlap'])
]

print(f"\nStarting set: {len(novel_conserved):,} translons")
print(f"  - Feature PhyloP >1.5")
print(f"  - No CDS overlap")

# Calculate max flank conservation for each translon
novel_conserved = novel_conserved.copy()
novel_conserved['max_flank'] = novel_conserved[['upstream_mean', 'downstream_mean']].max(axis=1)
novel_conserved['mean_flank'] = novel_conserved[['upstream_mean', 'downstream_mean']].mean(axis=1)

print(f"\n{'FLANKING CONSERVATION DISTRIBUTION':^80}")
print("-"*80)
print(f"  Mean of max_flank:  {novel_conserved['max_flank'].mean():6.3f}")
print(f"  Median of max_flank: {novel_conserved['max_flank'].median():6.3f}")
print(f"  Mean of mean_flank:  {novel_conserved['mean_flank'].mean():6.3f}")

# Try different thresholds for "low" flanking conservation
print(f"\n{'RESTRICTED CONSERVATION BY FLANK THRESHOLD':^80}")
print("-"*80)

thresholds = [
    (0.0, "Negative/neutral flanks (max_flank ≤ 0.0)"),
    (0.5, "Weak flanks (max_flank ≤ 0.5)"),
    (1.0, "Below strong (max_flank ≤ 1.0)"),
    (1.5, "Below feature level (max_flank ≤ 1.5)")
]

for threshold, label in thresholds:
    restricted = novel_conserved[novel_conserved['max_flank'] <= threshold]

    print(f"\n{label}: {len(restricted):,} translons ({100*len(restricted)/len(novel_conserved):.1f}%)")
    if len(restricted) > 0:
        print(f"  Mean feature PhyloP:     {restricted['feature_mean'].mean():6.3f}")
        print(f"  Mean max_flank PhyloP:   {restricted['max_flank'].mean():6.3f}")
        print(f"  Mean specificity:        {restricted['conservation_specificity'].mean():6.3f}")

# Use conservation_specificity metric (already calculated: feature - max(flanks))
print(f"\n{'USING CONSERVATION SPECIFICITY METRIC':^80}")
print("-"*80)
print("(feature_mean - max(upstream_mean, downstream_mean))")

spec_thresholds = [
    (2.0, "Very high specificity (>2.0)"),
    (1.5, "High specificity (>1.5)"),
    (1.0, "Moderate specificity (>1.0)"),
    (0.5, "Weak specificity (>0.5)")
]

for spec_thresh, label in spec_thresholds:
    restricted = novel_conserved[novel_conserved['conservation_specificity'] > spec_thresh]

    print(f"\n{label}: {len(restricted):,} translons ({100*len(restricted)/len(novel_conserved):.1f}%)")
    if len(restricted) > 0:
        print(f"  Mean feature PhyloP:     {restricted['feature_mean'].mean():6.3f}")
        print(f"  Mean max_flank PhyloP:   {restricted['max_flank'].mean():6.3f}")
        print(f"  Mean upstream PhyloP:    {restricted['upstream_mean'].mean():6.3f}")
        print(f"  Mean downstream PhyloP:  {restricted['downstream_mean'].mean():6.3f}")

# Key analysis: Split into two groups
print(f"\n{'KEY SPLIT: REGIONAL vs RESTRICTED CONSERVATION':^80}")
print("="*80)

# Define "conserved flanks" as max_flank > 1.0
regional = novel_conserved[novel_conserved['max_flank'] > 1.0]
restricted = novel_conserved[novel_conserved['max_flank'] <= 1.0]

print(f"\nOf the {len(novel_conserved):,} conserved translons without CDS overlap:")
print()
print(f"REGIONAL CONSERVATION (max_flank >1.0):")
print(f"  Count: {len(regional):,} ({100*len(regional)/len(novel_conserved):.1f}%)")
print(f"  - Feature sits in broadly conserved genomic region")
print(f"  - Conservation likely NOT specific to translon")
print(f"  Mean feature PhyloP:    {regional['feature_mean'].mean():6.3f}")
print(f"  Mean max_flank PhyloP:  {regional['max_flank'].mean():6.3f}")

print(f"\nRESTRICTED CONSERVATION (max_flank ≤1.0):")
print(f"  Count: {len(restricted):,} ({100*len(restricted)/len(novel_conserved):.1f}%)")
print(f"  - Conservation appears SPECIFIC to translon")
print(f"  - Flanking regions show weak/neutral conservation")
print(f"  Mean feature PhyloP:    {restricted['feature_mean'].mean():6.3f}")
print(f"  Mean max_flank PhyloP:  {restricted['max_flank'].mean():6.3f}")

# Show top candidates with restricted conservation
print(f"\n{'TOP 20 RESTRICTED CONSERVATION CANDIDATES':^80}")
print("-"*80)

restricted_sorted = restricted.sort_values('conservation_specificity', ascending=False)

cols = ['translon_id', 'chrom', 'start', 'end', 'exonic_length', 'blockCount',
        'feature_mean', 'upstream_mean', 'downstream_mean', 'conservation_specificity']

if len(restricted_sorted) > 0:
    print()
    print(restricted_sorted[cols].head(20).to_string(index=False))

    # Save
    output = 'notebooks/results/phylop/restricted_conservation_candidates.tsv'
    restricted_sorted.to_csv(output, sep='\t', index=False)
    print(f"\n\nFull list saved to: {output}")

# Summary for the email thread
print(f"\n\n{'='*80}")
print("SUMMARY FOR MANUSCRIPT RESPONSE")
print("="*80)

total = len(df)
conserved = (df['feature_mean'] > 1.5).sum()
conserved_no_cds = len(novel_conserved)
truly_restricted = len(restricted)

print(f"\nTotal translons analyzed: {total:,}")
print(f"Strongly conserved (PhyloP >1.5): {conserved:,} ({100*conserved/total:.1f}%)")
print(f"\nOf the {conserved:,} conserved translons:")
print(f"  - Overlap CDS: {conserved - conserved_no_cds:,} ({100*(conserved - conserved_no_cds)/conserved:.1f}%)")
print(f"    → Conservation explained by overlapping protein-coding gene")
print(f"\n  - No CDS overlap: {conserved_no_cds:,} ({100*conserved_no_cds/conserved:.1f}%)")
print(f"    Of these:")
print(f"    • Regional conservation (flanks >1.0): {len(regional):,} ({100*len(regional)/conserved_no_cds:.1f}%)")
print(f"      → Sit in broadly conserved genomic regions")
print(f"    • RESTRICTED conservation (flanks ≤1.0): {truly_restricted:,} ({100*truly_restricted/conserved_no_cds:.1f}%)")
print(f"      → Conservation appears specific to translon itself")

print(f"\n\nKEY FINDING:")
print(f"  Of {total:,} translons, only {truly_restricted:,} ({100*truly_restricted/total:.1f}%) show")
print(f"  conservation that is:")
print(f"    1. Strong (PhyloP >1.5)")
print(f"    2. Independent of known CDS")
print(f"    3. Specific to the feature (not regional)")
