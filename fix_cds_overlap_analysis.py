#!/usr/bin/env python3
"""
Fix the CDS overlap analysis - PyRanges uses -1 to mean NO overlap
"""

import pandas as pd
import numpy as np

# Load results
df = pd.read_csv('notebooks/results/phylop/translon_phylop_v2_with_cds_overlap.tsv', sep='\t')

print("="*80)
print("CORRECTED CDS OVERLAP ANALYSIS")
print("="*80)

# Fix the has_cds_overlap column
# -1 in overlapping_genes means NO overlap (PyRanges sentinel value)
df['has_cds_overlap_corrected'] = df['overlapping_genes'] != '-1'
df['same_strand_corrected'] = df['has_cds_overlap_corrected'] & df['same_strand']

print(f"\nTotal translons: {len(df):,}")

# Overall CDS overlap
with_cds = df['has_cds_overlap_corrected'].sum()
same_strand = df['same_strand_corrected'].sum()
no_cds = (~df['has_cds_overlap_corrected']).sum()

print(f"\n{'OVERALL CDS OVERLAP':^80}")
print("-"*80)
print(f"  Overlap CDS (any):       {with_cds:6,} ({100*with_cds/len(df):5.1f}%)")
print(f"    Same strand:           {same_strand:6,} ({100*same_strand/len(df):5.1f}%)")
print(f"    Antisense:             {with_cds - same_strand:6,} ({100*(with_cds - same_strand)/len(df):5.1f}%)")
print(f"  NO CDS overlap:          {no_cds:6,} ({100*no_cds/len(df):5.1f}%)")

# Conservation by CDS overlap status
print(f"\n{'PHYLOP BY CDS OVERLAP STATUS':^80}")
print("-"*80)

cds_overlap = df[df['has_cds_overlap_corrected']]
no_cds_overlap = df[~df['has_cds_overlap_corrected']]

print(f"\nWith CDS overlap (n={len(cds_overlap):,}):")
print(f"  Mean PhyloP:       {cds_overlap['feature_mean'].mean():6.3f}")
print(f"  Median PhyloP:     {cds_overlap['feature_median'].median():6.3f}")
print(f"  PhyloP >1.5:       {(cds_overlap['feature_mean'] > 1.5).sum():6,} ({100*(cds_overlap['feature_mean'] > 1.5).sum()/len(cds_overlap):5.1f}%)")
print(f"  PhyloP >0.5:       {(cds_overlap['feature_mean'] > 0.5).sum():6,} ({100*(cds_overlap['feature_mean'] > 0.5).sum()/len(cds_overlap):5.1f}%)")

print(f"\nWithout CDS overlap (n={len(no_cds_overlap):,}):")
print(f"  Mean PhyloP:       {no_cds_overlap['feature_mean'].mean():6.3f}")
print(f"  Median PhyloP:     {no_cds_overlap['feature_median'].median():6.3f}")
print(f"  PhyloP >1.5:       {(no_cds_overlap['feature_mean'] > 1.5).sum():6,} ({100*(no_cds_overlap['feature_mean'] > 1.5).sum()/len(no_cds_overlap):5.1f}%)")
print(f"  PhyloP >0.5:       {(no_cds_overlap['feature_mean'] > 0.5).sum():6,} ({100*(no_cds_overlap['feature_mean'] > 0.5).sum()/len(no_cds_overlap):5.1f}%)")

# CDS overlap by PhyloP threshold
print(f"\n{'CDS OVERLAP BY CONSERVATION LEVEL':^80}")
print("-"*80)

thresholds = [
    (2.0, "Very strong (>2.0)"),
    (1.5, "Strong (>1.5)"),
    (1.0, "Moderate+ (>1.0)"),
    (0.5, "Weak+ (>0.5)"),
    (0.0, "Any positive (>0.0)")
]

for threshold, label in thresholds:
    subset = df[df['feature_mean'] > threshold]

    if len(subset) > 0:
        with_cds_sub = subset['has_cds_overlap_corrected'].sum()
        without_cds_sub = (~subset['has_cds_overlap_corrected']).sum()

        print(f"\n{label}: {len(subset):,} translons")
        print(f"  With CDS overlap:      {with_cds_sub:6,} ({100*with_cds_sub/len(subset):5.1f}%)")
        print(f"  WITHOUT CDS overlap:   {without_cds_sub:6,} ({100*without_cds_sub/len(subset):5.1f}%)")

# High confidence novel candidates
print(f"\n{'NOVEL CONSERVED CANDIDATES (PhyloP >1.5, NO CDS overlap)':^80}")
print("-"*80)

novel_strong = df[
    (df['feature_mean'] > 1.5) &
    (~df['has_cds_overlap_corrected'])
].sort_values('feature_mean', ascending=False)

print(f"\nCount: {len(novel_strong):,}")

if len(novel_strong) > 0:
    print("\nTop 20 truly novel conserved translons:")
    cols = ['translon_id', 'chrom', 'start', 'end', 'exonic_length', 'blockCount',
            'feature_mean', 'upstream_mean', 'downstream_mean']
    print(novel_strong[cols].head(20).to_string(index=False))

    # Save
    output = 'notebooks/results/phylop/truly_novel_conserved_translons.tsv'
    novel_strong.to_csv(output, sep='\t', index=False)
    print(f"\nFull list saved to: {output}")

# Summary stats
print(f"\n{'SUMMARY STATISTICS':^80}")
print("="*80)

print(f"\nOverall (all 7,263 translons):")
print(f"  Mean PhyloP: {df['feature_mean'].mean():.3f}")
print(f"  CDS overlap: {100*with_cds/len(df):.1f}%")

print(f"\nConserved (PhyloP >1.5, n={(df['feature_mean'] > 1.5).sum():,}):")
conserved = df[df['feature_mean'] > 1.5]
print(f"  CDS overlap: {100*conserved['has_cds_overlap_corrected'].sum()/len(conserved):.1f}%")
print(f"  Novel (no CDS): {100*(~conserved['has_cds_overlap_corrected']).sum()/len(conserved):.1f}%")

print(f"\nKey Finding:")
print(f"  - {100*no_cds/len(df):.1f}% of translons do NOT overlap CDS")
print(f"  - Of conserved translons (>1.5), {100*len(novel_strong)/len(conserved):.1f}% are truly novel")
print(f"  - This represents {len(novel_strong):,} high-confidence novel candidates")

# Save corrected results
df_corrected = df.copy()
df_corrected['has_cds_overlap'] = df_corrected['has_cds_overlap_corrected']
df_corrected['same_strand'] = df_corrected['same_strand_corrected']
df_corrected = df_corrected.drop(columns=['has_cds_overlap_corrected', 'same_strand_corrected'])

output_corrected = 'notebooks/results/phylop/translon_phylop_v2_with_cds_overlap_CORRECTED.tsv'
df_corrected.to_csv(output_corrected, sep='\t', index=False)
print(f"\nCorrected full results saved to: {output_corrected}")
