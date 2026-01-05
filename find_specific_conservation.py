#!/usr/bin/env python3
import pandas as pd
import numpy as np

df = pd.read_csv('notebooks/results/phylop/translon_phylop_comprehensive.tsv', sep='\t')

print("="*80)
print("SPECIFICALLY CONSERVED TRANSLONS ANALYSIS")
print("High feature conservation + Low flanking conservation")
print("="*80)

for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset].copy()

    # Calculate specific conservation metrics
    subset['flank_max'] = subset[['upstream_mean', 'downstream_mean']].max(axis=1)
    subset['conservation_specificity'] = subset['feature_mean'] - subset['flank_max']

    print(f"\n{'='*80}")
    print(f"Dataset: PhyloP {dataset}")
    print(f"{'='*80}")

    # Define criteria for "specific conservation"
    # Strong feature conservation (>1.5) + weak flanking (<0.5)
    strong_specific = subset[
        (subset['feature_mean'] > 1.5) &
        (subset['flank_max'] < 0.5)
    ].sort_values('conservation_specificity', ascending=False)

    # Moderate feature conservation (>1.0) + weak flanking (<0.3)
    moderate_specific = subset[
        (subset['feature_mean'] > 1.0) &
        (subset['flank_max'] < 0.3)
    ].sort_values('conservation_specificity', ascending=False)

    # Any positive feature with negative flanks
    positive_vs_negative = subset[
        (subset['feature_mean'] > 0.5) &
        (subset['flank_max'] < 0)
    ].sort_values('conservation_specificity', ascending=False)

    print(f"\nCRITERIA 1: Strong feature (>1.5) + Weak flanks (<0.5)")
    print(f"  Count: {len(strong_specific):,}")
    if len(strong_specific) > 0:
        print(f"  Mean specificity: {strong_specific['conservation_specificity'].mean():.3f}")
        print(f"  Top 10:")
        top = strong_specific.head(10)[['translon_id', 'chrom', 'length', 'feature_mean',
                                         'upstream_mean', 'downstream_mean', 'conservation_specificity']]
        print(top.to_string(index=False))

    print(f"\nCRITERIA 2: Moderate feature (>1.0) + Weak flanks (<0.3)")
    print(f"  Count: {len(moderate_specific):,}")
    if len(moderate_specific) > 0:
        print(f"  Mean specificity: {moderate_specific['conservation_specificity'].mean():.3f}")

    print(f"\nCRITERIA 3: Positive feature (>0.5) + Negative flanks (<0)")
    print(f"  Count: {len(positive_vs_negative):,}")
    if len(positive_vs_negative) > 0:
        print(f"  Mean specificity: {positive_vs_negative['conservation_specificity'].mean():.3f}")
        print(f"  Top 10:")
        top = positive_vs_negative.head(10)[['translon_id', 'chrom', 'length', 'feature_mean',
                                              'upstream_mean', 'downstream_mean', 'conservation_specificity']]
        print(top.to_string(index=False))

print(f"\n\n{'='*80}")
print("CONSENSUS CANDIDATES: Specifically conserved in ALL THREE datasets")
print(f"{'='*80}\n")

# Find translons specifically conserved across all datasets
specific_threshold = 0.5  # Feature must be at least 0.5 higher than flanks

specific_translons = {}
for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset].copy()
    subset['flank_max'] = subset[['upstream_mean', 'downstream_mean']].max(axis=1)
    subset['conservation_specificity'] = subset['feature_mean'] - subset['flank_max']

    specific = subset[subset['conservation_specificity'] > specific_threshold]
    specific_translons[dataset] = set(specific['translon_id'])

# Find intersection
all_three = specific_translons['30way'] & specific_translons['100way'] & specific_translons['470way']

print(f"Translons with specificity > 0.5 in all three datasets: {len(all_three)}")

if len(all_three) > 0:
    # Get details for consensus candidates
    consensus_df = df[df['translon_id'].isin(all_three)].pivot_table(
        index='translon_id',
        columns='phylop_dataset',
        values=['feature_mean', 'upstream_mean', 'downstream_mean'],
        aggfunc='first'
    )

    # Calculate average across datasets
    avg_feature = df[df['translon_id'].isin(all_three)].groupby('translon_id')['feature_mean'].mean()
    avg_flank = df[df['translon_id'].isin(all_three)].groupby('translon_id')[['upstream_mean', 'downstream_mean']].mean().max(axis=1)
    avg_specificity = avg_feature - avg_flank

    consensus_summary = pd.DataFrame({
        'avg_feature_conservation': avg_feature,
        'avg_flank_conservation': avg_flank,
        'avg_specificity': avg_specificity
    }).sort_values('avg_specificity', ascending=False)

    print("\nTop 20 consensus candidates:")
    print(consensus_summary.head(20).to_string())

    # Export full list
    output_file = 'notebooks/results/phylop/specifically_conserved_consensus.tsv'
    consensus_summary.to_csv(output_file, sep='\t')
    print(f"\nFull list saved to: {output_file}")

print("\n" + "="*80)
print("DISTRIBUTION ANALYSIS")
print("="*80)

for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset].copy()
    subset['flank_max'] = subset[['upstream_mean', 'downstream_mean']].max(axis=1)
    subset['conservation_specificity'] = subset['feature_mean'] - subset['flank_max']

    print(f"\n{dataset}:")
    print(f"  Specificity > 1.0:  {(subset['conservation_specificity'] > 1.0).sum():6,} ({100*(subset['conservation_specificity'] > 1.0).sum()/len(subset):5.1f}%)")
    print(f"  Specificity > 0.5:  {(subset['conservation_specificity'] > 0.5).sum():6,} ({100*(subset['conservation_specificity'] > 0.5).sum()/len(subset):5.1f}%)")
    print(f"  Specificity > 0.0:  {(subset['conservation_specificity'] > 0.0).sum():6,} ({100*(subset['conservation_specificity'] > 0.0).sum()/len(subset):5.1f}%)")
    print(f"  Specificity < 0.0:  {(subset['conservation_specificity'] < 0.0).sum():6,} ({100*(subset['conservation_specificity'] < 0.0).sum()/len(subset):5.1f}%)")
