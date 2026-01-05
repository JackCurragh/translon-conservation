#!/usr/bin/env python3
import pandas as pd
import numpy as np

df = pd.read_csv('notebooks/results/phylop/translon_phylop_comprehensive.tsv', sep='\t')

print("="*80)
print("PHYLOP CONSERVATION ANALYSIS - COMPREHENSIVE REPORT")
print("="*80)

for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset]

    print(f"\n{'='*80}")
    print(f"Dataset: PhyloP {dataset} Vertebrate Alignment")
    print(f"{'='*80}")
    print(f"Translons analyzed: {len(subset):,}")

    print(f"\n{'OVERALL CONSERVATION':^80}")
    print("-"*80)
    print(f"  Mean PhyloP score:        {subset['feature_mean'].mean():7.3f} ± {subset['feature_mean'].std():.3f}")
    print(f"  Median PhyloP score:      {subset['feature_median'].median():7.3f}")
    print(f"  Range:                    [{subset['feature_mean'].min():.3f}, {subset['feature_mean'].max():.3f}]")
    print(f"  Fraction positive (>0):   {subset['feature_frac_positive'].mean():7.1%}")
    print(f"  Fraction conserved (>1.5): {subset['feature_frac_conserved_strong'].mean():7.1%}")

    strongly_conserved = (subset['feature_mean'] > 1.5).sum()
    conserved = ((subset['feature_mean'] > 0.5) & (subset['feature_mean'] <= 1.5)).sum()
    neutral = ((subset['feature_mean'] >= -0.5) & (subset['feature_mean'] <= 0.5)).sum()
    depleted = (subset['feature_mean'] < -0.5).sum()

    print(f"\n{'CONSERVATION DISTRIBUTION':^80}")
    print("-"*80)
    print(f"  {'Category':<30} {'Count':>10} {'Percentage':>15} {'PhyloP Range':>20}")
    print(f"  {'-'*30} {'-'*10} {'-'*15} {'-'*20}")
    print(f"  {'Strongly conserved':<30} {strongly_conserved:10,} {100*strongly_conserved/len(subset):14.1f}% {'> 1.5':>20}")
    print(f"  {'Conserved':<30} {conserved:10,} {100*conserved/len(subset):14.1f}% {'0.5 to 1.5':>20}")
    print(f"  {'Neutral':<30} {neutral:10,} {100*neutral/len(subset):14.1f}% {'-0.5 to 0.5':>20}")
    print(f"  {'Depleted':<30} {depleted:10,} {100*depleted/len(subset):14.1f}% {'< -0.5':>20}")

    codon_df = subset[subset['n_codons'] > 0]
    if len(codon_df) > 0:
        print(f"\n{'CODON POSITION ANALYSIS':^80}")
        print("-"*80)
        print(f"  Translons with codon data: {len(codon_df):,}")
        print(f"  Position 1 (first base):   {codon_df['codon_pos1_mean'].mean():7.3f} (std: {codon_df['codon_pos1_std'].mean():.3f})")
        print(f"  Position 2 (second base):  {codon_df['codon_pos2_mean'].mean():7.3f} (std: {codon_df['codon_pos2_std'].mean():.3f})")
        print(f"  Position 3 (wobble):       {codon_df['codon_pos3_mean'].mean():7.3f} (std: {codon_df['codon_pos3_std'].mean():.3f})")
        print(f"  Wobble depletion:          {codon_df['wobble_depletion_mean'].mean():7.3f}")
        print(f"\n  Conservation fractions:")
        print(f"    Position 1 (>1.5):       {codon_df['codon_pos1_frac_conserved'].mean():7.1%}")
        print(f"    Position 2 (>1.5):       {codon_df['codon_pos2_frac_conserved'].mean():7.1%}")
        print(f"    Position 3 (>1.5):       {codon_df['codon_pos3_frac_conserved'].mean():7.1%}")

    flanking_df = subset[subset['flanking_mean'].notna()]
    if len(flanking_df) > 0:
        print(f"\n{'FLANKING REGION CONTEXT':^80}")
        print("-"*80)
        print(f"  Feature vs flanking diff:       {flanking_df['feature_vs_flanking_diff'].mean():7.3f}")
        print(f"  Upstream mean:                  {flanking_df['upstream_mean'].mean():7.3f}")
        print(f"  Downstream mean:                {flanking_df['downstream_mean'].mean():7.3f}")
        print(f"  Specifically conserved (>0.5):  {flanking_df['specifically_conserved'].sum():6,} ({100*flanking_df['specifically_conserved'].sum()/len(flanking_df):5.1f}%)")
        print(f"  More conserved than flanks:     {flanking_df['more_conserved_than_flanks'].sum():6,} ({100*flanking_df['more_conserved_than_flanks'].sum()/len(flanking_df):5.1f}%)")

print(f"\n\n{'='*80}")
print("CROSS-DATASET COMPARISON")
print(f"{'='*80}\n")

comparison_metrics = [
    ('Mean conservation', 'feature_mean'),
    ('Median conservation', 'feature_median'),
    ('Fraction positive', 'feature_frac_positive'),
    ('Fraction strong (>1.5)', 'feature_frac_conserved_strong'),
    ('Wobble depletion', 'wobble_depletion_mean')
]

print(f"{'Metric':<30} {'30way':>12} {'100way':>12} {'470way':>12}")
print("-"*80)
for label, metric in comparison_metrics:
    values = []
    for dataset in ['30way', '100way', '470way']:
        subset = df[df['phylop_dataset'] == dataset]
        if metric == 'wobble_depletion_mean':
            subset = subset[subset['n_codons'] > 0]
        values.append(subset[metric].mean())
    print(f"{label:<30} {values[0]:12.3f} {values[1]:12.3f} {values[2]:12.3f}")

print("\n" + "="*80)
print("TOP 10 MOST CONSERVED TRANSLONS")
print("="*80 + "\n")

for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset]
    top10 = subset.nlargest(10, 'feature_mean')[['translon_id', 'chrom', 'length', 'feature_mean', 'feature_frac_conserved_strong']]
    print(f"\n{dataset}:")
    print(top10.to_string(index=False))
