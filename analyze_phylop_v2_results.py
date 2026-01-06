#!/usr/bin/env python3
import pandas as pd
import numpy as np

df = pd.read_csv('notebooks/results/phylop/translon_phylop_v2_comprehensive.tsv', sep='\t')

print("="*80)
print("PHYLOP V2 CONSERVATION ANALYSIS - TRANSCRIPT-AWARE FLANKS")
print("="*80)

for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset]

    print(f"\n{'='*80}")
    print(f"Dataset: PhyloP {dataset}")
    print(f"{'='*80}")
    print(f"Translons analyzed: {len(subset):,}")

    print(f"\n{'OVERALL CONSERVATION':^80}")
    print("-"*80)
    print(f"  Mean PhyloP:              {subset['feature_mean'].mean():7.3f} ± {subset['feature_mean'].std():.3f}")
    print(f"  Median PhyloP:            {subset['feature_median'].median():7.3f}")
    print(f"  Range:                    [{subset['feature_mean'].min():.3f}, {subset['feature_mean'].max():.3f}]")

    print(f"\n{'FLANKING REGION COMPARISON':^80}")
    print("-"*80)
    print(f"  Feature mean:             {subset['feature_mean'].mean():7.3f}")
    print(f"  Upstream mean:            {subset['upstream_mean'].mean():7.3f}")
    print(f"  Downstream mean:          {subset['downstream_mean'].mean():7.3f}")
    print(f"  Feature vs flanking diff: {subset['feature_vs_flanking_diff'].mean():7.3f}")

    print(f"\n{'CONSERVATION SPECIFICITY (feature - max(flanks))':^80}")
    print("-"*80)
    spec_gt_2 = (subset['conservation_specificity'] > 2.0).sum()
    spec_gt_1 = (subset['conservation_specificity'] > 1.0).sum()
    spec_gt_05 = (subset['conservation_specificity'] > 0.5).sum()
    spec_gt_0 = (subset['conservation_specificity'] > 0.0).sum()
    spec_lt_0 = (subset['conservation_specificity'] < 0.0).sum()

    print(f"  Specificity > 2.0:  {spec_gt_2:6,} ({100*spec_gt_2/len(subset):5.1f}%)")
    print(f"  Specificity > 1.0:  {spec_gt_1:6,} ({100*spec_gt_1/len(subset):5.1f}%)")
    print(f"  Specificity > 0.5:  {spec_gt_05:6,} ({100*spec_gt_05/len(subset):5.1f}%)")
    print(f"  Specificity > 0.0:  {spec_gt_0:6,} ({100*spec_gt_0/len(subset):5.1f}%)")
    print(f"  Specificity < 0.0:  {spec_lt_0:6,} ({100*spec_lt_0/len(subset):5.1f}%) [less conserved than flanks]")

    print(f"\n{'HIGH CONFIDENCE CANDIDATES':^80}")
    print("-"*80)

    # Strong feature + high specificity
    high_conf = subset[
        (subset['feature_mean'] > 1.5) &
        (subset['conservation_specificity'] > 1.0)
    ].sort_values('conservation_specificity', ascending=False)

    print(f"  Feature >1.5 AND specificity >1.0: {len(high_conf):,}")

    if len(high_conf) > 0:
        print(f"\n  Top 10:")
        top = high_conf.head(10)[['translon_id', 'chrom', 'exonic_length', 'blockCount',
                                   'feature_mean', 'upstream_mean', 'downstream_mean',
                                   'conservation_specificity']]
        print(top.to_string(index=False))

    # Any positive feature with negative flanks
    pos_vs_neg = subset[
        (subset['feature_mean'] > 0.5) &
        ((subset['upstream_mean'] < 0) | (subset['downstream_mean'] < 0))
    ].sort_values('conservation_specificity', ascending=False)

    print(f"\n  Positive feature (>0.5) + negative flank: {len(pos_vs_neg):,}")

print(f"\n\n{'='*80}")
print("CONSENSUS ACROSS ALL THREE DATASETS")
print(f"{'='*80}\n")

# Find translons with high specificity in ALL three datasets
specific_translons = {}
for dataset in ['30way', '100way', '470way']:
    subset = df[df['phylop_dataset'] == dataset]
    specific = subset[subset['conservation_specificity'] > 0.5]
    specific_translons[dataset] = set(specific['translon_id'])

all_three = specific_translons['30way'] & specific_translons['100way'] & specific_translons['470way']
print(f"Translons with specificity > 0.5 in ALL datasets: {len(all_three)}")

if len(all_three) > 0:
    consensus_df = df[df['translon_id'].isin(all_three)].pivot_table(
        index='translon_id',
        columns='phylop_dataset',
        values=['feature_mean', 'conservation_specificity'],
        aggfunc='first'
    )

    avg_feature = df[df['translon_id'].isin(all_three)].groupby('translon_id')['feature_mean'].mean()
    avg_spec = df[df['translon_id'].isin(all_three)].groupby('translon_id')['conservation_specificity'].mean()

    consensus_summary = pd.DataFrame({
        'avg_feature_conservation': avg_feature,
        'avg_specificity': avg_spec
    }).sort_values('avg_specificity', ascending=False)

    print("\nTop 20 consensus candidates:")
    print(consensus_summary.head(20).to_string())

    output_file = 'notebooks/results/phylop/specifically_conserved_v2_consensus.tsv'
    consensus_summary.to_csv(output_file, sep='\t')
    print(f"\nFull list saved to: {output_file}")

print(f"\n{'='*80}")
print("COMPARISON: V1 vs V2 (genomic flanks vs transcript-aware)")
print(f"{'='*80}\n")

# Load v1 results for comparison
try:
    df_v1 = pd.read_csv('notebooks/results/phylop/translon_phylop_comprehensive.tsv', sep='\t')

    print("Changes in flanking region means (should be HIGHER in v2 if hitting exons):")
    print(f"{'Dataset':<10} {'V1 Upstream':<15} {'V2 Upstream':<15} {'V1 Downstream':<15} {'V2 Downstream':<15}")
    print("-"*80)

    for dataset in ['30way', '100way', '470way']:
        v1_subset = df_v1[df_v1['phylop_dataset'] == dataset]
        v2_subset = df[df['phylop_dataset'] == dataset]

        print(f"{dataset:<10} {v1_subset['upstream_mean'].mean():15.3f} {v2_subset['upstream_mean'].mean():15.3f} "
              f"{v1_subset['downstream_mean'].mean():15.3f} {v2_subset['downstream_mean'].mean():15.3f}")

    print("\nNote: If V2 flanks are SIMILAR to V1, most translons may be single-exon")
    print("      or flanks genuinely fall in low-conservation regions adjacent to exons.")

except FileNotFoundError:
    print("V1 results not found for comparison")

print(f"\n{'='*80}")
print("EXON STRUCTURE ANALYSIS")
print(f"{'='*80}\n")

subset_470 = df[df['phylop_dataset'] == '470way']
print(f"Single-exon translons: {(subset_470['blockCount'] == 1).sum():,} ({100*(subset_470['blockCount'] == 1).sum()/len(subset_470):.1f}%)")
print(f"Multi-exon translons:  {(subset_470['blockCount'] > 1).sum():,} ({100*(subset_470['blockCount'] > 1).sum()/len(subset_470):.1f}%)")

print("\nSpecificity by exon structure (470way):")
single_exon = subset_470[subset_470['blockCount'] == 1]
multi_exon = subset_470[subset_470['blockCount'] > 1]

print(f"  Single-exon mean specificity: {single_exon['conservation_specificity'].mean():.3f}")
print(f"  Multi-exon mean specificity:  {multi_exon['conservation_specificity'].mean():.3f}")

print(f"\n  Single-exon with spec >0.5: {(single_exon['conservation_specificity'] > 0.5).sum():,} ({100*(single_exon['conservation_specificity'] > 0.5).sum()/len(single_exon):.1f}%)")
print(f"  Multi-exon with spec >0.5:  {(multi_exon['conservation_specificity'] > 0.5).sum():,} ({100*(multi_exon['conservation_specificity'] > 0.5).sum()/len(multi_exon):.1f}%)")
