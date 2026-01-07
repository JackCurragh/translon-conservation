#!/usr/bin/env python3
"""
Analyze v3 results with exonic flanks.

Full pipeline:
1. Add CDS overlap analysis
2. Identify restricted conservation candidates
3. Compare v2 vs v3 results
"""

import pandas as pd
import numpy as np
import pyranges as pr
import gzip
from pathlib import Path

print("="*80)
print("V3 PHYLOP ANALYSIS - EXONIC FLANKS")
print("="*80)

# ============================================================================
# Load v3 results
# ============================================================================

print("\n" + "="*80)
print("LOADING V3 RESULTS")
print("="*80)

v3_df = pd.read_csv('notebooks/results/phylop/translon_phylop_v3_exonic_flanks.tsv', sep='\t')
print(f"Loaded {len(v3_df):,} translons")

print(f"\nMean PhyloP scores:")
print(f"  Feature:    {v3_df['feature_mean'].mean():6.3f}")
print(f"  Upstream:   {v3_df['upstream_mean'].mean():6.3f}")
print(f"  Downstream: {v3_df['downstream_mean'].mean():6.3f}")
print(f"  Specificity: {v3_df['conservation_specificity'].mean():6.3f}")

# ============================================================================
# Add CDS overlap analysis
# ============================================================================

print("\n" + "="*80)
print("ADDING CDS OVERLAP ANALYSIS")
print("="*80)

gencode_gtf = Path('data/gencode.v46.annotation.gtf.gz')

if not gencode_gtf.exists():
    print(f"WARNING: {gencode_gtf} not found, cannot do CDS overlap analysis")
    v3_df['has_cds_overlap'] = True
    v3_df['overlapping_genes'] = 'unknown'
else:
    print(f"Extracting CDS features from {gencode_gtf}...")
    cds_records = []

    with gzip.open(gencode_gtf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]

            attrs = {}
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr:
                    parts = attr.split(' ', 1)
                    if len(parts) == 2:
                        key, val = parts
                        attrs[key] = val.strip('"')

            cds_records.append({
                'Chromosome': chrom,
                'Start': start,
                'End': end,
                'Strand': strand,
                'gene_name': attrs.get('gene_name', ''),
            })

    print(f"Found {len(cds_records):,} CDS features")

    # Create PyRanges
    cds_pr = pr.PyRanges(pd.DataFrame(cds_records))

    translon_records = []
    for _, row in v3_df.iterrows():
        translon_records.append({
            'Chromosome': row['chrom'],
            'Start': row['start'],
            'End': row['end'],
            'Strand': row['strand'],
            'translon_id': row['translon_id']
        })

    translon_pr = pr.PyRanges(pd.DataFrame(translon_records))

    print("Finding overlaps...")
    overlaps = translon_pr.join(cds_pr, how='left', suffix='_cds')

    # Process overlaps
    overlap_results = []
    for idx, row in v3_df.iterrows():
        translon_id = row['translon_id']

        translon_overlaps = overlaps.df[
            (overlaps.df['Chromosome'] == row['chrom']) &
            (overlaps.df['Start'] == row['start']) &
            (overlaps.df['End'] == row['end'])
        ]

        if len(translon_overlaps) > 0 and 'Start_cds' in translon_overlaps.columns:
            has_overlap = translon_overlaps['Start_cds'].notna().any()

            if has_overlap:
                overlapping_genes = translon_overlaps['gene_name'].dropna().unique()
                same_strand = (translon_overlaps['Strand'] == translon_overlaps['Strand_cds']).any()

                overlap_results.append({
                    'translon_id': translon_id,
                    'has_cds_overlap': True,
                    'same_strand': same_strand,
                    'overlapping_genes': ','.join(overlapping_genes[:5]) if len(overlapping_genes) > 0 else '-1'
                })
            else:
                overlap_results.append({
                    'translon_id': translon_id,
                    'has_cds_overlap': False,
                    'same_strand': False,
                    'overlapping_genes': '-1'
                })
        else:
            overlap_results.append({
                'translon_id': translon_id,
                'has_cds_overlap': False,
                'same_strand': False,
                'overlapping_genes': '-1'
            })

    overlap_df = pd.DataFrame(overlap_results)

    # Correct for PyRanges -1 sentinel
    overlap_df['has_cds_overlap'] = overlap_df['overlapping_genes'] != '-1'

    # Merge with v3_df
    v3_df = v3_df.merge(overlap_df, on='translon_id', how='left')

    print(f"\nCDS overlap results:")
    print(f"  With CDS overlap: {v3_df['has_cds_overlap'].sum():,} ({100*v3_df['has_cds_overlap'].sum()/len(v3_df):.1f}%)")
    print(f"  NO CDS overlap:   {(~v3_df['has_cds_overlap']).sum():,} ({100*(~v3_df['has_cds_overlap']).sum()/len(v3_df):.1f}%)")

# Save with CDS overlap
output_with_cds = 'notebooks/results/phylop/translon_phylop_v3_with_cds_overlap.tsv'
v3_df.to_csv(output_with_cds, sep='\t', index=False)
print(f"\nSaved: {output_with_cds}")

# ============================================================================
# Identify conserved translons without CDS overlap
# ============================================================================

print("\n" + "="*80)
print("CONSERVED TRANSLONS WITHOUT CDS OVERLAP")
print("="*80)

conserved = v3_df[v3_df['feature_mean'] > 1.5]
conserved_no_cds = conserved[~conserved['has_cds_overlap']]

print(f"\nTotal translons: {len(v3_df):,}")
print(f"Conserved (PhyloP >1.5): {len(conserved):,} ({100*len(conserved)/len(v3_df):.1f}%)")
print(f"  With CDS overlap: {conserved['has_cds_overlap'].sum():,}")
print(f"  Without CDS overlap: {len(conserved_no_cds):,}")

# ============================================================================
# Identify RESTRICTED conservation (low exonic flanks)
# ============================================================================

print("\n" + "="*80)
print("RESTRICTED CONSERVATION ANALYSIS (EXONIC FLANKS)")
print("="*80)

conserved_no_cds = conserved_no_cds.copy()
conserved_no_cds['max_flank'] = conserved_no_cds[['upstream_mean', 'downstream_mean']].max(axis=1)

print(f"\nStarting with {len(conserved_no_cds):,} conserved translons without CDS overlap")
print(f"Mean max_flank: {conserved_no_cds['max_flank'].mean():.3f}")

# Apply threshold
restricted = conserved_no_cds[conserved_no_cds['max_flank'] <= 1.0]
regional = conserved_no_cds[conserved_no_cds['max_flank'] > 1.0]

print(f"\nREGIONAL conservation (exonic flanks >1.0):")
print(f"  Count: {len(regional):,} ({100*len(regional)/len(conserved_no_cds):.1f}%)")
print(f"  Mean feature PhyloP: {regional['feature_mean'].mean():.3f}")
print(f"  Mean max_flank: {regional['max_flank'].mean():.3f}")

print(f"\nRESTRICTED conservation (exonic flanks ≤1.0):")
print(f"  Count: {len(restricted):,} ({100*len(restricted)/len(conserved_no_cds):.1f}%)")
print(f"  Mean feature PhyloP: {restricted['feature_mean'].mean():.3f}")
print(f"  Mean max_flank: {restricted['max_flank'].mean():.3f}")

# Save restricted candidates
restricted_sorted = restricted.sort_values('conservation_specificity', ascending=False)
output_restricted = 'notebooks/results/phylop/v3_restricted_conservation_candidates.tsv'
restricted_sorted.to_csv(output_restricted, sep='\t', index=False)
print(f"\nSaved: {output_restricted}")

print(f"\n{'='*80}")
print("V3 SUMMARY FOR MANUSCRIPT")
print(f"{'='*80}")

total = len(v3_df)
conserved_total = len(conserved)
conserved_with_cds = conserved['has_cds_overlap'].sum()
conserved_no_cds_count = len(conserved_no_cds)
regional_count = len(regional)
restricted_count = len(restricted)

print(f"\nTotal translons: {total:,}")
print(f"Strongly conserved (PhyloP >1.5): {conserved_total:,} ({100*conserved_total/total:.1f}%)")
print(f"\nOf the {conserved_total:,} conserved translons:")
print(f"  - Overlap CDS: {conserved_with_cds:,} ({100*conserved_with_cds/conserved_total:.1f}%)")
print(f"    → Conservation explained by overlapping protein-coding gene")
print(f"\n  - No CDS overlap: {conserved_no_cds_count:,} ({100*conserved_no_cds_count/conserved_total:.1f}%)")
print(f"    Of these:")
print(f"    • Regional conservation (exonic flanks >1.0): {regional_count:,} ({100*regional_count/conserved_no_cds_count:.1f}%)")
print(f"      → Sit in broadly conserved exonic regions")
print(f"    • RESTRICTED conservation (exonic flanks ≤1.0): {restricted_count:,} ({100*restricted_count/conserved_no_cds_count:.1f}%)")
print(f"      → Conservation appears specific to translon")

print(f"\n\nKEY FINDING (V3 with exonic flanks):")
print(f"  Of {total:,} translons, {restricted_count:,} ({100*restricted_count/total:.1f}%) show")
print(f"  conservation that is:")
print(f"    1. Strong (PhyloP >1.5)")
print(f"    2. Independent of known CDS")
print(f"    3. Specific to feature (exonic flanks ≤1.0)")

# ============================================================================
# Compare v2 vs v3
# ============================================================================

print("\n" + "="*80)
print("COMPARISON: V2 (genomic flanks) vs V3 (exonic flanks)")
print("="*80)

v2_df = pd.read_csv('notebooks/results/phylop/translon_phylop_v2_with_cds_overlap_CORRECTED.tsv', sep='\t')

# Merge for comparison
comparison = v3_df[['translon_id', 'upstream_mean', 'downstream_mean', 'conservation_specificity']].merge(
    v2_df[['translon_id', 'upstream_mean', 'downstream_mean', 'conservation_specificity']],
    on='translon_id',
    suffixes=('_v3', '_v2'),
    how='inner'
)

comparison['upstream_diff'] = comparison['upstream_mean_v3'] - comparison['upstream_mean_v2']
comparison['downstream_diff'] = comparison['downstream_mean_v3'] - comparison['downstream_mean_v2']

print(f"\nFlank conservation changes (v3 - v2):")
print(f"  Upstream mean diff:   {comparison['upstream_diff'].mean():+6.3f}")
print(f"  Downstream mean diff: {comparison['downstream_diff'].mean():+6.3f}")

print(f"\nInterpretation:")
if comparison['upstream_diff'].mean() > 0.1 or comparison['downstream_diff'].mean() > 0.1:
    print(f"  → V3 exonic flanks are MORE conserved than v2 genomic flanks")
    print(f"  → V2 was hitting introns (low conservation)")
    print(f"  → V3 correctly samples exonic sequence")
else:
    print(f"  → V3 and v2 flanks show similar conservation")
    print(f"  → Most translons may be in single-exon transcripts")
    print(f"  → Or genomic flanks were already exonic")

# Count changes in restricted classification
v2_restricted = pd.read_csv('notebooks/results/phylop/restricted_conservation_candidates.tsv', sep='\t')
v2_restricted_ids = set(v2_restricted['translon_id'])
v3_restricted_ids = set(restricted['translon_id'])

overlap = v2_restricted_ids & v3_restricted_ids
v2_only = v2_restricted_ids - v3_restricted_ids
v3_only = v3_restricted_ids - v2_restricted_ids

print(f"\n{'='*80}")
print("RESTRICTED CONSERVATION: V2 vs V3")
print(f"{'='*80}")

print(f"\nV2 (genomic flanks ≤1.0):  {len(v2_restricted_ids):,} translons")
print(f"V3 (exonic flanks ≤1.0):   {len(v3_restricted_ids):,} translons")
print(f"\nOverlap (in both):         {len(overlap):,} translons")
print(f"V2 only (lost in v3):      {len(v2_only):,} translons")
print(f"V3 only (new in v3):       {len(v3_only):,} translons")

if len(v2_only) > 0:
    print(f"\nTranslons lost in v3 (were v2 'restricted', now v3 'regional'):")
    print(f"  → Their genomic flanks were low (hitting introns)")
    print(f"  → Their exonic flanks are high (in conserved exonic regions)")

if len(v3_only) > 0:
    print(f"\nTranslons gained in v3 (were v2 'regional', now v3 'restricted'):")
    print(f"  → Their genomic flanks were high (hitting conserved exons)")
    print(f"  → Their exonic flanks are low (not in conserved exonic context)")

# Save comparison
comparison_output = 'notebooks/results/phylop/v2_v3_flank_comparison.tsv'
comparison.to_csv(comparison_output, sep='\t', index=False)
print(f"\nFull comparison saved: {comparison_output}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
