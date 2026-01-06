#!/usr/bin/env python3
"""
PhyloP Conservation Analysis v3 - Exonic Flanking Regions

Strategy:
1. Load translon exon structures from v2 comprehensive results (has blockSizes, blockStarts)
2. Merge with phase1_w_transcript.tsv to get transcript IDs
3. Load transcript exon structures from GENCODE GTF
4. Calculate exonic flanks by walking along transcript exons (skip introns)

Key improvement over v2:
- v2: Used genomic flanks (before first exon of translon, after last exon of translon)
- v3: Uses EXONIC flanks (walks along transcript exons, skips introns)
"""

import pandas as pd
import numpy as np
import pyBigWig
import gzip
from collections import defaultdict
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

# Input files
V2_RESULTS = 'notebooks/results/phylop/translon_phylop_v2_comprehensive.tsv'
TRANSCRIPT_ANNOTATIONS = 'phase1_w_transcript.tsv'
GENCODE_GTF = 'data/gencode.v46.annotation.gtf.gz'

# PhyloP bigWig files
PHYLOP_FILES = {
    '470way': 'data/phylop/hg38.phyloP470way.bw'
}

# Output
OUTPUT_DIR = Path('notebooks/results/phylop')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FLANK_SIZE = 150  # bases of EXONIC sequence to extract

print("="*80)
print("PHYLOP CONSERVATION ANALYSIS v3 - EXONIC FLANKING REGIONS")
print("="*80)
print(f"\nFlank size: {FLANK_SIZE}bp of exonic sequence (skips introns)")
print(f"Translon structures from: {V2_RESULTS}")
print(f"Transcript IDs from: {TRANSCRIPT_ANNOTATIONS}")
print(f"GENCODE GTF: {GENCODE_GTF}")

# ============================================================================
# Load transcript exon structures from GENCODE
# ============================================================================

print("\n" + "="*80)
print("LOADING TRANSCRIPT EXON STRUCTURES FROM GENCODE")
print("="*80)

transcript_exons = defaultdict(list)

print(f"Reading {GENCODE_GTF}...")
with gzip.open(GENCODE_GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        if fields[2] != 'exon':
            continue

        chrom = fields[0]
        start = int(fields[3]) - 1  # GTF is 1-based, convert to 0-based
        end = int(fields[4])
        strand = fields[6]

        # Parse attributes
        attrs = {}
        for attr in fields[8].split(';'):
            attr = attr.strip()
            if attr:
                parts = attr.split(' ', 1)
                if len(parts) == 2:
                    key, val = parts
                    attrs[key] = val.strip('"')

        transcript_id = attrs.get('transcript_id', '')
        if not transcript_id:
            continue

        # Store exon
        transcript_exons[transcript_id].append({
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand
        })

# Sort exons by start position for each transcript
for transcript_id in transcript_exons:
    transcript_exons[transcript_id].sort(key=lambda x: x['start'])

print(f"Loaded exons for {len(transcript_exons):,} transcripts")

# ============================================================================
# Load translon data
# ============================================================================

print("\n" + "="*80)
print("LOADING TRANSLON DATA")
print("="*80)

# Load v2 results to get exon structures (only need 470way)
print(f"Loading {V2_RESULTS}...")
v2_df = pd.read_csv(V2_RESULTS, sep='\t')
v2_df = v2_df[v2_df['phylop_dataset'] == '470way'].copy()
print(f"Loaded {len(v2_df):,} translons from v2 (470way only)")

# Load transcript annotations
print(f"Loading {TRANSCRIPT_ANNOTATIONS}...")
transcript_df = pd.read_csv(TRANSCRIPT_ANNOTATIONS, sep='\t')
print(f"Loaded {len(transcript_df):,} translon-transcript mappings")

# Merge
print("Merging datasets...")
translons = v2_df.merge(
    transcript_df[['orf_name', 'transcript']],
    left_on='translon_id',
    right_on='orf_name',
    how='left'
)

print(f"Merged: {len(translons):,} translons")
print(f"  With transcript ID: {translons['transcript'].notna().sum():,}")
print(f"  Without transcript ID: {translons['transcript'].isna().sum():,}")

# ============================================================================
# Helper functions
# ============================================================================

def parse_bed_blocks(chrom_start, block_count, block_sizes_str, block_starts_str):
    """
    Parse BED12 block notation into list of (start, end) tuples.
    """
    if pd.isna(block_sizes_str) or pd.isna(block_starts_str):
        return []

    try:
        block_sizes = [int(x) for x in str(block_sizes_str).rstrip(',').split(',')]
        block_starts = [int(x) for x in str(block_starts_str).rstrip(',').split(',')]

        blocks = []
        for i in range(int(block_count)):
            block_start = chrom_start + block_starts[i]
            block_end = block_start + block_sizes[i]
            blocks.append((block_start, block_end))

        return blocks
    except:
        return []


def get_exonic_flanks(translon_chrom, translon_blocks, translon_strand,
                      transcript_id, flank_size=150):
    """
    Get exonic flanking regions for a translon.

    Returns:
        upstream_coords: list of (start, end) tuples for upstream flank
        downstream_coords: list of (start, end) tuples for downstream flank
    """

    if not translon_blocks:
        return [], []

    translon_min = min(s for s, e in translon_blocks)
    translon_max = max(e for s, e in translon_blocks)

    # Get exons for this transcript
    exons = transcript_exons.get(transcript_id, [])

    if not exons or pd.isna(transcript_id):
        # No transcript annotation - fall back to genomic flanks
        if translon_strand == '+':
            upstream_coords = [(max(0, translon_min - flank_size), translon_min)]
            downstream_coords = [(translon_max, translon_max + flank_size)]
        else:
            upstream_coords = [(translon_max, translon_max + flank_size)]
            downstream_coords = [(max(0, translon_min - flank_size), translon_min)]

        return upstream_coords, downstream_coords

    # Separate exons into: before translon, overlapping translon, after translon
    before_exons = []
    after_exons = []
    translon_exons = []

    for exon in exons:
        if exon['end'] <= translon_min:
            before_exons.append(exon)
        elif exon['start'] >= translon_max:
            after_exons.append(exon)
        else:
            # Exon overlaps translon
            translon_exons.append(exon)

    strand = exons[0]['strand']

    # Collect upstream/downstream based on strand
    if strand == '+':
        upstream_coords = []
        downstream_coords = []

        # Upstream: walk backwards from translon
        remaining = flank_size

        # First, check if translon exons have bases before translon start
        for exon in reversed(translon_exons):
            if remaining <= 0:
                break
            if exon['start'] < translon_min:
                take_start = max(exon['start'], translon_min - remaining)
                take_end = translon_min
                upstream_coords.insert(0, (take_start, take_end))
                remaining -= (take_end - take_start)

        # Walk through before_exons (backwards)
        for exon in reversed(before_exons):
            if remaining <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining:
                upstream_coords.insert(0, (exon['start'], exon['end']))
                remaining -= exon_len
            else:
                upstream_coords.insert(0, (exon['end'] - remaining, exon['end']))
                remaining = 0

        # If still need more, extend genomically
        if remaining > 0 and exons:
            transcript_start = min(e['start'] for e in exons)
            upstream_coords.insert(0, (max(0, transcript_start - remaining), transcript_start))

        # Downstream: walk forwards from translon
        remaining = flank_size

        # First, check if translon exons have bases after translon end
        for exon in translon_exons:
            if remaining <= 0:
                break
            if exon['end'] > translon_max:
                take_start = translon_max
                take_end = min(exon['end'], translon_max + remaining)
                downstream_coords.append((take_start, take_end))
                remaining -= (take_end - take_start)

        # Walk through after_exons
        for exon in after_exons:
            if remaining <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining:
                downstream_coords.append((exon['start'], exon['end']))
                remaining -= exon_len
            else:
                downstream_coords.append((exon['start'], exon['start'] + remaining))
                remaining = 0

        # If still need more, extend genomically
        if remaining > 0 and exons:
            transcript_end = max(e['end'] for e in exons)
            downstream_coords.append((transcript_end, transcript_end + remaining))

    else:  # strand == '-'
        upstream_coords = []
        downstream_coords = []

        # For minus strand, "upstream" (5' of gene) is genomically downstream
        remaining = flank_size

        # Check if translon exons have bases after translon end (genomically)
        for exon in translon_exons:
            if remaining <= 0:
                break
            if exon['end'] > translon_max:
                take_start = translon_max
                take_end = min(exon['end'], translon_max + remaining)
                upstream_coords.append((take_start, take_end))
                remaining -= (take_end - take_start)

        # Walk through after_exons
        for exon in after_exons:
            if remaining <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining:
                upstream_coords.append((exon['start'], exon['end']))
                remaining -= exon_len
            else:
                upstream_coords.append((exon['start'], exon['start'] + remaining))
                remaining = 0

        # Extend genomically if needed
        if remaining > 0 and exons:
            transcript_end = max(e['end'] for e in exons)
            upstream_coords.append((transcript_end, transcript_end + remaining))

        # Downstream (3' of gene) = genomically upstream
        remaining = flank_size

        # Check if translon exons have bases before translon start
        for exon in reversed(translon_exons):
            if remaining <= 0:
                break
            if exon['start'] < translon_min:
                take_start = max(exon['start'], translon_min - remaining)
                take_end = translon_min
                downstream_coords.insert(0, (take_start, take_end))
                remaining -= (take_end - take_start)

        # Walk through before_exons (backwards)
        for exon in reversed(before_exons):
            if remaining <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining:
                downstream_coords.insert(0, (exon['start'], exon['end']))
                remaining -= exon_len
            else:
                downstream_coords.insert(0, (exon['end'] - remaining, exon['end']))
                remaining = 0

        # Extend genomically if needed
        if remaining > 0 and exons:
            transcript_start = min(e['start'] for e in exons)
            downstream_coords.insert(0, (max(0, transcript_start - remaining), transcript_start))

    return upstream_coords, downstream_coords


def extract_phylop_scores(bw, chrom, coord_list):
    """Extract PhyloP scores from multiple coordinate ranges."""
    all_scores = []

    for start, end in coord_list:
        scores = bw.values(chrom, start, end, numpy=True)
        if scores is not None:
            valid_scores = scores[~np.isnan(scores)]
            if len(valid_scores) > 0:
                all_scores.extend(valid_scores)

    return np.array(all_scores) if all_scores else np.array([])


# ============================================================================
# Process PhyloP data
# ============================================================================

print(f"\n{'='*80}")
print(f"Processing PhyloP 470way with exonic flanks")
print(f"{'='*80}")

phylop_path = PHYLOP_FILES['470way']

if not Path(phylop_path).exists():
    print(f"ERROR: {phylop_path} not found")
    exit(1)

print(f"Opening {phylop_path}")
bw = pyBigWig.open(phylop_path)

results = []

for idx, row in translons.iterrows():
    if idx % 1000 == 0:
        print(f"  Processing translon {idx+1:,}/{len(translons):,}")

    translon_id = row['translon_id']
    chrom = row['chrom']
    start = row['start']
    strand = row['strand']
    transcript_id = row.get('transcript', None)

    # Parse blocks
    block_count = row.get('blockCount', 1)
    block_sizes = row.get('blockSizes', '')
    block_starts = row.get('blockStarts', '')

    translon_blocks = parse_bed_blocks(start, block_count, block_sizes, block_starts)

    if not translon_blocks:
        # Single exon fallback
        translon_blocks = [(start, row['end'])]

    # Get exonic flanks
    upstream_coords, downstream_coords = get_exonic_flanks(
        chrom, translon_blocks, strand, transcript_id, FLANK_SIZE
    )

    # Extract PhyloP scores
    feature_scores = extract_phylop_scores(bw, chrom, translon_blocks)
    upstream_scores = extract_phylop_scores(bw, chrom, upstream_coords)
    downstream_scores = extract_phylop_scores(bw, chrom, downstream_coords)

    # Calculate statistics
    def calc_stats(scores):
        if len(scores) == 0:
            return {'n_bases': 0, 'mean': np.nan, 'median': np.nan, 'std': np.nan}
        return {
            'n_bases': len(scores),
            'mean': np.mean(scores),
            'median': np.median(scores),
            'std': np.std(scores)
        }

    feature_stats = calc_stats(feature_scores)
    upstream_stats = calc_stats(upstream_scores)
    downstream_stats = calc_stats(downstream_scores)

    # Derived metrics
    feature_mean = feature_stats['mean']
    upstream_mean = upstream_stats['mean']
    downstream_mean = downstream_stats['mean']

    if not np.isnan(feature_mean) and not np.isnan(upstream_mean) and not np.isnan(downstream_mean):
        conservation_specificity = feature_mean - max(upstream_mean, downstream_mean)
    else:
        conservation_specificity = np.nan

    # Store result
    result = {
        'translon_id': translon_id,
        'chrom': chrom,
        'start': start,
        'end': row['end'],
        'strand': strand,
        'transcript_id': transcript_id if pd.notna(transcript_id) else 'no_annotation',
        'exonic_length': row['exonic_length'],
        'blockCount': block_count,
        'flank_size': FLANK_SIZE,
        'flank_type': 'exonic',

        'feature_n_bases': feature_stats['n_bases'],
        'feature_mean': feature_stats['mean'],
        'feature_median': feature_stats['median'],
        'feature_std': feature_stats['std'],

        'upstream_n_bases': upstream_stats['n_bases'],
        'upstream_mean': upstream_stats['mean'],
        'upstream_median': upstream_stats['median'],
        'upstream_std': upstream_stats['std'],

        'downstream_n_bases': downstream_stats['n_bases'],
        'downstream_mean': downstream_stats['mean'],
        'downstream_median': downstream_stats['median'],
        'downstream_std': downstream_stats['std'],

        'conservation_specificity': conservation_specificity,
    }

    results.append(result)

bw.close()

# ============================================================================
# Save results
# ============================================================================

print(f"\n{'='*80}")
print("SAVING RESULTS")
print(f"{'='*80}")

df_results = pd.DataFrame(results)

output_file = OUTPUT_DIR / 'translon_phylop_v3_exonic_flanks.tsv'
df_results.to_csv(output_file, sep='\t', index=False)

print(f"\nResults saved to: {output_file}")
print(f"Total translons: {len(df_results):,}")
print(f"  With transcript annotation: {(df_results['transcript_id'] != 'no_annotation').sum():,}")
print(f"  Without transcript annotation: {(df_results['transcript_id'] == 'no_annotation').sum():,}")

# Quick summary
print(f"\n{'='*80}")
print("QUICK SUMMARY")
print(f"{'='*80}")
print(f"Mean feature PhyloP: {df_results['feature_mean'].mean():.3f}")
print(f"Mean upstream PhyloP: {df_results['upstream_mean'].mean():.3f}")
print(f"Mean downstream PhyloP: {df_results['downstream_mean'].mean():.3f}")
print(f"Mean conservation specificity: {df_results['conservation_specificity'].mean():.3f}")

print("\n" + "="*80)
print("COMPLETE")
print("="*80)
