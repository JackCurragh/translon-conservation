#!/usr/bin/env python3
"""
PhyloP Conservation Analysis v3 - Exonic Flanking Regions

Key improvement over v2:
- v2: Used genomic flanks (before first exon, after last exon)
- v3: Uses EXONIC flanks (walks along transcript exons, skips introns)

This ensures flanks don't have artificially low conservation from introns.
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
TRANSLON_FILE = 'phase1_w_transcript.tsv'
GENCODE_GTF = 'data/gencode.v46.annotation.gtf.gz'

# PhyloP bigWig files
PHYLOP_FILES = {
    '30way': 'data/phylop/hg38.phyloP30way.bw',
    '100way': 'data/phylop/hg38.phyloP100way.bw',
    '470way': 'data/phylop/hg38.phyloP470way.bw'
}

# Output
OUTPUT_DIR = Path('notebooks/results/phylop')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FLANK_SIZE = 150  # bases of EXONIC sequence to extract

print("="*80)
print("PHYLOP CONSERVATION ANALYSIS v3 - EXONIC FLANKING REGIONS")
print("="*80)
print(f"\nFlank size: {FLANK_SIZE}bp of exonic sequence")
print(f"Translon file: {TRANSLON_FILE}")
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
# Load translons with transcript annotations
# ============================================================================

print("\n" + "="*80)
print("LOADING TRANSLONS")
print("="*80)

translons = pd.read_csv(TRANSLON_FILE, sep='\t')
print(f"Loaded {len(translons):,} translons")

# Parse start/end coordinates (they're comma-separated for multi-exon)
def parse_coords(coord_str):
    """Parse comma-separated coordinate string"""
    return [int(x) for x in str(coord_str).split(',') if x.strip()]

translons['start_list'] = translons['starts'].apply(parse_coords)
translons['end_list'] = translons['ends'].apply(parse_coords)

# Calculate overall boundaries and exonic length
translons['translon_start'] = translons['start_list'].apply(min)
translons['translon_end'] = translons['end_list'].apply(max)
translons['exonic_length'] = translons['start_list'].apply(len)  # number of exons

# Rename columns to match output format
translons['chrom'] = translons['chrm'].apply(lambda x: f'chr{x}' if not str(x).startswith('chr') else str(x))
translons['translon_id'] = translons['orf_name']

print(f"Transcripts with annotations: {translons['transcript'].notna().sum():,}")

# ============================================================================
# Function: Get exonic flanking regions
# ============================================================================

def get_exonic_flanks(translon_chrom, translon_starts, translon_ends, translon_strand,
                       transcript_id, flank_size=150):
    """
    Get exonic flanking regions for a translon.

    Strategy:
    1. Get all exons for the transcript
    2. Find which exons contain the translon
    3. Walk along exonic sequence to collect flank_size bases:
       - Upstream: collect exonic bases before translon
       - Downstream: collect exonic bases after translon
    4. If we run out of transcript exons, extend genomically

    Returns:
        upstream_coords: list of (start, end) tuples for upstream flank
        downstream_coords: list of (start, end) tuples for downstream flank
    """

    # Get exons for this transcript
    exons = transcript_exons.get(transcript_id, [])

    if not exons:
        # No transcript annotation - fall back to genomic flanks
        translon_min = min(translon_starts)
        translon_max = max(translon_ends)

        if translon_strand == '+':
            upstream_coords = [(max(0, translon_min - flank_size), translon_min)]
            downstream_coords = [(translon_max, translon_max + flank_size)]
        else:
            # For minus strand, "upstream" is downstream genomically
            upstream_coords = [(translon_max, translon_max + flank_size)]
            downstream_coords = [(max(0, translon_min - flank_size), translon_min)]

        return upstream_coords, downstream_coords

    # Build complete exonic sequence coordinates for transcript
    # Already sorted by genomic position

    translon_min = min(translon_starts)
    translon_max = max(translon_ends)

    # Separate exons into: before translon, translon, after translon
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
        # Upstream = before translon genomically
        upstream_coords = []
        downstream_coords = []

        # Upstream: walk backwards from translon
        remaining_upstream = flank_size

        # First, get any bases in translon exons before translon start
        for exon in reversed(translon_exons):
            if remaining_upstream <= 0:
                break
            if exon['start'] < translon_min:
                # Part of this exon is before translon
                take_start = max(exon['start'], translon_min - remaining_upstream)
                take_end = translon_min
                upstream_coords.insert(0, (take_start, take_end))
                remaining_upstream -= (take_end - take_start)

        # Then walk through before_exons (backwards)
        for exon in reversed(before_exons):
            if remaining_upstream <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining_upstream:
                # Take whole exon
                upstream_coords.insert(0, (exon['start'], exon['end']))
                remaining_upstream -= exon_len
            else:
                # Take partial exon
                upstream_coords.insert(0, (exon['end'] - remaining_upstream, exon['end']))
                remaining_upstream = 0

        # If still need more, extend genomically before transcript
        if remaining_upstream > 0:
            transcript_start = min(e['start'] for e in exons)
            upstream_coords.insert(0, (max(0, transcript_start - remaining_upstream), transcript_start))

        # Downstream: walk forwards from translon
        remaining_downstream = flank_size

        # First, get any bases in translon exons after translon end
        for exon in translon_exons:
            if remaining_downstream <= 0:
                break
            if exon['end'] > translon_max:
                # Part of this exon is after translon
                take_start = translon_max
                take_end = min(exon['end'], translon_max + remaining_downstream)
                downstream_coords.append((take_start, take_end))
                remaining_downstream -= (take_end - take_start)

        # Then walk through after_exons
        for exon in after_exons:
            if remaining_downstream <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining_downstream:
                # Take whole exon
                downstream_coords.append((exon['start'], exon['end']))
                remaining_downstream -= exon_len
            else:
                # Take partial exon
                downstream_coords.append((exon['start'], exon['start'] + remaining_downstream))
                remaining_downstream = 0

        # If still need more, extend genomically after transcript
        if remaining_downstream > 0:
            transcript_end = max(e['end'] for e in exons)
            downstream_coords.append((transcript_end, transcript_end + remaining_downstream))

    else:  # strand == '-'
        # For minus strand, "upstream" (5' of translon) is genomically downstream
        upstream_coords = []
        downstream_coords = []

        # Upstream (5' of gene) = genomically downstream
        remaining_upstream = flank_size

        # First, get any bases in translon exons after translon end (genomically)
        for exon in translon_exons:
            if remaining_upstream <= 0:
                break
            if exon['end'] > translon_max:
                take_start = translon_max
                take_end = min(exon['end'], translon_max + remaining_upstream)
                upstream_coords.append((take_start, take_end))
                remaining_upstream -= (take_end - take_start)

        # Then walk through after_exons (genomically downstream = transcript upstream)
        for exon in after_exons:
            if remaining_upstream <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining_upstream:
                upstream_coords.append((exon['start'], exon['end']))
                remaining_upstream -= exon_len
            else:
                upstream_coords.append((exon['start'], exon['start'] + remaining_upstream))
                remaining_upstream = 0

        # If still need more, extend genomically
        if remaining_upstream > 0:
            transcript_end = max(e['end'] for e in exons)
            upstream_coords.append((transcript_end, transcript_end + remaining_upstream))

        # Downstream (3' of gene) = genomically upstream
        remaining_downstream = flank_size

        # First, get any bases in translon exons before translon start
        for exon in reversed(translon_exons):
            if remaining_downstream <= 0:
                break
            if exon['start'] < translon_min:
                take_start = max(exon['start'], translon_min - remaining_downstream)
                take_end = translon_min
                downstream_coords.insert(0, (take_start, take_end))
                remaining_downstream -= (take_end - take_start)

        # Then walk through before_exons (backwards)
        for exon in reversed(before_exons):
            if remaining_downstream <= 0:
                break
            exon_len = exon['end'] - exon['start']
            if exon_len <= remaining_downstream:
                downstream_coords.insert(0, (exon['start'], exon['end']))
                remaining_downstream -= exon_len
            else:
                downstream_coords.insert(0, (exon['end'] - remaining_downstream, exon['end']))
                remaining_downstream = 0

        # If still need more, extend genomically
        if remaining_downstream > 0:
            transcript_start = min(e['start'] for e in exons)
            downstream_coords.insert(0, (max(0, transcript_start - remaining_downstream), transcript_start))

    return upstream_coords, downstream_coords


def extract_phylop_scores(bw, chrom, coord_list):
    """
    Extract PhyloP scores from multiple coordinate ranges.

    Args:
        bw: pyBigWig object
        chrom: chromosome name
        coord_list: list of (start, end) tuples

    Returns:
        numpy array of scores
    """
    all_scores = []

    for start, end in coord_list:
        scores = bw.values(chrom, start, end, numpy=True)
        if scores is not None:
            # Filter out NaN values
            valid_scores = scores[~np.isnan(scores)]
            if len(valid_scores) > 0:
                all_scores.extend(valid_scores)

    return np.array(all_scores) if all_scores else np.array([])


# ============================================================================
# Process each PhyloP dataset
# ============================================================================

all_results = []

for dataset_name, phylop_path in PHYLOP_FILES.items():
    print(f"\n{'='*80}")
    print(f"Processing PhyloP {dataset_name}")
    print(f"{'='*80}")

    if not Path(phylop_path).exists():
        print(f"WARNING: {phylop_path} not found, skipping")
        continue

    print(f"Opening {phylop_path}")
    bw = pyBigWig.open(phylop_path)

    results = []

    for idx, row in translons.iterrows():
        if idx % 1000 == 0:
            print(f"  Processing translon {idx+1:,}/{len(translons):,}")

        translon_id = row['translon_id']
        chrom = row['chrom']
        strand = row['strand']
        transcript_id = row['transcript']

        # Get translon exon coordinates
        translon_starts = row['start_list']
        translon_ends = row['end_list']

        # Get exonic flanks
        upstream_coords, downstream_coords = get_exonic_flanks(
            chrom, translon_starts, translon_ends, strand,
            transcript_id, FLANK_SIZE
        )

        # Extract PhyloP scores
        # Feature scores
        feature_coords = list(zip(translon_starts, translon_ends))
        feature_scores = extract_phylop_scores(bw, chrom, feature_coords)

        # Flank scores
        upstream_scores = extract_phylop_scores(bw, chrom, upstream_coords)
        downstream_scores = extract_phylop_scores(bw, chrom, downstream_coords)

        # Calculate statistics
        def calc_stats(scores):
            if len(scores) == 0:
                return {
                    'n_bases': 0,
                    'mean': np.nan,
                    'median': np.nan,
                    'std': np.nan,
                    'min': np.nan,
                    'max': np.nan,
                    'q25': np.nan,
                    'q75': np.nan
                }
            return {
                'n_bases': len(scores),
                'mean': np.mean(scores),
                'median': np.median(scores),
                'std': np.std(scores),
                'min': np.min(scores),
                'max': np.max(scores),
                'q25': np.percentile(scores, 25),
                'q75': np.percentile(scores, 75)
            }

        feature_stats = calc_stats(feature_scores)
        upstream_stats = calc_stats(upstream_scores)
        downstream_stats = calc_stats(downstream_scores)

        # Calculate derived metrics
        feature_mean = feature_stats['mean']
        upstream_mean = upstream_stats['mean']
        downstream_mean = downstream_stats['mean']

        if not np.isnan(feature_mean) and not np.isnan(upstream_mean) and not np.isnan(downstream_mean):
            conservation_specificity = feature_mean - max(upstream_mean, downstream_mean)
            flanking_mean = (upstream_mean + downstream_mean) / 2
        else:
            conservation_specificity = np.nan
            flanking_mean = np.nan

        # Store result
        result = {
            'translon_id': translon_id,
            'chrom': chrom,
            'translon_start': row['translon_start'],
            'translon_end': row['translon_end'],
            'strand': strand,
            'transcript_id': transcript_id,
            'exonic_length': sum(e - s for s, e in feature_coords),
            'phylop_dataset': dataset_name,
            'flank_size': FLANK_SIZE,
            'flank_type': 'exonic',

            # Feature stats
            'feature_n_bases': feature_stats['n_bases'],
            'feature_mean': feature_stats['mean'],
            'feature_median': feature_stats['median'],
            'feature_std': feature_stats['std'],
            'feature_min': feature_stats['min'],
            'feature_max': feature_stats['max'],
            'feature_q25': feature_stats['q25'],
            'feature_q75': feature_stats['q75'],

            # Upstream stats
            'upstream_n_bases': upstream_stats['n_bases'],
            'upstream_mean': upstream_stats['mean'],
            'upstream_median': upstream_stats['median'],
            'upstream_std': upstream_stats['std'],
            'upstream_min': upstream_stats['min'],
            'upstream_max': upstream_stats['max'],
            'upstream_q25': upstream_stats['q25'],
            'upstream_q75': upstream_stats['q75'],

            # Downstream stats
            'downstream_n_bases': downstream_stats['n_bases'],
            'downstream_mean': downstream_stats['mean'],
            'downstream_median': downstream_stats['median'],
            'downstream_std': downstream_stats['std'],
            'downstream_min': downstream_stats['min'],
            'downstream_max': downstream_stats['max'],
            'downstream_q25': downstream_stats['q25'],
            'downstream_q75': downstream_stats['q75'],

            # Derived metrics
            'conservation_specificity': conservation_specificity,
            'flanking_mean': flanking_mean,
        }

        results.append(result)

    bw.close()

    all_results.extend(results)
    print(f"Completed {dataset_name}: {len(results):,} translons")

# ============================================================================
# Save results
# ============================================================================

print(f"\n{'='*80}")
print("SAVING RESULTS")
print(f"{'='*80}")

df_results = pd.DataFrame(all_results)

output_file = OUTPUT_DIR / 'translon_phylop_v3_exonic_flanks.tsv'
df_results.to_csv(output_file, sep='\t', index=False)

print(f"\nResults saved to: {output_file}")
print(f"Total rows: {len(df_results):,}")
print(f"Datasets: {df_results['phylop_dataset'].unique()}")

print("\n" + "="*80)
print("COMPLETE")
print("="*80)
