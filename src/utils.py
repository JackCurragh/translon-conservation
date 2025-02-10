"""
Utilities for genomic conservation analysis.
Shared functionality for analyzing conservation patterns across different genomic features.
"""

import polars as pl
import pyBigWig
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class ConservationTracker:
    """Handles loading and management of PhyloCSF and PhyloP tracks."""

    def __init__(self, phylocsf_dir: Path, phylop_path: Path, phylocsf_type: str = "Raw"):
        """
        Initialize conservation tracker.
        
        Args:
            phylocsf_dir: Directory containing PhyloCSF bigwig files
            phylop_path: Path to PhyloP bigwig file
            phylocsf_type: Type of PhyloCSF scores to use ("Raw" or "Smooth")
        """
        self.phylocsf_type = phylocsf_type
        self.phylocsf_tracks = self._load_phylocsf_tracks(phylocsf_dir)
        self.phylop_track = self._load_phylop_track(phylop_path)

    def _load_phylocsf_tracks(self, phylocsf_dir: Path) -> Dict[str, pyBigWig.pyBigWig]:
        """Load PhyloCSF bigwig files for all frames and strands."""
        tracks = {}
        for frame in [1, 2, 3]:
            for strand in ['+', '-']:
                sign = 'plus' if strand == '+' else 'minus'
                filename = f"PhyloCSF{self.phylocsf_type}_{sign}{frame}.bw"
                filepath = phylocsf_dir / filename
                if filepath.exists():
                    tracks[f"{strand}{frame}"] = pyBigWig.open(str(filepath))
        return tracks

    def _load_phylop_track(self, phylop_path: Path) -> pyBigWig.pyBigWig:
        """Load PhyloP bigwig file."""
        return pyBigWig.open(str(phylop_path))

    def get_conservation_scores(self,
                                chrom: str,
                                start: int,
                                end: int,
                                strand: str) -> Tuple[Dict[int, np.ndarray], np.ndarray]:
        """Get both PhyloCSF and PhyloP scores for a region."""
        # Get PhyloCSF scores for all frames
        phylocsf_scores = {}
        for frame in [1, 2, 3]:
            track_key = f"{strand}{frame}"
            if track_key in self.phylocsf_tracks:
                scores = self._get_scores(
                    self.phylocsf_tracks[track_key],
                    chrom, start, end
                )
                if scores is not None:
                    phylocsf_scores[frame] = scores

        # Get PhyloP scores
        phylop_scores = self._get_scores(self.phylop_track, chrom, start, end)

        return phylocsf_scores, phylop_scores

    def _get_scores(self, 
                    bw_file: pyBigWig.pyBigWig, 
                    chrom: str, 
                    start: int, 
                    end: int) -> Optional[np.ndarray]:
        """Get conservation scores for a specific genomic region."""
        try:
            # Ensure chromosome format
            chrom = f"chr{chrom}" if not chrom.startswith('chr') else chrom

            # Basic validation
            if chrom not in bw_file.chroms():
                return None

            chrom_length = bw_file.chroms(chrom)
            if start < 0 or end > chrom_length:
                return None

            # Get and process scores
            scores = bw_file.values(chrom, start, end)
            return np.array(scores)

        except Exception as e:
            print(f"Error accessing bigwig for {chrom}:{start}-{end}: {str(e)}")
            return None


class ConservationAnalyzer:
    """Analyzes conservation patterns in genomic regions."""

    def __init__(self, window_size: int = 30):
        self.window_size = window_size

    def calculate_metrics(self,
                          phylocsf_scores: Dict[int, np.ndarray],
                          phylop_scores: np.ndarray) -> Dict[str, float]:
        """Calculate comprehensive conservation metrics."""
        metrics = {}

        # Process PhyloCSF scores for each frame
        for frame, scores in phylocsf_scores.items():
            frame_metrics = self._analyze_score_array(scores, f"phylocsf_frame{frame}")
            metrics.update(frame_metrics)

        # Process PhyloP scores
        if phylop_scores is not None:
            phylop_metrics = self._analyze_score_array(phylop_scores, "phylop")
            metrics.update(phylop_metrics)

        return metrics

    def _analyze_score_array(self, scores: np.ndarray, prefix: str) -> Dict[str, float]:
        """Calculate metrics for a single score array."""
        if scores is None or len(scores) == 0:
            return {}

        metrics = {
            f"{prefix}_mean": float(np.nanmean(scores)),
            f"{prefix}_median": float(np.nanmedian(scores)),
            f"{prefix}_max": float(np.nanmax(scores)),
            f"{prefix}_q75": float(np.nanpercentile(scores, 75)),
            f"{prefix}_start": float(np.nanmean(scores[:15]))
        }

        # Window analysis
        window_scores = self._get_window_scores(scores)
        metrics[f"{prefix}_max_window"] = max(window_scores) if window_scores else np.nan

        # Positive run analysis
        metrics[f"{prefix}_longest_positive_run"] = self._get_longest_positive_run(scores)

        return metrics

    def _get_window_scores(self, scores: np.ndarray) -> List[float]:
        """Calculate scores for sliding windows."""
        if len(scores) < self.window_size:
            return [float(np.nanmean(scores))]

        window_scores = []
        for i in range(len(scores) - self.window_size + 1):
            window = scores[i:i+self.window_size]
            window_scores.append(float(np.nanmean(window)))

        return window_scores

    def _get_longest_positive_run(self, scores: np.ndarray) -> int:
        """Find the longest continuous run of positive scores."""
        current_run = max_run = 0
        for score in scores:
            if not np.isnan(score) and score > 0:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        return max_run


class BEDHandler:
    """Handles BED file operations for conservation analysis."""
    @staticmethod
    def load_bed12(path: Path) -> pl.DataFrame:
        """Load a BED12 file with proper schema."""
        schema = {
            "chrom": pl.Utf8,
            "chromStart": pl.Int64,
            "chromEnd": pl.Int64,
            "name": pl.Utf8,
            "score": pl.Int64,
            "strand": pl.Utf8,
            "thickStart": pl.Int64,
            "thickEnd": pl.Int64,
            "itemRgb": pl.Int64,
            "blockCount": pl.Int64,
            "blockSizes": pl.Utf8,
            "blockStarts": pl.Utf8
        }

        return pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            new_columns=list(schema.keys()),
            schema=schema
        )

    @staticmethod
    def write_scored_bed(df: pl.DataFrame, 
                         output_path: Path,
                         score_col: str = "conservation_score",
                         scale_factor: float = 1000):
        """Write a BED file with normalized scores."""
        bed_df = (df
                  .with_columns([
                      ((pl.col(score_col) - pl.col(score_col).min()) /
                       (pl.col(score_col).max() - pl.col(score_col).min()) * scale_factor)
                      .floor()
                      .cast(pl.Int64)
                      .clip(0, 1000)
                      .alias("score")
                  ])
                  .select([
                      "chrom",
                      "chromStart",
                      "chromEnd",
                      "name",
                      "score",
                      "strand"
                  ]))

        bed_df.write_csv(output_path, separator="\t", include_header=False)