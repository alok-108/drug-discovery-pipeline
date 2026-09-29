"""GWAS summary statistics parser and quality-control filter.

Provides standardized ingestion, formatting, genomic filtering, and
Manhattan/QQ preprocessing for complex trait genetics.
"""

from typing import Optional, Dict, Any, List
import pandas as pd
import numpy as np


class GWASParser:
    """Parser and quality control processor for GWAS summary statistics.

    Attributes:
        filepath: Path to the tab-delimited GWAS summary statistics file.
        df: Pandas DataFrame containing the ingested summary data.
    """

    REQUIRED_COLUMNS = ["SNP", "CHR", "BP", "REF", "ALT", "BETA", "SE", "PVALUE"]

    def __init__(self, filepath: str) -> None:
        """Initializes the parser with the GWAS summary file path.

        Args:
            filepath: Path to a TSV/CSV GWAS summary statistics file.
        """
        self.filepath = filepath
        self.df: pd.DataFrame = pd.DataFrame()

    def load_data(self) -> pd.DataFrame:
        """Loads and standardizes GWAS summary data.

        Returns:
            pd.DataFrame: Cleaned DataFrame with validated columns.

        Raises:
            ValueError: If mandatory columns are missing from the input file.
        """
        sep = "\t" if self.filepath.endswith(".tsv") else ","
        self.df = pd.read_csv(self.filepath, sep=sep)
        self._validate_schema()
        self._compute_derived_metrics()
        return self.df

    def _validate_schema(self) -> None:
        """Checks for the presence of mandatory GWAS columns."""
        missing = [col for col in self.REQUIRED_COLUMNS if col not in self.df.columns]
        if missing:
            raise ValueError(f"Missing required GWAS columns: {missing}")

    def _compute_derived_metrics(self) -> None:
        """Computes statistical metrics like -log10(p) and z-scores."""
        # Avoid log(0) by clipping minimum p-value
        p_clipped = np.clip(self.df["PVALUE"].astype(float), a_min=1e-300, a_max=1.0)
        self.df["LOG10_P"] = -np.log10(p_clipped)
        self.df["ZSCORE"] = self.df["BETA"] / np.clip(self.df["SE"].astype(float), a_min=1e-12, a_max=None)

    def filter_significant_variants(
        self,
        p_threshold: float = 5e-8,
        min_maf: Optional[float] = None
    ) -> pd.DataFrame:
        """Filters variants passing genome-wide statistical significance.

        Args:
            p_threshold: Significance cutoff (default: 5e-8).
            min_maf: Minimum minor allele frequency threshold (optional).

        Returns:
            pd.DataFrame: Filtered subset of lead and associated variants.
        """
        if self.df.empty:
            self.load_data()

        mask = self.df["PVALUE"] <= p_threshold
        if min_maf is not None and "AF" in self.df.columns:
            # Minor allele frequency is min(AF, 1 - AF)
            maf = np.minimum(self.df["AF"], 1.0 - self.df["AF"])
            mask = mask & (maf >= min_maf)

        return self.df[mask].copy().sort_values(by="PVALUE", ascending=True)

    def get_lead_snps_per_locus(self, window_bp: int = 500_000, p_threshold: float = 5e-8) -> pd.DataFrame:
        """Identifies lead independent variants across genomic windows (distance-based clumping).

        Args:
            window_bp: Distance window in base pairs for grouping loci.
            p_threshold: Significance threshold for selecting candidate lead SNPs.

        Returns:
            pd.DataFrame: Table of lead SNPs per genomic region.
        """
        sig_df = self.filter_significant_variants(p_threshold=p_threshold)
        if sig_df.empty:
            return pd.DataFrame()

        lead_snps: List[Dict[str, Any]] = []
        for chrom, chr_df in sig_df.groupby("CHR"):
            chr_df = chr_df.sort_values(by="PVALUE", ascending=True)
            assigned_positions: List[int] = []

            for _, row in chr_df.iterrows():
                pos = int(row["BP"])
                # Check if within window of an already picked lead SNP
                if any(abs(pos - p) < window_bp for p in assigned_positions):
                    continue
                assigned_positions.append(pos)
                lead_snps.append(row.to_dict())

        return pd.DataFrame(lead_snps)
