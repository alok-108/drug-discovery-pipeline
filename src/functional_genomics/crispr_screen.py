"""CRISPR-Cas9 and CRISPRi pooled screen analysis engine.

Processes guide RNA (sgRNA) read counts, performs size-factor normalization,
calculates log2 fold changes (LFC), and ranks target gene essentiality.
"""

from typing import Dict, List
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm


class CRISPRScreenAnalyzer:
    """Analyzes pooled CRISPR-Cas9 or CRISPRi depletion/enrichment screens.

    Calculates sgRNA-level log2 fold changes between treatment (T) and control (T0),
    and aggregates sgRNAs to rank genes using robust rank aggregation (RRA).
    """

    def __init__(self, pseudocount: float = 1.0) -> None:
        """Initializes CRISPR screen analyzer.

        Args:
            pseudocount: Small scalar added to count reads prior to log2 transform.
        """
        self.pseudocount = pseudocount

    def normalize_counts(self, count_df: pd.DataFrame, sample_cols: List[str]) -> pd.DataFrame:
        """Applies median-of-ratios (DESeq-style) size factor normalization.

        Args:
            count_df: DataFrame with sgRNA identifier, Gene, and raw count columns.
            sample_cols: List of column names representing experimental samples.

        Returns:
            pd.DataFrame: DataFrame with normalized counts added.
        """
        df_norm = count_df.copy()
        counts = count_df[sample_cols].values + self.pseudocount
        # Geometric mean across samples for each sgRNA
        log_geo_means = np.mean(np.log(counts), axis=1)
        geo_means = np.exp(log_geo_means)

        # Ratios to geometric mean
        ratios = counts / geo_means[:, np.newaxis]
        size_factors = np.median(ratios, axis=0)

        # Scale counts by size factors
        for idx, col in enumerate(sample_cols):
            df_norm[f"{col}_NORM"] = df_norm[col] / size_factors[idx]

        return df_norm

    def compute_sgrna_fold_change(
        self,
        count_df: pd.DataFrame,
        control_col: str,
        treatment_col: str
    ) -> pd.DataFrame:
        """Computes log2 fold change (LFC) between control and treated samples.

        Args:
            count_df: Normalized count DataFrame.
            control_col: Column name for baseline or control arm (e.g. Day 0).
            treatment_col: Column name for end-point treated arm (e.g. Day 21).

        Returns:
            pd.DataFrame: Table with added LOG2FC column.
        """
        df = count_df.copy()
        c_val = df[control_col] + self.pseudocount
        t_val = df[treatment_col] + self.pseudocount
        df["LOG2FC"] = np.log2(t_val / c_val)
        return df

    def aggregate_gene_essentiality(
        self,
        sgrna_df: pd.DataFrame,
        gene_col: str = "GENE",
        lfc_col: str = "LOG2FC"
    ) -> pd.DataFrame:
        """Aggregates sgRNAs per gene to identify essential/lethal vulnerabilities.

        Args:
            sgrna_df: DataFrame with individual sgRNA fold changes.
            gene_col: Column containing gene name.
            lfc_col: Column with log2 fold change values.

        Returns:
            pd.DataFrame: Gene-level essentiality scores, p-values, and FDR.
        """
        records: List[Dict[str, float]] = []
        all_lfc = sgrna_df[lfc_col].values
        # Empirical CDF for depletion (negative LFC = depletion)
        ranks = rankdata(all_lfc) / len(all_lfc)
        sgrna_df = sgrna_df.copy()
        sgrna_df["_RANK_SCORE"] = ranks

        for gene, group in sgrna_df.groupby(gene_col):
            gene_lfcs = group[lfc_col].values
            mean_lfc = float(np.mean(gene_lfcs))
            median_lfc = float(np.median(gene_lfcs))
            n_guides = len(gene_lfcs)

            # Modified rank test for depletion
            guide_ranks = group["_RANK_SCORE"].values
            z_score = (np.mean(guide_ranks) - 0.5) / (np.sqrt(1.0 / (12.0 * max(1, n_guides))))
            pval = float(norm.cdf(z_score))

            records.append({
                "GENE": gene,
                "NUM_GUIDES": n_guides,
                "MEAN_LOG2FC": mean_lfc,
                "MEDIAN_LOG2FC": median_lfc,
                "DEPLETION_PVALUE": pval,
            })

        out = pd.DataFrame(records)
        if out.empty:
            return out

        out = out.sort_values(by="DEPLETION_PVALUE", ascending=True).reset_index(drop=True)
        m = len(out)
        out["RANK"] = np.arange(1, m + 1)
        out["FDR"] = np.minimum.accumulate((out["DEPLETION_PVALUE"] * m / out["RANK"])[::-1])[::-1]
        out["FDR"] = np.clip(out["FDR"], 0.0, 1.0)
        out["IS_ESSENTIAL"] = (out["FDR"] < 0.05) & (out["MEAN_LOG2FC"] < -0.5)
        return out
