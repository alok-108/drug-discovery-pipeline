"""GenT and MuGenT gene-level association testing wrapper.

Implements gene-level aggregation from summary statistics using Cauchy
combination and multi-ancestry meta-analytic burden methods.
"""

from typing import Dict, List, Optional
import numpy as np
import pandas as pd
from scipy.stats import cauchy


class GenTWrapper:
    """Computes gene-based association tests and multi-ancestry aggregation (MuGenT).

    Aggregates variant-level GWAS summary statistics into gene-level p-values
    using the Cauchy Combination Test (CCT), which is robust to arbitrary
    linkage disequilibrium (LD) correlation structures.
    """

    def __init__(self, gene_window_kb: int = 50) -> None:
        """Initializes GenT wrapper.

        Args:
            gene_window_kb: Upstream/downstream window (in kb) around gene boundaries.
        """
        self.gene_window_kb = gene_window_kb

    def aggregate_gene_pvalues(
        self,
        gwas_df: pd.DataFrame,
        gene_col: str = "GENE",
        p_col: str = "PVALUE"
    ) -> pd.DataFrame:
        """Performs Cauchy combination test (GenT algorithm) per candidate gene.

        Args:
            gwas_df: DataFrame containing GWAS summary stats with gene annotations.
            gene_col: Column indicating gene name or symbol.
            p_col: Column containing single-variant p-values.

        Returns:
            pd.DataFrame: Table of genes ranked by combined p-value and FDR q-value.
        """
        if gene_col not in gwas_df.columns:
            raise KeyError(f"Gene annotation column '{gene_col}' not found in GWAS dataframe.")

        results: List[Dict[str, float]] = []

        for gene, group in gwas_df.groupby(gene_col):
            pvals = group[p_col].dropna().values
            if len(pvals) == 0:
                continue

            # Clip p-values for numerical stability in tan transform
            pvals_clipped = np.clip(pvals, 1e-15, 1.0 - 1e-15)
            # Cauchy test statistic: T = sum(tan((0.5 - p_i) * pi)) / k
            weights = np.ones_like(pvals_clipped) / len(pvals_clipped)
            t_cauchy = np.sum(weights * np.tan((0.5 - pvals_clipped) * np.pi))

            # Combined p-value from Cauchy survival function
            combined_p = cauchy.sf(t_cauchy)
            combined_p = float(np.clip(combined_p, 1e-300, 1.0))

            results.append({
                "GENE": gene,
                "NUM_VARIANTS": len(pvals),
                "MIN_PVALUE": float(np.min(pvals)),
                "COMBINED_PVALUE": combined_p,
                "CAUCHY_STAT": float(t_cauchy)
            })

        out_df = pd.DataFrame(results)
        if out_df.empty:
            return out_df

        # Sort and calculate Benjamini-Hochberg FDR
        out_df = out_df.sort_values(by="COMBINED_PVALUE", ascending=True).reset_index(drop=True)
        m = len(out_df)
        out_df["RANK"] = np.arange(1, m + 1)
        out_df["FDR_QVALUE"] = np.minimum.accumulate(
            (out_df["COMBINED_PVALUE"] * m / out_df["RANK"])[::-1]
        )[::-1]
        out_df["FDR_QVALUE"] = np.clip(out_df["FDR_QVALUE"], 0.0, 1.0)
        return out_df

    def run_mugent_multiancestry(
        self,
        ancestry_gwas_dict: Dict[str, pd.DataFrame],
        ancestry_weights: Optional[Dict[str, float]] = None
    ) -> pd.DataFrame:
        """Runs Multi-ancestry Gene-based Test (MuGenT) combining diverse cohorts.

        Args:
            ancestry_gwas_dict: Mapping of ancestry codes (EUR, AFR, EAS, SAS, AMR) to GWAS DataFrames.
            ancestry_weights: Optional weights per ancestry. Defaults to equal weighting.

        Returns:
            pd.DataFrame: Cross-ancestry harmonized gene prioritization table.
        """
        if not ancestry_gwas_dict:
            return pd.DataFrame()

        num_ancestries = len(ancestry_gwas_dict)
        weights = ancestry_weights or {k: 1.0 / num_ancestries for k in ancestry_gwas_dict}

        # Run single-ancestry GenT on each cohort
        gene_stats_by_ancestry: Dict[str, pd.DataFrame] = {}
        all_genes = set()

        for ancestry, df in ancestry_gwas_dict.items():
            res = self.aggregate_gene_pvalues(df)
            gene_stats_by_ancestry[ancestry] = res.set_index("GENE")
            all_genes.update(res["GENE"].tolist())

        meta_results = []
        for gene in all_genes:
            p_list = []
            w_list = []
            for anc, g_df in gene_stats_by_ancestry.items():
                if gene in g_df.index:
                    pval = g_df.loc[gene, "COMBINED_PVALUE"]
                    p_list.append(pval)
                    w_list.append(weights.get(anc, 1.0))

            if not p_list:
                continue

            p_arr = np.clip(np.array(p_list), 1e-15, 1.0 - 1e-15)
            w_arr = np.array(w_list) / np.sum(w_list)

            # Weighted Cauchy Combination
            t_stat = np.sum(w_arr * np.tan((0.5 - p_arr) * np.pi))
            meta_p = float(np.clip(cauchy.sf(t_stat), 1e-300, 1.0))

            meta_results.append({
                "GENE": gene,
                "NUM_ANCESTRIES": len(p_list),
                "MUGENT_PVALUE": meta_p,
                "CAUCHY_STAT": float(t_stat)
            })

        meta_df = pd.DataFrame(meta_results)
        if meta_df.empty:
            return meta_df

        meta_df = meta_df.sort_values(by="MUGENT_PVALUE", ascending=True).reset_index(drop=True)
        m = len(meta_df)
        meta_df["RANK"] = np.arange(1, m + 1)
        meta_df["FDR_QVALUE"] = np.minimum.accumulate(
            (meta_df["MUGENT_PVALUE"] * m / meta_df["RANK"])[::-1]
        )[::-1]
        return meta_df
