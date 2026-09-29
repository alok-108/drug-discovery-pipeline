"""GenT and MuGenT gene-level association testing wrapper.

Implements gene-level aggregation from summary statistics using Cauchy
combination and multi-ancestry meta-analytic burden methods.
"""

from typing import Dict, List, Optional, Any
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

    # Common canonical Ensembl mappings for top targets
    GENE_TO_ENSEMBL = {
        "TCF7L2": "ENSG00000148737",
        "SLC30A8": "ENSG00000129965",
        "KCNJ11": "ENSG00000187486",
        "PPARG": "ENSG00000132170",
        "GLP1R": "ENSG00000112164",
        "NTRK1": "ENSG00000198400",
        "APOE": "ENSG00000130203",
        "TREM2": "ENSG00000095970",
        "BIN1": "ENSG00000136717",
        "FLT3": "ENSG00000122025",
        "NPM1": "ENSG00000181163",
        "IDH1": "ENSG00000138413",
        "DNMT3A": "ENSG00000119772",
        "ADRB2": "ENSG00000169252",
        "TNF": "ENSG00000232810",
        "PTPN22": "ENSG00000134242",
        "SNCA": "ENSG00000145335",
        "MTNR1B": "ENSG00000169213",
        "KCNQ1": "ENSG00000053918",
        "FTO": "ENSG00000140718",
        "EGFR": "ENSG00000146648",
        "BRAF": "ENSG00000157764"
    }

    def aggregate_gene_pvalues(
        self,
        gwas_df: pd.DataFrame,
        gene_col: str = "GENE",
        p_col: str = "PVALUE",
        use_opentargets: bool = False,
        ot_client: Optional[Any] = None
    ) -> pd.DataFrame:
        """Performs Cauchy combination test (GenT algorithm) per candidate gene.

        Args:
            gwas_df: DataFrame containing GWAS summary stats with gene annotations.
            gene_col: Column indicating gene name or symbol.
            p_col: Column containing single-variant p-values.
            use_opentargets: If True, queries Open Targets to enrich genes with clinical evidence.
            ot_client: Optional OpenTargetsClient instance.

        Returns:
            pd.DataFrame: Table of genes ranked by combined p-value, FDR q-value, and optional OT evidence.
        """
        if gene_col not in gwas_df.columns:
            raise KeyError(f"Gene annotation column '{gene_col}' not found in GWAS dataframe.")

        results: List[Dict[str, Any]] = []

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

        if use_opentargets:
            out_df = self._enrich_opentargets_evidence(out_df, ot_client)

        return out_df

    def _enrich_opentargets_evidence(self, df: pd.DataFrame, client: Optional[Any] = None) -> pd.DataFrame:
        """Enriches ranked genes with Open Targets associations and drug count."""
        try:
            from target_discovery.opentargets_client import OpenTargetsClient
            ot = client or OpenTargetsClient()
        except Exception:
            return df

        enriched_df = df.copy()
        top_diseases = []
        top_scores = []
        known_drug_counts = []

        for gene in enriched_df["GENE"]:
            ensembl_id = self.GENE_TO_ENSEMBL.get(str(gene).upper(), "")
            if not ensembl_id:
                top_diseases.append("N/A")
                top_scores.append(0.0)
                known_drug_counts.append(0)
                continue

            try:
                assocs = ot.get_target_associations(ensembl_id, limit=3)
                drugs = ot.get_known_drugs(ensembl_id)
                if assocs:
                    top_diseases.append(assocs[0].get("disease_name", "Unknown"))
                    top_scores.append(assocs[0].get("score", 0.0))
                else:
                    top_diseases.append("None")
                    top_scores.append(0.0)
                known_drug_counts.append(len(drugs))
            except Exception:
                top_diseases.append("Error")
                top_scores.append(0.0)
                known_drug_counts.append(0)

        enriched_df["OT_TOP_DISEASE"] = top_diseases
        enriched_df["OT_ASSOCIATION_SCORE"] = top_scores
        enriched_df["OT_KNOWN_DRUGS"] = known_drug_counts
        return enriched_df

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
