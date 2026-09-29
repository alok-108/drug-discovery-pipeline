"""xQTL (eQTL + pQTL) integrator with Bayesian colocalization analysis.

Evaluates whether GWAS trait association and molecular quantitative trait
loci (eQTL / pQTL) share a single causal variant (coloc hypothesis H4).
"""

from typing import Dict, Any
import numpy as np
import pandas as pd


class XQTLIntegrator:
    """Integrates expression (eQTL) and protein (pQTL) summary statistics with GWAS.

    Implements approximate Bayesian colocalization testing following the Giambartolomei et al.
    framework, computing posterior probabilities for:
        H0: No association with either trait
        H1: Association with trait 1 (GWAS) only
        H2: Association with trait 2 (xQTL) only
        H3: Independent associations with both traits (different causal SNPs)
        H4: Shared causal variant for both GWAS and molecular QTL (causal target)
    """

    def __init__(
        self,
        prior_p1: float = 1e-4,
        prior_p2: float = 1e-4,
        prior_p12: float = 1e-5
    ) -> None:
        """Initializes the colocalization priors.

        Args:
            prior_p1: Prior probability that a variant is associated with trait 1.
            prior_p2: Prior probability that a variant is associated with trait 2.
            prior_p12: Prior probability of variant association with both traits.
        """
        self.p1 = prior_p1
        self.p2 = prior_p2
        self.p12 = prior_p12

    def compute_approx_bayes_factor(self, beta: np.ndarray, se: np.ndarray, prior_variance: float = 0.15**2) -> np.ndarray:
        """Calculates Wakefield's Approximate Bayes Factor (ABF) for single-variant associations.

        Args:
            beta: Effect size coefficients.
            se: Standard errors of effect sizes.
            prior_variance: Prior variance of the true effect size parameter W.

        Returns:
            np.ndarray: Approximate Bayes factor for each variant.
        """
        v = se**2
        z = beta / np.clip(se, 1e-12, None)
        r = prior_variance / (prior_variance + v)
        log_abf = 0.5 * (np.log(1.0 - r) + (r * (z**2)))
        return np.exp(np.clip(log_abf, -300, 300))

    def run_colocalization(
        self,
        gwas_sub_df: pd.DataFrame,
        qtl_sub_df: pd.DataFrame,
        target_gene: str
    ) -> Dict[str, Any]:
        """Runs Bayesian colocalization between GWAS and xQTL variants for a locus.

        Args:
            gwas_sub_df: Subset of GWAS variants for the locus (SNP, BETA, SE).
            qtl_sub_df: Subset of QTL variants for the target gene (SNP, BETA, SE).
            target_gene: Identifier of the tested candidate gene.

        Returns:
            Dict[str, Any]: Posterior probabilities (PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
                            and identification of candidate causal variants.
        """
        merged = pd.merge(
            gwas_sub_df[["SNP", "BETA", "SE"]],
            qtl_sub_df[["SNP", "BETA", "SE"]],
            on="SNP",
            suffixes=("_GWAS", "_QTL")
        ).dropna()

        if merged.empty:
            return {
                "gene": target_gene,
                "num_snps": 0,
                "PP_H0": 1.0,
                "PP_H1": 0.0,
                "PP_H2": 0.0,
                "PP_H3": 0.0,
                "PP_H4": 0.0,
                "is_colocalized": False
            }

        abf_gwas = self.compute_approx_bayes_factor(
            merged["BETA_GWAS"].values, merged["SE_GWAS"].values
        )
        abf_qtl = self.compute_approx_bayes_factor(
            merged["BETA_QTL"].values, merged["SE_QTL"].values
        )

        # Hypothesis weights
        w_h0 = 1.0
        w_h1 = np.sum(self.p1 * abf_gwas)
        w_h2 = np.sum(self.p2 * abf_qtl)
        w_h3 = np.sum(self.p1 * self.p2 * np.outer(abf_gwas, abf_qtl)) - np.sum(self.p1 * self.p2 * abf_gwas * abf_qtl)
        w_h4 = np.sum(self.p12 * abf_gwas * abf_qtl)

        total_weight = w_h0 + w_h1 + w_h2 + w_h3 + w_h4
        pp = [w_h0 / total_weight, w_h1 / total_weight, w_h2 / total_weight, w_h3 / total_weight, w_h4 / total_weight]

        top_snp_idx = int(np.argmax(abf_gwas * abf_qtl))
        lead_causal_snp = merged.iloc[top_snp_idx]["SNP"]

        return {
            "gene": target_gene,
            "num_snps": len(merged),
            "PP_H0": float(pp[0]),
            "PP_H1": float(pp[1]),
            "PP_H2": float(pp[2]),
            "PP_H3": float(pp[3]),
            "PP_H4": float(pp[4]),
            "lead_snp": lead_causal_snp,
            "is_colocalized": bool(pp[4] >= 0.75)
        }
