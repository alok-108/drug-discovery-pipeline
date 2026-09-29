"""Visualization utilities for genomic, functional, and structural outputs."""

from typing import Optional, List
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_manhattan(
    gwas_df: pd.DataFrame,
    output_path: Optional[str] = None,
    p_threshold: float = 5e-8
) -> plt.Figure:
    """Generates a multi-chromosome Manhattan plot from GWAS summary statistics.

    Args:
        gwas_df: DataFrame with CHR, BP, and PVALUE or LOG10_P.
        output_path: Optional path to save the generated figure.
        p_threshold: Genome-wide significance line threshold.

    Returns:
        plt.Figure: Matplotlib figure object.
    """
    df = gwas_df.copy()
    if "LOG10_P" not in df.columns:
        df["LOG10_P"] = -np.log10(np.clip(df["PVALUE"].astype(float), 1e-300, 1.0))

    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").fillna(99).astype(int)
    df = df.sort_values(by=["CHR", "BP"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, 5))
    colors = ["#1f77b4", "#aec7e8"]

    x_labels = []
    x_ticks = []
    current_pos = 0

    for idx, (chrom, group) in enumerate(df.groupby("CHR")):
        color = colors[idx % len(colors)]
        xs = current_pos + group["BP"].values
        ax.scatter(xs, group["LOG10_P"], color=color, s=25, alpha=0.85, edgecolors="none")
        x_ticks.append(current_pos + (group["BP"].values[-1] / 2.0))
        x_labels.append(str(chrom) if chrom != 99 else "Other")
        current_pos += int(group["BP"].values[-1]) + 10_000_000

    # Draw genome-wide significance threshold
    thresh_line = -np.log10(p_threshold)
    ax.axhline(thresh_line, color="#d62728", linestyle="--", linewidth=1.5, label=f"p = {p_threshold:.1e}")

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=9)
    ax.set_xlabel("Chromosome", fontsize=11, fontweight="bold")
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=11, fontweight="bold")
    ax.set_title("Genome-Wide Association Study (GWAS) Manhattan Plot", fontsize=13, fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(axis="y", linestyle=":", alpha=0.6)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300)
    return fig


def plot_volcano(
    crispr_df: pd.DataFrame,
    output_path: Optional[str] = None,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = -0.5
) -> plt.Figure:
    """Generates a volcano plot for CRISPR essentiality screen.

    Args:
        crispr_df: DataFrame with MEAN_LOG2FC and DEPLETION_PVALUE (or FDR).
        output_path: Optional path to save figure.
        fdr_threshold: Statistical significance line.
        lfc_threshold: Depletion effect size cutoff.

    Returns:
        plt.Figure: Matplotlib figure object.
    """
    df = crispr_df.copy()
    df["LOG10_P"] = -np.log10(np.clip(df["DEPLETION_PVALUE"], 1e-100, 1.0))

    fig, ax = plt.subplots(figsize=(8, 6))

    is_sig = (df["FDR"] < fdr_threshold) & (df["MEAN_LOG2FC"] < lfc_threshold)
    ax.scatter(df.loc[~is_sig, "MEAN_LOG2FC"], df.loc[~is_sig, "LOG10_P"], color="#7f7f7f", s=30, alpha=0.6, label="Non-essential")
    ax.scatter(df.loc[is_sig, "MEAN_LOG2FC"], df.loc[is_sig, "LOG10_P"], color="#e377c2", s=45, alpha=0.9, label="Essential Hits")

    ax.axvline(lfc_threshold, color="#17becf", linestyle="--", linewidth=1.2)
    ax.set_xlabel(r"Depletion Log$_2$ Fold Change", fontsize=11, fontweight="bold")
    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=11, fontweight="bold")
    ax.set_title("CRISPR-Cas9 Dependency Screen Volcano Plot", fontsize=13, fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300)
    return fig


def plot_plddt_profile(
    plddt_scores: List[float],
    target_name: str = "Target",
    output_path: Optional[str] = None
) -> plt.Figure:
    """Plots per-residue AlphaFold pLDDT structural confidence profile.

    Args:
        plddt_scores: List of pLDDT values (0-100) per amino acid position.
        target_name: Target protein symbol.
        output_path: Optional path to save figure.

    Returns:
        plt.Figure: Matplotlib figure object.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    res_indices = np.arange(1, len(plddt_scores) + 1)

    ax.plot(res_indices, plddt_scores, color="#2ca02c", linewidth=2.0, label="pLDDT")
    ax.axhline(90, color="#1f77b4", linestyle="--", alpha=0.7, label="Very High Confidence (>90)")
    ax.axhline(70, color="#ff7f0e", linestyle="--", alpha=0.7, label="Confident Fold (>70)")
    ax.axhline(50, color="#d62728", linestyle="--", alpha=0.7, label="Low Confidence (<50)")

    ax.set_ylim(0, 105)
    ax.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
    ax.set_ylabel("pLDDT Confidence", fontsize=11, fontweight="bold")
    ax.set_title(f"AlphaFold Confidence Profile: {target_name}", fontsize=13, fontweight="bold")
    ax.legend(loc="lower right")
    ax.grid(True, linestyle=":", alpha=0.5)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300)
    return fig
