"""Unit tests for Target Discovery module (GWASParser, GenTWrapper, XQTLIntegrator)."""

import os
import pytest
import pandas as pd

from target_discovery.gwas_parser import GWASParser
from target_discovery.gent_wrapper import GenTWrapper
from target_discovery.xqtl_integrator import XQTLIntegrator


@pytest.fixture
def sample_gwas_path():
    """Returns absolute path to sample GWAS data file."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(base_dir, "data", "sample", "example_gwas.tsv")
    assert os.path.exists(path), f"File {path} must exist"
    return path


def test_gwas_parser_load_and_filter(sample_gwas_path):
    """Verifies GWAS ingestion, schema validation, and p-value threshold filtering."""
    parser = GWASParser(sample_gwas_path)
    df = parser.load_data()
    assert not df.empty
    assert "LOG10_P" in df.columns
    assert "ZSCORE" in df.columns

    sig_df = parser.filter_significant_variants(p_threshold=5e-8)
    assert not sig_df.empty
    assert (sig_df["PVALUE"] <= 5e-8).all()


def test_gwas_parser_clumping(sample_gwas_path):
    """Verifies lead SNP genomic distance-based window clumping."""
    parser = GWASParser(sample_gwas_path)
    lead_df = parser.get_lead_snps_per_locus(window_bp=500_000, p_threshold=5e-8)
    assert not lead_df.empty
    assert "SNP" in lead_df.columns
    assert len(lead_df) <= len(parser.filter_significant_variants(5e-8))


def test_gent_gene_aggregation(sample_gwas_path):
    """Verifies GenT Cauchy combination test aggregating variant p-values to gene level."""
    parser = GWASParser(sample_gwas_path)
    gwas_df = parser.load_data()

    gent = GenTWrapper()
    gene_df = gent.aggregate_gene_pvalues(gwas_df)
    assert not gene_df.empty
    assert "COMBINED_PVALUE" in gene_df.columns
    assert "FDR_QVALUE" in gene_df.columns
    assert (gene_df["COMBINED_PVALUE"] >= 0.0).all()
    assert (gene_df["COMBINED_PVALUE"] <= 1.0).all()

    # Top target should have low p-value
    top_p = gene_df.iloc[0]["COMBINED_PVALUE"]
    assert top_p < 1e-10


def test_mugent_multiancestry():
    """Verifies MuGenT cross-ancestry meta-analysis combination across cohorts."""
    gent = GenTWrapper()
    df_eur = pd.DataFrame({
        "GENE": ["NTRK1", "GLP1R"],
        "PVALUE": [1e-12, 1e-6],
    })
    df_afr = pd.DataFrame({
        "GENE": ["NTRK1", "GLP1R"],
        "PVALUE": [1e-9, 1e-7],
    })

    meta_df = gent.run_mugent_multiancestry({"EUR": df_eur, "AFR": df_afr})
    assert not meta_df.empty
    assert "MUGENT_PVALUE" in meta_df.columns
    assert len(meta_df) == 2


def test_xqtl_colocalization():
    """Verifies Bayesian colocalization probability calculation (PP.H4 for causal sharing)."""
    integrator = XQTLIntegrator()
    gwas_sub = pd.DataFrame({
        "SNP": [f"rs{i}" for i in range(1, 11)],
        "BETA": [0.05, 0.10, 0.08, 0.95, 0.04, 0.02, 0.01, 0.03, 0.06, 0.02],
        "SE": [0.02] * 10
    })
    qtl_sub = pd.DataFrame({
        "SNP": [f"rs{i}" for i in range(1, 11)],
        "BETA": [0.03, 0.06, 0.05, 0.88, 0.02, 0.01, 0.04, 0.02, 0.05, 0.01],
        "SE": [0.02] * 10
    })

    res = integrator.run_colocalization(gwas_sub, qtl_sub, target_gene="TEST_GENE")
    assert res["num_snps"] == 10
    assert 0.0 <= res["PP_H4"] <= 1.0
    assert res["PP_H4"] > 0.5  # Strong colocalization at rs4
    assert res["lead_snp"] == "rs4"
