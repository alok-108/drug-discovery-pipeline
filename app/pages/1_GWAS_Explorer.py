"""Page 1: GWAS Explorer & Target Prioritization."""

import os
import sys
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.join(base_dir, "src")
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from target_discovery.gwas_parser import GWASParser
from target_discovery.gent_wrapper import GenTWrapper

st.set_page_config(page_title="GWAS Explorer", page_icon="🧬", layout="wide")
st.title("🧬 GWAS Explorer & Target Prioritization")

st.markdown("""
Upload custom GWAS summary statistics (TSV/CSV) or explore the built-in multi-trait benchmark dataset.
""")

sample_gwas = os.path.join(base_dir, "data", "sample", "example_gwas.tsv")
uploaded_file = st.file_uploader("Upload GWAS Summary Statistics (TSV or CSV)", type=["tsv", "csv", "txt"])

if uploaded_file is not None:
    df = pd.read_csv(uploaded_file, sep="\t" if uploaded_file.name.endswith(".tsv") else ",")
else:
    st.info("Using pre-loaded benchmark GWAS dataset (`data/sample/example_gwas.tsv`)")
    parser = GWASParser(sample_gwas)
    df = parser.load_data()

# Process metrics
if "LOG10_P" not in df.columns:
    df["LOG10_P"] = -np.log10(np.clip(df["PVALUE"].astype(float), 1e-300, 1.0))

p_threshold = st.slider("Genome-Wide Significance Threshold (-log10 P)", min_value=5.0, max_value=50.0, value=7.3, step=0.5)

st.subheader("Interactive Manhattan Plot")
df_plot = df.copy()
df_plot["CHR"] = pd.to_numeric(df_plot["CHR"], errors="coerce").fillna(99).astype(int)
df_plot = df_plot.sort_values(by=["CHR", "BP"]).reset_index(drop=True)

fig = px.scatter(
    df_plot,
    x="BP",
    y="LOG10_P",
    color="CHR",
    hover_data=["SNP", "GENE", "PVALUE", "BETA"],
    title="Genome-Wide Association Study (GWAS) Interactive Manhattan Plot",
    labels={"BP": "Base Pair Position", "LOG10_P": "-log10(P-value)", "CHR": "Chromosome"},
    template="plotly_dark",
    height=450
)
fig.add_hline(y=p_threshold, line_dash="dash", line_color="red", annotation_text=f"Threshold (y={p_threshold})")
st.plotly_chart(fig, use_container_width=True)

st.subheader("Gene-Level Prioritization (GenT Cauchy Combination)")
gent = GenTWrapper()
ranked_genes = gent.aggregate_gene_pvalues(df)

colA, colB = st.columns([3, 1])
with colB:
    enrich_btn = st.button("🚀 Enrich with Open Targets Evidence")

if enrich_btn:
    with st.spinner("Connecting to Open Targets GraphQL API (live query)..."):
        ranked_genes = gent.aggregate_gene_pvalues(df, use_opentargets=True)
    st.success("Successfully enriched target genes with Open Targets clinical associations!")

st.dataframe(
    ranked_genes.head(20),
    use_container_width=True
)
