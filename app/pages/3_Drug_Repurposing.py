"""Page 3: Clinical Evidence & Drug Repurposing."""

import os
import sys
import pandas as pd
import plotly.express as px
import streamlit as st

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.join(base_dir, "src")
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from app.utils.ot_api import fetch_known_drugs, fetch_target_associations, GENE_TO_ENSEMBL

st.set_page_config(page_title="Drug Repurposing", page_icon="💊", layout="wide")
st.title("💊 Clinical Evidence & Drug Repurposing")

st.markdown("""
Search candidate therapeutic targets to retrieve known approved and investigational drugs,
clinical trial phases, and disease association scores from the **Open Targets Platform**.
""")

col1, col2 = st.columns([1, 2])
with col1:
    gene_symbol = st.selectbox(
        "Select Target Gene:",
        list(GENE_TO_ENSEMBL.keys()),
        index=0
    )
    custom_gene = st.text_input("Or enter custom gene symbol:", value="")
    target = custom_gene.strip().upper() if custom_gene.strip() else gene_symbol

with col2:
    ensembl = GENE_TO_ENSEMBL.get(target, "ENSG00000146648")
    st.info(f"Target Symbol: **{target}** | Mapped Ensembl ID: `{ensembl}`")

with st.spinner(f"Querying Open Targets GraphQL API for {target}..."):
    drugs = fetch_known_drugs(target)
    associations = fetch_target_associations(target, limit=10)

st.subheader(f"Known Approved & Investigational Drugs Targeting {target}")

if drugs:
    drug_df = pd.DataFrame(drugs)
    st.dataframe(drug_df, use_container_width=True)
else:
    st.warning(f"No existing approved or clinical drugs found for {target}. This represents a true novel de novo target opportunity!")

st.markdown("---")
st.subheader(f"Top Disease Indications & Association Scores: {target}")

if associations:
    assoc_df = pd.DataFrame(associations)
    fig = px.bar(
        assoc_df,
        x="score",
        y="disease_name",
        orientation="h",
        color="score",
        color_continuous_scale="Viridis",
        title=f"Harmonized Disease Association Scores for {target} (Open Targets)",
        labels={"score": "Overall Association Score (0 to 1)", "disease_name": "Disease Indication"},
        template="plotly_dark",
        height=400
    )
    fig.update_layout(yaxis={"categoryorder": "total ascending"})
    st.plotly_chart(fig, use_container_width=True)
else:
    st.info("No disease associations available.")
