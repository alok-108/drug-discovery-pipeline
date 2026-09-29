"""AI-Driven Drug Discovery Pipeline: Streamlit Dashboard.

Interactive web application demonstrating:
1. GWAS locus exploration and Manhattan plotting.
2. 3D Protein structure visualization and pLDDT mapping.
3. Open Targets clinical evidence and drug repurposing queries.
"""

import streamlit as st

st.set_page_config(
    page_title="AI Drug Discovery Pipeline",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("💊 AI-Driven Drug Discovery Pipeline: From Genes to Medicines")
st.markdown("### Production-grade, open-source reference pipeline for computational therapeutics")

st.markdown("""
Welcome to the interactive companion dashboard for the **AI-Driven Drug Discovery Pipeline**! 
This web application unites statistical genetics, live clinical evidence from the **Open Targets Platform**, 
and 3D structural biology models into an intuitive exploration suite.
""")

col1, col2, col3 = st.columns(3)
with col1:
    st.info("**Stage 1 & 2: Target Discovery**\n\nParse GWAS summary statistics, run GenT Cauchy gene-based tests, and identify candidate targets.")
with col2:
    st.success("**Stage 3: Molecular Design**\n\nInspect 3D coordinates from the AlphaFold Database, evaluate pLDDT confidence scores, and locate druggable pockets.")
with col3:
    st.warning("**Stage 4 & 5: Repurposing & Translation**\n\nQuery known approved drugs, phase 1-4 clinical pipelines, and directed knowledge graphs for fast-track translation.")

st.markdown("---")
st.subheader("Interactive Module Sitemap")
st.markdown("""
Select a module from the **left sidebar navigation** to begin:
* 🧬 **Page 1: GWAS Explorer**: Upload summary statistics, render interactive multi-chromosome Manhattan plots, and enrich prioritized genes with live Open Targets evidence.
* 🔬 **Page 2: Protein Viewer**: Inspect 3D AlphaFold structures by UniProt ID or gene symbol, examine per-residue pLDDT profiles, and assess structural druggability.
* 💊 **Page 3: Drug Repurposing**: Discover approved and investigational molecules targeting candidate disease proteins, review clinical mechanisms of action, and compare evidence scores.
""")

st.markdown("---")
st.markdown("""
#### Free Tier Architecture
* **Open Targets Platform**: 100% Free GraphQL API (no API keys required).
* **AlphaFold DB**: Open access 3D coordinates hosted by EMBL-EBI.
* **Hosting**: Compatible with Streamlit Community Cloud (free deployment, 1GB RAM).
""")
