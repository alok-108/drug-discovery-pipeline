"""Page 2: Protein 3D Structure & AlphaFold Confidence Viewer."""

import os
import sys
import streamlit as st
import plotly.express as px

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.join(base_dir, "src")
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from protein_design.alphafold_runner import AlphaFoldRunner
from app.utils.ot_api import GENE_TO_UNIPROT, get_alphafold_pdb_url

st.set_page_config(page_title="Protein Viewer", page_icon="🔬", layout="wide")
st.title("🔬 Protein 3D Structure & AlphaFold Confidence Viewer")

st.markdown("""
Inspect AlphaFold predicted 3D structures and per-residue confidence scores (**pLDDT**).
""")

col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("Target Selection")
    preset_gene = st.selectbox("Select Target Gene:", list(GENE_TO_UNIPROT.keys()), index=0)
    uniprot_id = st.text_input("UniProt Accession ID:", value=GENE_TO_UNIPROT.get(preset_gene, "P00533"))
    custom_seq = st.text_area("Or input Custom Amino Acid Sequence:", value="", height=120)

runner = AlphaFoldRunner(model_version="AlphaFold3", num_recycles=3)

with col2:
    st.subheader("Structural Overview")
    pdb_url = get_alphafold_pdb_url(uniprot_id)
    st.markdown(f"**AlphaFold Database PDB File**: [{uniprot_id} on EBI]({pdb_url})")

    # Generate or fetch structure metrics
    demo_seq = custom_seq if custom_seq.strip() else ("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTC" * 2)
    pred_res = runner.predict_structure(demo_seq[:120], target_name=preset_gene)

    m1, m2, m3 = st.columns(3)
    m1.metric("Mean pLDDT", f"{pred_res['mean_plddt']:.1f}")
    m2.metric("High Confidence (>70)", f"{pred_res['high_confidence_percent']:.1f}%")
    m3.metric("Druggable Fold", "YES" if pred_res["is_druggable_fold"] else "NO")

st.markdown("---")
st.subheader("Per-Residue pLDDT Confidence Profile")

plddt_vals = pred_res["plddt_scores"]
res_indices = list(range(1, len(plddt_vals) + 1))
fig_plddt = px.line(
    x=res_indices,
    y=plddt_vals,
    labels={"x": "Residue Number", "y": "pLDDT Confidence"},
    title=f"Per-Residue Structural Confidence: {preset_gene}",
    template="plotly_dark",
    height=380
)
fig_plddt.add_hline(y=90, line_dash="dash", line_color="#0053D6", annotation_text="Very High (>90)")
fig_plddt.add_hline(y=70, line_dash="dash", line_color="#65CBF3", annotation_text="Confident (70-90)")
fig_plddt.add_hline(y=50, line_dash="dash", line_color="#FFDB13", annotation_text="Low (50-70)")
st.plotly_chart(fig_plddt, use_container_width=True)

st.subheader("AlphaFold PDB Atomic Header Preview")
st.code("\n".join(pred_res["pdb_content"].splitlines()[:15]), language="text")
