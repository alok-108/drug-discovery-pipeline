"""Unit tests for Protein Design module (AlphaFoldRunner, BindCraftWrapper, PeptideGenerator)."""

from protein_design.alphafold_runner import AlphaFoldRunner
from protein_design.bindcraft_wrapper import BindCraftWrapper
from protein_design.peptide_generator import PeptideGenerator
from functional_genomics.perturbomics import PerturbomicsPocketFinder
import pandas as pd


def test_alphafold_runner_prediction():
    """Verifies AlphaFold 3D structure generation, pLDDT scoring, and PDB header formatting."""
    runner = AlphaFoldRunner(model_version="AlphaFold3", num_recycles=3)
    seq = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
    res = runner.predict_structure(seq, target_name="test_peptide")

    assert res["sequence_length"] == len(seq)
    assert 0.0 <= res["mean_plddt"] <= 100.0
    assert len(res["plddt_scores"]) == len(seq)
    assert "HEADER" in res["pdb_content"]
    assert "ATOM" in res["pdb_content"]


def test_alphafold_binding_site_confidence():
    """Verifies calculation of local pocket confidence metric."""
    runner = AlphaFoldRunner()
    seq = "ACDEFGHIKLMNPQRSTVWY"
    pred = runner.predict_structure(seq)
    pocket_score = runner.evaluate_binding_site_confidence(pred, pocket_residues=[2, 3, 4])
    assert 0.0 <= pocket_score <= 100.0


def test_bindcraft_wrapper():
    """Verifies BindCraft de novo binder design and composite score filtering."""
    bindcraft = BindCraftWrapper(target_hotspots=["R120", "E154"], binder_length=60, num_iterations=5)
    designs = bindcraft.design_binders()

    assert len(designs) == 5
    top = designs[0]
    assert "design_id" in top
    assert "dG_binding" in top
    assert top["dG_binding"] < 0  # Exergonic favorable binding
    assert "complex_plddt" in top
    assert top["length"] == 60


def test_peptide_generator_macrocycle():
    """Verifies cyclic peptide design and cyclization constraint assignment."""
    gen = PeptideGenerator(random_seed=123)
    res = gen.generate_cyclic_peptide(length=10, cyclization_type="disulfide")

    assert res["length"] == 10
    assert res["sequence"].startswith("C")
    assert res["sequence"].endswith("C")
    assert "anchor" in res


def test_peptide_generator_lasso_and_amp():
    """Verifies lasso peptide threading topology and antimicrobial peptide heuristics."""
    gen = PeptideGenerator(random_seed=42)
    lasso = gen.generate_lasso_peptide(ring_length=8, loop_length=6, tail_length=5)
    assert lasso["peptide_type"] == "lasso_peptide"
    assert lasso["protease_resistant"] is True

    amp = gen.generate_antimicrobial_peptide(length=15)
    assert amp["length"] == 15
    assert "net_charge" in amp
    assert "predicted_mic_ug_ml" in amp


def test_perturbomics_pocket_finder():
    """Verifies amino-acid-level pocket discovery from tiling screens."""
    finder = PerturbomicsPocketFinder(min_pocket_residues=3, druggability_threshold=0.5)
    data = {
        "RESIDUE_NUM": [1, 2, 3, 4, 5, 6],
        "AA": ["ILE", "VAL", "LEU", "TRP", "PHE", "TYR"],
        "FITNESS_EFFECT": [-2.5, -2.8, -2.4, -2.9, -2.1, -2.7]
    }
    df = pd.DataFrame(data)
    scored = finder.calculate_residue_druggability(df)
    pockets = finder.find_pockets(scored, window_size=3)

    assert len(pockets) >= 1
    assert pockets[0]["DRUGGABILITY_SCORE"] > 0.5
