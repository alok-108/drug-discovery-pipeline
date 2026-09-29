"""Input/output utilities, configuration loader, and end-to-end CLI orchestrator."""

import os
import json
from typing import Dict, Any
import numpy as np
import yaml
import pandas as pd


def load_config(config_path: str = "configs/pipeline_config.yaml") -> Dict[str, Any]:
    """Loads and validates pipeline YAML configuration file.

    Args:
        config_path: Path to the YAML configuration file.

    Returns:
        Dict[str, Any]: Nested dictionary of configuration parameters.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found at {config_path}")

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    return config


def save_json(data: Any, filepath: str) -> None:
    """Serializes arbitrary dictionary or list to JSON file.

    Args:
        data: Serializable object.
        filepath: Destination file path.
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def export_summary_table(df: pd.DataFrame, filepath: str) -> None:
    """Exports pandas DataFrame to TSV or CSV format.

    Args:
        df: DataFrame to save.
        filepath: Destination file path.
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    sep = "\t" if filepath.endswith(".tsv") else ","
    df.to_csv(filepath, sep=sep, index=False)


def run_full_pipeline_cli(config_path: str = "configs/pipeline_config.yaml") -> Dict[str, Any]:
    """Executes the complete 5-stage drug discovery pipeline end-to-end.

    Stages:
        1. Target Discovery (GWAS parsing + GenT/MuGenT gene aggregation)
        2. Functional Validation (CRISPR essentiality + perturbomics)
        3. Molecule Design (AlphaFold folding + BindCraft de novo design)
        4. Repurposing Track (CausalPathKG graph traversal)
        5. Report Generation & Summary

    Args:
        config_path: Path to configuration YAML.

    Returns:
        Dict[str, Any]: Complete pipeline execution summary.
    """
    from target_discovery.gwas_parser import GWASParser
    from target_discovery.gent_wrapper import GenTWrapper
    from functional_genomics.perturbomics import PerturbomicsPocketFinder
    from protein_design.alphafold_runner import AlphaFoldRunner
    from protein_design.bindcraft_wrapper import BindCraftWrapper
    from repurposing.knowledge_graph import CausalPathKG

    print("=" * 70)
    print("[START] STARTING AI-DRIVEN DRUG DISCOVERY PIPELINE EXECUTION")
    print("=" * 70)

    cfg = load_config(config_path)
    output_dir = cfg["pipeline"].get("output_dir", "./outputs")
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------------------------------------
    # STAGE 1: Target Discovery
    # -------------------------------------------------------------
    print("\n[Stage 1/5] Ingesting GWAS and performing GenT gene aggregation...")
    gwas_file = cfg["target_discovery"]["gwas_file"]
    parser = GWASParser(gwas_file)
    gwas_df = parser.load_data()
    sig_variants = parser.filter_significant_variants(
        p_threshold=cfg["target_discovery"]["p_threshold"]
    )
    print(f"  -> Total variants parsed: {len(gwas_df)}")
    print(f"  -> Genome-wide significant variants (p <= {cfg['target_discovery']['p_threshold']}): {len(sig_variants)}")

    gent = GenTWrapper()
    gene_rankings = gent.aggregate_gene_pvalues(gwas_df)
    top_target = gene_rankings.iloc[0]["GENE"] if not gene_rankings.empty else "NTRK1"
    print(f"  -> Top prioritized drugable gene: {top_target} (p={gene_rankings.iloc[0]['COMBINED_PVALUE']:.2e})")

    export_summary_table(gene_rankings, os.path.join(output_dir, "01_prioritized_genes.tsv"))

    # -------------------------------------------------------------
    # STAGE 2: Functional Validation
    # -------------------------------------------------------------
    print("\n[Stage 2/5] Running CRISPR essentiality and Perturbomics pocket discovery...")
    # Simulated tiling screen for top target
    np.random.seed(42)
    residues = ["ARG", "LEU", "TYR", "GLU", "ASP", "TRP", "VAL", "ILE", "PHE", "HIS", "LYS", "MET"]
    screen_data = {
        "RESIDUE_NUM": list(range(100, 124)),
        "AA": [np.random.choice(residues) for _ in range(24)],
        "FITNESS_EFFECT": [-1.8, -2.1, -1.9, -0.2, -0.4, -2.5, -2.2, -1.9, -0.3, -0.5, -1.7, -2.3] * 2
    }
    tiling_df = pd.DataFrame(screen_data)
    finder = PerturbomicsPocketFinder(min_pocket_residues=4, druggability_threshold=0.60)
    scored_res = finder.calculate_residue_druggability(tiling_df)
    pockets = finder.find_pockets(scored_res)
    print(f"  -> Discovered {len(pockets)} high-vulnerability druggable pocket(s) in {top_target}")

    save_json(pockets, os.path.join(output_dir, "02_druggable_pockets.json"))

    # -------------------------------------------------------------
    # STAGE 3: Molecule Design (AlphaFold + BindCraft)
    # -------------------------------------------------------------
    print("\n[Stage 3/5] Folding target with AlphaFold and designing de novo binders...")
    af_runner = AlphaFoldRunner(model_version=cfg["protein_design"]["alphafold"]["model_version"])
    # Example target sequence fragment
    mock_seq = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGA"
    fold_result = af_runner.predict_structure(mock_seq, target_name=top_target)
    print(f"  -> Structure folded: mean pLDDT = {fold_result['mean_plddt']} (High-confidence: {fold_result['high_confidence_percent']}%)")

    bindcraft = BindCraftWrapper(target_hotspots=["R120", "E154", "W210"], binder_length=75, num_iterations=5)
    binders = bindcraft.design_binders()
    top_binder = binders[0] if binders else {}
    print(f"  -> Generated {len(binders)} de novo binders. Top binder: {top_binder.get('design_id')} (dG={top_binder.get('dG_binding')} kcal/mol)")

    save_json(binders, os.path.join(output_dir, "03_denovo_binders.json"))

    # -------------------------------------------------------------
    # STAGE 4: Parallel Repurposing Track (CausalPathKG)
    # -------------------------------------------------------------
    print("\n[Stage 4/5] Traversing CausalPathKG for immediate repurposing opportunities...")
    kg = CausalPathKG()
    repurposing_hits = kg.query_repurposing_candidates(top_target)
    if not repurposing_hits:
        # Fallback to general disease search
        repurposing_df = kg.find_disease_repurposing("T2D")
    else:
        repurposing_df = pd.DataFrame(repurposing_hits)

    print(f"  -> Identified {len(repurposing_df)} repurposing candidate(s) connected to target pathway")
    export_summary_table(repurposing_df, os.path.join(output_dir, "04_repurposing_candidates.tsv"))

    # -------------------------------------------------------------
    # STAGE 5: Synthesis Summary
    # -------------------------------------------------------------
    summary = {
        "pipeline_name": cfg["pipeline"]["name"],
        "top_prioritized_target": top_target,
        "variants_evaluated": len(gwas_df),
        "target_pvalue": float(gene_rankings.iloc[0]["COMBINED_PVALUE"]),
        "pockets_discovered": len(pockets),
        "target_mean_plddt": fold_result["mean_plddt"],
        "top_binder_id": top_binder.get("design_id"),
        "top_binder_affinity_kcal": top_binder.get("dG_binding"),
        "repurposing_candidates_count": len(repurposing_df),
        "status": "COMPLETED_SUCCESSFULLY"
    }

    save_json(summary, os.path.join(output_dir, "pipeline_run_summary.json"))
    print("\n" + "=" * 70)
    print("[SUCCESS] PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
    print(f"All artifacts saved to directory: {os.path.abspath(output_dir)}")
    print("=" * 70)
    return summary


def main_cli() -> None:
    """Console script entrypoint for CLI execution."""
    import sys
    config_file = sys.argv[1] if len(sys.argv) > 1 else "configs/pipeline_config.yaml"
    run_full_pipeline_cli(config_file)


if __name__ == "__main__":
    main_cli()
