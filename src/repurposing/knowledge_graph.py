"""CausalPathKG knowledge graph and drug repurposing engine.

Traverses multi-modal biological graphs linking GWAS candidate targets,
pathways, disease phenotypes, and FDA-approved compounds to discover repurposing candidates.
"""

from typing import Dict, List, Any
import networkx as nx
import pandas as pd


class CausalPathKG:
    """Knowledge Graph traversal engine for causal target-to-drug repurposing.

    Constructs a directed semantic graph connecting:
    - Disease nodes (e.g. Type 2 Diabetes, Alzheimer's, AML)
    - Target genes (e.g. TCF7L2, NTRK1, FLT3)
    - Biological pathways (e.g. Wnt signaling, Neurotrophin pathway, RTK-RAS)
    - Approved and investigational drugs (e.g. Metformin, Larotrectinib, Midostaurin)
    """

    def __init__(self) -> None:
        """Initializes empty DiGraph with standard biomedical ontology schema."""
        self.graph = nx.DiGraph()
        self._populate_reference_graph()

    def _populate_reference_graph(self) -> None:
        """Pre-populates canonical clinical targets and approved repurposing pairs."""
        edges = [
            # T2D Target relationships
            ("T2D", "TCF7L2", {"relation": "GWAS_ASSOCIATED", "weight": 0.95}),
            ("T2D", "GLP1R", {"relation": "GWAS_ASSOCIATED", "weight": 0.90}),
            ("T2D", "PPARG", {"relation": "GWAS_ASSOCIATED", "weight": 0.88}),
            ("GLP1R", "Semaglutide", {"relation": "AGONIST_DRUG", "weight": 0.99, "status": "Approved"}),
            ("PPARG", "Pioglitazone", {"relation": "AGONIST_DRUG", "weight": 0.95, "status": "Approved"}),
            ("TCF7L2", "Wnt_Signaling", {"relation": "PARTICIPATES_IN", "weight": 0.85}),
            ("Wnt_Signaling", "Pyrvinium", {"relation": "INHIBITOR_DRUG", "weight": 0.78, "status": "Repurposable"}),

            # Alzheimer's Target relationships
            ("Alzheimers_Disease", "NTRK1", {"relation": "GWAS_ASSOCIATED", "weight": 0.92}),
            ("Alzheimers_Disease", "TREM2", {"relation": "GWAS_ASSOCIATED", "weight": 0.89}),
            ("Alzheimers_Disease", "APOE", {"relation": "GWAS_ASSOCIATED", "weight": 0.99}),
            ("NTRK1", "Larotrectinib", {"relation": "INHIBITOR_DRUG", "weight": 0.94, "status": "Approved_Oncology"}),
            ("TREM2", "AL002", {"relation": "AGONIST_ANTIBODY", "weight": 0.85, "status": "Phase_2"}),

            # AML Target relationships
            ("AML", "FLT3", {"relation": "ONCOGENIC_DRIVER", "weight": 0.98}),
            ("AML", "IDH1", {"relation": "MUTATED_METABOLIC", "weight": 0.92}),
            ("FLT3", "Midostaurin", {"relation": "TARGETED_INHIBITOR", "weight": 0.97, "status": "Approved"}),
            ("FLT3", "Gilteritinib", {"relation": "TARGETED_INHIBITOR", "weight": 0.99, "status": "Approved"}),
            ("IDH1", "Ivosidenib", {"relation": "ALLOSTERIC_INHIBITOR", "weight": 0.96, "status": "Approved"}),
        ]

        for u, v, data in edges:
            self.graph.add_edge(u, v, **data)

    def query_repurposing_candidates(
        self,
        target_gene: str,
        max_depth: int = 2
    ) -> List[Dict[str, Any]]:
        """Finds drug candidates connected to a prioritized target gene within graph hops.

        Args:
            target_gene: Gene identifier (e.g. 'NTRK1', 'GLP1R').
            max_depth: Maximum network distance to traverse.

        Returns:
            List[Dict[str, Any]]: Repurposing candidates ranked by path confidence weight.
        """
        if target_gene not in self.graph:
            return []

        candidates = []

        for drug in self.graph.nodes():
            if drug == target_gene:
                continue

            if nx.has_path(self.graph, target_gene, drug):
                paths = list(nx.all_simple_paths(self.graph, target_gene, drug, cutoff=max_depth))
                for path in paths:
                    # Calculate cumulative path weight
                    weight_prod = 1.0
                    edge_labels = []
                    for i in range(len(path) - 1):
                        edge_data = self.graph.get_edge_data(path[i], path[i+1])
                        weight_prod *= edge_data.get("weight", 0.5)
                        edge_labels.append(edge_data.get("relation", "LINKED"))

                    dest_edge = self.graph.get_edge_data(path[-2], path[-1])
                    status = dest_edge.get("status", "Unknown") if dest_edge else "Unknown"

                    if "DRUG" in edge_labels[-1] or "INHIBITOR" in edge_labels[-1] or "AGONIST" in edge_labels[-1]:
                        candidates.append({
                            "target_gene": target_gene,
                            "drug_candidate": drug,
                            "path_length": len(path) - 1,
                            "causal_path": " -> ".join(path),
                            "mechanism": edge_labels[-1],
                            "regulatory_status": status,
                            "confidence_score": round(weight_prod, 4)
                        })

        # Rank by confidence score descending
        candidates.sort(key=lambda x: x["confidence_score"], reverse=True)
        return candidates

    def find_disease_repurposing(self, disease: str) -> pd.DataFrame:
        """Discovers all drug repurposing possibilities for a disease starting from genetics.

        Args:
            disease: Disease node identifier (e.g. 'T2D', 'Alzheimers_Disease', 'AML').

        Returns:
            pd.DataFrame: Complete ranked candidate table.
        """
        if disease not in self.graph:
            return pd.DataFrame()

        # Find direct genetic targets
        targets = [
            v for _, v, d in self.graph.out_edges(disease, data=True)
            if "ASSOCIATED" in d.get("relation", "") or "DRIVER" in d.get("relation", "")
        ]

        all_candidates = []
        for t in targets:
            cands = self.query_repurposing_candidates(t)
            all_candidates.extend(cands)

        df = pd.DataFrame(all_candidates)
        if not df.empty:
            df = df.drop_duplicates(subset=["drug_candidate"]).sort_values(
                by="confidence_score", ascending=False
            ).reset_index(drop=True)
        return df
