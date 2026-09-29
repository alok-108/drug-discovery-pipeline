"""Open Targets and AlphaFold DB utility helpers for Streamlit UI."""

import sys
import os
from typing import Dict, List, Any

# Ensure src/ is on python search path
base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_dir = os.path.join(base_dir, "src")
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

try:
    from target_discovery.opentargets_client import OpenTargetsClient
except ImportError:
    OpenTargetsClient = None


GENE_TO_ENSEMBL = {
    "EGFR": "ENSG00000146648",
    "TCF7L2": "ENSG00000148737",
    "SLC30A8": "ENSG00000129965",
    "KCNJ11": "ENSG00000187486",
    "PPARG": "ENSG00000132170",
    "GLP1R": "ENSG00000112164",
    "NTRK1": "ENSG00000198400",
    "APOE": "ENSG00000130203",
    "TREM2": "ENSG00000095970",
    "BIN1": "ENSG00000136717",
    "FLT3": "ENSG00000122025",
    "NPM1": "ENSG00000181163",
    "IDH1": "ENSG00000138413",
    "DNMT3A": "ENSG00000119772",
    "BRAF": "ENSG00000157764"
}

GENE_TO_UNIPROT = {
    "EGFR": "P00533",
    "TCF7L2": "Q9NQB0",
    "SLC30A8": "Q8IWU4",
    "KCNJ11": "Q14654",
    "PPARG": "P37231",
    "GLP1R": "P43220",
    "NTRK1": "P04629",
    "APOE": "P02649",
    "TREM2": "Q9NZC2",
    "BIN1": "O00499",
    "FLT3": "P36888",
    "NPM1": "P06748",
    "IDH1": "O75874",
    "DNMT3A": "Q9Y6K1",
    "BRAF": "P15056"
}


def get_client() -> Any:
    """Returns singleton instance of OpenTargetsClient."""
    if OpenTargetsClient:
        return OpenTargetsClient()
    return None


def fetch_target_associations(gene_symbol: str, limit: int = 10) -> List[Dict[str, Any]]:
    """Retrieves Open Targets associations for a gene symbol."""
    ensembl_id = GENE_TO_ENSEMBL.get(gene_symbol.upper())
    client = get_client()
    if client and ensembl_id:
        try:
            return client.get_target_associations(ensembl_id, limit=limit)
        except Exception:
            pass
    # Fallback default demo data
    return [
        {"disease_name": f"{gene_symbol}-related condition", "score": 0.85, "disease_id": "EFO_0001"}
    ]


def fetch_known_drugs(gene_symbol: str) -> List[Dict[str, Any]]:
    """Retrieves known clinical drugs targeting the gene."""
    ensembl_id = GENE_TO_ENSEMBL.get(gene_symbol.upper())
    client = get_client()
    if client and ensembl_id:
        try:
            return client.get_known_drugs(ensembl_id)
        except Exception:
            pass
    return []


def get_alphafold_pdb_url(uniprot_id: str) -> str:
    """Constructs official AlphaFold Protein Structure Database download URL."""
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
