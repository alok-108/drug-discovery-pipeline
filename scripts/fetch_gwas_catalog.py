"""NHGRI-EBI GWAS Catalog REST API client and data downloader.

Enables searching curated GWAS studies by disease trait, querying lead associations,
and retrieving open-access summary statistics without authentication.
"""

import os
import sys
import argparse
from typing import Dict, List, Any, Optional
import requests


class GWASCatalogClient:
    """Client for interacting with the NHGRI-EBI GWAS Catalog REST API.

    API Documentation: https://www.ebi.ac.uk/gwas/rest/api
    """

    BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"

    def __init__(self, session: Optional[requests.Session] = None) -> None:
        """Initializes GWAS Catalog client."""
        self.session = session or requests.Session()
        self.session.headers.update({
            "User-Agent": "drug-discovery-pipeline/0.2.0 (open-source-research)",
            "Accept": "application/json"
        })

    def list_studies(self, disease_trait: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Searches studies associated with a given disease trait or phenotype.

        Args:
            disease_trait: Disease keyword query (e.g. 'type 2 diabetes', 'Alzheimer').
            limit: Maximum number of study records to retrieve.

        Returns:
            List[Dict[str, Any]]: Parsed study accession metadata records.
        """
        endpoint = f"{self.BASE_URL}/studies/search/findByDiseaseTrait"
        params = {"diseaseTrait": disease_trait, "size": limit}
        try:
            res = self.session.get(endpoint, params=params, timeout=20)
            res.raise_for_status()
            data = res.json()
            embedded = data.get("_embedded", {}).get("studies", [])
            studies = []
            for s in embedded:
                studies.append({
                    "accession_id": s.get("accessionId"),
                    "title": s.get("studyTitle") or s.get("publicationInfo", {}).get("title"),
                    "pubmed_id": s.get("publicationInfo", {}).get("pubmedId"),
                    "author": s.get("publicationInfo", {}).get("author", {}).get("fullname"),
                    "trait": s.get("diseaseTrait", {}).get("trait"),
                    "full_pvalue_set": bool(s.get("fullPvalueSet")),
                    "ftp_path": s.get("summaryStatisticsFile") or f"ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{s.get('accessionId')}"
                })
            return studies
        except requests.RequestException as e:
            print(f"[GWAS Catalog] Search warning: {e}", file=sys.stderr)
            return []

    def fetch_associations(self, study_accession: str) -> List[Dict[str, Any]]:
        """Retrieves lead variant associations curated for a specific study accession.

        Args:
            study_accession: Catalog accession ID (e.g. 'GCST006867').

        Returns:
            List[Dict[str, Any]]: Variant-level association records.
        """
        endpoint = f"{self.BASE_URL}/studies/{study_accession}/associations"
        try:
            res = self.session.get(endpoint, timeout=20)
            res.raise_for_status()
            data = res.json()
            embedded = data.get("_embedded", {}).get("associations", [])
            results = []
            for item in embedded:
                loci = item.get("loci", [])
                snp = "Unknown"
                if loci and loci[0].get("strongestRiskAlleles"):
                    snp = loci[0]["strongestRiskAlleles"][0].get("riskAlleleName", "Unknown")
                results.append({
                    "study_accession": study_accession,
                    "lead_variant": snp,
                    "pvalue": item.get("pvalue"),
                    "beta": item.get("betaNum"),
                    "se": item.get("standardError")
                })
            return results
        except requests.RequestException as e:
            print(f"[GWAS Catalog] Association fetch error for {study_accession}: {e}", file=sys.stderr)
            return []

    def download_summary_stats(
        self,
        study_accession: str,
        out_dir: str = "data/real"
    ) -> str:
        """Downloads summary statistics metadata file for a study accession.

        Args:
            study_accession: Catalog accession ID.
            out_dir: Output base directory.

        Returns:
            str: Path to the written summary stats or metadata record.
        """
        target_dir = os.path.join(out_dir, study_accession)
        os.makedirs(target_dir, exist_ok=True)
        out_path = os.path.join(target_dir, f"{study_accession}_associations.json")

        assocs = self.fetch_associations(study_accession)
        import json
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump({
                "accession_id": study_accession,
                "associations_count": len(assocs),
                "associations": assocs
            }, f, indent=2)

        return out_path


def main() -> None:
    """CLI execution entrypoint."""
    parser = argparse.ArgumentParser(description="Query and download NHGRI-EBI GWAS Catalog studies.")
    parser.add_argument("--disease", type=str, default="type 2 diabetes", help="Disease trait query keyword")
    parser.add_argument("--limit", type=int, default=5, help="Maximum number of studies to fetch")
    parser.add_argument("--out", type=str, default="data/real", help="Destination directory")
    args = parser.parse_args()

    client = GWASCatalogClient()
    print(f"Searching GWAS Catalog for: '{args.disease}' (limit={args.limit})...")
    studies = client.list_studies(args.disease, limit=args.limit)
    print(f"Found {len(studies)} study records.")

    for s in studies:
        print(f" - [{s['accession_id']}] {s['trait']} ({s['author'] or 'Author N/A'})")
        file_path = client.download_summary_stats(s["accession_id"], out_dir=args.out)
        print(f"   Saved metadata to: {file_path}")


if __name__ == "__main__":
    main()
