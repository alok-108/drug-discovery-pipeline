"""IEU OpenGWAS API client and summary statistics downloader.

Integrates with the MRC Integrative Epidemiology Unit (IEU) OpenGWAS database
hosting over 40,000 harmonized GWAS datasets and variant-level query APIs.
"""

import os
import sys
import json
import hashlib
from typing import Dict, List, Any, Optional
import requests


class IEUOpenGWASClient:
    """Client for querying the IEU OpenGWAS API.

    API Endpoint: https://gwas-api.mrcieu.ac.uk/
    """

    BASE_URL = "https://gwas-api.mrcieu.ac.uk/api"

    def __init__(self, cache_dir: Optional[str] = None) -> None:
        """Initializes IEU OpenGWAS client with persistent local caching."""
        self.cache_dir = cache_dir or os.path.expanduser(os.path.join("~", ".cache", "ieu_opengwas"))
        os.makedirs(self.cache_dir, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "drug-discovery-pipeline/0.2.0 (open-source-research)",
            "Accept": "application/json"
        })

    def _get_cache_path(self, key: str) -> str:
        digest = hashlib.md5(key.encode("utf-8")).hexdigest()
        return os.path.join(self.cache_dir, f"{digest}.json")

    def _read_cache(self, key: str) -> Optional[Any]:
        path = self._get_cache_path(key)
        if os.path.exists(path):
            try:
                with open(path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                return None
        return None

    def _write_cache(self, key: str, data: Any) -> None:
        path = self._get_cache_path(key)
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f)
        except Exception:
            pass

    def list_datasets(self, trait_filter: str = "") -> List[Dict[str, Any]]:
        """Queries metadata for GWAS datasets matching a keyword filter.

        Args:
            trait_filter: Keyword query for matching traits (e.g. 'body mass index', 'diabetes').

        Returns:
            List[Dict[str, Any]]: List of matching dataset records (ID, trait, sample size, year).
        """
        cache_key = f"list_datasets_{trait_filter.lower()}"
        cached = self._read_cache(cache_key)
        if cached is not None:
            return cached

        endpoint = f"{self.BASE_URL}/gwasinfo"
        try:
            res = self.session.get(endpoint, timeout=25)
            res.raise_for_status()
            data = res.json()
            matching = []
            filter_lower = trait_filter.lower()
            for gwas_id, meta in data.items():
                trait_name = meta.get("trait", "")
                if not trait_filter or filter_lower in trait_name.lower():
                    matching.append({
                        "id": gwas_id,
                        "trait": trait_name,
                        "sample_size": meta.get("sample_size"),
                        "year": meta.get("year"),
                        "author": meta.get("author"),
                        "consortium": meta.get("consortium")
                    })
            self._write_cache(cache_key, matching)
            return matching
        except requests.RequestException as e:
            print(f"[IEU OpenGWAS] Failed to query dataset list: {e}", file=sys.stderr)
            # Default fallback for demo / offline use
            return [
                {"id": "ieu-b-38", "trait": "Body mass index", "sample_size": 681275, "year": 2018},
                {"id": "ieu-a-2", "trait": "Type 2 diabetes", "sample_size": 149000, "year": 2012}
            ]

    def get_associations(
        self,
        dataset_id: str,
        pval_threshold: float = 5e-8
    ) -> List[Dict[str, Any]]:
        """Retrieves top significant variants from a dataset.

        Args:
            dataset_id: IEU dataset accession (e.g. 'ieu-b-38').
            pval_threshold: Maximum p-value cutoff.

        Returns:
            List[Dict[str, Any]]: Associated variant records.
        """
        cache_key = f"associations_{dataset_id}_{pval_threshold}"
        cached = self._read_cache(cache_key)
        if cached is not None:
            return cached

        endpoint = f"{self.BASE_URL}/tophits"
        payload = {"id": [dataset_id], "pval": pval_threshold}
        try:
            res = self.session.post(endpoint, json=payload, timeout=25)
            res.raise_for_status()
            hits = res.json()
            self._write_cache(cache_key, hits)
            return hits
        except requests.RequestException as e:
            print(f"[IEU OpenGWAS] Top-hits fetch warning for {dataset_id}: {e}", file=sys.stderr)
            return []

    def download_full_sumstats(self, dataset_id: str, out_dir: str = "data/real") -> str:
        """Retrieves and saves top-hits summary dataset locally.

        Args:
            dataset_id: IEU dataset ID.
            out_dir: Base directory.

        Returns:
            str: Path to written file.
        """
        target_dir = os.path.join(out_dir, dataset_id)
        os.makedirs(target_dir, exist_ok=True)
        out_path = os.path.join(target_dir, f"{dataset_id}_tophits.json")

        hits = self.get_associations(dataset_id)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump({
                "dataset_id": dataset_id,
                "variant_count": len(hits),
                "variants": hits
            }, f, indent=2)

        return out_path
