"""Open Targets Platform GraphQL API client.

Provides live integration with the Open Targets Platform (https://platform.opentargets.org/)
for target-disease association scoring, known drugs, and functional annotations.
Includes caching, exponential backoff, and rate-limiting (max 10 req/s).
"""

import os
import json
import time
import hashlib
import functools
from typing import Dict, List, Any, Optional
import requests


class OpenTargetsClient:
    """Client for querying the Open Targets GraphQL API (v4).

    Attributes:
        graphql_endpoint: URL for the Open Targets GraphQL gateway.
        cache_dir: Local filesystem path for cached response JSONs.
        min_request_interval: Minimum time interval between API requests (10 req/sec => 0.1s).
    """

    DEFAULT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(
        self,
        endpoint: Optional[str] = None,
        cache_dir: Optional[str] = None,
        max_rate_per_sec: float = 10.0,
        user_agent: str = "drug-discovery-pipeline/0.2.0 (student-open-source)"
    ) -> None:
        """Initializes OpenTargetsClient with caching and rate-limiting."""
        self.graphql_endpoint = endpoint or self.DEFAULT_ENDPOINT
        self.cache_dir = cache_dir or os.path.expanduser(os.path.join("~", ".cache", "ot"))
        os.makedirs(self.cache_dir, exist_ok=True)
        self.min_request_interval = 1.0 / max_rate_per_sec
        self.last_request_time = 0.0
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": user_agent,
            "Content-Type": "application/json",
            "Accept": "application/json"
        })

    def _rate_limit(self) -> None:
        """Enforces client-side rate limiting to stay below 10 requests/sec."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.min_request_interval:
            time.sleep(self.min_request_interval - elapsed)
        self.last_request_time = time.time()

    def _get_cache_path(self, query: str, variables: Dict[str, Any]) -> str:
        """Computes a deterministic MD5 cache key and filesystem path."""
        payload = json.dumps({"query": query, "variables": variables}, sort_keys=True)
        digest = hashlib.md5(payload.encode("utf-8")).hexdigest()
        return os.path.join(self.cache_dir, f"{digest}.json")

    def _read_cache(self, query: str, variables: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """Reads payload from disk cache if present."""
        path = self._get_cache_path(query, variables)
        if os.path.exists(path):
            try:
                with open(path, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                return None
        return None

    def _write_cache(self, query: str, variables: Dict[str, Any], data: Dict[str, Any]) -> None:
        """Writes response payload to disk cache."""
        path = self._get_cache_path(query, variables)
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f)
        except Exception:
            pass

    def execute_query(
        self,
        query: str,
        variables: Optional[Dict[str, Any]] = None,
        max_retries: int = 4
    ) -> Dict[str, Any]:
        """Executes a GraphQL query against the Open Targets API with retries.

        Args:
            query: GraphQL query string.
            variables: Query variables dictionary.
            max_retries: Max attempts for handling rate limits (429) or transient 5xx errors.

        Returns:
            Dict[str, Any]: Parsed JSON data field.

        Raises:
            RuntimeError: If query execution fails after retries.
        """
        vars_dict = variables or {}
        cached = self._read_cache(query, vars_dict)
        if cached is not None:
            return cached

        self._rate_limit()
        payload = {"query": query, "variables": vars_dict}
        backoff = 1.0

        for attempt in range(max_retries):
            try:
                response = self.session.post(self.graphql_endpoint, json=payload, timeout=20)
                if response.status_code == 429:
                    time.sleep(backoff)
                    backoff *= 2.0
                    continue
                response.raise_for_status()
                res_json = response.json()
                if "errors" in res_json:
                    raise RuntimeError(f"GraphQL Errors: {res_json['errors']}")
                data = res_json.get("data", {})
                self._write_cache(query, vars_dict, data)
                return data
            except (requests.RequestException, RuntimeError) as e:
                if attempt == max_retries - 1:
                    raise RuntimeError(f"Failed to query Open Targets API: {e}") from e
                time.sleep(backoff)
                backoff *= 2.0

        return {}

    @functools.lru_cache(maxsize=128)
    def get_target_associations(self, ensembl_id: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Retrieves prioritized disease associations for a given target.

        Args:
            ensembl_id: Target Ensembl gene identifier (e.g. 'ENSG00000146648').
            limit: Number of top associations to retrieve.

        Returns:
            List[Dict[str, Any]]: List of dicts with disease_name, score, and id.
        """
        query = """
        query TargetAssociations($ensemblId: String!, $size: Int!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            associatedDiseases(page: {index: 0, size: $size}) {
              rows {
                disease {
                  id
                  name
                }
                score
                datatypeScores {
                  id
                  score
                }
              }
            }
          }
        }
        """
        data = self.execute_query(query, {"ensemblId": ensembl_id, "size": limit})
        target_node = data.get("target") or {}
        associations = target_node.get("associatedDiseases", {}).get("rows", [])
        results = []
        for row in associations:
            disease_info = row.get("disease", {})
            results.append({
                "disease_id": disease_info.get("id"),
                "disease_name": disease_info.get("name"),
                "score": round(float(row.get("score", 0.0)), 4),
                "datatype_scores": {
                    d.get("id"): round(float(d.get("score", 0.0)), 3)
                    for d in row.get("datatypeScores", [])
                }
            })
        return results

    @functools.lru_cache(maxsize=128)
    def get_known_drugs(self, ensembl_id: str) -> List[Dict[str, Any]]:
        """Retrieves approved and investigational drugs targeting the specified protein.

        Args:
            ensembl_id: Target Ensembl gene identifier.

        Returns:
            List[Dict[str, Any]]: List of annotated drug molecules and clinical phases.
        """
        query = """
        query KnownDrugs($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            approvedSymbol
            knownDrugs {
              uniqueDrugs
              rows {
                drug {
                  id
                  name
                  maximumClinicalTrialPhase
                  isApproved
                  hasBeenWithdrawn
                }
                mechanismOfAction
                drugType
              }
            }
          }
        }
        """
        data = self.execute_query(query, {"ensemblId": ensembl_id})
        target_node = data.get("target") or {}
        rows = target_node.get("knownDrugs", {}).get("rows", [])
        drugs = []
        for r in rows:
            d = r.get("drug", {})
            drugs.append({
                "drug_id": d.get("id"),
                "drug_name": d.get("name"),
                "phase": d.get("maximumClinicalTrialPhase"),
                "is_approved": bool(d.get("isApproved")),
                "is_withdrawn": bool(d.get("hasBeenWithdrawn")),
                "mechanism_of_action": r.get("mechanismOfAction"),
                "drug_type": r.get("drugType")
            })
        return drugs

    @functools.lru_cache(maxsize=128)
    def get_target_info(self, ensembl_id: str) -> Dict[str, Any]:
        """Fetches basic gene symbol, biotype, description, and pathway annotations.

        Args:
            ensembl_id: Target Ensembl gene identifier.

        Returns:
            Dict[str, Any]: Target metadata summary.
        """
        query = """
        query TargetInfo($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            approvedName
            biotype
            pathways {
              pathway
              topLevelTerm
            }
          }
        }
        """
        data = self.execute_query(query, {"ensemblId": ensembl_id})
        target_node = data.get("target") or {}
        pathways = [
            {"pathway": p.get("pathway"), "top_level": p.get("topLevelTerm")}
            for p in target_node.get("pathways", [])
        ]
        return {
            "ensembl_id": target_node.get("id", ensembl_id),
            "approved_symbol": target_node.get("approvedSymbol", ""),
            "approved_name": target_node.get("approvedName", ""),
            "biotype": target_node.get("biotype", ""),
            "pathways": pathways
        }
