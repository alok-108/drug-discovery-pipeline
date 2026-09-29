"""Unit tests for Open Targets live client and GenT enrichment (mocked via responses)."""

import os
import time
import responses
from target_discovery.opentargets_client import OpenTargetsClient
from target_discovery.gent_wrapper import GenTWrapper
import pandas as pd


MOCK_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"


def test_client_initialization(tmp_path):
    """Verifies client default properties, cache directory, and user-agent."""
    client = OpenTargetsClient(cache_dir=str(tmp_path), max_rate_per_sec=5.0)
    assert client.graphql_endpoint == MOCK_ENDPOINT
    assert os.path.exists(str(tmp_path))
    assert client.min_request_interval == 0.2
    assert "drug-discovery-pipeline" in client.session.headers["User-Agent"]


@responses.activate
def test_get_target_associations_mocked(tmp_path):
    """Verifies target association extraction with mock GraphQL response."""
    mock_payload = {
        "data": {
            "target": {
                "approvedSymbol": "EGFR",
                "associatedDiseases": {
                    "rows": [
                        {
                            "disease": {"id": "EFO_0000616", "name": "non-small cell lung carcinoma"},
                            "score": 0.892,
                            "datatypeScores": [{"id": "genetic_association", "score": 0.95}]
                        },
                        {
                            "disease": {"id": "EFO_0000511", "name": "glioblastoma"},
                            "score": 0.741,
                            "datatypeScores": [{"id": "genetic_association", "score": 0.82}]
                        }
                    ]
                }
            }
        }
    }
    responses.add(responses.POST, MOCK_ENDPOINT, json=mock_payload, status=200)

    client = OpenTargetsClient(cache_dir=str(tmp_path))
    results = client.get_target_associations("ENSG00000146648", limit=2)

    assert len(results) == 2
    assert results[0]["disease_name"] == "non-small cell lung carcinoma"
    assert results[0]["score"] == 0.892
    assert results[0]["datatype_scores"]["genetic_association"] == 0.95


@responses.activate
def test_get_known_drugs_mocked(tmp_path):
    """Verifies known drugs extraction including clinical phase and mechanism of action."""
    mock_payload = {
        "data": {
            "target": {
                "approvedSymbol": "EGFR",
                "knownDrugs": {
                    "uniqueDrugs": 1,
                    "rows": [
                        {
                            "drug": {
                                "id": "CHEMBL1743034",
                                "name": "Osimertinib",
                                "maximumClinicalTrialPhase": 4,
                                "isApproved": True,
                                "hasBeenWithdrawn": False
                            },
                            "mechanismOfAction": "Epidermal growth factor receptor inhibitor",
                            "drugType": "Small molecule"
                        }
                    ]
                }
            }
        }
    }
    responses.add(responses.POST, MOCK_ENDPOINT, json=mock_payload, status=200)

    client = OpenTargetsClient(cache_dir=str(tmp_path))
    drugs = client.get_known_drugs("ENSG00000146648")

    assert len(drugs) == 1
    assert drugs[0]["drug_name"] == "Osimertinib"
    assert drugs[0]["phase"] == 4
    assert drugs[0]["is_approved"] is True
    assert "inhibitor" in drugs[0]["mechanism_of_action"].lower()


@responses.activate
def test_get_target_info_mocked(tmp_path):
    """Verifies target metadata and biological pathway retrieval."""
    mock_payload = {
        "data": {
            "target": {
                "id": "ENSG00000146648",
                "approvedSymbol": "EGFR",
                "approvedName": "epidermal growth factor receptor",
                "biotype": "protein_coding",
                "pathways": [
                    {"pathway": "Signaling by EGFR", "topLevelTerm": "Signal Transduction"}
                ]
            }
        }
    }
    responses.add(responses.POST, MOCK_ENDPOINT, json=mock_payload, status=200)

    client = OpenTargetsClient(cache_dir=str(tmp_path))
    info = client.get_target_info("ENSG00000146648")

    assert info["approved_symbol"] == "EGFR"
    assert info["biotype"] == "protein_coding"
    assert len(info["pathways"]) == 1
    assert info["pathways"][0]["top_level"] == "Signal Transduction"


def test_rate_limiting_timing(tmp_path):
    """Verifies client throttles requests to satisfy max rate per second."""
    client = OpenTargetsClient(cache_dir=str(tmp_path), max_rate_per_sec=20.0)
    t0 = time.time()
    for _ in range(3):
        client._rate_limit()
    elapsed = time.time() - t0
    # 3 calls at 20 req/s => interval 0.05s => minimum ~0.10s
    assert elapsed >= 0.08


@responses.activate
def test_disk_caching(tmp_path):
    """Verifies repeated queries use local disk cache without network requests."""
    mock_payload = {
        "data": {
            "target": {
                "id": "ENSG00000146648",
                "approvedSymbol": "EGFR",
                "approvedName": "EGFR",
                "biotype": "protein_coding",
                "pathways": []
            }
        }
    }
    responses.add(responses.POST, MOCK_ENDPOINT, json=mock_payload, status=200)

    client = OpenTargetsClient(cache_dir=str(tmp_path))
    res1 = client.get_target_info("ENSG00000146648")
    assert res1["approved_symbol"] == "EGFR"
    assert len(responses.calls) == 1

    # Clear memory LRU cache to force disk cache path execution
    OpenTargetsClient.get_target_info.cache_clear()

    # Second call should read from disk cache
    res2 = client.get_target_info("ENSG00000146648")
    assert res2["approved_symbol"] == "EGFR"
    assert len(responses.calls) == 1  # No additional network call was made


@responses.activate
def test_retry_on_429_backoff(tmp_path):
    """Verifies automatic retry with exponential backoff on HTTP 429 rate limit."""
    mock_payload = {
        "data": {
            "target": {
                "approvedSymbol": "EGFR",
                "knownDrugs": {"uniqueDrugs": 0, "rows": []}
            }
        }
    }
    # First response returns 429 Too Many Requests, second succeeds
    responses.add(responses.POST, MOCK_ENDPOINT, status=429)
    responses.add(responses.POST, MOCK_ENDPOINT, json=mock_payload, status=200)

    client = OpenTargetsClient(cache_dir=str(tmp_path))
    drugs = client.get_known_drugs("ENSG00000146648")
    assert drugs == []
    assert len(responses.calls) == 2


@responses.activate
def test_gent_wrapper_enrichment_mocked(tmp_path):
    """Verifies GenTWrapper enriches output with Open Targets associations and drug metrics."""
    assoc_payload = {
        "data": {
            "target": {
                "approvedSymbol": "NTRK1",
                "associatedDiseases": {
                    "rows": [
                        {"disease": {"name": "Alzheimer disease"}, "score": 0.94}
                    ]
                }
            }
        }
    }
    drug_payload = {
        "data": {
            "target": {
                "approvedSymbol": "NTRK1",
                "knownDrugs": {
                    "rows": [
                        {"drug": {"name": "Larotrectinib", "maximumClinicalTrialPhase": 4}}
                    ]
                }
            }
        }
    }
    responses.add(responses.POST, MOCK_ENDPOINT, json=assoc_payload, status=200)
    responses.add(responses.POST, MOCK_ENDPOINT, json=drug_payload, status=200)

    client = OpenTargetsClient(cache_dir=str(tmp_path))
    gent = GenTWrapper()
    gwas_df = pd.DataFrame({
        "GENE": ["NTRK1"],
        "PVALUE": [1e-18]
    })

    enriched = gent.aggregate_gene_pvalues(gwas_df, use_opentargets=True, ot_client=client)
    assert "OT_TOP_DISEASE" in enriched.columns
    assert "OT_ASSOCIATION_SCORE" in enriched.columns
    assert "OT_KNOWN_DRUGS" in enriched.columns
    assert enriched.iloc[0]["OT_TOP_DISEASE"] == "Alzheimer disease"
    assert enriched.iloc[0]["OT_KNOWN_DRUGS"] == 1
