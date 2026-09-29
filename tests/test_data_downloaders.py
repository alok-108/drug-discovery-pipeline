"""Unit tests for GWAS Catalog and IEU OpenGWAS data downloaders (mocked)."""

import os
import sys
import responses

# Ensure scripts/ is importable
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
scripts_dir = os.path.join(base_dir, "scripts")
if scripts_dir not in sys.path:
    sys.path.insert(0, scripts_dir)

from fetch_gwas_catalog import GWASCatalogClient
from fetch_ieu_opengwas import IEUOpenGWASClient
from download_real_data import compute_sha256, write_checksum_file, download_real_gwas_data


@responses.activate
def test_gwas_catalog_list_studies_mocked():
    """Verifies GWAS Catalog study listing with mocked REST response."""
    endpoint = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait"
    mock_payload = {
        "_embedded": {
            "studies": [
                {
                    "accessionId": "GCST006867",
                    "studyTitle": "Genome-wide meta-analysis of type 2 diabetes",
                    "diseaseTrait": {"trait": "type 2 diabetes"},
                    "publicationInfo": {
                        "pubmedId": "30297969",
                        "author": {"fullname": "Mahajan A"}
                    },
                    "fullPvalueSet": True
                }
            ]
        }
    }
    responses.add(responses.GET, endpoint, json=mock_payload, status=200)

    client = GWASCatalogClient()
    studies = client.list_studies("type 2 diabetes", limit=1)

    assert len(studies) == 1
    assert studies[0]["accession_id"] == "GCST006867"
    assert studies[0]["author"] == "Mahajan A"
    assert studies[0]["full_pvalue_set"] is True


@responses.activate
def test_gwas_catalog_fetch_associations_mocked(tmp_path):
    """Verifies fetching lead associations for a GWAS Catalog study accession."""
    endpoint = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST006867/associations"
    mock_payload = {
        "_embedded": {
            "associations": [
                {
                    "pvalue": 1.2e-45,
                    "betaNum": 0.18,
                    "standardError": 0.02,
                    "loci": [
                        {"strongestRiskAlleles": [{"riskAlleleName": "rs7903146-T"}]}
                    ]
                }
            ]
        }
    }
    responses.add(responses.GET, endpoint, json=mock_payload, status=200)

    client = GWASCatalogClient()
    assocs = client.fetch_associations("GCST006867")
    assert len(assocs) == 1
    assert assocs[0]["lead_variant"] == "rs7903146-T"

    # Test file saving
    out_file = client.download_summary_stats("GCST006867", out_dir=str(tmp_path))
    assert os.path.exists(out_file)


@responses.activate
def test_ieu_opengwas_list_datasets_mocked(tmp_path):
    """Verifies IEU OpenGWAS dataset catalog listing."""
    endpoint = "https://gwas-api.mrcieu.ac.uk/api/gwasinfo"
    mock_payload = {
        "ieu-b-38": {
            "trait": "Body mass index",
            "sample_size": 681275,
            "year": 2018,
            "author": "Yengo L"
        },
        "ieu-a-2": {
            "trait": "Type 2 diabetes",
            "sample_size": 149000,
            "year": 2012,
            "author": "Morris AP"
        }
    }
    responses.add(responses.GET, endpoint, json=mock_payload, status=200)

    client = IEUOpenGWASClient(cache_dir=str(tmp_path))
    datasets = client.list_datasets(trait_filter="Body mass index")

    assert len(datasets) == 1
    assert datasets[0]["id"] == "ieu-b-38"
    assert datasets[0]["sample_size"] == 681275


@responses.activate
def test_ieu_opengwas_get_associations_mocked(tmp_path):
    """Verifies IEU OpenGWAS top-hits extraction and caching."""
    endpoint = "https://gwas-api.mrcieu.ac.uk/api/tophits"
    mock_payload = [
        {"rsid": "rs1558902", "p": 2.1e-120, "beta": 0.08, "se": 0.003}
    ]
    responses.add(responses.POST, endpoint, json=mock_payload, status=200)

    client = IEUOpenGWASClient(cache_dir=str(tmp_path))
    hits = client.get_associations("ieu-b-38", pval_threshold=5e-8)

    assert len(hits) == 1
    assert hits[0]["rsid"] == "rs1558902"

    saved = client.download_full_sumstats("ieu-b-38", out_dir=str(tmp_path))
    assert os.path.exists(saved)


def test_checksum_computation(tmp_path):
    """Verifies SHA-256 calculation and companion checksum file generation."""
    sample_file = tmp_path / "sample.tsv"
    sample_file.write_text("SNP\tPVALUE\nrs1\t1e-8\n", encoding="utf-8")

    sha = compute_sha256(str(sample_file))
    assert len(sha) == 64  # SHA-256 hex length

    chk_path = write_checksum_file(str(sample_file), sha)
    assert os.path.exists(chk_path)
    with open(chk_path, "r", encoding="utf-8") as f:
        content = f.read()
    assert sha in content
    assert "sample.tsv" in content


@responses.activate
def test_download_real_gwas_data_cli_mocked(tmp_path):
    """Verifies end-to-end execution of unified download_real_gwas_data CLI function."""
    endpoint_studies = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait"
    endpoint_assocs = "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST006867/associations"

    responses.add(responses.GET, endpoint_studies, json={
        "_embedded": {"studies": [{"accessionId": "GCST006867", "diseaseTrait": {"trait": "T2D"}}]}
    }, status=200)
    responses.add(responses.GET, endpoint_assocs, json={
        "_embedded": {"associations": [{"pvalue": 1e-10, "loci": []}]}
    }, status=200)

    download_real_gwas_data(
        source="gwas-catalog",
        disease_trait="T2D",
        max_studies=1,
        out_dir=str(tmp_path)
    )

    expected_file = os.path.join(str(tmp_path), "GCST006867", "GCST006867_associations.json")
    assert os.path.exists(expected_file)
    assert os.path.exists(f"{expected_file}.sha256")
