"""Unified Real GWAS Data Downloader CLI.

Downloads real public summary statistics and top-hit associations from
NHGRI-EBI GWAS Catalog or IEU OpenGWAS with resume support and SHA-256 checksums.
"""

import os
import hashlib
import argparse
from tqdm import tqdm

from fetch_gwas_catalog import GWASCatalogClient
from fetch_ieu_opengwas import IEUOpenGWASClient


def compute_sha256(filepath: str) -> str:
    """Computes SHA-256 cryptographic checksum of a file."""
    hasher = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def write_checksum_file(filepath: str, checksum: str) -> str:
    """Writes a companion .sha256 verification file."""
    checksum_path = f"{filepath}.sha256"
    with open(checksum_path, "w", encoding="utf-8") as f:
        f.write(f"{checksum}  {os.path.basename(filepath)}\n")
    return checksum_path


def download_real_gwas_data(
    source: str = "gwas-catalog",
    disease_trait: str = "type 2 diabetes",
    max_studies: int = 5,
    out_dir: str = "data/real"
) -> None:
    """Executes search and download with resume support and progress indicators.

    Args:
        source: 'gwas-catalog' or 'ieu-opengwas'.
        disease_trait: Disease keyword query.
        max_studies: Maximum number of study records to download.
        out_dir: Destination storage directory.
    """
    os.makedirs(out_dir, exist_ok=True)
    print("=" * 70)
    print(f"📥 REAL GWAS DATA INGESTION: Source={source} | Trait='{disease_trait}'")
    print("=" * 70)

    if source.lower() == "gwas-catalog":
        client = GWASCatalogClient()
        print(f"Searching GWAS Catalog for trait '{disease_trait}'...")
        studies = client.list_studies(disease_trait, limit=max_studies)
        if not studies:
            print("No matching studies found.")
            return

        print(f"Discovered {len(studies)} study record(s). Starting download:")
        for s in tqdm(studies, desc="Ingesting GWAS Catalog Studies"):
            acc_id = s["accession_id"]
            expected_json = os.path.join(out_dir, acc_id, f"{acc_id}_associations.json")
            if os.path.exists(expected_json):
                print(f"\n  [SKIPPED] {acc_id} already exists (Resume support).")
                continue

            saved_file = client.download_summary_stats(acc_id, out_dir=out_dir)
            sha = compute_sha256(saved_file)
            write_checksum_file(saved_file, sha)
            print(f"\n  [SAVED] {acc_id} -> {saved_file} (SHA256: {sha[:12]}...)")

    elif source.lower() == "ieu-opengwas":
        client = IEUOpenGWASClient()
        print(f"Searching IEU OpenGWAS for trait filter '{disease_trait}'...")
        datasets = client.list_datasets(trait_filter=disease_trait)
        datasets = datasets[:max_studies]
        if not datasets:
            print("No matching datasets found.")
            return

        print(f"Discovered {len(datasets)} dataset(s). Starting download:")
        for d in tqdm(datasets, desc="Ingesting IEU Datasets"):
            ds_id = d["id"]
            expected_json = os.path.join(out_dir, ds_id, f"{ds_id}_tophits.json")
            if os.path.exists(expected_json):
                print(f"\n  [SKIPPED] {ds_id} already exists (Resume support).")
                continue

            saved_file = client.download_full_sumstats(ds_id, out_dir=out_dir)
            sha = compute_sha256(saved_file)
            write_checksum_file(saved_file, sha)
            print(f"\n  [SAVED] {ds_id} -> {saved_file} (SHA256: {sha[:12]}...)")
    else:
        raise ValueError(f"Unknown source '{source}'. Choose 'gwas-catalog' or 'ieu-opengwas'.")

    print("\n" + "=" * 70)
    print(f"✅ Ingestion completed! Data stored at: {os.path.abspath(out_dir)}")
    print("=" * 70)


def main() -> None:
    """CLI parser entrypoint."""
    parser = argparse.ArgumentParser(description="Download real GWAS datasets from free public repositories.")
    parser.add_argument("--source", type=str, default="gwas-catalog", choices=["gwas-catalog", "ieu-opengwas"], help="Public repository source")
    parser.add_argument("--disease", type=str, default="type 2 diabetes", help="Disease trait query keyword")
    parser.add_argument("--max-studies", type=int, default=3, help="Maximum number of studies/datasets to fetch")
    parser.add_argument("--out", type=str, default="data/real", help="Destination storage directory")
    args = parser.parse_args()

    download_real_gwas_data(
        source=args.source,
        disease_trait=args.disease,
        max_studies=args.max_studies,
        out_dir=args.out
    )


if __name__ == "__main__":
    main()
