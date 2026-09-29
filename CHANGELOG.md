# Changelog

All notable changes to the **AI-Driven Drug Discovery Pipeline** are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.2.0] - 2026-09-29

### Added
- **Open Targets Platform Live Integration**:
  - `src/target_discovery/opentargets_client.py`: High-performance GraphQL client querying target-disease associations, approved clinical drugs, and biological pathways.
  - Multi-tier caching architecture (`functools.lru_cache` + persistent disk cache at `~/.cache/ot/`) and client-side rate-limiting ($\le 10\text{ req/s}$) with exponential backoff on HTTP 429.
  - `--use-opentargets` enrichment flag in `GenTWrapper.aggregate_gene_pvalues()`.
  - Comprehensive documentation in `docs/08_live_integrations.md`.
- **Quantum Chemistry Module (Qiskit Nature & VQE)**:
  - `src/quantum/molecular_hamiltonian.py`: Second-quantized electronic Hamiltonians, Coulomb nuclear repulsion, and Jordan-Wigner qubit transformations.
  - `src/quantum/qiskit_runner.py`: Variational Quantum Eigensolver (`QuantumMoleculeSolver`) solving $\text{H}_2$ and $\text{LiH}$ ground-state energies with chemical accuracy ($<1.0\text{ kcal/mol}$).
  - Support for zero-cost local CPU simulation (AerSimulator) and optional IBM Quantum cloud hardware.
  - Automated Potential Energy Surface (PES) generator (`run_all_demos()`) producing `docs/assets/h2_pes.png`.
  - Comprehensive documentation in `docs/09_quantum_chemistry.md`.
- **Production Containerization & GHCR Deployment**:
  - Multi-stage slim `Dockerfile` reducing container footprint by $>65\%$.
  - Security hardening with unprivileged non-root execution (`appuser`, UID 1000) and automated healthcheck.
  - `docker-compose.yml` orchestrating Streamlit dashboard (port 8501) and JupyterLab (port 8888).
  - GitHub Actions workflow (`.github/workflows/docker-publish.yml`) publishing multi-arch images (`linux/amd64`, `linux/arm64`) to `ghcr.io/alok-108/drug-discovery-pipeline`.
  - Comprehensive documentation in `docs/10_docker.md`.
- **Interactive Streamlit Web Dashboard**:
  - `app/streamlit_app.py`: Modern dark-theme web interface with pipeline overview and sitemap.
  - `app/pages/1_GWAS_Explorer.py`: File uploader, interactive Plotly Manhattan plots, GenT ranking, and live Open Targets enrichment.
  - `app/pages/2_Protein_Viewer.py`: AlphaFold 3D structure access, per-residue pLDDT profiles, and druggability indicators.
  - `app/pages/3_Drug_Repurposing.py`: Query known approved and investigational drugs, clinical trial phases, and disease association bar charts.
  - Deployment instructions for Streamlit Community Cloud (free zero-credit-card hosting) in `docs/11_streamlit_deploy.md`.
- **Real-World Genomic Data Ingestion CLI**:
  - `scripts/fetch_gwas_catalog.py`: Search and download curated studies from NHGRI-EBI GWAS Catalog REST API.
  - `scripts/fetch_ieu_opengwas.py`: Query MRC IEU OpenGWAS database with local disk caching.
  - `scripts/download_real_data.py`: Unified CLI downloader with `tqdm` progress bars, resume support, and companion SHA-256 checksums.
  - Documentation and FAIR open data guidelines in `docs/12_real_data.md`.
- **Testing & Quality Assurance**:
  - Added 19 new unit tests across Open Targets (`test_opentargets.py`), Quantum VQE (`test_qiskit_runner.py`), and real data downloaders (`test_data_downloaders.py`). Total test suite now stands at 30 passing tests.

---

## [0.1.0] - 2026-09-29

### Added
- Initial public release of the AI-Driven Drug Discovery Pipeline.
- Complete 5-stage documentation guides (`docs/01_target_discovery.md` through `docs/07_tool_comparison.md`, `faq.md`).
- Core computational genetics modules (`GWASParser`, `GenTWrapper`, `XQTLIntegrator`).
- Functional genomics and perturbomics pocket discovery engines.
- Molecular design interfaces (`AlphaFoldRunner`, `BindCraftWrapper`, `PeptideGenerator`).
- CausalPathKG knowledge graph traversal engine.
- Interactive Jupyter tutorial notebooks 01 through 04.
