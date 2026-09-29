# Contributing to AI-Driven Drug Discovery Pipeline

Thank you for your interest in contributing to the **AI-Driven Drug Discovery Pipeline**! We welcome researchers, computational biologists, pharmacologists, bioinformaticians, and software engineers from all backgrounds.

## Code of Conduct

We are committed to providing a welcoming, inclusive, and harassment-free environment for all participants. Please treat all contributors with respect and professionalism.

## Getting Started

1. **Fork the Repository**: Click the `Fork` button on GitHub.
2. **Clone your fork**:
   ```bash
   git clone https://github.com/<your-username>/drug-discovery-pipeline.git
   cd drug-discovery-pipeline
   ```
3. **Set up a development environment**:
   ```bash
   conda env create -f environment.yml
   conda activate drug-discovery
   pip install -e ".[dev]"
   ```
4. **Create a topic branch**:
   ```bash
   git checkout -b feat/add-new-target-module
   ```

## Development Guidelines

### Coding Standards
- **Style**: Adhere to PEP8 conventions. Run `ruff check .` before pushing.
- **Type Annotations**: All public functions and classes must include Python type hints (`typing` module).
- **Documentation**: Use Google-style docstrings for every module, class, and method.
- **Testing**: Add unit tests in the `tests/` directory for any new logic. Verify that all tests pass:
  ```bash
  pytest tests/ -v
  ```

### Conventional Commits
Please format commit messages following [Conventional Commits](https://www.conventionalcommits.org/):
- `feat: add multi-ancestry meta-analysis in MuGenT wrapper`
- `fix: resolve zero-division edge case in p-value log transformation`
- `docs: update documentation for BindCraft de novo design parameters`
- `test: add unit tests for perturbomics pocket scoring`

### Pull Request Process
1. Ensure all CI checks (pytest, linting) pass.
2. Update relevant documentation in `docs/` or `README.md` if introducing new features.
3. Submit a Pull Request targeting the `main` branch with a clear description of the problem solved and test evidence.
4. Maintainers will review your PR and provide constructive feedback.

## Questions & Discussions
For bugs and feature requests, please open an Issue on GitHub. For scientific discussions and method proposals, feel free to start a GitHub Discussion!
