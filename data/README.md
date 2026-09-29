# Data Directory

This directory contains benchmark and demonstration datasets for the AI-Driven Drug Discovery Pipeline.

## Directory Structure

```
data/
├── README.md              # Dataset descriptions, schemas, and sources
└── sample/
    └── example_gwas.tsv   # Mock GWAS summary statistics across cardiometabolic and neurological targets
```

## Schema: `example_gwas.tsv`

| Column | Type | Description |
|--------|------|-------------|
| `SNP` | String | Variant identifier (dbSNP rsID or chr:pos:ref:alt) |
| `CHR` | Integer/String | Chromosome (1-22, X, Y) |
| `BP` | Integer | Base-pair position (GRCh38 coordinate system) |
| `REF` | String | Reference allele |
| `ALT` | String | Effect/alternative allele |
| `BETA` | Float | Effect size from additive genetic association model |
| `SE` | Float | Standard error of the beta estimate |
| `PVALUE` | Float | Association p-value from regression test |
| `AF` | Float | Effect allele frequency in the studied population |
| `GENE` | String | Proximal or fine-mapped candidate causal gene |

## Recommended External Public Datasets for Production Runs

1. **GWAS Catalog**: [https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/) - Harmonized summary statistics.
2. **eQTLGen Consortium**: [https://www.eqtlgen.org/](https://www.eqtlgen.org/) - Blood cis/trans eQTLs for >30,000 individuals.
3. **UK Biobank Pharma Proteomics Project (UKB-PPP)**: [https://www.ukbiobank.ac.uk/](https://www.ukbiobank.ac.uk/) - Plasma pQTLs for ~3,000 proteins in 50,000+ individuals.
4. **DepMap / Project Achilles**: [https://depmap.org/portal/](https://depmap.org/portal/) - Genome-scale CRISPR-Cas9 essentiality and dependency screens across >1,000 cancer cell lines.
5. **AlphaFold Protein Structure Database**: [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/) - High-accuracy 3D coordinates for human proteome.
