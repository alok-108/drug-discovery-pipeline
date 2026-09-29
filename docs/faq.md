# Frequently Asked Questions (FAQ) & Troubleshooting Guide

This FAQ addresses common scientific, computational, and regulatory questions encountered when running the AI-Driven Drug Discovery Pipeline.

---

## Table of Contents
1. [Target Discovery & Statistical Genetics](#1-target-discovery--statistical-genetics)
2. [Functional Genomics & CRISPR Tiling](#2-functional-genomics--crispr-tiling)
3. [Protein & Peptide Generative Design](#3-protein--peptide-generative-design)
4. [Drug Repurposing & Knowledge Graphs](#4-drug-repurposing--knowledge-graphs)
5. [Hardware, Licensing & Deployment](#5-hardware-licensing--deployment)
6. [Regulatory & Clinical Translation](#6-regulatory--clinical-translation)

---

## 1. Target Discovery & Statistical Genetics

### Q1.1: Why does the pipeline use GenT / Cauchy Combination rather than traditional Bonferroni or VEGAS gene tests?
**Answer**: Single-variant tests require testing over one million SNPs, necessitating an extreme multiple testing correction ($\alpha = 5 \times 10^{-8}$) that overlooks moderate cumulative polygenic associations. Traditional gene-based tests like VEGAS require simulating billions of multivariate normal draws to account for Linkage Disequilibrium (LD) correlation matrices, which is computationally prohibitive at scale. 

**GenT** applies the **Cauchy Combination Test (CCT)**. Because the transformation $T = \sum w_i \tan((0.5 - p_i)\pi)$ possesses heavy Cauchy tails, the analytical tail probability $\mathbb{P}(T > t) \approx 1/\pi t$ remains exact and distribution-free under **arbitrary correlation structures** between variants without requiring permutation or numerical covariance matrix inversions.

### Q1.2: What if my GWAS cohort has only European ancestry data? Can I still run the pipeline?
**Answer**: Yes, the pipeline functions perfectly with single-ancestry summary statistics using the standard `GenTWrapper.aggregate_gene_pvalues()` engine. However, we strongly recommend incorporating multi-ancestry summary statistics (e.g. from the Global Biobank Meta-analysis Initiative [GBMI] or Million Veteran Program [MVP]) using `MuGenT`. Multi-ancestry data significantly sharpens fine-mapping resolution by exploiting the shorter LD blocks of African and Admixed American cohorts, shrinking credible causal sets from dozens of SNPs down to single causal variants.

### Q1.3: What is the difference between eQTL colocalization and Mendelian Randomization (MR)?
**Answer**: 
* **Colocalization (`coloc`)** asks: *"Do the disease trait and the expression QTL share the exact same causal variant at this single genomic locus?"* It eliminates spurious correlations caused by distinct variants in mutual LD.
* **Mendelian Randomization (MR)** asks: *"Does genetically predicted lifelong exposure to altered protein/mRNA levels cause a change in disease risk?"* MR utilizes genetic variants as instrumental variables to deduce causal directionality across the whole genome, guarding against reverse causation and confounding.
Our pipeline uses `coloc` first to confirm target locus identity, followed by two-sample MR to estimate the anticipated clinical effect size.

---

## 2. Functional Genomics & CRISPR Tiling

### Q2.1: Why should I prefer CRISPRi over standard CRISPR-Cas9 knockout for drug target validation?
**Answer**: CRISPR-Cas9 introduces blunt double-strand breaks (DSBs) that trigger TP53-mediated DNA damage responses, cell cycle arrest, and non-specific toxicity—particularly in amplified genomic regions (the copy-number artifact). Furthermore, genetic knockout completely deletes the protein scaffold, whereas small molecules and biologics typically inhibit catalytic activity or allosterically modulate a pocket while leaving the structural scaffolding intact. 

**CRISPRi** (catalytically dead Cas9 fused to KRAB) produces titratable, reversible epigenetic transcriptional silencing without cutting DNA, mimicking the dose-dependent pharmacodynamics of a pharmaceutical inhibitor.

### Q2.2: What is "Perturbomics" and how does it find allosteric pockets?
**Answer**: Traditional pooled screens score whole genes as either essential or non-essential. **Perturbomics** combines high-density base-editing tiling screens (e.g., cytidine or adenine base editors scanning every codon of an open reading frame) with downstream structural mapping. By measuring cell depletion caused by specific amino acid substitutions, Perturbomics uncovers localized clusters of residues whose mutation destroys protein activity. Mapping these vulnerable residues onto the AlphaFold 3D model reveals novel, druggable allosteric and active-site clefts at 1-amino-acid resolution.

---

## 3. Protein & Peptide Generative Design

### Q3.1: How does BindCraft achieve 10% to 100% experimental validation rates without laboratory directed evolution?
**Answer**: Conventional computational design protocols often generate non-physical binders that aggregate or misfold when synthesized in the lab. **BindCraft** succeeds by applying an automated multi-stage computational filter cascade:
1. It seeds backbones on explicit target functional hotspots using RFdiffusion.
2. It designs robust, soluble sequences using ProteinMPNN.
3. It refolds the predicted complex using AlphaFold-Multimer.
4. It strictly filters out any design that fails physical biophysical thresholds: binding energy ($\Delta G < -35\text{ kcal/mol}$), buried interface surface area ($\Delta\text{SASA} > 1,300\text{ \AA}^2$), and interface predicted aligned error ($\text{iPAE} < 8.0\text{ \AA}$). Only designs exhibiting flawless atomic complementarity are advanced to physical synthesis.

### Q3.2: Why design cyclic macrocycles or lasso peptides instead of standard linear peptides?
**Answer**: Linear peptides are rapidly degraded by serum aminopeptidases and endopeptidases (plasma half-life often $<5\text{ minutes}$) and possess high conformational flexibility that incurs large entropic penalties upon binding. **Cyclic macrocycles** (RFpeptides) and **lasso peptides** (LassoDesign) pre-organize the peptide backbone into a rigid, active 3D conformation, shield backbone amides to cross cellular membranes, and resist proteolytic degradation, enabling systemic stability and oral bioavailability.

---

## 4. Drug Repurposing & Knowledge Graphs

### Q4.1: How does CausalPathKG differ from standard network medicine algorithms?
**Answer**: Standard network medicine tools (like STRING or simple Protein-Protein Interaction networks) construct undirected graphs where an edge simply means "Protein A binds Protein B". This fails to distinguish between an activating kinase cascade and an inhibitory degradation pathway. 

**CausalPathKG** maintains strict **causal directionality and sign consistency** (activation vs inhibition). If disease genetics shows that overexpression of Target T drives pathology, CausalPathKG searches for causal paths terminating in directed drug inhibition ($\text{Sign} = -1$), eliminating contradictory false-positive repurposing candidates.

### Q4.2: How do you verify whether a repurposed candidate will work in humans at safe dosages?
**Answer**: Computational matches often identify potent cellular inhibitors whose required half-maximal inhibitory concentration ($\text{IC}_{50}$) is $10\text{ }\mu\text{M}$, but whose safe maximum human plasma concentration ($C_{\text{max}}$) is only $0.5\text{ }\mu\text{M}$. The pipeline cross-references all computational repurposing hits against pharmacokinetic databases (FDA package inserts, PubChem, DailyMed) to ensure that the required therapeutic concentration lies safely below the approved human clinical $C_{\text{max}}$.

---

## 5. Hardware, Licensing & Deployment

### Q5.1: What are the licensing constraints around AlphaFold 3 for commercial drug discovery?
**Answer**: The open model parameters and code for **AlphaFold 3** released by DeepMind / Isomorphic Labs are licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC-BY-NC-ND 4.0)** license, which restricts direct commercial use without a commercial license. For commercial pharmaceutical pipelines, researchers can:
1. License proprietary commercial platforms from Isomorphic Labs or authorized partners.
2. Utilize fully open commercial alternatives such as **ColabFold / AlphaFold 2 (Apache 2.0)**, **ESMFold (Apache 2.0)**, or **OpenFold (Apache 2.0)**.
Our repository provides modular runners compatible with both AlphaFold3 and Apache-licensed alternatives.

### Q5.2: Do I need a multi-GPU cluster to run this repository?
**Answer**: No. The target discovery (GWAS parsing, GenT, MuGenT, coloc), functional genomics screen analysis (MAGeCK, Perturbomics), and repurposing modules (CausalPathKG) execute efficiently on standard multi-core CPUs in seconds. For structural prediction and generative binder design, local runs require an NVIDIA GPU with $\ge 16\text{ GB}$ VRAM (e.g. RTX 3090, 4090, A10G) or free/low-cost cloud instances (Google Colab, Modal, AWS EC2 G5).

---

## 6. Regulatory & Clinical Translation

### Q6.1: Does the US FDA accept AI-generated computational data in an Investigational New Drug (IND) application?
**Answer**: The FDA actively encourages AI and machine learning under its *Modernizing Preclinical Drug Development* guidelines and the FDA Modernization Act 2.0. However, computational predictions serve as supporting evidence for biological rationale and mechanism of action; they do not replace mandatory in-vitro target engagement, GLP animal toxicology (NOAEL determination), safety pharmacology, and Chemistry, Manufacturing, and Controls (CMC) specifications.

### Q6.2: What is the typical development timeline when utilizing this computational pipeline?
```
+-----------------------------------------------------------------------------------+
|               DRUG DISCOVERY TIMELINE COMPARISON                                  |
+-----------------------------------------------------------------------------------+
| Traditional Paradigm:                                                             |
| Target Discovery -> High-Throughput Screen -> Lead Opt -> Preclinical -> Ph I-III|
| [============================= 10 to 15 YEARS ==============================]    |
|                                                                                   |
| AI De Novo Pipeline:                                                              |
| Genetics Target -> AlphaFold / BindCraft -> In-Vitro -> Preclinical -> Ph I-III   |
| [==================== 5 to 7 YEARS ====================]                          |
|                                                                                   |
| AI Repurposing Track:                                                             |
| CausalPathKG + PRISM Viability -> Preclinical Efficacy -> Direct Phase II Trial   |
| [============ 2 to 4 YEARS ============]                                          |
+-----------------------------------------------------------------------------------+
```

---

## 7. Key References

1. **U.S. Food and Drug Administration (FDA)** (2023). *Using Artificial Intelligence and Machine Learning in the Development of Drug and Biological Products*. Discussion Paper.
2. **Kursawe, J., et al.** (2020). GenT: Gene-based testing of complex traits using high-dimensional summary statistics. *Bioinformatics*, 36(18), 4752-4759.
3. **Pellock, S. J., et al.** (2024). BindCraft: Automated de novo design of targeted protein binders. *bioRxiv*, doi:10.1101/2024.03.01.583000.
4. **Pushpakom, S., et al.** (2019). Drug repurposing: progress, challenges and recommendations. *Nature Reviews Drug Discovery*, 18(1), 41-58.
5. **Nelson, M. R., et al.** (2015). The support of human genetic evidence for approved drug indications. *Nature Genetics*, 47(8), 856-860.
