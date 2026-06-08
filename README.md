# Control Group Selection Affects Biomarker Discovery in Colorectal Cancer

A Python-based exploratory computational biology project investigating how the choice of normal control tissue affects machine learning-based biomarker discovery in colon adenocarcinoma (COAD).

This project was developed as an independent learning project to practice Python, GitHub, Claude Code, public biomedical data analysis, and statistical methods used in research-oriented engineering work.

---

## Research Question

Does changing the definition of "normal" tissue change which genes are selected as colorectal cancer biomarker candidates?

More specifically, this project compares TCGA COAD tumor samples against several different normal reference groups:

1. TCGA tumor-adjacent normal tissue (NAT)
2. GTEx colon transverse tissue
3. GTEx colon sigmoid tissue
4. Combined GTEx colon tissue

The goal is to test whether machine learning-selected biomarker candidates remain stable when the normal control group changes.

---

## Motivation

Many cancer biomarker discovery studies compare tumor tissue against a normal reference group. However, "normal" can be defined in multiple ways.

For colorectal cancer, possible normal controls include:

- Tumor-adjacent normal tissue from cancer patients
- Healthy colon tissue from non-cancer donors
- Colon tissue from different anatomical regions, such as transverse and sigmoid colon

These groups are not necessarily equivalent. Tumor-adjacent normal tissue may already be molecularly affected by the tumor environment, while GTEx healthy tissue comes from a different data source with potential batch effects.

This project explores whether these choices change the genes that appear important for cancer classification.

---

## Project Summary

The analysis pipeline includes:

1. Data extraction from public TCGA, GTEx, and GEO datasets
2. Quality control using PCA
3. Batch-effect evaluation and correction using neuroCombat
4. Differential expression analysis
5. Machine learning-based feature selection
6. Biomarker stability comparison across control groups
7. Biological validation using literature review, pathway analysis, survival analysis, and an external GEO cohort

---

## Data Sources

| Dataset | Description | Sample Count |
|---|---:|---:|
| TCGA COAD Tumor | Colon adenocarcinoma tumor tissue | 288 |
| TCGA COAD NAT | Tumor-adjacent normal tissue | 41 |
| GTEx Colon Transverse | Non-cancer transverse colon tissue | 167 |
| GTEx Colon Sigmoid | Non-cancer sigmoid colon tissue | 141 |
| GSE156451 | Independent paired colorectal cancer cohort | 72 paired samples |

Main TCGA and GTEx data were obtained from UCSC Xena. External validation was performed using GEO dataset GSE156451.

Raw and processed expression matrices are not included in this repository because of file size and data redistribution limitations. See `data/README.md` for data source notes.

---

## Methods Overview

### 1. Quality Control

PCA was used to visualize sample-level structure across tumor, NAT, and GTEx groups.

The first PCA plots showed clear separation between TCGA-derived samples and GTEx-derived samples. This suggested a strong batch effect caused by differences in data source, sample collection, or processing pipeline.

To evaluate this further, PCA was performed in several ways:

- All groups together
- TCGA tumor vs TCGA NAT only
- NAT vs GTEx normal tissues only
- Before and after ComBat correction

GTEx-containing comparisons were batch-corrected using `neuroCombat`.

Important note: TCGA tumor vs TCGA NAT comparisons used the original expression data because both groups came from the same source. GTEx-based comparisons used batch-corrected GTEx data to reduce technical differences between TCGA and GTEx.

---

### 2. Differential Expression Analysis

Differential expression analysis was performed separately for each tumor-control comparison.

Comparisons:

1. Tumor vs TCGA NAT
2. Tumor vs GTEx All
3. Tumor vs GTEx Colon Transverse
4. Tumor vs GTEx Colon Sigmoid

Statistical approach:

- Welch's t-test
- Benjamini-Hochberg FDR correction
- DEG threshold: `|log2FC| > 1` and `FDR < 0.05`

DEG results were visualized using volcano plots.

---

### 3. Machine Learning Feature Selection

For each tumor-control comparison, the top differentially expressed genes were used as candidate features for machine learning.

Models used:

- LASSO logistic regression
- Random Forest

The purpose was not to build a clinical diagnostic model, but to identify which genes were selected as important features under different definitions of normal tissue.

Because initial models produced AUC = 1.000, an additional corrected/nested validation step was added to reduce the risk of feature-selection leakage. The corrected models still showed very high classification performance, suggesting that the tumor-normal signal in this dataset is strong. However, this result should still be interpreted cautiously because the analysis uses highly separable public transcriptomic data.

---

### 4. Biomarker Stability Analysis

Biomarker candidate lists from the four tumor-control comparisons were compared.

The project evaluated:

- Whether the same genes appeared across all control groups
- Whether NAT-based candidates overlapped with GTEx-based candidates
- Whether GTEx transverse and sigmoid comparisons produced similar candidates
- How many comparisons each candidate appeared in

A stability score was assigned to each gene based on the number of comparisons in which it appeared.

---

### 5. Biological Validation

Candidate genes were further evaluated using:

- Gene symbol conversion
- Literature review
- Pathway analysis
- Survival analysis
- External validation in GSE156451

Two initial candidates were removed after further review because they appeared likely to be artifacts or weak candidates:

- A retired or invalid Ensembl ID
- A pseudogene-like candidate with inconsistent behavior across NAT and GTEx comparisons

The final retained biomarker candidates were:

- `CDH3`
- `CLDN1`
- `ETV4`
- `KRT80`

---

## Key Results

### DEG Counts

| Comparison | DEG Count |
|---|---:|
| Tumor vs NAT | 11,623 |
| Tumor vs GTEx All | 17,562 |
| Tumor vs GTEx Transverse | 16,516 |
| Tumor vs GTEx Sigmoid | 19,298 |

GTEx-based comparisons produced more differentially expressed genes than the NAT-based comparison. This suggests that the definition of normal control tissue can strongly affect downstream gene selection.

---

### Biomarker Stability

The biomarker overlap analysis showed that candidate genes were not fully stable across control definitions.

Main observations:

- No gene appeared across all four final candidate lists under the strict consensus criteria.
- NAT-based candidates overlapped weakly with GTEx-based candidates.
- GTEx-based comparisons showed stronger overlap with each other than with NAT.

This supports the central idea of the project: biomarker discovery results can change depending on how "normal" tissue is defined.

---

### Final Candidate Genes

| Gene | General Biological Role | Interpretation |
|---|---|---|
| CDH3 | Cell adhesion | Previously associated with colorectal cancer and epithelial tumor biology |
| CLDN1 | Tight junction / epithelial barrier | Strongly linked to colorectal cancer progression in prior studies |
| ETV4 | Transcription factor | Associated with tumor progression and cancer-related transcriptional regulation |
| KRT80 | Cytoskeleton / epithelial structure | Reported in cancer invasion and metastasis-related contexts |

These genes were retained because they were supported by stability analysis, expression patterns, and biological literature review.

---

### External Validation

The four final candidates were tested in an independent paired colorectal cancer dataset, GSE156451.

| Gene | log2FC | p-value |
|---|---:|---:|
| CDH3 | +3.16 | 1.1 × 10^-36 |
| CLDN1 | +3.35 | 7.8 × 10^-35 |
| ETV4 | +4.08 | 2.8 × 10^-44 |
| KRT80 | +2.61 | 1.5 × 10^-31 |

All four genes showed increased expression in tumor samples in the external dataset.

---

## Repository Structure

```text
COAD_control_stability_project/
├── README.md
├── work_log.md
├── notes.md
├── project_plan.md
├── data/
│   └── README.md
├── scripts/
│   ├── extract_coad_matrix.py
│   ├── 02_qc_pca.py
│   ├── 03_qc_pca_tcga_only.py
│   ├── 04_qc_pca_nat_vs_gtex.py
│   ├── 05_combat_correction.py
│   ├── deg_analysis.py
│   ├── ml_analysis.py
│   ├── ml_corrected.py
│   ├── stability_analysis.py
│   ├── biological_validation.py
│   └── external_validation.py
└── results/
    ├── qc/
    ├── deg/
    ├── ml/
    ├── ml_corrected/
    ├── stability/
    ├── biological_validation/
    ├── external_validation/
    └── report_Final.pages
Scripts
Script	Purpose
extract_coad_matrix.py	Extracts COAD-relevant expression samples from the larger public dataset
02_qc_pca.py	Runs PCA for overall quality control
03_qc_pca_tcga_only.py	Runs PCA for TCGA tumor vs NAT samples only
04_qc_pca_nat_vs_gtex.py	Runs PCA for NAT vs GTEx normal tissues
05_combat_correction.py	Applies neuroCombat batch correction
deg_analysis.py	Performs differential expression analysis and volcano plotting
ml_analysis.py	Runs initial ML feature selection
ml_corrected.py	Runs corrected validation to reduce leakage risk
stability_analysis.py	Compares biomarker candidates across control groups
biological_validation.py	Performs gene annotation, literature-based review, and biological validation
external_validation.py	Tests final candidates in an independent GEO dataset
Tools and Libraries
Main tools used:
Python
pandas
numpy
scipy
scikit-learn
matplotlib
seaborn
neuroCombat
lifelines
mygene
Claude Code and ChatGPT were used as coding and learning assistants throughout the project.
Important Limitations
This project is exploratory and should not be interpreted as a clinical diagnostic study.
Main limitations:
The project uses public retrospective datasets rather than newly generated experimental data.
TCGA and GTEx come from different data sources, so batch effects are a major concern.
ComBat correction reduces technical effects but cannot guarantee perfect removal of all dataset-specific bias.
The NAT group is much smaller than the GTEx groups.
The machine learning models classify highly separable tumor-normal transcriptomic samples and may overestimate real-world clinical performance.
Survival analysis was limited and did not include full multivariate clinical adjustment.
Final biomarker candidates require further validation in larger independent cohorts and experimental settings.
Learning Notes
This project was also a learning exercise.
At the beginning, many of the core concepts were new to me, including:
RNA-seq expression matrices
TCGA and GTEx data structure
PCA
Batch effects
ComBat correction
Differential expression analysis
Volcano plots
LASSO
Random Forest
AUC
Biomarker stability
The work log documents the learning process, including mistakes, confusion, AI-assisted coding, and the steps I took to understand the results.
See work_log.md for a detailed record of time spent and what I learned.
Use of AI Tools
This project used AI tools openly and intentionally.
Claude Code and ChatGPT were used to:
Help write and debug Python code
Explain unfamiliar statistical and biological concepts
Assist in organizing the analysis pipeline
Help interpret output files and plots
Draft and revise documentation
My role was to define the project direction, choose the comparison framework, run and review the workflow, ask questions about results, document what I learned, identify limitations, and organize the final project for GitHub.
This repository is therefore both a computational biology project and a record of learning how to use software and AI tools for research-oriented data analysis.
Author
Minhyuck Lee (Finn)
Valdosta State University
June 2026
