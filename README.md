# COAD_control_stability_project
A Python-based exploratory research project analyzing how normal control selection affects machine learning-based biomarker discovery in colon adenocarcinoma.

# Control Group Selection Affects Biomarker Discovery in Colorectal Cancer

An independent computational biology project investigating how the choice of normal reference tissue affects machine learning-based biomarker discovery in colorectal adenocarcinoma (COAD).

---

## Research Question

Does changing the definition of "normal" change which genes a machine learning model selects as cancer biomarker candidates?

---

## Background

Most studies that identify cancer biomarkers using public gene expression data treat the choice of normal control group as a fixed decision. This project systematically evaluates whether that choice — NAT (tumor-adjacent normal tissue) vs. healthy GTEx tissue — meaningfully changes the biomarker candidates selected by machine learning models.

---

## Data Sources

| Dataset | Description | Samples |
|---|---|---|
| TCGA COAD Tumor | Colorectal adenocarcinoma tumor tissue | 288 |
| TCGA COAD NAT | Tumor-adjacent normal tissue (same patients) | 41 |
| GTEx Colon Transverse | Healthy colon tissue (transverse) | 167 |
| GTEx Colon Sigmoid | Healthy colon tissue (sigmoid) | 141 |
| GSE156451 | Independent CRC cohort (Wuhan University) | 72 paired |

Data downloaded from [UCSC Xena](https://xenabrowser.net/) and [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/).

> Raw data files are not included in this repository due to size and licensing. See links above to download.

---

## Analysis Pipeline

```
1. Quality Control (QC)
   → PCA to visualize batch effects between TCGA and GTEx
   → neuroCombat batch correction applied to NAT + GTEx only

2. Differential Expression Analysis (DEG)
   → Welch's t-test + Benjamini-Hochberg FDR correction
   → Thresholds: |log2FC| > 1, FDR < 0.05
   → Performed separately for each of 4 control groups

3. Machine Learning
   → Feature pool: top 500 DEGs by FDR
   → LASSO (feature selection) + Random Forest (importance ranking)
   → Nested cross-validation to rule out feature-selection leakage

4. Biomarker Stability Analysis
   → Venn diagrams comparing candidate lists across 4 control groups
   → Stability scoring (1–4 points per gene)
   → Final candidates: genes with score ≥ 3

5. Biological Validation
   → Gene name conversion via mygene.info
   → PubMed literature review
   → Pathway analysis
   → Survival analysis (Kaplan-Meier)

6. External Validation
   → Independent replication in GSE156451 (n=72 paired CRC patients)
```

---

## Key Results

### DEG Counts by Control Group

| Comparison | DEGs |
|---|---|
| Tumor vs NAT | 11,623 |
| Tumor vs GTEx All | 17,562 |
| Tumor vs Transverse | 16,516 |
| Tumor vs Sigmoid | 19,298 |

### Biomarker Stability

- 4-way overlap across all control groups: **0 genes**
- NAT vs GTEx overlap: **1–2 genes**
- GTEx groups overlap: **5–8 genes**

→ Biomarker candidates changed substantially depending on how "normal" was defined.

### Final Biomarker Candidates (Stability Score ≥ 3)

| Gene | Role | Prior Literature |
|---|---|---|
| CDH3 | Cell adhesion | Established CRC biomarker |
| CLDN1 | Epithelial barrier function | Confirmed by meta-analysis |
| ETV4 | Transcription factor | Key role in CRC progression |
| KRT80 | Cytoskeleton | Promotes invasion and metastasis |

### External Validation (GSE156451)

All four candidates independently replicated in a Chinese CRC cohort (n=72 paired samples):

| Gene | log₂FC | p-value |
|---|---|---|
| CDH3 | +3.16 | 1.1×10⁻³⁶ |
| CLDN1 | +3.35 | 7.8×10⁻³⁵ |
| ETV4 | +4.08 | 2.8×10⁻⁴⁴ |
| KRT80 | +2.61 | 1.5×10⁻³¹ |

### Leakage Validation

Nested cross-validation confirmed AUC = 1.000 across all comparisons (delta = 0.000), ruling out feature-selection leakage.

---

## Tools and Libraries

```
Python 3.14
pandas, numpy, scipy
scikit-learn (LASSO, Random Forest, cross-validation)
neuroCombat (batch correction)
matplotlib, seaborn (visualization)
lifelines (survival analysis)
mygene (gene ID conversion)
```

---

## Repository Structure

```
COAD-biomarker-stability/
├── README.md
├── WORK_LOG.md
├── data/
│   └── README.md         # Data sources and download instructions
├── scripts/
│   ├── 01_qc.py
│   ├── 02_deg.py
│   ├── 03_ml.py
│   ├── 04_stability.py
│   ├── 05_biological_validation.py
│   └── 06_external_validation.py
└── results/
    ├── qc/
    ├── deg/
    ├── ml/
    ├── ml_corrected/
    ├── stability/
    ├── biological_validation/
    └── external_validation/
```

---

## Notes

This project was completed independently as a first research experience. Claude Code (Anthropic) was used as a coding assistant throughout the analysis. All research questions, design decisions, and interpretations were made by the author.

---

## Author

Minhyuck Lee (Finn)  
Valdosta State University
June 2026
