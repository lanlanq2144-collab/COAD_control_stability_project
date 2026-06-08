# Data

This project uses three publicly available datasets. Due to file size and licensing constraints, raw and processed data files are **not included in this repository**.

Download each dataset from the sources below and place them in the appropriate folders before running the analysis scripts.

---

## Data Sources

### TCGA COAD & GTEx Expression Data

- **Source:** UCSC Xena Browser
- **URL:** https://xenabrowser.net/
- **Dataset:** TCGA + GTEx combined RSEM TPM (hg38) — *"TCGA TARGET GTEx"* cohort
- **Files needed:**
  - `TcgaTargetGtex_rsem_gene_tpm` — expression matrix (log₂(TPM + 0.001))
  - `TcgaTargetGTEX_phenotype.txt` — sample phenotype/metadata
- **Place in:** `data_raw/`

### GSE156451 (Independent CRC Cohort)

- **Source:** NCBI GEO
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156451
- **Description:** RNA-seq (FPKM) from 72 paired colorectal cancer tumor and normal tissue samples (Wuhan University)
- **Files needed:** Download `GSE156451_RAW.tar` from the Supplementary Files section
- **Place in:** `data_raw/GSE156451/`

---

## Folder Structure (after download)

```
data_raw/
├── RSEM tpm (n=19,131)              ← UCSC Xena expression matrix
├── TCGA GTEX main categories.txt    ← sample group labels
├── TCGA TARGET GTEX selected phenotypes.txt
└── GSE156451/
    ├── GSE156451_RAW.tar
    └── GSM4731674_T1-RNA.txt.gz     ← extracted per-sample files
        ...

data_processed/
├── expression/                      ← filtered & batch-corrected matrices
├── metadata/                        ← sample labels used in analysis
├── qc/                              ← PCA scores and variance tables
└── sample_ids/                      ← curated sample ID lists

metadata/                            ← additional sample metadata
```

---

## Notes

- `data_raw/`, `data_processed/`, and `metadata/` are excluded from version control via `.gitignore`
- All preprocessing steps are documented in `scripts/` and reproducible from the raw downloads above
- TCGA data usage is subject to the [TCGA Data Use Certification](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/using-tcga/understanding-tcga-data-policies)
- GTEx data usage is subject to the [GTEx data use agreement](https://gtexportal.org/home/datasets)
