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

## Quick setup (helper scripts)

Direct download + group definitions are scripted so the pipeline is reproducible
from a clean clone:

```bash
python tools/fetch_data.py        # downloads the 3 files into data_raw/ (1.3 GB matrix + phenotype + GEO tar)
python tools/build_sample_ids.py  # writes data_processed/sample_ids/*.txt from the phenotype
```

Direct URLs (UCSC Xena Toil hub / NCBI GEO):

- `https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz`
- `https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156451/suppl/GSE156451_RAW.tar`

### Sample selection (the group definitions)

The four groups are defined by these phenotype filters (verified to reproduce the
exact group sizes used in the analysis):

| Group | n | Filter (`TcgaTargetGTEX_phenotype.txt`) |
|---|---:|---|
| TCGA_COAD_Tumor | 288 | `detailed_category == "Colon Adenocarcinoma"` & `_sample_type == "Primary Tumor"` |
| TCGA_COAD_NAT | 41 | `detailed_category == "Colon Adenocarcinoma"` & `_sample_type == "Solid Tissue Normal"` |
| GTEx_Colon_Transverse | 167 | `_study == "GTEX"` & `primary disease or tissue == "Colon - Transverse"` |
| GTEx_Colon_Sigmoid | 141 | `_study == "GTEX"` & `primary disease or tissue == "Colon - Sigmoid"` |

The resulting sample-ID lists are version-controlled in [`data/sample_ids/`](sample_ids/)
for provenance. `extract_coad_matrix.py` reads them from `data_processed/sample_ids/`,
so copy them there (or just run `tools/build_sample_ids.py`, which writes them directly).

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
