import pandas as pd
from pathlib import Path
import gzip

# ============================================================
# Project setting
# ============================================================

PROJECT_DIR = Path.home() / "Desktop" / "COAD_control_stability_project"

RAW_DIR = PROJECT_DIR / "data_raw"
SAMPLE_ID_DIR = PROJECT_DIR / "data_processed" / "sample_ids"
EXPRESSION_OUT_DIR = PROJECT_DIR / "data_processed" / "expression"
METADATA_OUT_DIR = PROJECT_DIR / "data_processed" / "metadata"

EXPRESSION_OUT_DIR.mkdir(parents=True, exist_ok=True)
METADATA_OUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Helper functions
# ============================================================

def is_gzip_file(path: Path) -> bool:
    """Check whether a file is gzip-compressed."""
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def read_first_line(path: Path) -> str:
    """Read the first line of a normal or gzip-compressed text file."""
    if is_gzip_file(path):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            return f.readline().strip()
    else:
        with open(path, "rt", encoding="utf-8", errors="replace") as f:
            return f.readline().strip()


def read_sample_ids(path: Path) -> list[str]:
    """Read sample IDs from a one-column text file."""
    if not path.exists():
        raise FileNotFoundError(f"Missing sample ID file: {path}")

    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


# ============================================================
# 1. Load sample ID files
# ============================================================

sample_files = {
    "TCGA_COAD_Tumor": SAMPLE_ID_DIR / "TCGA_COAD_Tumor_sample_ids.txt",
    "TCGA_COAD_NAT": SAMPLE_ID_DIR / "TCGA_COAD_NAT_sample_ids.txt",
    "GTEx_Colon_Sigmoid": SAMPLE_ID_DIR / "GTEx_Colon_Sigmoid_sample_ids.txt",
    "GTEx_Colon_Transverse": SAMPLE_ID_DIR / "GTEx_Colon_Transverse_sample_ids.txt",
}

sample_to_group = {}

for group_name, file_path in sample_files.items():
    ids = read_sample_ids(file_path)
    print(f"{group_name}: {len(ids)} samples")

    for sample_id in ids:
        sample_to_group[sample_id] = group_name

selected_sample_ids = list(sample_to_group.keys())

print("\nTotal selected samples:", len(selected_sample_ids))


# ============================================================
# 2. Automatically find the RSEM TPM file
# RSEM = RNA-Seq by Expectation-Maximization
# TPM = transcripts per million
# ============================================================

if not RAW_DIR.exists():
    raise FileNotFoundError(f"data_raw folder not found: {RAW_DIR}")

raw_files = [p for p in RAW_DIR.iterdir() if p.is_file()]

if not raw_files:
    raise FileNotFoundError(f"No files found in data_raw folder: {RAW_DIR}")

best_file = None
best_match_count = -1
best_header = None

print("\nSearching for RSEM TPM expression file in data_raw...")

for file_path in raw_files:
    try:
        first_line = read_first_line(file_path)
        columns = first_line.split("\t")

        if len(columns) < 10:
            continue

        header_samples = set(columns[1:])
        match_count = sum(1 for s in selected_sample_ids if s in header_samples)

        print(f"Candidate: {file_path.name} | matched samples: {match_count}")

        if match_count > best_match_count:
            best_match_count = match_count
            best_file = file_path
            best_header = columns

    except Exception as e:
        print(f"Skipped {file_path.name}: {e}")

if best_file is None or best_match_count <= 0:
    print("\nFiles found in data_raw:")
    for p in raw_files:
        print("-", p.name)

    raise ValueError(
        "Could not find an expression matrix file containing the selected sample IDs. "
        "Check whether the RSEM TPM file is inside data_raw."
    )

RSEM_TPM_PATH = best_file
compression = "gzip" if is_gzip_file(RSEM_TPM_PATH) else None

print("\nSelected RSEM TPM file:")
print(RSEM_TPM_PATH)
print("Compression:", compression)
print("Matched selected samples:", best_match_count)


# ============================================================
# 3. Check available samples in RSEM TPM file
# ============================================================

all_columns = best_header
gene_col = all_columns[0]
available_samples = set(all_columns[1:])

found_samples = [s for s in selected_sample_ids if s in available_samples]
missing_samples = [s for s in selected_sample_ids if s not in available_samples]

print("\nFirst column treated as gene column:", gene_col)
print("Found samples:", len(found_samples))
print("Missing samples:", len(missing_samples))

missing_df = pd.DataFrame({"missing_sample_id": missing_samples})
missing_df.to_csv(METADATA_OUT_DIR / "missing_sample_ids.tsv", sep="\t", index=False)


# ============================================================
# 4. Extract selected columns from RSEM TPM file
# ============================================================

usecols = [gene_col] + found_samples

print("\nReading selected expression columns...")
expr_genes_by_samples = pd.read_csv(
    RSEM_TPM_PATH,
    sep="\t",
    usecols=usecols,
    compression=compression
)

expr_genes_by_samples = expr_genes_by_samples.rename(columns={gene_col: "gene_id"})

genes_by_samples_path = EXPRESSION_OUT_DIR / "coad_expression_tpm_genes_by_samples.tsv"
expr_genes_by_samples.to_csv(genes_by_samples_path, sep="\t", index=False)

print("Saved genes x samples matrix:")
print(genes_by_samples_path)


# ============================================================
# 5. Convert to machine-learning format: samples x genes
# ============================================================

print("\nTransposing matrix to samples x genes format...")
expr_samples_by_genes = expr_genes_by_samples.set_index("gene_id").T
expr_samples_by_genes.index.name = "sample"

samples_by_genes_path = EXPRESSION_OUT_DIR / "coad_expression_tpm_samples_by_genes.tsv"
expr_samples_by_genes.to_csv(samples_by_genes_path, sep="\t")

print("Saved samples x genes matrix:")
print(samples_by_genes_path)


# ============================================================
# 6. Create sample label file
# ============================================================

labels = []

for sample_id in found_samples:
    group = sample_to_group[sample_id]

    if group == "TCGA_COAD_Tumor":
        tumor_vs_normal_label = "Tumor"
    else:
        tumor_vs_normal_label = "Normal"

    labels.append({
        "sample": sample_id,
        "group": group,
        "tumor_vs_normal_label": tumor_vs_normal_label
    })

label_df = pd.DataFrame(labels)

label_path = METADATA_OUT_DIR / "coad_sample_labels.tsv"
label_df.to_csv(label_path, sep="\t", index=False)

print("\nSaved sample labels:")
print(label_path)


# ============================================================
# 7. Save summary
# ============================================================

summary_df = label_df["group"].value_counts().reset_index()
summary_df.columns = ["group", "n_samples"]

summary_path = METADATA_OUT_DIR / "coad_extracted_sample_summary.tsv"
summary_df.to_csv(summary_path, sep="\t", index=False)

print("\nSample summary:")
print(summary_df)

print("\nDone.")