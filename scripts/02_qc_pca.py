import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# ============================================================
# Project paths
# ============================================================

PROJECT_DIR = Path.home() / "Desktop" / "COAD_control_stability_project"

EXPR_PATH = PROJECT_DIR / "data_processed" / "expression" / "coad_expression_tpm_samples_by_genes.tsv"
LABEL_PATH = PROJECT_DIR / "data_processed" / "metadata" / "coad_sample_labels.tsv"

QC_OUT_DIR = PROJECT_DIR / "data_processed" / "qc"
FIG_OUT_DIR = PROJECT_DIR / "results" / "qc"

QC_OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. Load expression matrix and labels
# ============================================================

print("Loading expression matrix...")
expr = pd.read_csv(EXPR_PATH, sep="\t", index_col=0)

print("Loading sample labels...")
labels = pd.read_csv(LABEL_PATH, sep="\t")

print("Expression matrix shape before alignment:", expr.shape)
print("Labels shape:", labels.shape)

# Make sure expression values are numeric
expr = expr.apply(pd.to_numeric, errors="coerce")

# Keep only samples that exist in both files
common_samples = expr.index.intersection(labels["sample"])

expr = expr.loc[common_samples]
labels = labels.set_index("sample").loc[common_samples].reset_index()

print("Expression matrix shape after alignment:", expr.shape)
print("Number of matched labels:", labels.shape[0])

# ============================================================
# 2. Decide whether to log-transform
# ============================================================

min_value = np.nanmin(expr.values)
max_value = np.nanmax(expr.values)

print("Expression value range:")
print("Min:", min_value)
print("Max:", max_value)

if min_value < 0:
    print("Negative values detected. Treating expression matrix as already log-transformed.")
    expr_processed = expr.copy()
else:
    print("No negative values detected. Applying log2(TPM + 1).")
    expr_processed = np.log2(expr + 1)

# Replace infinite values with NaN, then handle missing values
expr_processed = expr_processed.replace([np.inf, -np.inf], np.nan)

n_missing_before = expr_processed.isna().sum().sum()
print("Total missing values before imputation:", n_missing_before)

# Drop genes that are completely missing
expr_processed = expr_processed.dropna(axis=1, how="all")

# Fill remaining missing values with each gene's median
expr_processed = expr_processed.fillna(expr_processed.median(axis=0))

# If any missing values remain, fill with 0
expr_processed = expr_processed.fillna(0)

n_missing_after = expr_processed.isna().sum().sum()
print("Total missing values after imputation:", n_missing_after)

# ============================================================
# 3. Select variable genes
# ============================================================

print("Selecting variable genes...")

gene_variances = expr_processed.var(axis=0)

# Remove genes with zero variance
variable_genes = gene_variances[gene_variances > 0].index
expr_variable = expr_processed[variable_genes]

N_TOP_VARIABLE_GENES = 3000
top_genes = expr_variable.var(axis=0).sort_values(ascending=False).head(N_TOP_VARIABLE_GENES).index
expr_top = expr_variable[top_genes]

print("Expression matrix after filtering:", expr_top.shape)

# Save selected genes
pd.DataFrame({"gene_id": top_genes}).to_csv(
    QC_OUT_DIR / "pca_top_variable_genes.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 4. Standardize genes
# ============================================================

print("Scaling data...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr_top)

# Final NaN check
if np.isnan(X_scaled).sum() > 0:
    raise ValueError("NaN values still remain after preprocessing.")

# ============================================================
# 5. Run PCA
# ============================================================

print("Running PCA...")
pca = PCA(n_components=10)
pca_scores = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame(
    pca_scores[:, :5],
    index=expr_top.index,
    columns=["PC1", "PC2", "PC3", "PC4", "PC5"]
)

pca_df.index.name = "sample"
pca_df = pca_df.reset_index()
pca_df = pca_df.merge(labels, on="sample", how="left")

explained = pd.DataFrame({
    "principal_component": [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
    "explained_variance_ratio": pca.explained_variance_ratio_
})

pca_df.to_csv(QC_OUT_DIR / "pca_scores.tsv", sep="\t", index=False)
explained.to_csv(QC_OUT_DIR / "pca_explained_variance.tsv", sep="\t", index=False)

print("Saved PCA scores and explained variance.")

# ============================================================
# 6. Plot PCA
# ============================================================

def plot_pca(color_column, output_name, title):
    plt.figure(figsize=(8, 6))

    groups = pca_df[color_column].dropna().unique()

    for group in groups:
        subset = pca_df[pca_df[color_column] == group]
        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            label=group,
            alpha=0.75,
            s=35
        )

    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    plt.title(title)
    plt.legend(fontsize=8)
    plt.tight_layout()

    out_path = FIG_OUT_DIR / output_name
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Saved plot: {out_path}")

plot_pca(
    color_column="group",
    output_name="pca_by_group.png",
    title="PCA of COAD Tumor, NAT, and GTEx Colon Samples"
)

plot_pca(
    color_column="tumor_vs_normal_label",
    output_name="pca_tumor_vs_normal.png",
    title="PCA: Tumor vs Normal"
)

# ============================================================
# 7. Print summary
# ============================================================

print("\nExplained variance:")
print(explained.head())

print("\nGroup counts:")
print(pca_df["group"].value_counts())

print("\nDone.")