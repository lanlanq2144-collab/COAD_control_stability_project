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

expr = expr.apply(pd.to_numeric, errors="coerce")

# ============================================================
# 2. Keep only TCGA COAD Tumor and TCGA COAD NAT
# ============================================================

tcga_groups = ["TCGA_COAD_Tumor", "TCGA_COAD_NAT"]

labels_tcga = labels[labels["group"].isin(tcga_groups)].copy()

common_samples = expr.index.intersection(labels_tcga["sample"])

expr_tcga = expr.loc[common_samples]
labels_tcga = labels_tcga.set_index("sample").loc[common_samples].reset_index()

print("TCGA-only expression matrix shape:", expr_tcga.shape)
print("TCGA-only labels shape:", labels_tcga.shape)

print("\nGroup counts:")
print(labels_tcga["group"].value_counts())

# ============================================================
# 3. Check whether values are already log-transformed
# ============================================================

min_value = np.nanmin(expr_tcga.values)
max_value = np.nanmax(expr_tcga.values)

print("\nExpression value range:")
print("Min:", min_value)
print("Max:", max_value)

if min_value < 0:
    print("Negative values detected. Treating expression matrix as already log-transformed.")
    expr_processed = expr_tcga.copy()
else:
    print("No negative values detected. Applying log2(TPM + 1).")
    expr_processed = np.log2(expr_tcga + 1)

# Handle missing values
expr_processed = expr_processed.replace([np.inf, -np.inf], np.nan)

print("Missing values before imputation:", expr_processed.isna().sum().sum())

expr_processed = expr_processed.dropna(axis=1, how="all")
expr_processed = expr_processed.fillna(expr_processed.median(axis=0))
expr_processed = expr_processed.fillna(0)

print("Missing values after imputation:", expr_processed.isna().sum().sum())

# ============================================================
# 4. Select top variable genes
# ============================================================

gene_variances = expr_processed.var(axis=0)

# Remove zero-variance genes
variable_genes = gene_variances[gene_variances > 0].index
expr_variable = expr_processed[variable_genes]

N_TOP_VARIABLE_GENES = 3000

top_genes = (
    expr_variable
    .var(axis=0)
    .sort_values(ascending=False)
    .head(N_TOP_VARIABLE_GENES)
    .index
)

expr_top = expr_variable[top_genes]

print("\nExpression matrix after filtering:", expr_top.shape)

pd.DataFrame({"gene_id": top_genes}).to_csv(
    QC_OUT_DIR / "tcga_only_pca_top_variable_genes.tsv",
    sep="\t",
    index=False
)

# ============================================================
# 5. Standardize and run PCA
# ============================================================

print("Scaling data...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(expr_top)

if np.isnan(X_scaled).sum() > 0:
    raise ValueError("NaN values still remain after preprocessing.")

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
pca_df = pca_df.merge(labels_tcga, on="sample", how="left")

explained = pd.DataFrame({
    "principal_component": [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
    "explained_variance_ratio": pca.explained_variance_ratio_
})

pca_df.to_csv(QC_OUT_DIR / "tcga_only_pca_scores.tsv", sep="\t", index=False)
explained.to_csv(QC_OUT_DIR / "tcga_only_pca_explained_variance.tsv", sep="\t", index=False)

# ============================================================
# 6. Plot PCA
# ============================================================

plt.figure(figsize=(8, 6))

for group in pca_df["group"].unique():
    subset = pca_df[pca_df["group"] == group]
    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        alpha=0.75,
        s=40
    )

pc1_var = pca.explained_variance_ratio_[0] * 100
pc2_var = pca.explained_variance_ratio_[1] * 100

plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
plt.title("PCA of TCGA COAD Tumor vs TCGA NAT Only")
plt.legend(fontsize=9)
plt.tight_layout()

out_path = FIG_OUT_DIR / "pca_tcga_tumor_vs_nat_only.png"
plt.savefig(out_path, dpi=300)
plt.close()

print("\nSaved plot:")
print(out_path)

print("\nExplained variance:")
print(explained.head())

print("\nDone.")