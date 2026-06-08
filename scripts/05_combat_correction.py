import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from neuroCombat import neuroCombat
import os

# 1. 경로 설정
BASE_DIR = os.path.expanduser("~/Desktop/COAD_control_stability_project")
EXPR_FILE = os.path.join(BASE_DIR, "data_processed/expression/coad_expression_tpm_genes_by_samples.tsv")
META_FILE = os.path.join(BASE_DIR, "data_processed/metadata/coad_sample_labels.tsv")
OUT_DIR   = os.path.join(BASE_DIR, "results/qc")
EXPR_OUT  = os.path.join(BASE_DIR, "data_processed/expression")
os.makedirs(OUT_DIR, exist_ok=True)

# 2. 데이터 로드
print("Loading expression matrix...")
df = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)

print("Loading sample labels...")
meta = pd.read_csv(META_FILE, sep="\t", index_col=0)

# 3. GTEx + NAT만 필터링
keep_groups = ["TCGA_COAD_NAT", "GTEx_Colon_Transverse", "GTEx_Colon_Sigmoid"]
meta_sub = meta[meta["group"].isin(keep_groups)]
common_samples = meta_sub.index.intersection(df.columns)
df_sub = df[common_samples]
meta_sub = meta_sub.loc[common_samples]

print(f"샘플 수: {df_sub.shape[1]}")
print(meta_sub["group"].value_counts())

# 4. 배치 변수 (TCGA=0, GTEx=1)
batch = meta_sub["group"].apply(
    lambda x: 0 if x == "TCGA_COAD_NAT" else 1
).values

# 5. neuroCombat 보정
print("ComBat 보정 중...")
covars = {"batch": batch}
result = neuroCombat(
    dat=df_sub.values,
    covars=pd.DataFrame(covars),
    batch_col="batch"
)
df_corrected = pd.DataFrame(
    result["data"],
    index=df_sub.index,
    columns=df_sub.columns
)
print("보정 완료!")

# 6. 저장
out_path = os.path.join(EXPR_OUT, "coad_expression_combat_corrected.tsv")
df_corrected.to_csv(out_path, sep="\t")
print(f"저장 완료: {out_path}")

# 7. PCA 함수
def run_pca_plot(df_input, meta_input, title, filename):
    X = df_input.T
    gene_var = X.var(axis=0)
    top_genes = gene_var.nlargest(3000).index
    X = X[top_genes]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=2)
    result = pca.fit_transform(X_scaled)
    colors = {
        "TCGA_COAD_NAT":         "red",
        "GTEx_Colon_Transverse": "blue",
        "GTEx_Colon_Sigmoid":    "green"
    }
    plt.figure(figsize=(8, 6))
    for group in keep_groups:
        idx = meta_input["group"] == group
        plt.scatter(
            result[idx.values, 0],
            result[idx.values, 1],
            c=colors[group],
            label=group,
            alpha=0.6,
            s=40
        )
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    out = os.path.join(OUT_DIR, filename)
    plt.savefig(out, dpi=150)
    plt.show()
    print(f"저장: {out}")

# 8. 보정 전후 PCA
print("보정 전 PCA...")
run_pca_plot(df_sub, meta_sub, "Before ComBat", "pca_before_combat.png")

print("보정 후 PCA...")
run_pca_plot(df_corrected, meta_sub, "After ComBat", "pca_after_combat.png")