import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import os

# 1. 경로 설정
BASE_DIR = os.path.expanduser("~/Desktop/COAD_control_stability_project")
EXPR_FILE = os.path.join(BASE_DIR, "data_processed/expression/coad_expression_tpm_genes_by_samples.tsv")
META_FILE = os.path.join(BASE_DIR, "data_processed/metadata/coad_sample_labels.tsv")
OUT_DIR   = os.path.join(BASE_DIR, "results/qc")
os.makedirs(OUT_DIR, exist_ok=True)

# 2. 데이터 로드
print("Loading expression matrix...")
df = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)  # 행=유전자, 열=샘플

print("Loading sample labels...")
meta = pd.read_csv(META_FILE, sep="\t", index_col=0)

# 3. NAT + GTEx만 필터링 (Tumor 제외)
keep_groups = ["TCGA_COAD_NAT", "GTEx_Colon_Transverse", "GTEx_Colon_Sigmoid"]
meta_sub = meta[meta["group"].isin(keep_groups)]

# 4. 해당 샘플만 expression에서 추출 후 전치 (행=샘플, 열=유전자)
common_samples = meta_sub.index.intersection(df.columns)
X = df[common_samples].T

print(f"샘플 수: {X.shape[0]}, 유전자 수: {X.shape[1]}")
print(meta_sub["group"].value_counts())

# 5. log2 변환 (이미 돼있으면 스킵)
if X.values.max() > 50:
    print("log2(TPM+1) 변환 적용...")
    X = np.log2(X + 1)

# 6. 분산 상위 3000 유전자만 사용
gene_var = X.var(axis=0)
top_genes = gene_var.nlargest(3000).index
X = X[top_genes]

# 7. 스케일링 + PCA
print("PCA 실행 중...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(X_scaled)

# 8. 시각화
colors = {
    "TCGA_COAD_NAT":       "red",
    "GTEx_Colon_Transverse": "blue",
    "GTEx_Colon_Sigmoid":    "green"
}

plt.figure(figsize=(8, 6))
for group in keep_groups:
    idx = meta_sub.loc[common_samples, "group"] == group
    plt.scatter(
        pca_result[idx.values, 0],
        pca_result[idx.values, 1],
        c=colors[group],
        label=group,
        alpha=0.6,
        s=40
    )

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.title("PCA: NAT vs GTEx (Tumor 제외)")
plt.legend()
plt.tight_layout()

out_path = os.path.join(OUT_DIR, "pca_nat_vs_gtex.png")
plt.savefig(out_path, dpi=150)
plt.show()
print(f"저장 완료: {out_path}")