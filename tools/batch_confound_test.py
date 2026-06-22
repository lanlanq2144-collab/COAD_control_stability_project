"""
Batch-confound test for the COAD control-selection project.

His claim: Tumor-vs-GTEx classifiers hit AUC=1.0 even after nested-CV leakage
correction, therefore the signal is "genuinely strong" biology.

Counter-test: classify two tissues that are BOTH normal colon — TCGA tumor-
adjacent normal (NAT) vs GTEx normal. If a model still separates them ~perfectly,
the perfect Tumor-vs-GTEx AUC is partly data-source (batch), not cancer biology.
Nested CV cannot fix this — it only fixes feature-selection leakage.

Three setups, all NORMAL-vs-NORMAL:
  1. NAT(original)  vs GTEx(ComBat)   <- mirrors his Tumor(orig)-vs-GTEx(ComBat) asymmetry  [KEY]
  2. NAT(original)  vs GTEx(original) <- uncorrected baseline (expect trivially separable)
  3. NAT(ComBat)    vs GTEx(ComBat)   <- both inside the corrected space

Feature selection is unsupervised (top-variance genes) so it is leakage-free;
AUC is 5-fold stratified cross-validated.
"""
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score

BASE = "C:/coad_run"
EXPR_ORIG   = f"{BASE}/data_processed/expression/coad_expression_tpm_genes_by_samples.tsv"
EXPR_COMBAT = f"{BASE}/data_processed/expression/coad_expression_combat_corrected.tsv"
LABELS      = f"{BASE}/data_processed/metadata/coad_sample_labels.tsv"

N_TOP_VAR = 2000
CV = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

lab = pd.read_csv(LABELS, sep="\t", index_col=0)["group"].to_dict()
orig = pd.read_csv(EXPR_ORIG, sep="\t", index_col=0)
comb = pd.read_csv(EXPR_COMBAT, sep="\t", index_col=0)
genes = orig.index.intersection(comb.index)
orig, comb = orig.loc[genes], comb.loc[genes]

def cols(df, groups):
    return [c for c in df.columns if lab.get(c) in groups]

GTEX = ["GTEx_Colon_Transverse", "GTEx_Colon_Sigmoid"]
nat_orig   = cols(orig, ["TCGA_COAD_NAT"])
gtex_orig  = cols(orig, GTEX)
gtex_comb  = cols(comb, GTEX)
nat_comb   = cols(comb, ["TCGA_COAD_NAT"])

def auc(Xdf_pos, Xdf_neg, pos_cols, neg_cols):
    # build samples x genes, label pos=1 (NAT), neg=0 (GTEx); top-variance genes
    X = pd.concat([Xdf_pos[pos_cols], Xdf_neg[neg_cols]], axis=1).T
    y = np.array([1]*len(pos_cols) + [0]*len(neg_cols))
    top = X.var(axis=0).nlargest(N_TOP_VAR).index
    Xs = StandardScaler().fit_transform(X[top].values)
    out = {}
    for name, clf in [
        ("LogReg", LogisticRegression(class_weight="balanced", max_iter=2000, random_state=42)),
        ("RandomForest", RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                                 random_state=42, n_jobs=-1)),
    ]:
        p = cross_val_predict(clf, Xs, y, cv=CV, method="predict_proba")[:, 1]
        out[name] = roc_auc_score(y, p)
    return out, len(pos_cols), len(neg_cols)

print("="*74)
print("NORMAL-vs-NORMAL separability  (NAT=TCGA normal, GTEx=normal)  5-fold CV AUC")
print("="*74)
tests = [
    ("1. NAT(orig)   vs GTEx(ComBat)  [mirrors his cross-source setup]", orig, comb, nat_orig, gtex_comb),
    ("2. NAT(orig)   vs GTEx(orig)    [uncorrected baseline]",           orig, orig, nat_orig, gtex_orig),
    ("3. NAT(ComBat) vs GTEx(ComBat)  [both corrected]",                 comb, comb, nat_comb, gtex_comb),
]
for label, dp, dn, pc, nc in tests:
    res, npos, nneg = auc(dp, dn, pc, nc)
    print(f"\n{label}")
    print(f"   n: {npos} NAT vs {nneg} GTEx")
    for m, a in res.items():
        print(f"   {m:14s} AUC = {a:.4f}")

print("\n" + "="*74)
print("Reading: both groups are normal colon, so biology predicts AUC ~0.5-0.7.")
print("AUC ~1.0 (esp. setup 1) => the classifier is reading DATA SOURCE, not cancer.")
print("His Tumor-vs-GTEx AUC=1.0 therefore conflates source + biology; nested CV")
print("does not address this. The clean comparison is Tumor vs TCGA-NAT (same source).")
