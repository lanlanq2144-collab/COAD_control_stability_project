"""
ML analysis for TCGA COAD DEG results
Two modes:
  A) Classification: LASSO-CV + RF-CV → ROC-AUC
  B) Biomarker discovery: LASSO regularization path + RF importance
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

# ── paths ─────────────────────────────────────────────────────────────────────
BASE        = "/Users/georgiatechgogofinn/Desktop/COAD_control_stability_project"
TPM_FILE    = f"{BASE}/data_processed/expression/coad_expression_tpm_genes_by_samples.tsv"
COMBAT_FILE = f"{BASE}/data_processed/expression/coad_expression_combat_corrected.tsv"
LABELS_FILE = f"{BASE}/data_processed/metadata/coad_sample_labels.tsv"
DEG_DIR     = f"{BASE}/results/deg"
OUT_DIR     = f"{BASE}/results/ml"
os.makedirs(OUT_DIR, exist_ok=True)

# ── hyperparams ───────────────────────────────────────────────────────────────
N_TOP_DEG        = 500    # top significant DEGs as feature pool
N_CV             = 5      # CV folds
N_TREES          = 300    # RF trees
RF_TOP_N         = 50     # RF top-N genes for consensus
# For biomarker discovery, fix C_BIO so LASSO selects ~20-60 genes
C_BIO            = 0.5    # lenient regularization → richer gene list

# ── load data ─────────────────────────────────────────────────────────────────
print("Loading data...")
labels    = pd.read_csv(LABELS_FILE, sep="\t", index_col=0)
label_map = labels["group"].to_dict()

expr_orig   = pd.read_csv(TPM_FILE,    sep="\t", index_col=0)
expr_combat = pd.read_csv(COMBAT_FILE, sep="\t", index_col=0)

common_genes = expr_orig.index.intersection(expr_combat.index)
expr_orig    = expr_orig.loc[common_genes]
expr_combat  = expr_combat.loc[common_genes]

tumor_cols  = [c for c in expr_orig.columns if label_map.get(c) == "TCGA_COAD_Tumor"]
expr_merged = pd.concat([expr_orig[tumor_cols], expr_combat], axis=1)

print(f"  Original : {expr_orig.shape}  ComBat : {expr_combat.shape}")


def get_cols(expr_df, groups):
    if isinstance(groups, str):
        groups = [groups]
    return [c for c in expr_df.columns if label_map.get(c) in groups]


# ── comparisons ───────────────────────────────────────────────────────────────
COMPARISONS = [
    dict(name="Tumor_vs_TCGA_NAT",
         expr=expr_orig,  tumor="TCGA_COAD_Tumor",  control="TCGA_COAD_NAT",
         deg_file=f"{DEG_DIR}/Tumor_vs_TCGA_NAT_DEG.tsv"),
    dict(name="Tumor_vs_GTEx_All",
         expr=expr_merged, tumor="TCGA_COAD_Tumor",
         control=["GTEx_Colon_Transverse", "GTEx_Colon_Sigmoid"],
         deg_file=f"{DEG_DIR}/Tumor_vs_GTEx_All_DEG.tsv"),
    dict(name="Tumor_vs_GTEx_Colon_Transverse",
         expr=expr_merged, tumor="TCGA_COAD_Tumor",  control="GTEx_Colon_Transverse",
         deg_file=f"{DEG_DIR}/Tumor_vs_GTEx_Colon_Transverse_DEG.tsv"),
    dict(name="Tumor_vs_GTEx_Colon_Sigmoid",
         expr=expr_merged, tumor="TCGA_COAD_Tumor",  control="GTEx_Colon_Sigmoid",
         deg_file=f"{DEG_DIR}/Tumor_vs_GTEx_Colon_Sigmoid_DEG.tsv"),
]

summary_rows   = []
lasso_bio_sets = {}   # biomarker-mode LASSO genes per comparison
rf_top_sets    = {}   # RF top-N genes per comparison

# ══ main loop ════════════════════════════════════════════════════════════════
for comp in COMPARISONS:
    name = comp["name"]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    # ── feature pool: top N significant DEGs ──────────────────────────────────
    deg     = pd.read_csv(comp["deg_file"], sep="\t")
    deg_sig = deg[(deg["padj"] < 0.05) & (deg["log2FC"].abs() > 1)].copy()
    top_ids = (deg_sig.nsmallest(N_TOP_DEG, "padj")["gene_id"]
                      .tolist())
    top_ids = [g for g in top_ids if g in comp["expr"].index]
    print(f"  Feature pool : {len(top_ids)} DEGs")

    # ── X, y ─────────────────────────────────────────────────────────────────
    expr    = comp["expr"]
    t_cols  = get_cols(expr, comp["tumor"])
    c_cols  = get_cols(expr, comp["control"])
    all_col = t_cols + c_cols
    y       = np.array([1]*len(t_cols) + [0]*len(c_cols))

    X_raw   = expr.loc[top_ids, all_col].T.values.astype(float)
    scaler  = StandardScaler()
    X       = scaler.fit_transform(X_raw)
    cv      = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=42)
    print(f"  Samples : {len(t_cols)} Tumor + {len(c_cols)} Control  |  X {X.shape}")

    # ══════════════════════════════════════════════════════════════════════════
    # A. CLASSIFICATION — LASSO-CV optimal C
    # ══════════════════════════════════════════════════════════════════════════
    print("  [A] LASSO-CV (classification)...")
    lasso_cv = LogisticRegressionCV(
        penalty="l1", solver="liblinear",
        Cs=np.logspace(-3, 2, 30),
        cv=cv, scoring="roc_auc",
        class_weight="balanced", max_iter=2000, random_state=42,
    )
    lasso_cv.fit(X, y)
    best_C  = float(lasso_cv.C_[0])
    coef_cv = lasso_cv.coef_[0]
    nz      = np.where(coef_cv != 0)[0]
    lasso_cv_genes = [top_ids[i] for i in nz]
    print(f"    best C={best_C:.5f} | genes selected={len(lasso_cv_genes)}")

    # CV probabilities for ROC
    lasso_prob = cross_val_predict(
        LogisticRegression(penalty="l1", solver="liblinear", C=best_C,
                           class_weight="balanced", max_iter=2000, random_state=42),
        X, y, cv=cv, method="predict_proba")[:, 1]
    lasso_auc = roc_auc_score(y, lasso_prob)

    # ══════════════════════════════════════════════════════════════════════════
    # B. BIOMARKER DISCOVERY — LASSO at C_BIO (less regularized)
    # ══════════════════════════════════════════════════════════════════════════
    print(f"  [B] LASSO biomarker (C={C_BIO})...")
    lasso_bio = LogisticRegression(
        penalty="l1", solver="liblinear", C=C_BIO,
        class_weight="balanced", max_iter=2000, random_state=42,
    )
    lasso_bio.fit(X, y)
    coef_bio = lasso_bio.coef_[0]
    nz_bio   = np.where(coef_bio != 0)[0]
    bio_genes = [top_ids[i] for i in nz_bio]
    bio_coefs = coef_bio[nz_bio]
    print(f"    genes selected={len(bio_genes)}")

    lasso_bio_df = (pd.DataFrame({
        "gene_id"   : bio_genes,
        "lasso_coef": bio_coefs,
        "abs_coef"  : np.abs(bio_coefs),
    }).sort_values("abs_coef", ascending=False))
    lasso_bio_df.to_csv(f"{OUT_DIR}/{name}_lasso_genes.tsv", sep="\t", index=False)
    lasso_bio_sets[name] = set(bio_genes)

    # ── LASSO regularization path plot ───────────────────────────────────────
    Cs_path = np.logspace(-3, 1, 60)
    n_selected = []
    for c in Cs_path:
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=c,
                                 class_weight="balanced", max_iter=2000)
        clf.fit(X, y)
        n_selected.append(np.sum(clf.coef_[0] != 0))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogx(Cs_path, n_selected, lw=2, color="#27AE60")
    ax.axvline(best_C, ls="--", color="#E74C3C", lw=1.5,
               label=f"CV-optimal C={best_C:.4f}")
    ax.axvline(C_BIO,  ls="--", color="#2980B9", lw=1.5,
               label=f"Biomarker C={C_BIO}")
    ax.set_xlabel("C  (1 / regularization strength)", fontsize=11)
    ax.set_ylabel("# genes selected by LASSO", fontsize=11)
    ax.set_title(f"LASSO Regularization Path\n{name.replace('_',' ')}", fontsize=11)
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/{name}_lasso_path.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════════════════
    # C. RANDOM FOREST
    # ══════════════════════════════════════════════════════════════════════════
    print("  [C] Random Forest...")
    rf = RandomForestClassifier(n_estimators=N_TREES, class_weight="balanced",
                                random_state=42, n_jobs=-1)
    rf.fit(X, y)
    rf_prob = cross_val_predict(
        RandomForestClassifier(n_estimators=N_TREES, class_weight="balanced",
                               random_state=42, n_jobs=-1),
        X, y, cv=cv, method="predict_proba")[:, 1]
    rf_auc = roc_auc_score(y, rf_prob)

    rf_df = (pd.DataFrame({"gene_id": top_ids,
                            "rf_importance": rf.feature_importances_})
               .sort_values("rf_importance", ascending=False))
    rf_df.to_csv(f"{OUT_DIR}/{name}_rf_importance.tsv", sep="\t", index=False)
    rf_top_sets[name] = set(rf_df.head(RF_TOP_N)["gene_id"])

    print(f"    LASSO AUC (CV-optimal) = {lasso_auc:.4f}")
    print(f"    RF    AUC              = {rf_auc:.4f}")

    # ── ROC curve ─────────────────────────────────────────────────────────────
    fpr_l, tpr_l, _ = roc_curve(y, lasso_prob)
    fpr_r, tpr_r, _ = roc_curve(y, rf_prob)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(fpr_l, tpr_l, lw=2.2, color="#E74C3C",
            label=f"LASSO  (AUC = {lasso_auc:.3f}, {len(lasso_cv_genes)} genes)")
    ax.plot(fpr_r, tpr_r, lw=2.2, color="#2980B9",
            label=f"Random Forest (AUC = {rf_auc:.3f})")
    ax.plot([0, 1], [0, 1], lw=1, ls="--", color="gray", alpha=0.6)
    ax.fill_between(fpr_l, tpr_l, alpha=0.07, color="#E74C3C")
    ax.fill_between(fpr_r, tpr_r, alpha=0.07, color="#2980B9")
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title(f"ROC Curve — {name.replace('_',' ')}", fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.set_xlim([-0.01, 1.01])
    ax.set_ylim([-0.01, 1.02])
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/{name}_roc_curve.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── RF top-20 importance bar ───────────────────────────────────────────────
    top20 = rf_df.head(20)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(range(len(top20)), top20["rf_importance"].values[::-1],
            color="#2980B9", alpha=0.8)
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels([g.split(".")[0] for g in top20["gene_id"].values[::-1]],
                       fontsize=9)
    ax.set_xlabel("RF Feature Importance (Mean Decrease Impurity)", fontsize=11)
    ax.set_title(f"Top 20 RF Genes — {name.replace('_',' ')}", fontsize=11,
                 fontweight="bold")
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/{name}_rf_top20.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── per-comparison biomarker candidates ───────────────────────────────────
    consensus  = lasso_bio_sets[name] & rf_top_sets[name]
    print(f"  Consensus biomarkers (LASSO_bio ∩ RF top-{RF_TOP_N}): {len(consensus)}")

    deg_idx   = deg_sig.set_index("gene_id")
    cand_rows = []
    for gene in consensus:
        lrow = lasso_bio_df[lasso_bio_df["gene_id"] == gene]
        rrow = rf_df[rf_df["gene_id"] == gene]
        cand_rows.append({
            "gene_id"       : gene,
            "gene_id_clean" : gene.split(".")[0],
            "log2FC"        : deg_idx.loc[gene, "log2FC"] if gene in deg_idx.index else np.nan,
            "padj"          : deg_idx.loc[gene, "padj"]   if gene in deg_idx.index else np.nan,
            "lasso_coef"    : lrow["lasso_coef"].values[0] if len(lrow) else 0,
            "rf_importance" : rrow["rf_importance"].values[0] if len(rrow) else 0,
            "direction"     : ("Up" if deg_idx.loc[gene, "log2FC"] > 0
                               else "Down") if gene in deg_idx.index else "NA",
        })
    cand_df = (pd.DataFrame(cand_rows)
                 .sort_values("rf_importance", ascending=False))
    cand_df.to_csv(f"{OUT_DIR}/{name}_biomarker_candidates.tsv", sep="\t", index=False)

    summary_rows.append(dict(
        comparison           = name,
        n_tumor              = len(t_cols),
        n_control            = len(c_cols),
        n_input_genes        = len(top_ids),
        n_lasso_cv_genes     = len(lasso_cv_genes),
        n_lasso_bio_genes    = len(bio_genes),
        n_consensus          = len(consensus),
        lasso_auc            = round(lasso_auc, 4),
        rf_auc               = round(rf_auc,    4),
    ))

# ══ cross-comparison stable biomarkers ═══════════════════════════════════════
print("\n" + "="*60)
print("  Cross-comparison stable biomarkers (LASSO biomarker mode)")
print("="*60)

gene_freq = {}
for n, genes in lasso_bio_sets.items():
    for g in genes:
        gene_freq[g] = gene_freq.get(g, 0) + 1

for thr in [2, 3, 4]:
    cnt = sum(1 for c in gene_freq.values() if c >= thr)
    print(f"  Selected in ≥ {thr} comparisons: {cnt} genes")

# Build stable table
deg_nat = pd.read_csv(f"{DEG_DIR}/Tumor_vs_TCGA_NAT_DEG.tsv", sep="\t")
deg_nat = deg_nat.set_index("gene_id")

stable_rows = []
for gene, cnt in sorted(gene_freq.items(), key=lambda x: -x[1]):
    if cnt < 2:
        continue
    row = {"gene_id": gene,
           "gene_id_clean": gene.split(".")[0],
           "n_comparisons": cnt}
    if gene in deg_nat.index:
        row["log2FC_vs_NAT"] = deg_nat.loc[gene, "log2FC"]
        row["padj_vs_NAT"]   = deg_nat.loc[gene, "padj"]
    for cname in [c["name"] for c in COMPARISONS]:
        try:
            rdf = pd.read_csv(f"{OUT_DIR}/{cname}_rf_importance.tsv",
                              sep="\t", index_col=0)
            row[f"rf_{cname}"] = rdf.loc[gene, "rf_importance"] if gene in rdf.index else np.nan
        except Exception:
            pass
    stable_rows.append(row)

stable_df = pd.DataFrame(stable_rows)
stable_df.to_csv(f"{OUT_DIR}/stable_biomarker_candidates.tsv", sep="\t", index=False)
print(f"\n  Saved → {OUT_DIR}/stable_biomarker_candidates.tsv")

# ── summary ───────────────────────────────────────────────────────────────────
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{OUT_DIR}/ml_summary.tsv", sep="\t", index=False)

print("\n" + "="*60)
print("ML SUMMARY")
print("="*60)
print(summary_df.to_string(index=False))

# top stable candidates
if len(stable_df) > 0:
    print(f"\nTop stable biomarker candidates (≥3 comparisons):")
    show = stable_df[stable_df["n_comparisons"] >= 3][
        ["gene_id_clean", "n_comparisons", "log2FC_vs_NAT", "padj_vs_NAT"]
    ].head(20)
    if len(show):
        print(show.to_string(index=False))
    else:
        print("  (none selected in ≥3 comparisons — see ≥2 comparison candidates)")
        show2 = stable_df[["gene_id_clean", "n_comparisons",
                            "log2FC_vs_NAT", "padj_vs_NAT"]].head(20)
        print(show2.to_string(index=False))

print(f"\nAll results saved to: {OUT_DIR}")
