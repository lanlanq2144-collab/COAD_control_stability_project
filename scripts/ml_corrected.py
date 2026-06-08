"""
Nested CV ML pipeline — fixes feature-selection leakage.

Problem with original ml_analysis.py:
  - Top 500 DEGs were selected using ALL samples (tumor + control)
  - THEN cross-validation was applied to those pre-selected features
  - Test-set samples already influenced which features were chosen → data leakage
  - AUC = 1.0 could be partially inflated by this

Corrected approach (this script):
  Outer loop  (5-fold):  split whole dataset → train / test
    Inner loop (3-fold):  on TRAIN ONLY
        1. Recompute DEGs (t-test + BH) on train samples
        2. Select top 500 significant DEGs
        3. Fit LASSO-CV (hyperparameter search inside inner loop)
        4. Fit Random Forest
    Evaluate LASSO & RF on held-out TEST fold (unseen)
  → Aggregate test-fold probabilities → ROC-AUC

Biomarker discovery (same-data LASSO at C=0.5, labeled clearly):
  For gene lists, we still run full-data LASSO after corrected DEG selection.
  The AUC shown for this is the nested-CV AUC (honest).
  Gene selection is noted as exploratory / not for AUC claims.

Results saved to results/ml_corrected/
Comparison table saved to results/ml_corrected/comparison_leaked_vs_corrected.tsv
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve

# ── paths ─────────────────────────────────────────────────────────────────────
BASE        = "/Users/georgiatechgogofinn/Desktop/COAD_control_stability_project"
TPM_FILE    = f"{BASE}/data_processed/expression/coad_expression_tpm_genes_by_samples.tsv"
COMBAT_FILE = f"{BASE}/data_processed/expression/coad_expression_combat_corrected.tsv"
LABELS_FILE = f"{BASE}/data_processed/metadata/coad_sample_labels.tsv"
ML_OLD_DIR  = f"{BASE}/results/ml"          # original (leaked) results
OUT_DIR     = f"{BASE}/results/ml_corrected"
os.makedirs(OUT_DIR, exist_ok=True)

# ── settings ──────────────────────────────────────────────────────────────────
N_TOP_DEG   = 500
N_OUTER     = 5     # outer CV folds (honest AUC)
N_INNER     = 3     # inner CV folds (LASSO hyperparameter search)
N_TREES     = 300
C_BIO       = 0.5   # for exploratory biomarker gene list
DEG_LOG2FC  = 1.0
DEG_FDR     = 0.05

# ── load data ─────────────────────────────────────────────────────────────────
print("Loading data...")
labels      = pd.read_csv(LABELS_FILE, sep="\t", index_col=0)
label_map   = labels["group"].to_dict()

expr_orig   = pd.read_csv(TPM_FILE,    sep="\t", index_col=0)
expr_combat = pd.read_csv(COMBAT_FILE, sep="\t", index_col=0)

common_genes = expr_orig.index.intersection(expr_combat.index)
expr_orig    = expr_orig.loc[common_genes]
expr_combat  = expr_combat.loc[common_genes]

tumor_cols   = [c for c in expr_orig.columns if label_map.get(c) == "TCGA_COAD_Tumor"]
expr_merged  = pd.concat([expr_orig[tumor_cols], expr_combat], axis=1)

print(f"  Genes: {len(common_genes):,}  |  Tumor samples: {len(tumor_cols)}")


def get_cols(expr_df, groups):
    if isinstance(groups, str): groups = [groups]
    return [c for c in expr_df.columns if label_map.get(c) in groups]


# ── comparison definitions ────────────────────────────────────────────────────
COMPARISONS = [
    dict(name="Tumor_vs_TCGA_NAT",
         expr=expr_orig,   tumor="TCGA_COAD_Tumor",
         control="TCGA_COAD_NAT"),
    dict(name="Tumor_vs_GTEx_All",
         expr=expr_merged, tumor="TCGA_COAD_Tumor",
         control=["GTEx_Colon_Transverse", "GTEx_Colon_Sigmoid"]),
    dict(name="Tumor_vs_GTEx_Colon_Transverse",
         expr=expr_merged, tumor="TCGA_COAD_Tumor",
         control="GTEx_Colon_Transverse"),
    dict(name="Tumor_vs_GTEx_Colon_Sigmoid",
         expr=expr_merged, tumor="TCGA_COAD_Tumor",
         control="GTEx_Colon_Sigmoid"),
]


# ── DEG function (applied to train-fold samples only) ─────────────────────────
def compute_deg_on_train(expr_mat, tumor_idx, control_idx,
                          n_top=N_TOP_DEG, lfc_thr=DEG_LOG2FC, fdr_thr=DEG_FDR):
    """
    Run t-test + BH FDR on training samples only.
    Returns array of gene indices (top n_top by padj).
    """
    T = expr_mat[:, tumor_idx]    # genes × train-tumor
    C = expr_mat[:, control_idx]  # genes × train-control

    _, pval  = stats.ttest_ind(T, C, axis=1, equal_var=False)
    log2fc   = T.mean(axis=1) - C.mean(axis=1)

    _, padj, _, _ = multipletests(pval, method="fdr_bh")

    sig_mask = (padj < fdr_thr) & (np.abs(log2fc) > lfc_thr)
    # rank significant genes by padj; fill non-sig with 1.0
    ranked   = np.where(sig_mask, padj, 1.0)
    top_idx  = np.argsort(ranked)[:n_top]
    return top_idx


# ═════════════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═════════════════════════════════════════════════════════════════════════════
summary_rows = []

for comp in COMPARISONS:
    name = comp["name"]
    print(f"\n{'='*60}")
    print(f"  {name}")
    print(f"{'='*60}")

    expr      = comp["expr"]
    t_cols    = get_cols(expr, comp["tumor"])
    c_cols    = get_cols(expr, comp["control"])
    all_cols  = t_cols + c_cols
    y_all     = np.array([1]*len(t_cols) + [0]*len(c_cols))

    # Full expression matrix: numpy array  (genes × samples)
    X_full_raw = expr[all_cols].values.astype(float)  # genes × samples
    n_tumor    = len(t_cols)
    n_control  = len(c_cols)
    n_samples  = len(all_cols)
    print(f"  Samples: {n_tumor} Tumor + {n_control} Control = {n_samples}")

    # ── nested CV ─────────────────────────────────────────────────────────────
    outer_cv = StratifiedKFold(n_splits=N_OUTER, shuffle=True, random_state=42)

    lasso_probs = np.zeros(n_samples)
    rf_probs    = np.zeros(n_samples)
    lasso_n_genes_folds = []

    print(f"  Running nested CV ({N_OUTER}-outer × {N_INNER}-inner)...")

    for fold, (train_idx, test_idx) in enumerate(
            outer_cv.split(np.zeros(n_samples), y_all)):

        y_train = y_all[train_idx]
        y_test  = y_all[test_idx]

        # indices of tumor/control within train fold
        tumor_train_mask   = train_idx < n_tumor
        control_train_mask = train_idx >= n_tumor
        tumor_train_local  = np.where(tumor_train_mask)[0]
        control_train_local= np.where(control_train_mask)[0]

        # ── STEP 1: DEG on train only ──────────────────────────────────────
        top_gene_idx = compute_deg_on_train(
            X_full_raw,
            train_idx[tumor_train_mask],    # absolute column indices
            train_idx[control_train_mask],
        )

        # ── STEP 2: build feature matrices ────────────────────────────────
        X_train_raw = X_full_raw[top_gene_idx][:, train_idx].T  # samples × genes
        X_test_raw  = X_full_raw[top_gene_idx][:, test_idx].T

        scaler      = StandardScaler()
        X_train     = scaler.fit_transform(X_train_raw)
        X_test      = scaler.transform(X_test_raw)

        # ── STEP 3a: LASSO with inner CV (hyperparameter search on train) ──
        inner_cv   = StratifiedKFold(n_splits=N_INNER, shuffle=True, random_state=42)
        lasso_cv   = LogisticRegressionCV(
            penalty="l1", solver="liblinear",
            Cs=np.logspace(-3, 2, 20),
            cv=inner_cv, scoring="roc_auc",
            class_weight="balanced", max_iter=2000, random_state=42,
        )
        lasso_cv.fit(X_train, y_train)
        n_sel = np.sum(lasso_cv.coef_[0] != 0)
        lasso_n_genes_folds.append(n_sel)

        # ── STEP 3b: RF trained on train fold ─────────────────────────────
        rf = RandomForestClassifier(n_estimators=N_TREES, class_weight="balanced",
                                    random_state=42, n_jobs=-1)
        rf.fit(X_train, y_train)

        # ── STEP 4: predict on TEST fold (never seen during training) ──────
        lasso_probs[test_idx] = lasso_cv.predict_proba(X_test)[:, 1]
        rf_probs[test_idx]    = rf.predict_proba(X_test)[:, 1]

        print(f"    fold {fold+1}: LASSO selected {n_sel} genes")

    # ── aggregate AUC across folds ────────────────────────────────────────────
    lasso_auc_corr = roc_auc_score(y_all, lasso_probs)
    rf_auc_corr    = roc_auc_score(y_all, rf_probs)
    print(f"  Corrected LASSO AUC = {lasso_auc_corr:.4f}")
    print(f"  Corrected RF    AUC = {rf_auc_corr:.4f}")
    print(f"  Avg LASSO genes/fold: {np.mean(lasso_n_genes_folds):.1f}")

    # ── ROC curve ─────────────────────────────────────────────────────────────
    fpr_l, tpr_l, _ = roc_curve(y_all, lasso_probs)
    fpr_r, tpr_r, _ = roc_curve(y_all, rf_probs)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(fpr_l, tpr_l, lw=2.2, color="#E74C3C",
            label=f"LASSO  (AUC = {lasso_auc_corr:.3f})")
    ax.plot(fpr_r, tpr_r, lw=2.2, color="#2980B9",
            label=f"Random Forest (AUC = {rf_auc_corr:.3f})")
    ax.plot([0, 1], [0, 1], lw=1, ls="--", color="gray", alpha=0.6)
    ax.fill_between(fpr_l, tpr_l, alpha=0.07, color="#E74C3C")
    ax.fill_between(fpr_r, tpr_r, alpha=0.07, color="#2980B9")
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title(
        f"ROC Curve (Corrected Nested CV)\n{name.replace('_', ' ')}",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=10)
    ax.set_xlim([-0.01, 1.01]); ax.set_ylim([-0.01, 1.02])
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/{name}_roc_corrected.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ── exploratory biomarker genes (full data, labeled as exploratory) ───────
    # Recompute DEGs on ALL samples (for gene list only, not for AUC claims)
    top_gene_idx_full = compute_deg_on_train(
        X_full_raw,
        np.arange(n_tumor),               # all tumor columns
        np.arange(n_tumor, n_samples),    # all control columns
    )
    X_bio_raw = X_full_raw[top_gene_idx_full].T
    scaler_bio = StandardScaler()
    X_bio = scaler_bio.fit_transform(X_bio_raw)

    lasso_bio = LogisticRegression(
        penalty="l1", solver="liblinear", C=C_BIO,
        class_weight="balanced", max_iter=2000, random_state=42,
    )
    lasso_bio.fit(X_bio, y_all)
    coef_bio = lasso_bio.coef_[0]
    nz_bio   = np.where(coef_bio != 0)[0]
    bio_genes_idx   = top_gene_idx_full[nz_bio]
    bio_gene_ids    = expr.index[bio_genes_idx].tolist()
    bio_gene_coefs  = coef_bio[nz_bio]

    bio_df = (pd.DataFrame({
        "gene_id":    bio_gene_ids,
        "lasso_coef": bio_gene_coefs,
        "abs_coef":   np.abs(bio_gene_coefs),
    }).sort_values("abs_coef", ascending=False))
    bio_df.to_csv(f"{OUT_DIR}/{name}_lasso_genes_corrected.tsv",
                  sep="\t", index=False)
    print(f"  Exploratory biomarker genes (full-data LASSO C=0.5): {len(bio_df)}")

    # ── read old leaked AUC for comparison ────────────────────────────────────
    old_summary_path = f"{ML_OLD_DIR}/ml_summary.tsv"
    lasso_auc_old = rf_auc_old = np.nan
    if os.path.exists(old_summary_path):
        old_sum = pd.read_csv(old_summary_path, sep="\t")
        row = old_sum[old_sum["comparison"] == name]
        if len(row):
            lasso_auc_old = row["lasso_auc"].values[0]
            rf_auc_old    = row["rf_auc"].values[0]

    summary_rows.append({
        "comparison"         : name,
        "n_tumor"            : n_tumor,
        "n_control"          : n_control,
        "lasso_auc_leaked"   : round(lasso_auc_old,    4),
        "rf_auc_leaked"      : round(rf_auc_old,       4),
        "lasso_auc_corrected": round(lasso_auc_corr,   4),
        "rf_auc_corrected"   : round(rf_auc_corr,      4),
        "lasso_auc_delta"    : round(lasso_auc_old - lasso_auc_corr, 4)
                               if not np.isnan(lasso_auc_old) else np.nan,
        "rf_auc_delta"       : round(rf_auc_old    - rf_auc_corr,    4)
                               if not np.isnan(rf_auc_old)    else np.nan,
        "avg_lasso_genes_per_fold": round(np.mean(lasso_n_genes_folds), 1),
        "bio_lasso_genes"    : len(bio_df),
    })

# ── summary table ─────────────────────────────────────────────────────────────
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f"{OUT_DIR}/comparison_leaked_vs_corrected.tsv",
                  sep="\t", index=False)

# ── comparison bar chart ───────────────────────────────────────────────────────
print("\nDrawing comparison plot...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
comp_names = [r["comparison"].replace("Tumor_vs_", "vs ").replace("_", " ")
              for r in summary_rows]
x = np.arange(len(comp_names))
w = 0.35

for ax, model, col_old, col_new in [
    (axes[0], "LASSO", "lasso_auc_leaked", "lasso_auc_corrected"),
    (axes[1], "Random Forest", "rf_auc_leaked", "rf_auc_corrected"),
]:
    old_vals = summary_df[col_old].values
    new_vals = summary_df[col_new].values

    bars_old = ax.bar(x - w/2, old_vals, w, label="Leaked (original)",
                      color="#E74C3C", alpha=0.8, edgecolor="white")
    bars_new = ax.bar(x + w/2, new_vals, w, label="Corrected (nested CV)",
                      color="#2980B9", alpha=0.8, edgecolor="white")

    for bar, val in zip(bars_old, old_vals):
        ax.text(bar.get_x() + bar.get_width()/2,
                min(val + 0.003, 1.005),
                f"{val:.3f}", ha="center", va="bottom", fontsize=8,
                color="#C0392B", fontweight="bold")
    for bar, val in zip(bars_new, new_vals):
        ax.text(bar.get_x() + bar.get_width()/2,
                min(val + 0.003, 1.005),
                f"{val:.3f}", ha="center", va="bottom", fontsize=8,
                color="#1A6A99", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(comp_names, rotation=18, ha="right", fontsize=9)
    ax.set_ylabel("Cross-validated ROC-AUC", fontsize=11)
    ax.set_ylim(0.5, 1.08)
    ax.set_title(f"{model}: Leaked vs. Corrected AUC", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.axhline(1.0, ls="--", color="gray", lw=0.8, alpha=0.5)

plt.suptitle("Effect of Feature-Selection Leakage on Reported AUC",
             fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/comparison_leaked_vs_corrected.png",
            dpi=150, bbox_inches="tight")
plt.close(fig)

# ── print summary ─────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("LEAKAGE CORRECTION SUMMARY")
print("="*70)
print(summary_df[[
    "comparison", "lasso_auc_leaked", "lasso_auc_corrected",
    "lasso_auc_delta", "rf_auc_leaked", "rf_auc_corrected", "rf_auc_delta",
]].to_string(index=False))
print(f"\nAll results saved to: {OUT_DIR}")
print("\nInterpretation guide:")
print("  delta > 0  → original AUC was inflated by leakage")
print("  delta ≈ 0  → leakage had minimal effect (signal is genuinely strong)")
print("  delta < 0  → should not happen (corrected AUC higher would be unusual)")
