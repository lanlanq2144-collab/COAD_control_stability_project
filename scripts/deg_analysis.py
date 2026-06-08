"""
DEG analysis for TCGA COAD project
4 comparisons: Tumor vs TCGA_NAT, Tumor vs GTEx (all/Transverse/Sigmoid)
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

# ── paths ────────────────────────────────────────────────────────────────────
BASE = "/Users/georgiatechgogofinn/Desktop/COAD_control_stability_project"
TPM_FILE     = f"{BASE}/data_processed/expression/coad_expression_tpm_genes_by_samples.tsv"
COMBAT_FILE  = f"{BASE}/data_processed/expression/coad_expression_combat_corrected.tsv"
LABELS_FILE  = f"{BASE}/data_processed/metadata/coad_sample_labels.tsv"
OUT_DIR      = f"{BASE}/results/deg"
os.makedirs(OUT_DIR, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
print("Loading data...")
labels = pd.read_csv(LABELS_FILE, sep="\t", index_col=0)

expr_orig   = pd.read_csv(TPM_FILE,    sep="\t", index_col=0)
expr_combat = pd.read_csv(COMBAT_FILE, sep="\t", index_col=0)

print(f"  Original  : {expr_orig.shape[1]} samples × {expr_orig.shape[0]} genes")
print(f"  ComBat    : {expr_combat.shape[1]} samples × {expr_combat.shape[0]} genes")
print(f"  Labels    : {labels.shape[0]} samples")

# helper: map sample → group
label_map = labels["group"].to_dict()


def get_samples(expr_df, group):
    cols = [c for c in expr_df.columns if label_map.get(c) == group]
    return cols


def run_deg(name, expr_df, tumor_cols, control_cols,
            log2fc_thr=1.0, padj_thr=0.05):
    """t-test + BH FDR + volcano plot."""
    print(f"\n[{name}]  Tumor n={len(tumor_cols)}  Control n={len(control_cols)}")

    tumor   = expr_df[tumor_cols].values   # genes × samples
    control = expr_df[control_cols].values

    # log2FC (data is already log2-transformed → difference of means)
    log2fc = tumor.mean(axis=1) - control.mean(axis=1)

    # Welch t-test gene by gene
    t_stat, pval = stats.ttest_ind(tumor, control, axis=1, equal_var=False)

    # BH FDR correction
    _, padj, _, _ = multipletests(pval, method="fdr_bh")

    # build result table
    result = pd.DataFrame({
        "gene_id" : expr_df.index,
        "log2FC"  : log2fc,
        "pvalue"  : pval,
        "padj"    : padj,
        "mean_tumor"   : tumor.mean(axis=1),
        "mean_control" : control.mean(axis=1),
    })
    result = result.sort_values("padj")

    sig_up   = ((result["padj"] < padj_thr) & (result["log2FC"] >  log2fc_thr)).sum()
    sig_down = ((result["padj"] < padj_thr) & (result["log2FC"] < -log2fc_thr)).sum()
    print(f"  Up  (log2FC>{log2fc_thr}, padj<{padj_thr}): {sig_up}")
    print(f"  Down(log2FC<-{log2fc_thr}, padj<{padj_thr}): {sig_down}")

    # save table
    safe_name = name.replace(" ", "_").replace("/", "_")
    csv_path  = f"{OUT_DIR}/{safe_name}_DEG.tsv"
    result.to_csv(csv_path, sep="\t", index=False)
    print(f"  Saved → {csv_path}")

    # ── volcano plot ─────────────────────────────────────────────────────────
    neg_log10_padj = -np.log10(result["padj"].clip(lower=1e-300))

    up_mask   = (result["padj"] < padj_thr) & (result["log2FC"] >  log2fc_thr)
    down_mask = (result["padj"] < padj_thr) & (result["log2FC"] < -log2fc_thr)
    ns_mask   = ~(up_mask | down_mask)

    fig, ax = plt.subplots(figsize=(8, 7))

    ax.scatter(result.loc[ns_mask,   "log2FC"], neg_log10_padj[ns_mask],
               s=3, color="#AAAAAA", alpha=0.4, rasterized=True, label="NS")
    ax.scatter(result.loc[up_mask,   "log2FC"], neg_log10_padj[up_mask],
               s=6, color="#E74C3C", alpha=0.7, rasterized=True, label=f"Up ({sig_up})")
    ax.scatter(result.loc[down_mask, "log2FC"], neg_log10_padj[down_mask],
               s=6, color="#2980B9", alpha=0.7, rasterized=True, label=f"Down ({sig_down})")

    # threshold lines
    ax.axhline(-np.log10(padj_thr), color="black", lw=0.8, ls="--", alpha=0.6)
    ax.axvline( log2fc_thr,  color="black", lw=0.8, ls="--", alpha=0.6)
    ax.axvline(-log2fc_thr,  color="black", lw=0.8, ls="--", alpha=0.6)

    # label top 15 genes by significance
    top = result[up_mask | down_mask].head(15)
    for _, row in top.iterrows():
        ax.annotate(
            row["gene_id"].split(".")[0],
            xy=(row["log2FC"], -np.log10(row["padj"])),
            xytext=(3, 2), textcoords="offset points",
            fontsize=5, color="black",
        )

    ax.set_xlabel("log₂ Fold Change (Tumor / Control)", fontsize=12)
    ax.set_ylabel("-log₁₀ adjusted p-value (BH FDR)", fontsize=12)
    ax.set_title(f"Volcano Plot — {name}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, markerscale=2)

    plt.tight_layout()
    png_path = f"{OUT_DIR}/{safe_name}_volcano.png"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {png_path}")

    return result


# ════════════════════════════════════════════════════════════════════════════
# 1. Tumor vs TCGA_COAD_NAT  (original data)
# ════════════════════════════════════════════════════════════════════════════
tumor_orig_cols = get_samples(expr_orig, "TCGA_COAD_Tumor")
nat_cols        = get_samples(expr_orig, "TCGA_COAD_NAT")

deg1 = run_deg(
    "Tumor_vs_TCGA_NAT",
    expr_orig,
    tumor_orig_cols,
    nat_cols,
)

# ════════════════════════════════════════════════════════════════════════════
# 2–4. Tumor vs GTEx  (ComBat-corrected data)
#   Tumor samples are NOT in the ComBat file → taken from original expr,
#   then concatenated with ComBat-corrected GTEx columns.
# ════════════════════════════════════════════════════════════════════════════
# Tumor columns that exist in the original file
tumor_orig_cols_set = set(tumor_orig_cols)

# Build a merged expression matrix:
#   genes × (tumor_from_orig ∪ gtex_from_combat)
# Both files share the same gene index; align on gene_id.
common_genes = expr_orig.index.intersection(expr_combat.index)
print(f"\nCommon genes between orig & combat: {len(common_genes)}")

expr_tumor_aligned  = expr_orig.loc[common_genes, tumor_orig_cols]
expr_combat_aligned = expr_combat.loc[common_genes]

# Merged data for GTEx comparisons
expr_merged = pd.concat([expr_tumor_aligned, expr_combat_aligned], axis=1)

tumor_merged_cols = tumor_orig_cols  # same sample IDs, now in merged df

# 2. Tumor vs all GTEx
gtex_all_cols = get_samples(expr_combat_aligned, "GTEx_Colon_Transverse") + \
                get_samples(expr_combat_aligned, "GTEx_Colon_Sigmoid")

deg2 = run_deg(
    "Tumor_vs_GTEx_All",
    expr_merged,
    tumor_merged_cols,
    gtex_all_cols,
)

# 3. Tumor vs GTEx_Colon_Transverse
gtex_trans_cols = get_samples(expr_combat_aligned, "GTEx_Colon_Transverse")

deg3 = run_deg(
    "Tumor_vs_GTEx_Colon_Transverse",
    expr_merged,
    tumor_merged_cols,
    gtex_trans_cols,
)

# 4. Tumor vs GTEx_Colon_Sigmoid
gtex_sig_cols = get_samples(expr_combat_aligned, "GTEx_Colon_Sigmoid")

deg4 = run_deg(
    "Tumor_vs_GTEx_Colon_Sigmoid",
    expr_merged,
    tumor_merged_cols,
    gtex_sig_cols,
)

# ── summary ──────────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("SUMMARY (|log2FC|>1, padj<0.05)")
print("="*60)
for name, df in [
    ("Tumor vs TCGA_NAT",             deg1),
    ("Tumor vs GTEx All",             deg2),
    ("Tumor vs GTEx_Colon_Transverse",deg3),
    ("Tumor vs GTEx_Colon_Sigmoid",   deg4),
]:
    up   = ((df["padj"] < 0.05) & (df["log2FC"] >  1)).sum()
    down = ((df["padj"] < 0.05) & (df["log2FC"] < -1)).sum()
    print(f"  {name:<35}  Up:{up:>5}  Down:{down:>5}  Total:{up+down:>6}")

print(f"\nAll results saved to: {OUT_DIR}")
