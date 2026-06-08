"""
External validation — GSE156451 (CRC RNA-seq, 72 paired Tumor/Normal, FPKM)
Validates CDH3, CLDN1, ETV4, KRT80 expression difference in independent cohort.
"""

import os, gzip, glob, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# ── paths ─────────────────────────────────────────────────────────────────────
BASE     = "/Users/georgiatechgogofinn/Desktop/COAD_control_stability_project"
DATA_DIR = f"{BASE}/data_raw/GSE156451"
OUT_DIR  = f"{BASE}/results/external_validation"
os.makedirs(OUT_DIR, exist_ok=True)

TARGET_GENES = ["CDH3", "CLDN1", "ETV4", "KRT80"]

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Read all per-sample FPKM files
# ══════════════════════════════════════════════════════════════════════════════
print("="*60)
print("STEP 1: Parsing FPKM files")
print("="*60)

files = sorted(glob.glob(f"{DATA_DIR}/GSM*RNA.txt.gz"))
print(f"  Found {len(files)} sample files")

records = []
for fpath in files:
    fname   = os.path.basename(fpath)                       # GSM4731674_T1-RNA.txt.gz
    gsm_id  = fname.split("_")[0]
    label   = fname.split("_")[1].replace("-RNA.txt.gz","") # T1, T2, N1, N2 …
    group   = "Tumor" if label.startswith("T") else "Normal"

    gene_vals = {}
    with gzip.open(fpath, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2: continue
            gene = parts[0]
            if gene in TARGET_GENES:
                try:
                    gene_vals[gene] = float(parts[1])
                except ValueError:
                    pass

    row = {"gsm_id": gsm_id, "sample": label, "group": group}
    row.update(gene_vals)
    records.append(row)

df = pd.DataFrame(records)
df = df.dropna(subset=TARGET_GENES, how="all")

print(f"\n  Parsed: {len(df)} samples")
print(f"  Tumor : {sum(df['group']=='Tumor')}")
print(f"  Normal: {sum(df['group']=='Normal')}")
print(f"\n  Preview:")
print(df[["sample","group"] + TARGET_GENES].head(6).to_string(index=False))

# ── check for paired samples ──────────────────────────────────────────────────
# sample labels: T1…T72, N1…N72  → extract patient number
df["patient_num"] = df["sample"].str.extract(r"(\d+)").astype(int)
tumor_patients  = set(df[df["group"]=="Tumor"]["patient_num"])
normal_patients = set(df[df["group"]=="Normal"]["patient_num"])
paired_patients = tumor_patients & normal_patients
print(f"\n  Paired patients (T+N both present): {len(paired_patients)}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Log2 transform (FPKM → log2(FPKM + 1))
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("STEP 2: log2(FPKM + 1) transformation")
print("="*60)

for g in TARGET_GENES:
    df[f"{g}_log2"] = np.log2(df[g].clip(lower=0) + 1)

print("  FPKM range before transform:")
for g in TARGET_GENES:
    print(f"    {g:6s}: min={df[g].min():.2f}  max={df[g].max():.2f}  "
          f"mean={df[g].mean():.2f}")

# ── save metadata ─────────────────────────────────────────────────────────────
df.to_csv(f"{OUT_DIR}/GSE156451_expression.tsv", sep="\t", index=False)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 & 4 & 5: Statistics + Boxplots
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("STEP 3-5: Statistics + Boxplots")
print("="*60)

tumor_df  = df[df["group"] == "Tumor"]
normal_df = df[df["group"] == "Normal"]

stats_rows = []
for gene in TARGET_GENES:
    col = f"{gene}_log2"
    t_vals = tumor_df[col].dropna().values
    n_vals = normal_df[col].dropna().values

    t_stat, pval = stats.ttest_ind(t_vals, n_vals, equal_var=False)
    log2fc    = t_vals.mean() - n_vals.mean()
    direction = "Up" if log2fc > 0 else "Down"
    sig = ("***" if pval < 0.001 else "**" if pval < 0.01
           else "*" if pval < 0.05 else "ns")

    # paired t-test for paired samples
    paired_df = df[df["patient_num"].isin(paired_patients)].copy()
    t_paired  = paired_df[paired_df["group"]=="Tumor"][col].values
    n_paired  = paired_df[paired_df["group"]=="Normal"][col].values
    if len(t_paired) == len(n_paired) and len(t_paired) > 3:
        _, pval_paired = stats.ttest_rel(t_paired, n_paired)
    else:
        pval_paired = np.nan

    print(f"  {gene:6s}: log2FC={log2fc:+.3f}  p={pval:.3e} {sig}"
          f"  paired-p={pval_paired:.3e}"
          f"  (T n={len(t_vals)}, N n={len(n_vals)})")

    stats_rows.append({
        "gene"        : gene,
        "n_tumor"     : len(t_vals),
        "n_normal"    : len(n_vals),
        "mean_tumor"  : round(t_vals.mean(), 4),
        "mean_normal" : round(n_vals.mean(), 4),
        "log2FC"      : round(log2fc, 4),
        "pvalue"      : pval,
        "pvalue_paired": pval_paired,
        "direction"   : direction,
        "significance": sig,
    })

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(f"{OUT_DIR}/expression_stats.tsv", sep="\t", index=False)

# ── individual boxplots ───────────────────────────────────────────────────────
for row in stats_rows:
    gene    = row["gene"]
    col     = f"{gene}_log2"
    t_vals  = tumor_df[col].dropna().values
    n_vals  = normal_df[col].dropna().values
    pval    = row["pvalue"]
    log2fc  = row["log2FC"]
    sig     = row["significance"]

    fig, ax = plt.subplots(figsize=(5, 6.5))
    bp = ax.boxplot(
        [n_vals, t_vals],
        labels=["Normal", "Tumor"],
        patch_artist=True,
        medianprops=dict(color="black", linewidth=2.5),
        whiskerprops=dict(linewidth=1.5),
        capprops=dict(linewidth=1.5),
        flierprops=dict(marker="o", markersize=3,
                        markerfacecolor="gray", alpha=0.4),
        widths=0.55,
    )
    for patch, color in zip(bp["boxes"], ["#2980B9", "#E74C3C"]):
        patch.set_facecolor(color); patch.set_alpha(0.72)

    np.random.seed(42)
    for i, (vals, color) in enumerate(zip([n_vals, t_vals],
                                           ["#2980B9", "#E74C3C"]), 1):
        jx = np.random.uniform(-0.18, 0.18, len(vals))
        ax.scatter(np.full(len(vals), i) + jx, vals,
                   color=color, alpha=0.40, s=20, zorder=3)

    # significance bracket
    y_max = max(np.percentile(t_vals, 97), np.percentile(n_vals, 97))
    y_top = y_max + (y_max - min(n_vals.min(), t_vals.min())) * 0.08
    ax.plot([1, 1, 2, 2], [y_max, y_top, y_top, y_max],
            color="black", lw=1.3)
    ax.text(1.5, y_top + (y_top - y_max) * 0.1, sig,
            ha="center", va="bottom", fontsize=14, fontweight="bold")

    p_str = (f"p = {pval:.4f}" if pval >= 0.0001 else "p < 0.0001")
    ax.set_title(f"{gene}  |  GSE156451 (CRC)\n"
                 f"log₂FC = {log2fc:+.3f}   {p_str}  {sig}",
                 fontsize=11, fontweight="bold")
    ax.set_ylabel("log₂(FPKM + 1)", fontsize=11)
    ax.set_xlabel(f"Normal (n={len(n_vals)})  vs  Tumor (n={len(t_vals)})",
                  fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0.3, 2.7)
    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/boxplot_{gene}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved boxplot_{gene}.png")

# ── 4-panel combined boxplot ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(18, 6.5), sharey=False)

for ax, row in zip(axes, stats_rows):
    gene   = row["gene"]
    col    = f"{gene}_log2"
    t_vals = tumor_df[col].dropna().values
    n_vals = normal_df[col].dropna().values
    pval   = row["pvalue"]
    log2fc = row["log2FC"]
    sig    = row["significance"]

    bp = ax.boxplot(
        [n_vals, t_vals], labels=["Normal", "Tumor"],
        patch_artist=True,
        medianprops=dict(color="black", linewidth=2.5),
        widths=0.52,
        flierprops=dict(marker="o", markersize=3,
                        markerfacecolor="gray", alpha=0.35),
    )
    for patch, color in zip(bp["boxes"], ["#2980B9", "#E74C3C"]):
        patch.set_facecolor(color); patch.set_alpha(0.72)

    np.random.seed(42)
    for i, (vals, color) in enumerate(zip([n_vals, t_vals],
                                           ["#2980B9","#E74C3C"]), 1):
        jx = np.random.uniform(-0.16, 0.16, len(vals))
        ax.scatter(np.full(len(vals), i) + jx, vals,
                   color=color, alpha=0.35, s=16, zorder=3)

    y_max = max(np.percentile(t_vals, 97), np.percentile(n_vals, 97))
    rng   = y_max - min(n_vals.min(), t_vals.min())
    y_top = y_max + rng * 0.08
    ax.plot([1,1,2,2],[y_max, y_top, y_top, y_max], color="black", lw=1.2)
    ax.text(1.5, y_top + rng*0.01, sig,
            ha="center", va="bottom", fontsize=13, fontweight="bold")

    p_str = f"p={pval:.3e}"
    ax.set_title(f"**{gene}**\nlog₂FC={log2fc:+.2f}\n{p_str}  {sig}",
                 fontsize=10, fontweight="bold")
    ax.set_ylabel("log₂(FPKM+1)" if ax == axes[0] else "", fontsize=10)
    ax.set_xlabel(f"N n={len(n_vals)} | T n={len(t_vals)}", fontsize=9)
    ax.spines[["top","right"]].set_visible(False)
    ax.set_xlim(0.3, 2.7)

plt.suptitle("External Validation — GSE156451 (CRC, 72 patients)\n"
             "Tumor vs Paired Normal, log₂(FPKM+1)",
             fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/boxplot_all_genes.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\n  Saved boxplot_all_genes.png")

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("EXTERNAL VALIDATION SUMMARY — GSE156451")
print("="*60)
print(f"Dataset: Genome-wide enhancer landscape in CRC")
print(f"Samples: {len(tumor_df)} Tumor  |  {len(normal_df)} Normal  "
      f"|  {len(paired_patients)} paired")
print()
print(stats_df[["gene","n_tumor","n_normal","log2FC",
                "pvalue","pvalue_paired","significance","direction"]
              ].to_string(index=False))
print(f"\nAll results saved to: {OUT_DIR}")
