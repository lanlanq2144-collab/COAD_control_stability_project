"""
4-step biomarker stability analysis
1. 4-set Venn diagram
2. Stability scores
3. Genes by stability score
4. Final candidates (score >= 3)
"""

import os
import numpy as np
import pandas as pd
from itertools import combinations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib_venn import venn2, venn2_circles
import warnings
warnings.filterwarnings("ignore")

# ── paths ─────────────────────────────────────────────────────────────────────
BASE    = "/Users/georgiatechgogofinn/Desktop/COAD_control_stability_project"
ML_DIR  = f"{BASE}/results/ml"
DEG_DIR = f"{BASE}/results/deg"
OUT_DIR = f"{BASE}/results/stability"
os.makedirs(OUT_DIR, exist_ok=True)

# ── short labels ──────────────────────────────────────────────────────────────
COMP_INFO = {
    "NAT"       : ("Tumor_vs_TCGA_NAT",             "#E74C3C"),
    "GTEx_All"  : ("Tumor_vs_GTEx_All",             "#2980B9"),
    "Transverse": ("Tumor_vs_GTEx_Colon_Transverse", "#27AE60"),
    "Sigmoid"   : ("Tumor_vs_GTEx_Colon_Sigmoid",   "#8E44AD"),
}

# ── 1. Load biomarker sets ─────────────────────────────────────────────────────
print("Loading biomarker candidate lists...")
gene_sets = {}
gene_details = {}   # gene_id_clean → row from any file

for label, (fname, color) in COMP_INFO.items():
    path = f"{ML_DIR}/{fname}_biomarker_candidates.tsv"
    df   = pd.read_csv(path, sep="\t")
    # use versioned gene_id as key for uniqueness
    ids  = set(df["gene_id"].tolist())
    gene_sets[label] = ids
    for _, row in df.iterrows():
        gid = row["gene_id"]
        if gid not in gene_details:
            gene_details[gid] = row.to_dict()
    print(f"  {label:12s}: {len(ids):3d} genes")

all_genes = set().union(*gene_sets.values())
print(f"\n  Total unique genes across all comparisons: {len(all_genes)}")

# ── 2. Stability score matrix ─────────────────────────────────────────────────
print("\nBuilding stability score matrix...")
labels_order = list(COMP_INFO.keys())
colors_order = [COMP_INFO[l][1] for l in labels_order]

records = []
for gene in sorted(all_genes):
    membership = {l: (gene in gene_sets[l]) for l in labels_order}
    score      = sum(membership.values())
    det        = gene_details.get(gene, {})
    records.append({
        "gene_id"      : gene,
        "gene_id_clean": gene.split(".")[0],
        "stability_score": score,
        **{f"in_{l}": membership[l] for l in labels_order},
        "log2FC"       : det.get("log2FC", np.nan),
        "padj"         : det.get("padj",   np.nan),
        "direction"    : det.get("direction", "NA"),
        "lasso_coef"   : det.get("lasso_coef", np.nan),
        "rf_importance": det.get("rf_importance", np.nan),
    })

score_df = pd.DataFrame(records).sort_values(
    ["stability_score", "rf_importance"], ascending=[False, False]
)
score_df.to_csv(f"{OUT_DIR}/stability_scores_all.tsv", sep="\t", index=False)
print(f"  Saved stability_scores_all.tsv ({len(score_df)} genes)")

# ── count per score ───────────────────────────────────────────────────────────
score_counts = score_df["stability_score"].value_counts().sort_index(ascending=False)
print("\n  Stability score distribution:")
for sc, cnt in score_counts.items():
    print(f"    Score {sc}/4 : {cnt} genes")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: 4-set Venn diagram (custom ellipse layout)
# ══════════════════════════════════════════════════════════════════════════════
print("\nDrawing 4-set Venn diagram...")

def venn4_custom(sets, labels, colors, title, out_path):
    """
    4-set Venn with 4 overlapping ellipses + region counts.
    Layout: 4 ellipses arranged in a 2×2 offset grid.
    """
    # ellipse centers & radii (axes units)
    params = [
        (-0.45,  0.25, 0.72, 0.45, -30),   # NAT       top-left
        ( 0.45,  0.25, 0.72, 0.45,  30),   # GTEx_All  top-right
        (-0.25, -0.25, 0.72, 0.45, -15),   # Transverse bottom-left
        ( 0.25, -0.25, 0.72, 0.45,  15),   # Sigmoid    bottom-right
    ]

    # precompute all 2^4 - 1 region memberships
    n  = len(sets)
    # for each non-empty subset of {0,1,2,3}: count exclusive members
    region_data = {}
    for mask in range(1, 1 << n):
        in_idx  = [i for i in range(n) if mask & (1 << i)]
        out_idx = [i for i in range(n) if not (mask & (1 << i))]
        members = sets[in_idx[0]].copy()
        for i in in_idx[1:]:
            members &= sets[i]
        for i in out_idx:
            members -= sets[i]
        region_data[mask] = len(members)

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.1, 1.1)
    ax.set_aspect("equal")
    ax.axis("off")

    # draw ellipses
    for i, (cx, cy, w, h, angle) in enumerate(params):
        ell = Ellipse((cx, cy), w, h, angle=angle,
                      facecolor=colors[i], edgecolor="white",
                      linewidth=2, alpha=0.35)
        ax.add_patch(ell)

    # label positions (outside ellipses)
    label_pos = [(-1.15, 0.65), (1.15, 0.65), (-1.0, -0.75), (1.0, -0.75)]
    for i, (lx, ly) in enumerate(label_pos):
        ax.text(lx, ly, labels[i], ha="center", va="center",
                fontsize=11, fontweight="bold", color=colors[i])
        ax.annotate("", xy=label_pos[i],
                    xytext=(params[i][0], params[i][1]),
                    arrowprops=dict(arrowstyle="-", color=colors[i],
                                   lw=1.2, alpha=0.6))

    # region count positions (approximate centers of each region)
    # single sets
    region_text = {
        0b0001: (-0.95,  0.45),   # NAT only
        0b0010: ( 0.95,  0.45),   # GTEx_All only
        0b0100: (-0.75, -0.55),   # Transverse only
        0b1000: ( 0.75, -0.55),   # Sigmoid only
        # pairs
        0b0011: ( 0.00,  0.60),   # NAT ∩ GTEx_All
        0b0101: (-0.72,  0.00),   # NAT ∩ Transverse
        0b1001: (-0.25,  0.20),   # NAT ∩ Sigmoid
        0b0110: ( 0.25,  0.20),   # GTEx_All ∩ Transverse
        0b1010: ( 0.72,  0.00),   # GTEx_All ∩ Sigmoid
        0b1100: ( 0.00, -0.60),   # Transverse ∩ Sigmoid
        # triples
        0b0111: (-0.20,  0.35),   # NAT ∩ GTEx_All ∩ Transverse
        0b1011: ( 0.20,  0.35),   # NAT ∩ GTEx_All ∩ Sigmoid
        0b1101: (-0.20, -0.35),   # NAT ∩ Transverse ∩ Sigmoid
        0b1110: ( 0.20, -0.35),   # GTEx_All ∩ Transverse ∩ Sigmoid
        # all four
        0b1111: ( 0.00,  0.00),
    }
    for mask, (tx, ty) in region_text.items():
        cnt = region_data.get(mask, 0)
        if cnt > 0:
            ax.text(tx, ty, str(cnt), ha="center", va="center",
                    fontsize=10, fontweight="bold",
                    color="black" if mask != 0b1111 else "white",
                    bbox=dict(boxstyle="round,pad=0.15",
                              facecolor="white" if mask != 0b1111 else "#333333",
                              alpha=0.8 if mask != 0b1111 else 0.9,
                              edgecolor="none"))

    ax.set_title(title, fontsize=14, fontweight="bold", pad=16)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")

venn4_custom(
    sets   = [gene_sets[l] for l in labels_order],
    labels = labels_order,
    colors = colors_order,
    title  = "Biomarker Overlap Across 4 Comparison Groups\n(Tumor vs NAT / GTEx_All / Transverse / Sigmoid)",
    out_path = f"{OUT_DIR}/venn4_biomarker_overlap.png",
)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Stability score bar chart (score distribution)
# ══════════════════════════════════════════════════════════════════════════════
print("Drawing stability score distribution...")

score_colors = {4: "#C0392B", 3: "#E67E22", 2: "#F1C40F", 1: "#BDC3C7"}
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: count per score
ax = axes[0]
scores = [4, 3, 2, 1]
counts = [score_counts.get(s, 0) for s in scores]
bars = ax.bar([f"Score {s}\n({n} genes)" for s, n in zip(scores, counts)],
              counts, color=[score_colors[s] for s in scores],
              edgecolor="white", linewidth=1.5)
for bar, cnt in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.15,
            str(cnt), ha="center", va="bottom", fontsize=12, fontweight="bold")
ax.set_ylabel("Number of genes", fontsize=12)
ax.set_title("Stability Score Distribution", fontsize=13, fontweight="bold")
ax.set_ylim(0, max(counts) * 1.25)
ax.spines[["top", "right"]].set_visible(False)

# Right: horizontal gene list for score 3+4
ax2 = axes[1]
top_genes = score_df[score_df["stability_score"] >= 3].copy()
top_genes = top_genes.sort_values(["stability_score", "log2FC"],
                                   ascending=[True, False])
y_pos  = range(len(top_genes))
colors_bar = [score_colors[s] for s in top_genes["stability_score"]]
ax2.barh(list(y_pos), top_genes["log2FC"].fillna(0).values,
         color=colors_bar, edgecolor="white", linewidth=0.8, alpha=0.85)
ax2.set_yticks(list(y_pos))
ax2.set_yticklabels(top_genes["gene_id_clean"].values, fontsize=8)
ax2.axvline(0, color="black", lw=0.8)
ax2.set_xlabel("log₂FC (Tumor vs NAT)", fontsize=11)
ax2.set_title("Stable Biomarker Candidates\n(Stability Score ≥ 3)", fontsize=12, fontweight="bold")
# legend
patches = [mpatches.Patch(color=score_colors[4], label="Score 4 (all comparisons)"),
           mpatches.Patch(color=score_colors[3], label="Score 3")]
ax2.legend(handles=patches, fontsize=9, loc="lower right")
ax2.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
fig.savefig(f"{OUT_DIR}/stability_score_distribution.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {OUT_DIR}/stability_score_distribution.png")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Membership heatmap  (genes × comparisons)
# ══════════════════════════════════════════════════════════════════════════════
print("Drawing membership heatmap...")

heatmap_genes = score_df[score_df["stability_score"] >= 2].copy()
heatmap_genes = heatmap_genes.sort_values(
    ["stability_score", "log2FC"], ascending=[False, False]
)

mat = heatmap_genes[[f"in_{l}" for l in labels_order]].values.astype(float)
gnames = heatmap_genes["gene_id_clean"].values
scores_h = heatmap_genes["stability_score"].values
lfc_h    = heatmap_genes["log2FC"].fillna(0).values
dir_h    = heatmap_genes["direction"].values

n_genes = len(gnames)
fig, (ax_heat, ax_lfc, ax_score) = plt.subplots(
    1, 3, figsize=(13, max(6, n_genes * 0.35)),
    gridspec_kw={"width_ratios": [4, 1, 0.6]}
)

# heatmap
cmap_custom = matplotlib.colors.ListedColormap(["#ECF0F1", "#2C3E50"])
im = ax_heat.imshow(mat, cmap=cmap_custom, aspect="auto",
                    vmin=0, vmax=1, interpolation="nearest")
ax_heat.set_xticks(range(len(labels_order)))
ax_heat.set_xticklabels(labels_order, fontsize=10, rotation=30, ha="right")
ax_heat.set_yticks(range(n_genes))
ax_heat.set_yticklabels(gnames, fontsize=8)
ax_heat.set_title("Comparison Membership", fontsize=11, fontweight="bold")
for i in range(n_genes):
    for j in range(len(labels_order)):
        ax_heat.text(j, i, "✓" if mat[i, j] else "",
                     ha="center", va="center", fontsize=10,
                     color="white" if mat[i, j] else "#BDC3C7")

# log2FC bar
dir_colors = ["#E74C3C" if d == "Up" else "#2980B9" for d in dir_h]
ax_lfc.barh(range(n_genes), lfc_h, color=dir_colors, alpha=0.8)
ax_lfc.axvline(0, color="black", lw=0.6)
ax_lfc.set_yticks([])
ax_lfc.set_xlabel("log₂FC", fontsize=9)
ax_lfc.set_title("log₂FC\n(vs NAT)", fontsize=10, fontweight="bold")
ax_lfc.spines[["top", "right"]].set_visible(False)

# stability score dots
for i, sc in enumerate(scores_h):
    ax_score.scatter(0.5, i, s=180, color=score_colors[sc],
                     edgecolors="white", linewidth=1, zorder=3)
    ax_score.text(0.5, i, str(sc), ha="center", va="center",
                  fontsize=7, color="white", fontweight="bold")
ax_score.set_xlim(0, 1)
ax_score.set_ylim(-0.5, n_genes - 0.5)
ax_score.set_yticks([])
ax_score.set_xticks([])
ax_score.set_title("Score", fontsize=10, fontweight="bold")
ax_score.spines[:].set_visible(False)

plt.suptitle("Biomarker Stability Heatmap (Score ≥ 2)", fontsize=13,
             fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/stability_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {OUT_DIR}/stability_heatmap.png")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Pairwise Venn diagrams (6 pairs)
# ══════════════════════════════════════════════════════════════════════════════
print("Drawing pairwise Venn diagrams...")

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()
pairs = list(combinations(labels_order, 2))

for ax, (la, lb) in zip(axes, pairs):
    sa, sb = gene_sets[la], gene_sets[lb]
    only_a = len(sa - sb)
    only_b = len(sb - sa)
    both   = len(sa & sb)
    v = venn2(subsets=(only_a, only_b, both),
              set_labels=(la, lb), ax=ax)
    ca, cb = COMP_INFO[la][1], COMP_INFO[lb][1]
    if v.get_patch_by_id("10"): v.get_patch_by_id("10").set_facecolor(ca)
    if v.get_patch_by_id("01"): v.get_patch_by_id("01").set_facecolor(cb)
    if v.get_patch_by_id("11"):
        v.get_patch_by_id("11").set_facecolor("#7F8C8D")
    for patch_id in ["10", "01", "11"]:
        p = v.get_patch_by_id(patch_id)
        if p: p.set_alpha(0.6)
    ax.set_title(f"{la}  ∩  {lb}\n(shared: {both})", fontsize=10)

plt.suptitle("Pairwise Biomarker Overlap", fontsize=14,
             fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/venn_pairwise.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {OUT_DIR}/venn_pairwise.png")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Genes by stability score (separate TSVs)
# ══════════════════════════════════════════════════════════════════════════════
print("\nSaving genes by stability score...")
for sc in [4, 3, 2, 1]:
    subset = score_df[score_df["stability_score"] == sc].copy()
    subset.to_csv(f"{OUT_DIR}/stability_score_{sc}.tsv", sep="\t", index=False)
    print(f"  Score {sc}: {len(subset)} genes  → stability_score_{sc}.tsv")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Final biomarker candidates (score >= 3)
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating final biomarker candidates (score ≥ 3)...")

final = score_df[score_df["stability_score"] >= 3].copy()

# Add DEG stats from all 4 comparisons
for label, (fname, _) in COMP_INFO.items():
    deg_path = f"{DEG_DIR}/{fname}_DEG.tsv"
    deg      = pd.read_csv(deg_path, sep="\t").set_index("gene_id")
    for gid in final["gene_id"]:
        if gid in deg.index:
            final.loc[final["gene_id"] == gid, f"log2FC_{label}"] = \
                deg.loc[gid, "log2FC"]
            final.loc[final["gene_id"] == gid, f"padj_{label}"] = \
                deg.loc[gid, "padj"]

col_order = (["gene_id", "gene_id_clean", "stability_score", "direction"]
             + [f"in_{l}" for l in labels_order]
             + [f"log2FC_{l}" for l in labels_order]
             + [f"padj_{l}" for l in labels_order]
             + ["lasso_coef", "rf_importance"])
col_order = [c for c in col_order if c in final.columns]
final = final[col_order].sort_values(
    ["stability_score", "rf_importance"], ascending=[False, False]
)
final.to_csv(f"{OUT_DIR}/final_biomarker_candidates.tsv", sep="\t", index=False)

# ── Final candidates figure ───────────────────────────────────────────────────
print("Drawing final candidates plot...")

fig, axes = plt.subplots(1, 2, figsize=(14, max(5, len(final) * 0.55 + 1)))

# Left: log2FC across comparisons (dot plot)
ax_dot = axes[0]
lfc_cols = [f"log2FC_{l}" for l in labels_order if f"log2FC_{l}" in final.columns]
gene_names = final["gene_id_clean"].values
n_f = len(gene_names)

for gi, (_, row) in enumerate(final.iterrows()):
    for ci, col in enumerate(lfc_cols):
        val = row.get(col, np.nan)
        if pd.notna(val):
            color = "#E74C3C" if val > 0 else "#2980B9"
            ax_dot.scatter(ci, gi, s=abs(val) * 20 + 30,
                           color=color, alpha=0.75, edgecolors="white", lw=0.5)
            ax_dot.text(ci, gi, f"{val:.1f}", ha="center", va="center",
                        fontsize=6, color="white", fontweight="bold")

ax_dot.set_xticks(range(len(lfc_cols)))
ax_dot.set_xticklabels(labels_order, fontsize=9, rotation=20, ha="right")
ax_dot.set_yticks(range(n_f))
ax_dot.set_yticklabels(gene_names, fontsize=9)
ax_dot.set_title("log₂FC per Comparison\n(size ∝ |log₂FC|,  red=Up  blue=Down)",
                 fontsize=10, fontweight="bold")
ax_dot.spines[["top", "right"]].set_visible(False)
ax_dot.set_xlim(-0.6, len(lfc_cols) - 0.4)
ax_dot.set_ylim(-0.6, n_f - 0.4)
ax_dot.invert_yaxis()

# Right: membership + score
ax_r = axes[1]
for gi, (_, row) in enumerate(final.iterrows()):
    sc = row["stability_score"]
    ax_r.barh(gi, sc, color=score_colors[sc], alpha=0.8, edgecolor="white")
    for ci, lab in enumerate(labels_order):
        mark = "■" if row[f"in_{lab}"] else "□"
        ax_r.text(4.15 + ci * 0.45, gi, mark, ha="center", va="center",
                  fontsize=10,
                  color=COMP_INFO[lab][1] if row[f"in_{lab}"] else "#BDC3C7")

ax_r.set_yticks(range(n_f))
ax_r.set_yticklabels(gene_names, fontsize=9)
ax_r.set_xlabel("Stability Score", fontsize=10)
ax_r.set_xlim(0, 4 + len(labels_order) * 0.45 + 0.3)
ax_r.set_xticks([1, 2, 3, 4])
ax_r.set_title("Stability Score & Membership\n(■ = present  in comparison)",
               fontsize=10, fontweight="bold")
ax_r.invert_yaxis()
# column headers for membership squares
for ci, lab in enumerate(labels_order):
    ax_r.text(4.15 + ci * 0.45, -0.8, lab, ha="center", va="bottom",
              fontsize=7, color=COMP_INFO[lab][1], rotation=30)
ax_r.spines[["top", "right"]].set_visible(False)

plt.suptitle("Final Biomarker Candidates (Stability Score ≥ 3)",
             fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/final_biomarker_candidates.png", dpi=150,
            bbox_inches="tight")
plt.close(fig)

# ── print summary ─────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("STABILITY ANALYSIS SUMMARY")
print("="*60)
print(f"\nTotal unique biomarker genes: {len(all_genes)}")
print()
for sc in [4, 3, 2, 1]:
    sub = score_df[score_df["stability_score"] == sc]
    print(f"Score {sc}/4 ({len(sub)} genes):")
    for _, row in sub.iterrows():
        lfc = f"{row['log2FC']:+.2f}" if pd.notna(row["log2FC"]) else "  NA "
        print(f"  {row['gene_id_clean']:<22}  log2FC={lfc}  {row['direction']}")
    print()

print("="*60)
print("FINAL CANDIDATES (score ≥ 3)")
print("="*60)
display_cols = ["gene_id_clean", "stability_score", "direction",
                "log2FC_NAT", "log2FC_GTEx_All"]
display_cols = [c for c in display_cols if c in final.columns]
print(final[display_cols].to_string(index=False))
print(f"\nAll results saved to: {OUT_DIR}")
