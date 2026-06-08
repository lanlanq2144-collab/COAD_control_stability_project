"""
Biological validation of final biomarker candidates
1. ENSG ID → gene symbol (mygene.info)
2. Pathway enrichment (g:Profiler REST API)
3. Survival analysis (GDC clinical data + KM curves)
"""

import os, json, time, warnings
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
warnings.filterwarnings("ignore")

# ── paths ─────────────────────────────────────────────────────────────────────
BASE       = "/Users/georgiatechgogofinn/Desktop/COAD_control_stability_project"
FINAL_BM   = f"{BASE}/results/stability/final_biomarker_candidates.tsv"
EXPR_FILE  = f"{BASE}/data_processed/expression/coad_expression_tpm_genes_by_samples.tsv"
LABEL_FILE = f"{BASE}/data_processed/metadata/coad_sample_labels.tsv"
OUT_DIR    = f"{BASE}/results/biological_validation"
os.makedirs(OUT_DIR, exist_ok=True)

# ── load biomarker list ───────────────────────────────────────────────────────
bm     = pd.read_csv(FINAL_BM, sep="\t")
labels = pd.read_csv(LABEL_FILE, sep="\t", index_col=0)
label_map = labels["group"].to_dict()

print("Final biomarker candidates:")
for _, r in bm.iterrows():
    print(f"  {r['gene_id_clean']}  log2FC(NAT)={r['log2FC_NAT']:+.2f}  {r['direction']}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Gene symbol mapping  (mygene.info)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("STEP 1: Gene name mapping (mygene.info)")
print("="*60)

def query_mygene(ensg_id):
    try:
        r = requests.get(
            "https://mygene.info/v3/query",
            params={"q": f"ensembl.gene:{ensg_id}",
                    "fields": "symbol,name,summary",
                    "species": "human", "size": 1},
            timeout=10,
        )
        hits = r.json().get("hits", [])
        return hits[0] if hits else {}
    except Exception as e:
        print(f"  mygene error ({ensg_id}): {e}")
        return {}

gene_info = {}
for _, row in bm.iterrows():
    ec = row["gene_id_clean"]
    print(f"  {ec} ...", end=" ", flush=True)
    hit = query_mygene(ec)
    gene_info[ec] = {
        "gene_id"      : row["gene_id"],
        "ensg_clean"   : ec,
        "symbol"       : hit.get("symbol", "unknown"),
        "name"         : hit.get("name", ""),
        "summary"      : (hit.get("summary") or "")[:300],
        "log2FC_NAT"   : row.get("log2FC_NAT", np.nan),
        "log2FC_GTEx"  : row.get("log2FC_GTEx_All", np.nan),
        "direction"    : row.get("direction", ""),
        "stability"    : row.get("stability_score", ""),
    }
    print(hit.get("symbol", "unknown"))
    time.sleep(0.25)

gene_df = pd.DataFrame(list(gene_info.values()))
gene_df.to_csv(f"{OUT_DIR}/gene_name_mapping.tsv", sep="\t", index=False)

print("\nMapping result:")
for _, r in gene_df.iterrows():
    print(f"  {r['ensg_clean']} → {r['symbol']:<12} | {r['name'][:55]}")

# ── gene info summary table figure ────────────────────────────────────────────
print("\nDrawing gene info table...")
fig, ax = plt.subplots(figsize=(16, max(4, len(gene_df) * 1.1 + 1.5)))
ax.axis("off")

tbl_data = []
for _, r in gene_df.iterrows():
    lfc_nat  = f"{r['log2FC_NAT']:+.2f}"  if pd.notna(r["log2FC_NAT"])  else "NA"
    lfc_gtex = f"{r['log2FC_GTEx']:+.2f}" if pd.notna(r["log2FC_GTEx"]) else "NA"
    summ = (r["summary"][:90] + "…") if len(r["summary"]) > 90 else r["summary"]
    tbl_data.append([r["symbol"], r["name"][:38],
                     lfc_nat, lfc_gtex, r["direction"], summ])

tbl = ax.table(
    cellText=tbl_data,
    colLabels=["Symbol", "Gene Name", "log2FC\n(vs NAT)",
               "log2FC\n(vs GTEx)", "Direction", "Summary"],
    cellLoc="left", loc="center",
    colWidths=[0.07, 0.16, 0.07, 0.07, 0.07, 0.45],
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(8)
tbl.scale(1, 2.5)
for j in range(6):
    tbl[0, j].set_facecolor("#2C3E50")
    tbl[0, j].set_text_props(color="white", fontweight="bold")
for i, r in enumerate(gene_df.itertuples()):
    bg = "#F8F9FA" if i % 2 == 0 else "#FFFFFF"
    for j in range(6):
        tbl[i+1, j].set_facecolor(bg)
    tbl[i+1, 4].set_facecolor("#FDECEA" if r.direction == "Up" else "#EBF5FB")
ax.set_title("Final Biomarker Candidates — Gene Information",
             fontsize=13, fontweight="bold", pad=18)
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/gene_info_table.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {OUT_DIR}/gene_info_table.png")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Pathway enrichment  (g:Profiler)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("STEP 2: Pathway enrichment (g:Profiler)")
print("="*60)

symbols = [v["symbol"] for v in gene_info.values() if v["symbol"] != "unknown"]
print(f"  Querying with symbols: {symbols}")

def run_gprofiler(syms):
    try:
        r = requests.post(
            "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
            json={
                "organism": "hsapiens",
                "query"   : syms,
                "sources" : ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"],
                "significance_threshold_method": "fdr",
                "user_threshold": 0.05,
                "no_evidences": False,
            },
            timeout=30,
        )
        return r.json()
    except Exception as e:
        print(f"  g:Profiler error: {e}")
        return {}

gp = run_gprofiler(symbols)
pw_rows = []
if gp.get("result"):
    for item in gp["result"]:
        pw_rows.append({
            "source"    : item.get("source", ""),
            "term_id"   : item.get("native", ""),
            "term_name" : item.get("name", ""),
            "p_value"   : item.get("p_value", 1.0),
            "term_size" : item.get("term_size", 0),
            "intersection_size": item.get("intersection_size", 0),
            "genes"     : ",".join([
                                    (g if isinstance(g, str) else g[0] if g else "")
                                    for g in (item.get("intersections") or [])
                                    if g
                                ]),
        })
    print(f"  Found {len(pw_rows)} significant terms")
else:
    print("  No significant pathways (API unavailable or gene set too small)")

pw_df = pd.DataFrame(pw_rows) if pw_rows else pd.DataFrame(
    columns=["source","term_id","term_name","p_value","term_size",
             "intersection_size","genes"])
pw_df.to_csv(f"{OUT_DIR}/pathway_enrichment.tsv", sep="\t", index=False)

# pathway dot plot
source_colors = {
    "GO:BP": "#3498DB", "GO:MF": "#E74C3C", "GO:CC": "#2ECC71",
    "KEGG" : "#F39C12", "REAC" : "#9B59B6", "WP"   : "#1ABC9C",
}
top_pw = pw_df[pw_df["p_value"] < 0.05].sort_values("p_value").head(25)

fig, ax = plt.subplots(figsize=(11, max(5, len(top_pw) * 0.42 + 2)))
if len(top_pw) > 0:
    neg_lp = -np.log10(top_pw["p_value"].clip(1e-50))
    cols   = [source_colors.get(s, "#7F8C8D") for s in top_pw["source"]]
    sizes  = [max(30, is_ * 60) for is_ in top_pw["intersection_size"]]
    ypos   = range(len(top_pw))
    ax.scatter(neg_lp, list(ypos), s=sizes, c=cols, alpha=0.85, zorder=3)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(
        [f"{row['term_name'][:58]}  [{row['source']}]"
         for _, row in top_pw.iterrows()], fontsize=8)
    ax.axvline(-np.log10(0.05), ls="--", color="gray", lw=1, alpha=0.5)
    ax.set_xlabel("-log₁₀(p-value, FDR)", fontsize=11)
    patches = [mpatches.Patch(color=c, label=s)
               for s, c in source_colors.items()
               if s in top_pw["source"].values]
    ax.legend(handles=patches, title="Database", fontsize=8,
              title_fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
else:
    ax.text(0.5, 0.5,
            "No significant pathways found\n"
            "(n=6 genes; enrichment power is low with small gene sets)\n\n"
            "Individual gene functions are shown in gene_info_table.png",
            ha="center", va="center", fontsize=11,
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.6",
                      facecolor="#FFF9C4", alpha=0.9))
    ax.axis("off")

ax.set_title("Pathway Enrichment — Biomarker Candidates (g:Profiler, FDR<0.05)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
fig.savefig(f"{OUT_DIR}/pathway_dotplot.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {OUT_DIR}/pathway_dotplot.png")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Survival analysis
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("STEP 3: Survival analysis (TCGA-COAD, GDC API)")
print("="*60)

# ── 3a. Download clinical data ────────────────────────────────────────────────
def fetch_gdc_clinical(project="TCGA-COAD", size=600):
    try:
        r = requests.get(
            "https://api.gdc.cancer.gov/cases",
            params={
                "filters": json.dumps({
                    "op": "in",
                    "content": {"field": "project.project_id", "value": [project]},
                }),
                "fields": ("submitter_id,"
                            "demographic.vital_status,"
                            "demographic.days_to_death,"
                            "diagnoses.days_to_last_follow_up"),
                "expand": "demographic,diagnoses",
                "format": "JSON",
                "size"  : str(size),
            },
            timeout=30,
        )
        return r.json().get("data", {}).get("hits", [])
    except Exception as e:
        print(f"  GDC API error: {e}")
        return []

print("  Fetching TCGA-COAD clinical data...")
hits = fetch_gdc_clinical()
print(f"  Retrieved {len(hits)} cases")

clinical_rows = []
for case in hits:
    pid  = case.get("submitter_id", "")
    demo = case.get("demographic") or {}
    diag = (case.get("diagnoses") or [{}])[0]

    vital  = (demo.get("vital_status") or "").lower()
    days_d = demo.get("days_to_death")
    days_f = diag.get("days_to_last_follow_up")

    if vital == "dead" and days_d is not None and float(days_d) >= 0:
        os_time, os_event = float(days_d), 1
    elif days_f is not None and float(days_f) >= 0:
        os_time, os_event = float(days_f), 0
    else:
        continue

    clinical_rows.append({"patient_id": pid,
                           "os_days"   : os_time,
                           "os_event"  : os_event})

clinical_df = pd.DataFrame(clinical_rows)
print(f"  Usable OS records: {len(clinical_df)}")

if len(clinical_df) < 20:
    print("  Insufficient clinical data — skipping survival analysis")
else:
    clinical_df.to_csv(f"{OUT_DIR}/tcga_coad_clinical.tsv", sep="\t", index=False)

    # ── 3b. Load expression matrix (tumor only) ───────────────────────────────
    print("  Loading expression matrix (tumor samples)...")
    expr = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)
    tumor_cols  = [c for c in expr.columns if label_map.get(c) == "TCGA_COAD_Tumor"]
    expr_tumor  = expr[tumor_cols]
    # patient barcode = first 12 chars of TCGA sample ID
    sample_to_patient = {s: s[:12] for s in tumor_cols}

    # ── 3c. Match samples to clinical ────────────────────────────────────────
    clin_idx = clinical_df.set_index("patient_id")
    match_rows = []
    for samp, pid in sample_to_patient.items():
        if pid in clin_idx.index:
            match_rows.append({
                "sample"    : samp,
                "patient_id": pid,
                "os_days"   : clin_idx.loc[pid, "os_days"],
                "os_event"  : clin_idx.loc[pid, "os_event"],
            })
    match_df = pd.DataFrame(match_rows)
    print(f"  Matched {len(match_df)} tumor samples with OS data")

    # ── 3d. KM per gene ───────────────────────────────────────────────────────
    surv_summary = []

    for _, bm_row in bm.iterrows():
        gene_id = bm_row["gene_id"]
        ec      = bm_row["gene_id_clean"]
        symbol  = gene_info.get(ec, {}).get("symbol", ec)
        direction = bm_row.get("direction", "")

        if gene_id not in expr_tumor.index:
            print(f"  {symbol}: gene not in expression matrix — skip")
            continue

        surv = match_df.copy()
        surv["expr"] = surv["sample"].map(expr_tumor.loc[gene_id].to_dict())
        surv = surv.dropna(subset=["expr"])

        if len(surv) < 20:
            print(f"  {symbol}: too few matched samples — skip")
            continue

        median_val = surv["expr"].median()
        surv["group"] = np.where(surv["expr"] > median_val, "High", "Low")

        high = surv[surv["group"] == "High"]
        low  = surv[surv["group"] == "Low"]

        lr   = logrank_test(high["os_days"], low["os_days"],
                            event_observed_A=high["os_event"],
                            event_observed_B=low["os_event"])
        pval = lr.p_value
        sig  = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

        # KM plot
        fig, ax = plt.subplots(figsize=(8, 6))
        kmf_h = KaplanMeierFitter()
        kmf_l = KaplanMeierFitter()
        kmf_h.fit(high["os_days"], event_observed=high["os_event"],
                  label=f"High expression (n={len(high)})")
        kmf_l.fit(low["os_days"],  event_observed=low["os_event"],
                  label=f"Low expression  (n={len(low)})")
        kmf_h.plot_survival_function(ax=ax, color="#E74C3C",
                                     ci_show=True, ci_alpha=0.12, lw=2)
        kmf_l.plot_survival_function(ax=ax, color="#2980B9",
                                     ci_show=True, ci_alpha=0.12, lw=2)
        ax.set_xlabel("Overall Survival (days)", fontsize=12)
        ax.set_ylabel("Survival Probability", fontsize=12)
        p_str = f"p = {pval:.4f}" if pval >= 0.0001 else "p < 0.0001"
        ax.set_title(
            f"Kaplan-Meier Survival Curve\n"
            f"{symbol}  ({ec})   {direction}-regulated\n"
            f"Log-rank {p_str}  {sig}",
            fontsize=11, fontweight="bold",
        )
        ax.set_ylim(-0.05, 1.05)
        ax.text(0.97, 0.97,
                f"Split at median\n(log₂TPM = {median_val:.2f})",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9, color="gray")
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        out_km = f"{OUT_DIR}/survival_km_{symbol}.png"
        fig.savefig(out_km, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  {symbol:<12} {p_str}  {sig}")

        surv_summary.append({
            "gene_id"    : gene_id,
            "symbol"     : symbol,
            "direction"  : direction,
            "n_high"     : len(high),
            "n_low"      : len(low),
            "logrank_p"  : round(pval, 6),
            "significant": pval < 0.05,
        })

    if surv_summary:
        surv_df = pd.DataFrame(surv_summary)
        surv_df.to_csv(f"{OUT_DIR}/survival_summary.tsv", sep="\t", index=False)

        # ── survival summary bar ──────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(7, max(4, len(surv_df) * 0.7 + 1.5)))
        neg_lp  = -np.log10(np.clip(surv_df["logrank_p"].values, 1e-10, 1.0))
        bcolors = ["#E74C3C" if p < 0.05 else "#BDC3C7"
                   for p in surv_df["logrank_p"]]
        bars = ax.barh(surv_df["symbol"].values, neg_lp,
                       color=bcolors, alpha=0.85, edgecolor="white")
        ax.axvline(-np.log10(0.05), ls="--", color="#E74C3C",
                   lw=1.5, alpha=0.7, label="p = 0.05")
        for bar, p in zip(bars, surv_df["logrank_p"]):
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                    sig, va="center", fontsize=11,
                    color="#C0392B" if p < 0.05 else "gray")
        ax.set_xlabel("-log₁₀(log-rank p-value)", fontsize=11)
        ax.set_title("Survival Association per Gene\n(High vs Low expression, TCGA-COAD tumor)",
                     fontsize=12, fontweight="bold")
        ax.legend(fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        fig.savefig(f"{OUT_DIR}/survival_summary.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved → {OUT_DIR}/survival_summary.png")

# ── print final summary ───────────────────────────────────────────────────────
print("\n" + "="*60)
print("BIOLOGICAL VALIDATION SUMMARY")
print("="*60)

print("\n[ Gene Names ]")
print(gene_df[["ensg_clean", "symbol", "name", "direction",
               "log2FC_NAT"]].to_string(index=False))

if pw_rows:
    top5 = pw_df.head(5)
    print(f"\n[ Top Pathways (g:Profiler, n={len(pw_df)}) ]")
    for _, r in top5.iterrows():
        print(f"  {r['source']:<8} {r['term_name'][:55]}  p={r['p_value']:.2e}")

if "surv_summary" in dir() and surv_summary:
    print("\n[ Survival Analysis ]")
    print(surv_df[["symbol", "direction", "logrank_p",
                   "significant"]].to_string(index=False))

print(f"\nAll results saved to: {OUT_DIR}")
