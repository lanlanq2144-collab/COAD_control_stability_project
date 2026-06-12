"""
Build the four sample-ID lists his pipeline expects but does NOT ship:
    data_processed/sample_ids/TCGA_COAD_Tumor_sample_ids.txt      (288)
    data_processed/sample_ids/TCGA_COAD_NAT_sample_ids.txt        ( 41)
    data_processed/sample_ids/GTEx_Colon_Transverse_sample_ids.txt (167)
    data_processed/sample_ids/GTEx_Colon_Sigmoid_sample_ids.txt    (141)

Derived from TcgaTargetGTEX_phenotype.txt.gz (UCSC Xena). Filters were verified
against his README group sizes. Run AFTER tools/fetch_data.py.

Usage (from repo root):
    python tools/build_sample_ids.py
    python tools/build_sample_ids.py --data-dir D:/coad_data
"""
import argparse, gzip, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# (output filename, filter predicate on the phenotype row dict)
GROUPS = {
    "TCGA_COAD_Tumor": lambda r: r["detailed_category"] == "Colon Adenocarcinoma"
                                 and r["_sample_type"] == "Primary Tumor",
    "TCGA_COAD_NAT":   lambda r: r["detailed_category"] == "Colon Adenocarcinoma"
                                 and r["_sample_type"] == "Solid Tissue Normal",
    "GTEx_Colon_Transverse": lambda r: r["_study"] == "GTEX"
                                 and r["primary disease or tissue"] == "Colon - Transverse",
    "GTEx_Colon_Sigmoid":    lambda r: r["_study"] == "GTEX"
                                 and r["primary disease or tissue"] == "Colon - Sigmoid",
}
EXPECTED = {"TCGA_COAD_Tumor": 288, "TCGA_COAD_NAT": 41,
            "GTEx_Colon_Transverse": 167, "GTEx_Colon_Sigmoid": 141}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(REPO / "data_raw"))
    args = ap.parse_args()

    pheno = Path(args.data_dir) / "TcgaTargetGTEX_phenotype.txt.gz"
    if not pheno.exists():
        sys.exit(f"phenotype file not found: {pheno}\nRun tools/fetch_data.py first.")

    out_dir = REPO / "data_processed" / "sample_ids"
    out_dir.mkdir(parents=True, exist_ok=True)

    # read phenotype (tab-separated, latin-1 is safe for this file)
    with gzip.open(pheno, "rt", encoding="latin-1") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]

    print(f"phenotype rows: {len(rows)}\n")
    ok = True
    for group, pred in GROUPS.items():
        ids = [r["sample"] for r in rows if pred(r)]
        (out_dir / f"{group}_sample_ids.txt").write_text("\n".join(ids) + "\n")
        flag = "OK" if len(ids) == EXPECTED[group] else f"!! expected {EXPECTED[group]}"
        ok &= len(ids) == EXPECTED[group]
        print(f"  {group:<24} {len(ids):>4}  {flag}")

    print(f"\nWrote 4 lists to {out_dir}")
    print("Match his README group sizes:" , "YES" if ok else "NO (Xena version drift?)")
    print("\nNext: run notebooks/01_extract_coad_matrix.ipynb (or scripts/extract_coad_matrix.py)")


if __name__ == "__main__":
    sys.exit(main())
