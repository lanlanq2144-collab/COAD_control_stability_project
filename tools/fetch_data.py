"""
Download the raw data this project needs, into <repo>/data_raw/.

  * TcgaTargetGtex_rsem_gene_tpm.gz    (UCSC Xena Toil hub) ~1.23 GB
  * TcgaTargetGTEX_phenotype.txt.gz    (UCSC Xena Toil hub) ~130 KB
  * GSE156451_RAW.tar                  (NCBI GEO)           ~42 MB  -> extracted

Usage (from repo root):
    python tools/fetch_data.py                 # -> <repo>/data_raw
    python tools/fetch_data.py --data-dir D:/coad_data   # custom location

NOTE (Box): this repo lives inside a Box-synced folder. The default target
(<repo>/data_raw) is git-ignored but Box will still upload ~1.3 GB. To avoid
that, pass --data-dir pointing OUTSIDE Box (e.g. C:/data/coad), or mark
data_raw as online-only / exclude it from Box sync.
"""
import argparse, sys, tarfile
from pathlib import Path
import requests

REPO = Path(__file__).resolve().parent.parent

FILES = {
    "TcgaTargetGtex_rsem_gene_tpm.gz":
        "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "TcgaTargetGTEX_phenotype.txt.gz":
        "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz",
    "GSE156451_RAW.tar":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156451/suppl/GSE156451_RAW.tar",
}


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  skip (exists): {dest.name}")
        return
    print(f"  downloading {dest.name} ...")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        got = 0
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                got += len(chunk)
                if total:
                    pct = got * 100 // total
                    print(f"\r    {got>>20} / {total>>20} MB ({pct}%)", end="", flush=True)
        print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=str(REPO / "data_raw"),
                    help="target folder for data_raw (default: <repo>/data_raw)")
    args = ap.parse_args()

    raw = Path(args.data_dir)
    raw.mkdir(parents=True, exist_ok=True)
    print(f"data_raw -> {raw}\n")

    for name, url in FILES.items():
        download(url, raw / name)

    # extract the GEO tar into data_raw/GSE156451/
    tar_path = raw / "GSE156451_RAW.tar"
    geo_dir = raw / "GSE156451"
    if tar_path.exists():
        geo_dir.mkdir(exist_ok=True)
        print(f"\n  extracting GSE156451_RAW.tar -> {geo_dir} ...")
        with tarfile.open(tar_path) as t:
            t.extractall(geo_dir)
        n = len(list(geo_dir.glob("GSM*")))
        print(f"  extracted {n} GSM files")

    print("\nDone. If you used a custom --data-dir, point the pipeline at it "
          "(see tools/build_sample_ids.py and the notebook paths).")


if __name__ == "__main__":
    sys.exit(main())
