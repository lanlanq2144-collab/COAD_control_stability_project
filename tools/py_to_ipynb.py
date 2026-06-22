"""
Convert the project's analysis scripts (scripts/*.py) into Jupyter notebooks
(notebooks/*.ipynb), preserving code verbatim but giving each script a sensible
cell structure:

  * the module docstring + a short provenance note become the intro markdown,
  * each comment "banner" section (# ===== / # STEP n / # 1. ...) becomes a
    markdown heading attached to the code it introduces (no more orphaned
    one-line comment cells, no decoration-only cells),
  * each top-level function/class gets its own code cell,
  * hard-coded "~/Desktop/COAD_control_stability_project" paths are rewritten to
    a portable REPO_ROOT so the notebooks run from any clone.

No code lines are dropped or reordered: only comment decoration is lifted into
markdown, so running every cell top-to-bottom is identical to running the .py.

The original .py files are left untouched. Run from the repo root:
    python tools/py_to_ipynb.py
"""

import re
from pathlib import Path
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "scripts"
OUT = REPO / "notebooks"
OUT.mkdir(exist_ok=True)

ORDER = [
    "extract_coad_matrix.py", "02_qc_pca.py", "03_qc_pca_tcga_only.py",
    "04_qc_pca_nat_vs_gtex.py", "05_combat_correction.py", "deg_analysis.py",
    "ml_analysis.py", "ml_corrected.py", "stability_analysis.py",
    "biological_validation.py", "external_validation.py",
]

TITLES = {
    "extract_coad_matrix.py":   "Step 1 — Build the COAD expression matrix",
    "02_qc_pca.py":             "Step 2 — QC: PCA across all groups",
    "03_qc_pca_tcga_only.py":   "Step 3 — QC: PCA within TCGA (Tumor vs NAT)",
    "04_qc_pca_nat_vs_gtex.py": "Step 4 — QC: PCA NAT vs GTEx (batch check)",
    "05_combat_correction.py":  "Step 5 — ComBat batch correction",
    "deg_analysis.py":          "Step 6 — Differential expression (4 comparisons)",
    "ml_analysis.py":           "Step 7 — ML feature selection (LASSO + RF)",
    "ml_corrected.py":          "Step 8 — Nested-CV leakage correction",
    "stability_analysis.py":    "Step 9 — Biomarker stability across controls",
    "biological_validation.py": "Step 10 — Biological validation (symbols, pathways, survival)",
    "external_validation.py":   "Step 11 — External validation (GSE156451)",
}

PREAMBLE = '''\
# --- portable repo-root resolution (added during .py -> .ipynb conversion) ---
from pathlib import Path as _Path
def _find_repo_root():
    for _p in [_Path.cwd(), *_Path.cwd().parents]:
        if (_p / "scripts").is_dir() and (_p / "README.md").is_file():
            return _p
    return _Path.cwd()
REPO_ROOT = _find_repo_root()
print("REPO_ROOT =", REPO_ROOT)
'''

PATH_SUBS = [
    (re.compile(r'BASE_DIR\s*=\s*os\.path\.expanduser\([^)]*\)'), 'BASE_DIR = str(REPO_ROOT)'),
    (re.compile(r'BASE\s*=\s*["\'][^"\']*COAD_control_stability_project["\']'), 'BASE = str(REPO_ROOT)'),
    (re.compile(r'PROJECT_DIR\s*=\s*Path\.home\(\)[^\n]*'), 'PROJECT_DIR = REPO_ROOT'),
]

# chars that count as pure decoration in a banner comment
_DECO = set("#=-─═_* ")


def rewrite_paths(text):
    for pat, repl in PATH_SUBS:
        text = pat.sub(repl, text)
    return text


def is_col0_comment(ln):
    return ln[:1] not in (" ", "\t") and ln.lstrip().startswith("#")


def is_decoration(ln):
    body = ln.lstrip()[1:]  # drop the leading '#'
    return body.strip() == "" or set(body) <= _DECO


def comment_text(ln):
    """Human text of a comment line, decoration stripped."""
    return ln.lstrip()[1:].strip(" #=-─═_*")


def split_sections(lines):
    """
    Yield (heading:str|None, code:list[str]) sections.
    A heading is the human text of the comment-banner that introduces the code.
    Top-level def/class and the post-function return to col-0 also start cells.
    """
    sections, i, n = [], 0, len(lines)
    while i < n:
        # 1. gather a leading comment/blank header run
        header = []
        while i < n and (is_col0_comment(lines[i]) or lines[i].strip() == ""):
            header.append(lines[i]); i += 1
        heading = " — ".join(
            t for t in (comment_text(h) for h in header
                        if is_col0_comment(h) and not is_decoration(h)) if t
        ) or None

        # 2. gather the code body for this section
        code, had_def = [], False
        while i < n:
            ln = lines[i]
            if is_col0_comment(ln):                       # next header begins
                break
            col0_code = ln[:1] not in (" ", "\t") and ln.strip() != ""
            if col0_code and re.match(r'(def|class)\s', ln):
                if code:                                  # new top-level def
                    break
                had_def = True
            elif col0_code and had_def:                   # returned to top level after a def
                break
            code.append(ln); i += 1

        if heading or any(s.strip() for s in code):
            sections.append((heading, code))
    return sections


def extract_docstring(code):
    """If the code block opens with a module docstring, pop it off."""
    j = 0
    while j < len(code) and code[j].strip() == "":
        j += 1
    if j < len(code) and code[j].lstrip()[:3] in ('"""', "'''"):
        q = code[j].lstrip()[:3]
        # single-line docstring
        if code[j].count(q) >= 2 and len(code[j].strip()) > 3:
            doc = [code[j].strip().strip(q)]
            return doc, code[j + 1:]
        body = [code[j].split(q, 1)[1]]
        k = j + 1
        while k < len(code) and q not in code[k]:
            body.append(code[k]); k += 1
        if k < len(code):
            body.append(code[k].split(q, 1)[0])
        return [b for b in body], code[k + 1:]
    return None, code


def convert(fname):
    src = rewrite_paths((SCRIPTS / fname).read_text(encoding="utf-8"))
    sections = split_sections(src.splitlines())

    nb = new_notebook()
    intro = (f"# {TITLES.get(fname, fname)}\n\n"
             f"Auto-converted from `scripts/{fname}` — code preserved verbatim; "
             f"section banners lifted into headings; paths made portable "
             f"(`REPO_ROOT`).")
    # fold a leading module docstring into the intro
    if sections:
        doc, rest = extract_docstring(sections[0][1])
        if doc is not None:
            text = "\n".join(doc).strip()
            if text:
                intro += "\n\n---\n\n" + text
            sections[0] = (sections[0][0], rest)
    nb.cells.append(new_markdown_cell(intro))
    nb.cells.append(new_code_cell(PREAMBLE))

    for heading, code in sections:
        if heading:
            nb.cells.append(new_markdown_cell(f"### {heading}"))
        body = "\n".join(code).strip("\n")
        if body.strip():
            nb.cells.append(new_code_cell(body))

    stem = re.sub(r"^\d+_", "", Path(fname).stem)
    out_name = f"{ORDER.index(fname)+1:02d}_{stem}.ipynb"
    nbformat.write(nb, OUT / out_name)
    n_md = sum(1 for c in nb.cells if c.cell_type == "markdown")
    n_code = sum(1 for c in nb.cells if c.cell_type == "code")
    return out_name, n_md, n_code


if __name__ == "__main__":
    print(f"Repo: {REPO}\nWriting notebooks to: {OUT}\n")
    for f in ORDER:
        if (SCRIPTS / f).exists():
            name, n_md, n_code = convert(f)
            print(f"  {f:<28} -> notebooks/{name}  ({n_md} md + {n_code} code)")
        else:
            print(f"  {f:<28} -- MISSING, skipped")
    print("\nDone.")
