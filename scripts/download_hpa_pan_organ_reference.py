"""Download official Human Protein Atlas normal pan-organ single-cell transcriptomics reference.

Used as the healthy pan-organ baseline tau_0 reference (81 human single-cell types across
normal human organs) in Fig. 3, with full SHA-256 provenance logging.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import urllib.request
import zipfile
import pandas as pd

DEST_DIR = pathlib.Path(r"F:\OmniNiche_Universal\data_vault\raw_single_cell\hpa_pan_organ_normal")
URL = "https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip"


def main() -> None:
    DEST_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = DEST_DIR / "rna_single_cell_type.tsv.zip"
    if not zip_path.exists():
        req = urllib.request.Request(URL, headers={"User-Agent": "OmniNiche/2.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            zip_path.write_bytes(resp.read())

    sha256 = hashlib.sha256(zip_path.read_bytes()).hexdigest()
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(DEST_DIR)
        extracted = zf.namelist()

    tsv_path = DEST_DIR / extracted[0]
    df = pd.read_csv(tsv_path, sep="\t")
    n_genes = int(df["Gene name"].nunique())
    n_cell_types = int(df["Cell type"].nunique())
    manifest = {
        "source_url": URL,
        "zip_file": str(zip_path),
        "zip_bytes": zip_path.stat().st_size,
        "zip_sha256": sha256,
        "extracted_tsv": str(tsv_path),
        "tsv_sha256": hashlib.sha256(tsv_path.read_bytes()).hexdigest(),
        "n_rows": int(len(df)),
        "n_unique_genes": n_genes,
        "n_unique_cell_types": n_cell_types,
        "columns": list(df.columns),
        "sample_cell_types": sorted(df["Cell type"].unique().tolist())[:20],
    }
    (DEST_DIR / "manifest_hpa_single_cell.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"HPA Pan-Organ Normal Single-Cell Reference Ready: {n_genes} genes x {n_cell_types} cell types | SHA256={sha256[:16]}...")


if __name__ == "__main__":
    main()
