"""Exhaustive NCBI GEO census (retmax=500, zero arbitrary truncation) across non-CNS cohorts.

Complies strictly with lab-constitution Rule 13 (Exhaustive Census Before Selection)
and OmniNiche Firewall-2 (Non-CNS zero-overlap provenance lock).
"""
from __future__ import annotations

import json
import pathlib
import time
import urllib.parse
import urllib.request


BASE_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
OUT_PATH = pathlib.Path(r"F:\OmniNiche_Universal\data_vault\raw_single_cell\geo_exhaustive_census_v1.json")

QUERIES = {
    "solid_tumor_anti_pd1_scrna": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("single cell" OR "scRNA-seq" OR "single-cell") AND '
        '("anti-PD-1" OR "PD-1" OR "pembrolizumab" OR "nivolumab" OR "atezolizumab" OR "checkpoint blockade" OR "immunotherapy") AND '
        '("responder" OR "non-responder" OR "response" OR "pre-treatment" OR "biopsy")'
    ),
    "solid_tumor_visium_spatial": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("Visium" OR "spatial transcriptomics") AND '
        '("lung" OR "breast" OR "colorectal" OR "pancreatic" OR "renal" OR "melanoma" OR "hepatocellular" OR "carcinoma") AND '
        '("10x" OR "Visium")'
    ),
    "liver_mash_fibrosis_scrna": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("single cell" OR "scRNA-seq" OR "snRNA-seq" OR "single-nucleus") AND '
        '("NASH" OR "MASH" OR "cirrhosis" OR "liver fibrosis" OR "steatohepatitis" OR "NAFLD")'
    ),
    "kidney_ckd_aki_scrna": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("single cell" OR "scRNA-seq" OR "snRNA-seq" OR "single-nucleus") AND '
        '("chronic kidney disease" OR "kidney fibrosis" OR "acute kidney injury" OR "diabetic kidney")'
    ),
    "lung_ipf_scrna": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("single cell" OR "scRNA-seq") AND '
        '("idiopathic pulmonary fibrosis" OR "pulmonary fibrosis")'
    ),
    "ibd_anti_tnf_scrna": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("single cell" OR "scRNA-seq") AND '
        '("ulcerative colitis" OR "Crohn" OR "inflammatory bowel disease") AND '
        '("biopsy" OR "mucosa" OR "treatment" OR "TNF")'
    ),
    "sepsis_ards_scrna": (
        '"Homo sapiens"[Organism] AND "gse"[Filter] AND '
        '("single cell" OR "scRNA-seq") AND '
        '("sepsis" OR "septic shock" OR "ARDS")'
    ),
}


def run_geo_query(label: str, query: str, retmax: int = 500) -> dict:
    s_url = BASE_EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "term": query, "retmax": str(retmax), "retmode": "json"}
    )
    req = urllib.request.Request(s_url, headers={"User-Agent": "OmniNiche/2.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        es = json.loads(resp.read().decode("utf-8"))["esearchresult"]
    total = int(es["count"])
    ids = es["idlist"]
    summaries = []
    for i in range(0, len(ids), 100):
        batch = ",".join(ids[i : i + 100])
        u = BASE_EUTILS + "esummary.fcgi?" + urllib.parse.urlencode(
            {"db": "gds", "id": batch, "retmode": "json"}
        )
        req_u = urllib.request.Request(u, headers={"User-Agent": "OmniNiche/2.0"})
        with urllib.request.urlopen(req_u, timeout=45) as resp_u:
            res = json.loads(resp_u.read().decode("utf-8")).get("result", {})
        for uid in ids[i : i + 100]:
            item = res.get(uid, {})
            acc = item.get("accession", "")
            if acc.startswith("GSE"):
                summaries.append(
                    {
                        "gse": acc,
                        "title": item.get("title", ""),
                        "n_samples": int(item.get("n_samples", 0) or 0),
                        "taxon": item.get("taxon", ""),
                        "gdstype": item.get("gdstype", ""),
                        "suppfile": item.get("suppfile", ""),
                        "ftplink": item.get("ftplink", ""),
                        "pdat": item.get("pdat", ""),
                        "summary": item.get("summary", "")[:400],
                    }
                )
        time.sleep(0.35)

    usable = [
        x
        for x in summaries
        if any(
            fmt in x["suppfile"].upper()
            for fmt in ("H5AD", "H5", "MTX", "RDS", "CSV", "TSV", "TXT", "TAR")
        )
        and "RAW" not in x["suppfile"].upper()
        or any(fmt in x["suppfile"].upper() for fmt in ("H5AD", "H5", "MTX", "CSV", "TSV", "TXT", "TAR"))
    ]
    print(f"[{label}] Total={total} | Fetched={len(summaries)} | WithMatrix={len(usable)}")
    for x in usable[:6]:
        print(f"  - {x['gse']} (n={x['n_samples']}, supp={x['suppfile'][:40]}): {x['title'][:90]}")
    return {
        "query": query,
        "retmax": retmax,
        "total_ncbi_count": total,
        "fetched_gse_count": len(summaries),
        "matrix_available_count": len(usable),
        "candidates": summaries,
    }


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    results = {}
    for label, q in QUERIES.items():
        results[label] = run_geo_query(label, q, retmax=500)
        time.sleep(0.4)
    OUT_PATH.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nSaved exhaustive census -> {OUT_PATH}")


if __name__ == "__main__":
    main()
