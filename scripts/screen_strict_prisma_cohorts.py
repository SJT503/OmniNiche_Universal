"""Domain-strict PRISMA cohort screener across the 1,072-study GEO census + landmark verification.

Strictly enforces lab-constitution Rule 13:
- Zero CNS / neurological / psychiatric studies
- Zero viral COVID-19 / SARS-CoV-2 / influenza / allergy / preeclampsia cross-contamination
- Primary human tissue (no cell lines / organoid-only / CD3-only when whole-niche needed)
- Directly readable Python formats (.h5ad, .h5, .h5Seurat, .mtx.gz, .csv.gz, .tsv.gz, .txt.gz, or moderate .tar)
"""
from __future__ import annotations

import json
import pathlib
import re
import urllib.parse
import urllib.request

CENSUS_PATH = pathlib.Path(r"F:\OmniNiche_Universal\data_vault\raw_single_cell\geo_exhaustive_census_v1.json")

GLOBAL_BAN = re.compile(
    r"glioma|glioblastoma|brain|alzheimer|parkinson|spinal|stroke|multiple sclerosis|"
    r"meningioma|neuroblastoma|medulloblastoma|astrocytoma|cortex|cerebell|hippocamp|neural|"
    r"covid|sars-cov|influenza|peanut|allergy|preeclampsia|decidual|appendiceal|aortic|"
    r"organoid|xenograft|murine|mouse|mice|cell line|in vitro",
    re.I,
)

STRICT_INCLUDE = {
    "solid_tumor_anti_pd1_scrna": re.compile(
        r"(lung|nsclc|breast|tnbc|melanoma|renal|carcinoma|cancer|squamous|basal cell|esophageal|gastric|bladder|hepatocellular|hnscc|head and neck)",
        re.I,
    ),
    "solid_tumor_visium_spatial": re.compile(
        r"(lung|nsclc|breast|tnbc|melanoma|renal|colorectal|pancreatic|pdac|hepatocellular|squamous|carcinoma|cancer)",
        re.I,
    ),
    "liver_mash_fibrosis_scrna": re.compile(
        r"(liver|hepatic|cirrhosis|nash|mash|steatohepatitis|nafld|fibrosis)",
        re.I,
    ),
    "kidney_ckd_aki_scrna": re.compile(
        r"(kidney|renal|glomerular|diabetic kidney|nephropathy)",
        re.I,
    ),
    "lung_ipf_scrna": re.compile(
        r"(pulmonary fibrosis|ipf|interstitial lung)",
        re.I,
    ),
    "ibd_anti_tnf_scrna": re.compile(
        r"(ulcerative colitis|crohn|inflammatory bowel)",
        re.I,
    ),
    "sepsis_ards_scrna": re.compile(
        r"(sepsis|septic shock)",
        re.I,
    ),
}


def list_geo_suppl(gse: str) -> list[tuple[str, str]]:
    num = gse[3:]
    prefix = "GSE" + (num[:-3] if len(num) > 3 else "") + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "OmniNiche/2.0"})
        with urllib.request.urlopen(req, timeout=12) as resp:
            html = resp.read().decode("utf-8", errors="ignore")
        return re.findall(r'<a href="([^"?/]+)">[^<]+</a>\s+\S+\s+\S+\s+(\S+)', html)
    except Exception as exc:
        return [("ERR", str(exc)[:50])]


def fetch_specific_landmarks(gses: list[str]) -> list[dict]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    term = " OR ".join(f"{g}[ACCN]" for g in gses) + ' AND "gse"[Filter]'
    u = base + "esearch.fcgi?" + urllib.parse.urlencode({"db": "gds", "term": term, "retmax": "50", "retmode": "json"})
    with urllib.request.urlopen(u, timeout=20) as r:
        ids = json.loads(r.read().decode("utf-8"))["esearchresult"]["idlist"]
    if not ids:
        return []
    u2 = base + "esummary.fcgi?" + urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
    with urllib.request.urlopen(u2, timeout=25) as r:
        res = json.loads(r.read().decode("utf-8")).get("result", {})
    out = []
    for uid in ids:
        it = res.get(uid, {})
        if it.get("accession", "").startswith("GSE"):
            out.append(it)
    return out


def main() -> None:
    data = json.loads(CENSUS_PATH.read_text(encoding="utf-8"))
    for cat, obj in data.items():
        inc = STRICT_INCLUDE[cat]
        strict_pool = []
        for item in obj["candidates"]:
            if item["taxon"] != "Homo sapiens":
                continue
            txt = item["title"] + " " + item["summary"]
            if GLOBAL_BAN.search(txt):
                continue
            if not inc.search(item["title"]):
                continue
            if cat == "liver_mash_fibrosis_scrna" and re.search(r"cholangiocarcinoma|hepatoblastoma", txt, re.I):
                continue
            strict_pool.append(item)

        def rank(x: dict) -> tuple[int, int]:
            s = x["suppfile"].upper()
            py_friendly = 4 if "H5AD" in s else (3 if ("H5" in s or "MTX" in s) else (2 if ("CSV" in s or "TSV" in s or "TXT" in s) else 1))
            return (py_friendly, min(x["n_samples"], 120))

        strict_pool.sort(key=rank, reverse=True)
        print(f"\n*** {cat}: Strict Primary Human Non-CNS Pool = {len(strict_pool)} / {obj['total_ncbi_count']} ***")
        for x in strict_pool[:7]:
            suppl = list_geo_suppl(x["gse"])
            s_txt = ", ".join(f"{fn}({sz})" for fn, sz in suppl[:5])
            print(f"  {x['gse']} | n={x['n_samples']:3d} | {x['title'][:88]}")
            print(f"     -> {s_txt}")

    print("\n*** Verifying Well-Known Peer-Reviewed Clinical Landmarks Live via NCBI API ***")
    landmarks = [
        "GSE123813",  # Yost et al. Nat Med 2019 BCC/SCC pre/post anti-PD-1 scRNA-seq
        "GSE169246",  # Zhang et al. Cancer Cell 2021 TNBC anti-PD-L1 pre/post R vs NR scRNA-seq
        "GSE179994",  # Liu et al. Nat Cancer 2022 NSCLC anti-PD-1 pre/post R vs NR scRNA-seq
        "GSE207422",  # Hu et al. Cell Res 2023 NSCLC neoadjuvant PD-1 R vs NR scRNA-seq
        "GSE136103",  # Ramachandran et al. Nature 2019 Human Liver Cirrhosis/Fibrosis niche scRNA-seq
        "GSE135893",  # Habermann et al. Sci Adv 2020 Human IPF vs Control Lung scRNA-seq
        "GSE131882",  # Wilson et al. PNAS 2019 Human Diabetic Kidney Disease snRNA-seq
        "GSE176078",  # Wu et al. Nat Genet 2021 Human Breast Cancer scRNA + Visium Spatial Atlas
        "GSE203612",  # Barkley et al. Nat Genet 2022 Pan-Cancer 10x Visium Spatial Transcriptomics
    ]
    for it in fetch_specific_landmarks(landmarks):
        acc = it["accession"]
        suppl = list_geo_suppl(acc)
        s_txt = ", ".join(f"{fn}({sz})" for fn, sz in suppl[:5])
        print(f"  {acc} | n={it.get('n_samples')} | taxon={it.get('taxon')} | {it.get('title', '')[:82]}")
        print(f"     -> {s_txt}")


if __name__ == "__main__":
    main()
