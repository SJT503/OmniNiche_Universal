"""PRISMA filtering and live GEO HTTPS supplementary directory inspection for Stage II cohorts."""
from __future__ import annotations

import json
import pathlib
import re
import urllib.request

CENSUS_PATH = pathlib.Path(r"F:\OmniNiche_Universal\data_vault\raw_single_cell\geo_exhaustive_census_v1.json")
CNS_BAN = re.compile(
    r"glioma|glioblastoma|brain|alzheimer|parkinson|spinal|stroke|multiple sclerosis|"
    r"meningioma|neuroblastoma|medulloblastoma|astrocytoma|cortex|cerebell|hippocamp|neural",
    re.I,
)


def list_geo_suppl(gse: str) -> list[tuple[str, str]]:
    num = gse[3:]
    prefix = "GSE" + (num[:-3] if len(num) > 3 else "") + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "OmniNiche/2.0"})
        with urllib.request.urlopen(req, timeout=12) as resp:
            html = resp.read().decode("utf-8", errors="ignore")
        # Parse <a href="filename">...</a> date size
        rows = re.findall(r'<a href="([^"?/]+)">[^<]+</a>\s+\S+\s+\S+\s+(\S+)', html)
        return rows
    except Exception as exc:
        return [("ERR", str(exc)[:60])]


def main() -> None:
    data = json.loads(CENSUS_PATH.read_text(encoding="utf-8"))
    for cat, obj in data.items():
        print(f"\n================ {cat} (NCBI Total={obj['total_ncbi_count']}) ================")
        cands = []
        for item in obj["candidates"]:
            txt = item["title"] + " " + item["summary"]
            if CNS_BAN.search(txt):
                continue
            if item["taxon"] != "Homo sapiens":
                continue
            supp = item["suppfile"].upper()
            if not any(k in supp for k in ("H5AD", "H5", "MTX", "CSV", "TSV", "TXT", "TAR")):
                continue
            if item["n_samples"] < 6 and "H5AD" not in supp:
                continue
            cands.append(item)

        def score(x: dict) -> tuple[int, int]:
            s = x["suppfile"].upper()
            # Prefer directly loadable H5AD / H5 / MTX / TSV / CSV over giant TAR-only if possible
            has_h5ad = 4 if "H5AD" in s else (3 if ("MTX" in s or "H5" in s) else (2 if ("CSV" in s or "TSV" in s or "TXT" in s) else 1))
            n_ok = 1 if 8 <= x["n_samples"] <= 150 else 0
            return (has_h5ad, n_ok * 100 + min(x["n_samples"], 100))

        cands.sort(key=score, reverse=True)
        print(f"PRISMA Human Non-CNS Matrix-Bearing Pool: {len(cands)}")
        for x in cands[:8]:
            files = list_geo_suppl(x["gse"])
            file_str = ", ".join(f"{fn}({sz})" for fn, sz in files[:4])
            print(f"  {x['gse']} | n={x['n_samples']:3d} | {x['pdat'][:4]} | {x['title'][:82]}")
            print(f"     suppl: {file_str}")


if __name__ == "__main__":
    main()
