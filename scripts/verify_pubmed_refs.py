"""Live NCBI PubMed E-utilities verification for all 30 OmniNiche Stage II Flagship references."""

from __future__ import annotations

import json
import pathlib
import time
import urllib.parse
import urllib.request


# 19 already-confirmed exact primary PMIDs + 11 targeted FirstAuthor[1AU] + Year[DP] + Title keyword queries
CONFIRMED_PMIDS = {
    "sun_2022_attrition": "35865092",
    "nelson_2015_genetics": "26121088",
    "buniello_2025_opentargets": "39657122",
    "kryuchkova_2017_tau": "26891983",
    "swanson_2025_virtual_lab": "40730228",
    "ramachandran_2019_gse136103": "31597160",
    "wang_2023_gse212837": "36599008",
    "fabre_2023_spp1_mac": "37027478",
    "anstee_2024_aurora_cvc": "37061109",
    "habermann_2020_gse135893": "32832598",
    "raghu_2024_pamrevlumab": "38762797",
    "mariathasan_2018_tgfb_pd1": "29443960",
    "kieffer_2020_mycaf": "32434947",
    "gouin_2021_gse171351": "34385456",
    "barkley_2022_gse203612": "35931863",
    "yost_2019_gse123813": "31359002",
    "hu_2023_gse207422": "36869384",
    "wang_2008_landscape": "18719111",
    "raue_2009_profile_lik": "19505944",
}

TARGETED_QUERIES = [
    ("wong_2019_success", "Wong CH[1AU] AND Lo AW[AU] AND Biostatistics[TA] AND 2019[DP] NOT Corrigendum[TI]"),
    ("minikel_2024_nature", "Minikel EV[1AU] AND Nature[TA] AND 2024[DP] AND clinical[TI]"),
    ("karlsson_2021_hpa_sc", "Karlsson M[1AU] AND Uhlen M[AU] AND 2021[DP] AND single-cell[TI]"),
    ("bunne_2024_virtual_cell", "Bunne C[1AU] AND Cell[TA] AND 2024[DP] AND virtual[TI]"),
    ("williams_2015_ctgov_stop", "Williams RJ[1AU] AND Zarin DA[AU] AND 2015[DP] AND Terminated[TI]"),
    ("harrison_2024_resmetirom", "Harrison SA[1AU] AND N Engl J Med[TA] AND 2024[DP] AND Resmetirom[TI] NOT Reply[TI]"),
    ("harrison_2020_selonsertib", "Harrison SA[1AU] AND 2020[DP] AND Selonsertib[TI]"),
    ("richeldi_2014_nintedanib", "Richeldi L[1AU] AND N Engl J Med[TA] AND 2014[DP] AND nintedanib[TI]"),
    ("bakris_2020_finerenone", "Bakris GL[1AU] AND N Engl J Med[TA] AND 2020[DP] AND Finerenone[TI]"),
    ("dezeeuw_2013_bardoxolone", "de Zeeuw D[1AU] AND N Engl J Med[TA] AND 2013[DP] AND Bardoxolone[TI]"),
    ("qiu_2021_gse167363_sepsis", "Qiu X[1AU] AND 2021[DP] AND sepsis[TI] AND single-cell[TIAB]"),
]


def fetch_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "OmniNiche-Stage2-RefVerify/1.0"})
    with urllib.request.urlopen(req, timeout=25) as resp:
        return json.loads(resp.read().decode("utf-8"))


def main() -> None:
    out_path = pathlib.Path(r"F:\OmniNiche_Universal\outputs\stage2_tables\verified_pubmed_references.json")
    found_pmids: dict[str, str] = dict(CONFIRMED_PMIDS)
    for key, q in TARGETED_QUERIES:
        enc = urllib.parse.quote(q)
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=3&term={enc}"
        try:
            data = fetch_json(url)
            idlist = data.get("esearchresult", {}).get("idlist", [])
            if idlist:
                found_pmids[key] = idlist[0]
                print(f"[OK] {key} -> PMID {idlist[0]}")
            else:
                print(f"[MISS] {key} -> NO HIT for query: {q}")
        except Exception as exc:
            print(f"[ERR] {key} -> {exc}")
        time.sleep(0.35)

    all_ids = ",".join(found_pmids.values())
    sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id={all_ids}"
    summaries = fetch_json(sum_url).get("result", {})

    verified_records = []
    for key, pmid in found_pmids.items():
        rec = summaries.get(pmid, {})
        authors = [a.get("name", "") for a in rec.get("authors", [])]
        doi = ""
        for aid in rec.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
        entry = {
            "key": key,
            "pmid": pmid,
            "authors": authors,
            "first_author": authors[0] if authors else "",
            "title": rec.get("title", ""),
            "journal": rec.get("source", ""),
            "fulljournalname": rec.get("fulljournalname", ""),
            "pubdate": rec.get("pubdate", ""),
            "year": rec.get("pubdate", "")[:4],
            "volume": rec.get("volume", ""),
            "issue": rec.get("issue", ""),
            "pages": rec.get("pages", ""),
            "doi": doi,
        }
        verified_records.append(entry)

    out_path.write_text(json.dumps(verified_records, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nSaved {len(verified_records)} live-verified PubMed records to {out_path}")


if __name__ == "__main__":
    main()
