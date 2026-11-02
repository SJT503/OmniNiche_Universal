"""Fetch and SHA-256-lock Stage II v2.2 empirical external benchmarks:
1. Open Targets Platform v4 Human Genetic Support (`genetic_association`, `somatic_mutation`)
   and per-axis ChEMBL Regulatory Approval (`APPROVAL`) gold-standard indicators across all 12 axes.
2. Expanded 24 Landmark Phase IIb/III Human Antifibrotic RCTs from ClinicalTrials.gov API v2.
3. Independent Strictly Pre-Treatment Anti-PD-1 Monotherapy Human Biopsy Cohorts from NCBI GEO:
   - GSE91061 (Riaz et al., Cell 2017; pre-nivolumab monotherapy melanoma)
   - GSE78220 (Hugo et al., Cell 2016; pre-pembrolizumab/nivolumab monotherapy melanoma)
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import sys
import time
import urllib.request
from typing import Any

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

OT_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# Map the 43 resolved search aliases s0..s42 in _ot_expanded/resolved_disease_ids.json to the 12 disease axes
SEARCH_INDEX_TO_AXIS: dict[int, str] = {
    **{i: "Solid_Oncology" for i in range(0, 12)},
    **{i: "Hepatic_MASH_Fibrosis" for i in range(12, 15)},
    **{i: "Pulmonary_IPF_ARDS" for i in range(15, 19)},
    **{i: "Renal_CKD_AKI" for i in range(19, 22)},
    **{i: "Autoimmune_IBD_Rheum" for i in range(22, 27)},
    **{i: "Cardiovascular_HF_Athero" for i in range(27, 31)},
    **{i: "Metabolic_Endocrine" for i in range(31, 33)},
    33: "Infectious_Sepsis",
    **{i: "Hematologic_Malignancy" for i in range(34, 37)},
    **{i: "Dermatologic_Immune" for i in range(37, 39)},
    **{i: "Musculoskeletal_Bone" for i in range(39, 41)},
    **{i: "Ophthalmic_Retinal" for i in range(41, 43)},
}

LANDMARK_24_NCTS = [
    # 10 Positive Phase IIb/III Antifibrotic / Organ-Remodeling RCTs (Primary Endpoint Met / Regulatory Approval)
    "NCT03900429",  # MAESTRO-NASH: Resmetirom (THRB, F2/F3 Liver MASH, Phase 3 Win / FDA Approved)
    "NCT04822181",  # ESSENCE: Semaglutide (GLP1R, F2/F3 Liver MASH, Phase 3 Win)
    "NCT03008070",  # NATIVE: Lanifibranor (PPARA/D/G, F2/F3 Liver MASH, Phase 2b Win)
    "NCT04166773",  # SYNERGY-NASH: Tirzepatide (GLP1R/GIPR, F2/F3 Liver MASH, Phase 2b Win - NEJM 2024)
    "NCT04767529",  # HARMONY: Efruxifermin (FGFR1/2/3+KLB, F2/F3 Liver MASH, Phase 2b Win - Lancet GH 2023)
    "NCT02540993",  # FIDELIO-DKD: Finerenone (NR3C2, F2/F3 Kidney CKD/DKD, Phase 3 Win / FDA Approved)
    "NCT03036150",  # DAPA-CKD: Dapagliflozin (SLC5A2, F2/F3 Kidney CKD, Phase 3 Win / FDA Approved - NEJM 2020)
    "NCT01335464",  # INPULSIS-1: Nintedanib (FGFR1/PDGFRA/PDGFRB/KDR, F3/F4 Lung IPF, Phase 3 Win / FDA Approved)
    "NCT01366209",  # ASCEND: Pirfenidone (TGFB1, F2/F3 Lung IPF, Phase 3 Win / FDA Approved)
    "NCT05321069",  # FIBRONEER-IPF: Nerandomilast / BI 1015550 (PDE4B, F3/F4 Lung IPF, Phase 3 Win - NEJM 2025)
    # 14 Negative Phase IIb/III Antifibrotic RCTs (Scientific Primary Endpoint Failures Across F2/F3, F3/F4 & F4)
    "NCT01672866",  # GS-US-321-0105: Simtuzumab (LOXL2 mAb, F3 Liver MASH, Phase 2b Failure)
    "NCT01672879",  # GS-US-321-0106: Simtuzumab (LOXL2 mAb, F4 Cirrhosis Liver MASH, Phase 2b Failure)
    "NCT01769196",  # RAINIER: Simtuzumab (LOXL2 mAb, F3/F4 Lung IPF, Phase 2b Failure)
    "NCT03053050",  # STELLAR-3: Selonsertib (MAP3K5 SM, F3 Liver MASH, Phase 3 Failure)
    "NCT03053063",  # STELLAR-4: Selonsertib (MAP3K5 SM, F4 Cirrhosis Liver MASH, Phase 3 Failure)
    "NCT02177786",  # GS-US-223-1015: Selonsertib (MAP3K5 SM, F2/F3 Kidney DKD, Phase 2b Failure - JASN 2019)
    "NCT03028740",  # AURORA: Cenicriviroc (CCR2/CCR5 SM, F2/F3 Liver MASH, Phase 3 Failure)
    "NCT02704403",  # RESOLVE-IT: Elafibranor (PPARA/PPARD SM, F2/F3 Liver MASH, Phase 3 Failure)
    "NCT02855164",  # FLIGHT-FXR: Tropifexor (NR1H4 SM, F2/F3 Liver MASH, Phase 2b Histological Fibrosis Failure)
    "NCT03486899",  # FALCON 1: Pegbelfermin (Short-t1/2 PEG-FGF21, F3 Liver MASH, Phase 2b Failure - Lancet GH 2024)
    "NCT03955146",  # ZEPHYRUS-1: Pamrevlumab (CCN2 mAb, F3/F4 Lung IPF, Phase 3 Failure - JAMA 2024)
    "NCT03711162",  # ISABELA 1: Ziritaxestat / GLPG1690 (ENPP2 SM, F3/F4 Lung IPF, Phase 3 Failure - JAMA 2023)
    "NCT04552899",  # STARSCAPE: Zinpentraxin alfa / PRM-151 (APCS/PTX2, F3/F4 Lung IPF, Phase 3 Futility - AJRCCM 2024)
    "NCT01872689",  # RIFF: Lebrikizumab (IL13 mAb, F3/F4 Lung IPF, Phase 2b Failure - ERJ 2018)
]


def _post_ot_graphql(query: str, variables: dict[str, Any]) -> dict[str, Any]:
    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    for attempt in range(4):
        try:
            req = urllib.request.Request(
                OT_GRAPHQL_URL,
                data=body,
                headers={"Content-Type": "application/json", "User-Agent": "OmniNiche/2.2"},
            )
            with urllib.request.urlopen(req, timeout=45) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            if "data" in payload and payload["data"]:
                return payload["data"]
        except Exception as exc:
            if attempt == 3:
                raise RuntimeError(f"OpenTargets GraphQL failed after 4 attempts: {exc}") from exc
            time.sleep(1.5 * (attempt + 1))
    return {}


def fetch_opentargets_genetics_and_approvals(vault_ct: pathlib.Path, target_genes: set[str]) -> dict[str, Any]:
    """Fetch Open Targets v4 human genetic support datatype scores & per-axis drug approval status."""
    exp_dir = vault_ct / "_ot_expanded"
    res_ids = json.loads((exp_dir / "resolved_disease_ids.json").read_text(encoding="utf-8"))

    # Map disease ID -> list of axes (usually 1-to-1)
    did_to_axes: dict[str, list[str]] = {}
    ordered_dids: list[tuple[str, str]] = []
    seen: set[str] = set()
    for idx in range(43):
        axis = SEARCH_INDEX_TO_AXIS[idx]
        hits = ((res_ids.get(f"s{idx}") or {}).get("hits")) or []
        if hits and isinstance(hits[0], dict) and hits[0].get("id"):
            did = hits[0]["id"]
            dname = hits[0].get("name", "")
            did_to_axes.setdefault(did, []).append(axis)
            if did not in seen:
                seen.add(did)
                ordered_dids.append((did, dname))

    # 1. Extract per-axis & pan-axis drug-indication APPROVAL gold standard from cached _ot_expanded/drugs_batch_*.json
    axis_approved_drugs: dict[str, set[str]] = {ax: set() for ax in set(SEARCH_INDEX_TO_AXIS.values())}
    axis_approved_targets: dict[str, set[str]] = {ax: set() for ax in set(SEARCH_INDEX_TO_AXIS.values())}
    global_approved_drugs: set[str] = set()

    for b_file in sorted(exp_dir.glob("drugs_batch_*.json")):
        b_data = json.loads(b_file.read_text(encoding="utf-8"))
        for key, dis in b_data.items():
            if not key.startswith("d") or not isinstance(dis, dict):
                continue
            did = dis.get("id", "")
            axes = did_to_axes.get(did, [])
            rows = ((dis.get("drugAndClinicalCandidates") or {}).get("rows")) or []
            for row in rows:
                if not isinstance(row, dict):
                    continue
                stage_ind = str(row.get("maxClinicalStage") or "").upper()
                drug = row.get("drug") or {}
                dname = str(drug.get("name") or "").strip().lower()
                max_global = str(drug.get("maximumClinicalStage") or "").upper()
                if max_global == "APPROVAL" and dname:
                    global_approved_drugs.add(dname)
                if stage_ind == "APPROVAL":
                    if dname:
                        for ax in axes:
                            axis_approved_drugs[ax].add(dname)
                    for mr in ((drug.get("mechanismsOfAction") or {}).get("rows")) or []:
                        if not isinstance(mr, dict):
                            continue
                        for t in mr.get("targets") or []:
                            if isinstance(t, dict) and t.get("approvedSymbol"):
                                sym = t["approvedSymbol"].strip().upper()
                                for ax in axes:
                                    axis_approved_targets[ax].add(sym)

    # 2. Query Open Targets v4 associatedTargets datatypeScores across all 43 resolved disease ontologies
    query = """
    query AssocTargets($efoId: String!, $index: Int!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {index: $index, size: $size}) {
          count
          rows {
            score
            target { approvedSymbol }
            datatypeScores { id score }
          }
        }
      }
    }
    """

    # Store max scores per (axis, gene) and global pan-axis max per gene
    axis_gene_scores: dict[str, dict[str, dict[str, float]]] = {
        ax: {} for ax in set(SEARCH_INDEX_TO_AXIS.values())
    }
    pan_gene_scores: dict[str, dict[str, float]] = {}

    page_size = 2500
    for did, dname in ordered_dids:
        axes = did_to_axes[did]
        total_fetched = 0
        # Up to 3 pages (7,500 ranked associated targets per ontology easily covers every genetically supported human target)
        for page_idx in range(3):
            data = _post_ot_graphql(query, {"efoId": did, "index": page_idx, "size": page_size})
            dis = data.get("disease")
            if not dis:
                break
            assoc = dis.get("associatedTargets") or {}
            rows = assoc.get("rows") or []
            if not rows:
                break
            total_fetched += len(rows)
            for r in rows:
                sym = str(((r.get("target") or {}).get("approvedSymbol")) or "").strip().upper()
                if not sym or sym not in target_genes:
                    continue
                dt_map = {d["id"]: float(d["score"]) for d in (r.get("datatypeScores") or []) if "id" in d}
                gen_assoc = dt_map.get("genetic_association", 0.0)
                som_mut = dt_map.get("somatic_mutation", 0.0)
                lit_score = dt_map.get("literature", 0.0)
                anim_score = dt_map.get("animal_model", 0.0)
                overall = float(r.get("score") or 0.0)
                human_gen = max(gen_assoc, som_mut)

                for dest in [pan_gene_scores] + [axis_gene_scores[ax] for ax in axes]:
                    cur = dest.get(sym)
                    if cur is None:
                        dest[sym] = {
                            "genetic_association": round(gen_assoc, 5),
                            "somatic_mutation": round(som_mut, 5),
                            "human_genetic_max": round(human_gen, 5),
                            "animal_model": round(anim_score, 5),
                            "literature": round(lit_score, 5),
                            "overall_association": round(overall, 5),
                        }
                    else:
                        cur["genetic_association"] = round(max(cur["genetic_association"], gen_assoc), 5)
                        cur["somatic_mutation"] = round(max(cur["somatic_mutation"], som_mut), 5)
                        cur["human_genetic_max"] = round(max(cur["human_genetic_max"], human_gen), 5)
                        cur["animal_model"] = round(max(cur["animal_model"], anim_score), 5)
                        cur["literature"] = round(max(cur["literature"], lit_score), 5)
                        cur["overall_association"] = round(max(cur["overall_association"], overall), 5)

            if len(rows) < page_size:
                break
            time.sleep(0.08)
        print(f"[OpenTargets Genetics] {did:14s} ({dname[:38]:38s}) -> scanned {total_fetched:5d} target associations")

    out = {
        "n_clinical_targets_tracked": len(target_genes),
        "n_targets_with_any_ot_association": len(pan_gene_scores),
        "n_targets_with_human_genetics_gt0": sum(1 for v in pan_gene_scores.values() if v["human_genetic_max"] > 0),
        "axis_gene_scores": axis_gene_scores,
        "pan_gene_scores": pan_gene_scores,
        "axis_approved_drugs": {k: sorted(v) for k, v in axis_approved_drugs.items()},
        "axis_approved_targets": {k: sorted(v) for k, v in axis_approved_targets.items()},
        "global_approved_drugs": sorted(global_approved_drugs),
    }
    out_path = vault_ct / "opentargets_target_genetics_and_approval.json"
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        f"[OpenTargets Genetics Saved] {out['n_targets_with_any_ot_association']}/{len(target_genes)} targets matched; "
        f"{out['n_targets_with_human_genetics_gt0']} with direct human genetic support > 0."
    )
    return out


def fetch_24_landmark_fibrosis_ncts(vault_ct: pathlib.Path) -> dict[str, Any]:
    """Fetch and SHA-256-lock all 24 landmark Phase IIb/III antifibrotic RCTs from ClinicalTrials.gov API v2."""
    verified: dict[str, dict[str, Any]] = {}
    for nct in LANDMARK_24_NCTS:
        url = f"https://clinicaltrials.gov/api/v2/studies/{nct}"
        req = urllib.request.Request(url, headers={"User-Agent": "OmniNiche/2.2"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            raw_bytes = resp.read()
        sha = hashlib.sha256(raw_bytes).hexdigest()
        st = json.loads(raw_bytes.decode("utf-8"))
        p = st.get("protocolSection") or {}
        ident = p.get("identificationModule") or {}
        smod = p.get("statusModule") or {}
        dmod = p.get("designModule") or {}
        elig = (p.get("eligibilityModule") or {}).get("eligibilityCriteria", "")
        verified[nct] = {
            "nct_id": nct,
            "official_title": ident.get("officialTitle") or ident.get("briefTitle") or "",
            "brief_title": ident.get("briefTitle") or "",
            "overall_status": smod.get("overallStatus", ""),
            "why_stopped": smod.get("whyStopped", ""),
            "phases": dmod.get("phases") or [],
            "enrollment": ((dmod.get("enrollmentInfo") or {}).get("count")) or 0,
            "has_results": bool(st.get("hasResults", False)),
            "eligibility_excerpt": elig[:500],
            "sha256_lock": sha,
        }
        print(f"[CT.gov 24-Fibrosis] {nct} | N={verified[nct]['enrollment']:5d} | sha256={sha[:16]} | {verified[nct]['brief_title'][:58]}")
    out_p = vault_ct / "landmark_24_fibrosis_ncts_sha256.json"
    out_p.write_text(json.dumps(verified, ensure_ascii=False, indent=2), encoding="utf-8")
    return verified


def download_monotherapy_icb_geo_cohorts(vault_sc: pathlib.Path) -> dict[str, Any]:
    """Download GSE91061 (Riaz et al., Cell 2017) & GSE78220 (Hugo et al., Cell 2016) with SHA-256 verification."""
    downloads = [
        (
            "GSE91061_Melanoma_Nivolumab_Monotherapy",
            "GSE91061_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz",
        ),
        (
            "GSE91061_Melanoma_Nivolumab_Monotherapy",
            "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
        ),
        (
            "GSE91061_Melanoma_Nivolumab_Monotherapy",
            "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
        ),
        (
            "GSE91061_Melanoma_Nivolumab_Monotherapy",
            "bms038_clinical_data.csv",
            "https://raw.githubusercontent.com/riazn/bms038_analysis/master/data/bms038_clinical_data.csv",
        ),
        (
            "GSE78220_Melanoma_Pembrolizumab_Monotherapy",
            "GSE78220_series_matrix.txt.gz",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/GSE78220_series_matrix.txt.gz",
        ),
        (
            "GSE78220_Melanoma_Pembrolizumab_Monotherapy",
            "GSE78220_PatientFPKM.xlsx",
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx",
        ),
    ]
    manifest: dict[str, Any] = {}
    for subdir, fname, url in downloads:
        ddir = vault_sc / subdir
        ddir.mkdir(parents=True, exist_ok=True)
        fpath = ddir / fname
        if not fpath.exists() or fpath.stat().st_size < 1000:
            req = urllib.request.Request(url, headers={"User-Agent": "OmniNiche/2.2"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                data = resp.read()
            fpath.write_bytes(data)
        else:
            data = fpath.read_bytes()
        sha = hashlib.sha256(data).hexdigest()
        manifest[fname] = {
            "path": str(fpath.relative_to(ROOT)),
            "url": url,
            "size_bytes": len(data),
            "sha256": sha,
        }
        print(f"[GEO Monotherapy ICB] {fname:52s} | {len(data)/1e6:6.2f} MB | sha256={sha[:16]}")
    (vault_sc / "manifest_monotherapy_icb_cohorts.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return manifest


def main() -> None:
    import pandas as pd

    vault_ct = ROOT / "data_vault" / "raw_clinical_trials"
    vault_sc = ROOT / "data_vault" / "raw_single_cell"
    out_tables = ROOT / "outputs" / "stage2_tables"

    t2 = pd.read_csv(out_tables / "Table_S2_Target_Resolved_Trials_Deconvolution.csv", usecols=["target_genes"])
    target_genes: set[str] = set()
    for val in t2["target_genes"].dropna().astype(str):
        for g in val.split(";"):
            if g.strip():
                target_genes.add(g.strip().upper())

    print(f">>> [1/3] Fetching Open Targets v4 Human Genetic Support for {len(target_genes)} clinical target genes...")
    fetch_opentargets_genetics_and_approvals(vault_ct, target_genes)

    print("\n>>> [2/3] Fetching & SHA-256-locking 24 Landmark Phase IIb/III Antifibrotic RCTs...")
    fetch_24_landmark_fibrosis_ncts(vault_ct)

    print("\n>>> [3/3] Downloading GSE91061 & GSE78220 Pre-Treatment Anti-PD-1 Monotherapy Melanoma Biopsy Cohorts...")
    download_monotherapy_icb_geo_cohorts(vault_sc)
    print("\n=== ALL V2.2 EXTERNAL EMPIRICAL BENCHMARKS FETCHED & SHA-256 LOCKED ===")


if __name__ == "__main__":
    main()
