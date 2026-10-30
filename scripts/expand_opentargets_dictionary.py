"""Expand OpenTargets Drug-Target Dictionary across all 12 Non-CNS Disease Axes via live MONDO/EFO resolution.

Routes 100% of queries through the official `query_opentargets.py` wrapper (`uv run`) using batched
GraphQL aliases for fast, polite-http-compliant execution.
"""
from __future__ import annotations

import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.stage2.macro_clinical_census import (
    OPENTARGETS_SCRIPT,
    PLACEBO_COMPARATOR_STOPWORDS,
    normalize_intervention_drug_name,
)

DISEASE_SEARCH_TERMS = [
    "neoplasm",
    "non-small cell lung carcinoma",
    "breast carcinoma",
    "colorectal carcinoma",
    "pancreatic carcinoma",
    "melanoma",
    "renal cell carcinoma",
    "hepatocellular carcinoma",
    "gastric carcinoma",
    "ovarian carcinoma",
    "prostate carcinoma",
    "urinary bladder carcinoma",
    "metabolic dysfunction-associated steatohepatitis",
    "liver cirrhosis",
    "hepatic fibrosis",
    "idiopathic pulmonary fibrosis",
    "chronic obstructive pulmonary disease",
    "asthma",
    "acute respiratory distress syndrome",
    "chronic kidney disease",
    "diabetic kidney disease",
    "acute kidney injury",
    "ulcerative colitis",
    "Crohn disease",
    "rheumatoid arthritis",
    "systemic lupus erythematosus",
    "systemic scleroderma",
    "heart failure",
    "myocardial infarction",
    "atherosclerosis",
    "pulmonary arterial hypertension",
    "type 2 diabetes mellitus",
    "obesity",
    "sepsis",
    "acute myeloid leukemia",
    "multiple myeloma",
    "lymphoma",
    "psoriasis",
    "atopic dermatitis",
    "osteoarthritis",
    "osteoporosis",
    "age-related macular degeneration",
    "diabetic retinopathy",
]


def run_ot_custom_query(query: str, out_file: pathlib.Path, limit: int = 15000) -> dict:
    cmd = [
        "uv",
        "run",
        str(OPENTARGETS_SCRIPT),
        "--output",
        str(out_file),
        "--limit",
        str(limit),
        "custom-query",
        query,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"OpenTargets wrapper failed: {res.stderr[-400:]}")
    return json.loads(out_file.read_text(encoding="utf-8"))


def main() -> None:
    vault_ct = ROOT / "data_vault" / "raw_clinical_trials"
    tmp_dir = vault_ct / "_ot_expanded"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Resolve exact live OpenTargets v26 MONDO/EFO IDs via batched search aliases
    alias_lines = []
    for idx, term in enumerate(DISEASE_SEARCH_TERMS):
        safe_term = term.replace('"', "")
        alias_lines.append(
            f'  s{idx}: search(queryString: "{safe_term}", entityNames: ["disease"], page: {{index: 0, size: 2}}) {{ hits {{ id name }} }}'
        )
    search_q = "query resolveDiseaseIds {\n" + "\n".join(alias_lines) + "\n}"
    search_data = run_ot_custom_query(search_q, tmp_dir / "resolved_disease_ids.json", limit=50)

    disease_ids: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    for idx in range(len(DISEASE_SEARCH_TERMS)):
        hits = ((search_data.get(f"s{idx}") or {}).get("hits")) or []
        for h in hits[:1]:
            if isinstance(h, dict) and h.get("id") and h["id"] not in seen_ids:
                seen_ids.add(h["id"])
                disease_ids.append((h["id"], h.get("name", "")))

    print(f"[OpenTargets] Resolved {len(disease_ids)} live non-CNS MONDO/EFO disease ontologies.")

    # Step 2: Fetch full drugAndClinicalCandidates + synonyms + mechanismsOfAction in batches of 6
    ot_map_path = vault_ct / "opentargets_drug_target_map.json"
    lookup: dict[str, dict] = {}
    if ot_map_path.exists():
        try:
            lookup = json.loads(ot_map_path.read_text(encoding="utf-8"))
        except Exception:
            lookup = {}

    batch_size = 5
    for b_start in range(0, len(disease_ids), batch_size):
        batch = disease_ids[b_start : b_start + batch_size]
        parts = []
        for j, (did, _) in enumerate(batch):
            parts.append(
                f"""  d{j}: disease(efoId: "{did}") {{
    id
    name
    drugAndClinicalCandidates {{
      count
      rows {{
        maxClinicalStage
        drug {{
          id
          name
          synonyms {{ label }}
          drugType
          maximumClinicalStage
          mechanismsOfAction {{
            rows {{
              mechanismOfAction
              actionType
              targets {{ id approvedSymbol }}
            }}
          }}
        }}
      }}
    }}
  }}"""
            )
        q_batch = "query fetchDrugsBatch {\n" + "\n".join(parts) + "\n}"
        b_out = tmp_dir / f"drugs_batch_{b_start}.json"
        b_data = run_ot_custom_query(q_batch, b_out, limit=15000)

        for j in range(len(batch)):
            dis = b_data.get(f"d{j}") or {}
            rows = ((dis.get("drugAndClinicalCandidates") or {}).get("rows")) or []
            print(f"  + {dis.get('id', batch[j][0]):15s} ({dis.get('name', '')[:42]:42s}): {len(rows):5d} drug rows")
            for row in rows:
                if not isinstance(row, dict) or "_truncated" in row:
                    continue
                drug = row.get("drug") or {}
                dname = (drug.get("name") or "").strip().lower()
                if not dname:
                    continue
                moa_rows = ((drug.get("mechanismsOfAction") or {}).get("rows")) or []
                targets: list[str] = []
                moas: list[str] = []
                action_types: list[str] = []
                for mr in moa_rows:
                    if not isinstance(mr, dict) or "_truncated" in mr:
                        continue
                    if mr.get("mechanismOfAction"):
                        moas.append(mr["mechanismOfAction"])
                    if mr.get("actionType"):
                        action_types.append(mr["actionType"])
                    for t in mr.get("targets") or []:
                        if isinstance(t, dict) and t.get("approvedSymbol"):
                            sym = t["approvedSymbol"].strip()
                            if sym not in targets:
                                targets.append(sym)
                if not targets:
                    continue

                entry = {
                    "chembl_id": drug.get("id", ""),
                    "canonical_name": dname,
                    "drug_type": drug.get("drugType", "Unknown"),
                    "max_phase": drug.get("maximumClinicalStage", ""),
                    "targets": targets,
                    "moa": "; ".join(dict.fromkeys(moas)),
                    "action_types": list(dict.fromkeys(action_types)),
                }
                norm_d = normalize_intervention_drug_name(dname)
                if norm_d:
                    lookup[norm_d] = entry
                lookup[dname] = entry
                for syn in drug.get("synonyms") or []:
                    if isinstance(syn, dict) and syn.get("label"):
                        sn = normalize_intervention_drug_name(syn["label"])
                        if len(sn) >= 4 and sn not in PLACEBO_COMPARATOR_STOPWORDS and sn not in lookup:
                            lookup[sn] = entry

    ot_map_path.write_text(json.dumps(lookup, ensure_ascii=False, indent=2), encoding="utf-8")
    unique_drugs = len({v["chembl_id"] for v in lookup.values() if v.get("chembl_id")})
    unique_genes = len({g for v in lookup.values() for g in v.get("targets", [])})
    print(
        f"\n[OpenTargets Expanded Dictionary Complete] "
        f"Normalized Lookup Keys={len(lookup):,} | Unique ChEMBL Drugs={unique_drugs:,} | Unique Target Genes={unique_genes:,}"
    )


if __name__ == "__main__":
    main()
