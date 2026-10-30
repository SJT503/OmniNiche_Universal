"""Stage II Task 1: Macro-Atlas Non-CNS ClinicalTrials.gov v2 Census & OpenTargets Drug-Target Engine.

Harvests >= 30,000 verbatim Phase II-IV non-CNS interventional drug/biologic clinical trials
from ClinicalTrials.gov API v2 across 12 major non-CNS disease axes, locks every study with
SHA-256 (Firewall-1), strictly blocks any CNS/neurological/psychiatric studies (Firewall-2),
deconvolves scientific vs administrative trial terminations, and maps interventions to official
human target gene symbols via OpenTargets Platform (ChEMBL MoA).
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import subprocess
import time
import urllib.parse
import urllib.request
from typing import Any

import pandas as pd

CT_API_V2_URL = "https://clinicaltrials.gov/api/v2/studies"
OPENTARGETS_SCRIPT = pathlib.Path(
    r"C:\Users\sheng\.gemini\config\plugins\science\skills\opentargets_database\scripts\query_opentargets.py"
)

# Strict CNS/Neurological/Psychiatric exclusion regex (0% overlap with nBME-26-6122)
CNS_EXCLUSION_REGEX = re.compile(
    r"\b("
    r"alzheimer|parkinson|huntington|amyotrophic lateral sclerosis|als\b|multiple sclerosis|"
    r"traumatic brain|tbi\b|spinal cord injur|stroke|cerebral ischemia|intracerebral hemorrhage|"
    r"subarachnoid hemorrhage|epilepsy|seizure|glioma|glioblastoma|gbm\b|astrocytoma|"
    r"medulloblastoma|meningioma|neuroblastoma|brain metastas|brain tumor|cns lymphoma|"
    r"schizophrenia|bipolar disorder|major depressive|depression|autism|adhd|migraine|"
    r"dementia|frontotemporal|lewy body|myasthenia gravis|guillain-barre|neuropathy|"
    r"encephalitis|meningitis|narcolepsy|insomnia|anxiety disorder|post-traumatic stress|ptsd"
    r")",
    re.IGNORECASE,
)

EFFICACY_FUTILITY_REGEX = re.compile(
    r"(lack of efficacy|futility|insufficient efficacy|no efficacy|primary endpoint not met|"
    r"did not meet primary|failed to meet|unfavorable risk-benefit|benefit-risk|"
    r"interim analysis.*futil|no clinical benefit|suboptimal efficacy|low response rate|"
    r"inferior efficacy|disease progression)",
    re.IGNORECASE,
)

SAFETY_TOXICITY_REGEX = re.compile(
    r"(toxicity|adverse event|safety|hepatotox|cardiotox|nephrotox|pneumonitis|"
    r"dose-limiting|dlt\b|mortality imbalance|increased mortality|death|fatal|"
    r"tolerability|intolerance|clinical hold|fda hold|sae\b|serious adverse)",
    re.IGNORECASE,
)

ADMIN_ENROLLMENT_REGEX = re.compile(
    r"(accrual|recruitment|enrollment|enrolment|slow recruiting|low enrollment|"
    r"business|strategic|commercial|sponsor decision|portfolio|reprioritiz|"
    r"funding|financial|budget|resource|covid-19|pandemic|supply|manufacturing|"
    r"drug availability|investigator left|pi relocated|site closure|administrative|"
    r"protocol redesign|replaced by another study|regulatory delay)",
    re.IGNORECASE,
)

PLACEBO_COMPARATOR_STOPWORDS = {
    "placebo",
    "standard of care",
    "soc",
    "normal saline",
    "saline",
    "vehicle",
    "sham",
    "observation",
    "best supportive care",
    "bsc",
    "matching placebo",
    "oral placebo",
    "iv placebo",
    "usual care",
    "control",
    "no intervention",
}

SALT_AND_FORMULATION_SUFFIXES = (
    " hydrochloride",
    " dihydrochloride",
    " mesylate",
    " esylate",
    " tosylate",
    " maleate",
    " fumarate",
    " tartrate",
    " sulfate",
    " sodium",
    " potassium",
    " calcium",
    " acetate",
    " succinate",
    " phosphate",
    " citrate",
    " hydrobromide",
    " monohydrate",
    " dihydrate",
    " tablet",
    " tablets",
    " capsule",
    " capsules",
    " injection",
    " infusion",
    " solution",
    " oral",
    " intravenous",
    " subcutaneous",
)

# 12 Non-CNS Disease Axes with ClinicalTrials.gov condition queries and OpenTargets EFO/MONDO IDs
NON_CNS_12_AXES: list[dict[str, Any]] = [
    {
        "axis_id": "Solid_Oncology",
        "ct_cond": '(lung cancer OR breast cancer OR colorectal cancer OR pancreatic cancer OR melanoma OR renal cell carcinoma OR hepatocellular carcinoma OR gastric cancer OR ovarian cancer OR prostate cancer OR bladder cancer OR esophageal cancer OR head and neck squamous)',
        "ot_ids": ["EFO_0003060", "EFO_0000305", "EFO_0005842", "EFO_0002618", "EFO_0000756", "EFO_0000681", "EFO_0000182", "EFO_0000292"],
    },
    {
        "axis_id": "Hepatic_MASH_Fibrosis",
        "ct_cond": '(steatohepatitis OR NASH OR MASH OR NAFLD OR MASLD OR liver fibrosis OR liver cirrhosis OR primary biliary cholangitis OR primary sclerosing cholangitis)',
        "ot_ids": ["MONDO_0007027", "EFO_0003095", "EFO_0001422", "EFO_1001486", "EFO_0004264"],
    },
    {
        "axis_id": "Pulmonary_IPF_ARDS",
        "ct_cond": '(idiopathic pulmonary fibrosis OR pulmonary fibrosis OR interstitial lung disease OR acute respiratory distress syndrome OR COPD OR asthma)',
        "ot_ids": ["EFO_0000768", "EFO_0004244", "EFO_0000272", "EFO_0000341", "EFO_0000270"],
    },
    {
        "axis_id": "Renal_CKD_AKI",
        "ct_cond": '(chronic kidney disease OR diabetic nephropathy OR diabetic kidney disease OR acute kidney injury OR glomerulonephritis OR IgA nephropathy OR focal segmental glomerulosclerosis)',
        "ot_ids": ["EFO_0003884", "EFO_0000401", "EFO_0003086", "EFO_0004194", "MONDO_0005301"],
    },
    {
        "axis_id": "Autoimmune_IBD_Rheum",
        "ct_cond": '(ulcerative colitis OR Crohn disease OR inflammatory bowel disease OR rheumatoid arthritis OR systemic lupus erythematosus OR systemic sclerosis OR ankylosing spondylitis)',
        "ot_ids": ["EFO_0000729", "EFO_0000384", "EFO_0003767", "EFO_0000685", "EFO_0002690", "EFO_0000717"],
    },
    {
        "axis_id": "Cardiovascular_HF_Athero",
        "ct_cond": '(heart failure OR myocardial infarction OR atherosclerosis OR coronary artery disease OR pulmonary arterial hypertension OR cardiomyopathy)',
        "ot_ids": ["EFO_0003144", "EFO_0000612", "EFO_0003914", "EFO_0001645", "EFO_0001361"],
    },
    {
        "axis_id": "Metabolic_Endocrine",
        "ct_cond": '(type 2 diabetes OR obesity OR dyslipidemia OR hypercholesterolemia OR metabolic syndrome)',
        "ot_ids": ["EFO_0001360", "EFO_0001073", "EFO_0004530"],
    },
    {
        "axis_id": "Infectious_Sepsis",
        "ct_cond": '(sepsis OR septic shock OR bacteremia OR bacterial pneumonia OR systemic fungal infection)',
        "ot_ids": ["EFO_0005208", "EFO_0003106"],
    },
    {
        "axis_id": "Hematologic_Malignancy",
        "ct_cond": '(acute myeloid leukemia OR chronic lymphocytic leukemia OR multiple myeloma OR diffuse large B-cell lymphoma OR myelofibrosis)',
        "ot_ids": ["EFO_0000222", "EFO_0000095", "EFO_0001378", "EFO_0000403"],
    },
    {
        "axis_id": "Dermatologic_Immune",
        "ct_cond": '(psoriasis OR atopic dermatitis OR hidradenitis suppurativa OR vitiligo)',
        "ot_ids": ["EFO_0000676", "EFO_0000274", "EFO_1000967"],
    },
    {
        "axis_id": "Musculoskeletal_Bone",
        "ct_cond": '(osteoarthritis OR osteoporosis OR psoriatic arthritis OR gout)',
        "ot_ids": ["EFO_0002506", "EFO_0003882", "EFO_0003770"],
    },
    {
        "axis_id": "Ophthalmic_Retinal",
        "ct_cond": '(macular degeneration OR diabetic retinopathy OR diabetic macular edema OR non-infectious uveitis)',
        "ot_ids": ["EFO_0001365", "EFO_0003773"],
    },
]


def is_cns_excluded_trial(conditions: list[str], title: str) -> bool:
    """Return True if any condition or title hits CNS/neurological/psychiatric terms."""
    blob = " | ".join(conditions) + " | " + (title or "")
    return bool(CNS_EXCLUSION_REGEX.search(blob))


def classify_trial_stop_deconvolution(overall_status: str, why_stopped: str) -> str:
    """Deconvolve terminated/withdrawn trials into scientific efficacy, safety, or administrative reasons."""
    status = (overall_status or "").upper().strip()
    if status == "COMPLETED":
        return "Completed"
    reason = (why_stopped or "").strip()
    if not reason:
        return "Unreported_Stop_Reason"
    # Check scientific safety first if explicit toxicity/SAE caused termination
    has_eff = bool(EFFICACY_FUTILITY_REGEX.search(reason))
    has_saf = bool(SAFETY_TOXICITY_REGEX.search(reason))
    has_adm = bool(ADMIN_ENROLLMENT_REGEX.search(reason))
    if has_eff and not has_saf:
        return "Scientific_Efficacy_Futility"
    if has_saf and not has_eff:
        return "Scientific_Safety_Toxicity"
    if has_eff and has_saf:
        return "Scientific_Efficacy_Futility"
    if has_adm:
        return "Administrative_Commercial_Enrollment"
    return "Other_Unclassified_Stop"


def normalize_intervention_drug_name(raw_name: str) -> str:
    """Normalize clinical trial intervention string to canonical lowercase drug stem."""
    s = (raw_name or "").strip().lower()
    s = re.sub(r"^(drug|biological|biologic|combination|treatment)\s*:\s*", "", s)
    # Remove parenthetical research codes or brand explanations
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub(r"\[[^\]]*\]", " ", s)
    # Strip dosage numbers and units onward (e.g. '200 mg iv every 3 weeks', '150mg')
    s = re.sub(r"\b\d+(\.\d+)?\s*(mg|mcg|ug|g|ml|iu|units|%|mg/kg|mg/m2)\b.*$", " ", s)
    s = re.sub(r"\b(iv|po|sc|sq|im|qd|bid|tid|q2w|q3w|q4w|every\s+\d+\s+\w+)\b.*$", " ", s)
    s = " ".join(s.split())
    for suf in SALT_AND_FORMULATION_SUFFIXES:
        if s.endswith(suf):
            s = s[: -len(suf)].strip()
    s = re.sub(r"[^a-z0-9\- ]+", " ", s)
    return " ".join(s.split())


def resolve_intervention_targets(
    interventions: list[str],
    drug_target_lookup: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Resolve raw clinical trial interventions to ChEMBL/OpenTargets drug names and target gene symbols."""
    matched_drugs: list[str] = []
    target_genes: list[str] = []
    modalities: list[str] = []
    moas: list[str] = []

    for raw in interventions:
        norm = normalize_intervention_drug_name(raw)
        if not norm or norm in PLACEBO_COMPARATOR_STOPWORDS or "placebo" in norm:
            continue
        # Try exact normalized match, then first-token match (for multi-word drug + regimen descriptions)
        candidates = [norm]
        parts = re.split(r"\s+(?:\+|plus|and|with|/)\s+", norm)
        candidates.extend(p.strip() for p in parts if len(p.strip()) >= 4)
        first_tok = norm.split()[0] if norm.split() else ""
        if len(first_tok) >= 5:
            candidates.append(first_tok)

        for cand in candidates:
            info = drug_target_lookup.get(cand)
            if info and info.get("targets"):
                canonical = info.get("canonical_name", cand)
                if canonical not in matched_drugs:
                    matched_drugs.append(canonical)
                for g in info["targets"]:
                    if g and g not in target_genes:
                        target_genes.append(g)
                if info.get("drug_type") and info["drug_type"] not in modalities:
                    modalities.append(info["drug_type"])
                if info.get("moa") and info["moa"] not in moas:
                    moas.append(info["moa"])
                break

    return {
        "matched_drugs": matched_drugs,
        "target_genes": target_genes,
        "drug_modalities": modalities,
        "mechanisms_of_action": moas,
    }


def build_opentargets_drug_target_dictionary(out_json_path: pathlib.Path) -> dict[str, dict[str, Any]]:
    """Query OpenTargets Platform via official skill wrapper across all 12 non-CNS disease ontologies."""
    out_json_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_json_path.parent / "_ot_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    query = """query getDrugsWithMoA($efoId: String!) {
      disease(efoId: $efoId) {
        id
        name
        drugAndClinicalCandidates {
          count
          rows {
            maxClinicalStage
            drug {
              id
              name
              synonyms { label }
              drugType
              maximumClinicalStage
              mechanismsOfAction {
                rows {
                  mechanismOfAction
                  actionType
                  targets { id approvedSymbol }
                }
              }
            }
          }
        }
      }
    }"""

    all_efo_ids: list[str] = []
    for ax in NON_CNS_12_AXES:
        for eid in ax["ot_ids"]:
            if eid not in all_efo_ids:
                all_efo_ids.append(eid)
    # Also include broad non-CNS root ontologies in OpenTargets to maximize drug-target recall
    broad_roots = ["EFO_0000616", "EFO_0000408", "EFO_0000319", "EFO_0000540", "EFO_0000589"]
    for br in broad_roots:
        if br not in all_efo_ids:
            all_efo_ids.append(br)

    lookup: dict[str, dict[str, Any]] = {}
    queried_ok = 0

    for eid in all_efo_ids:
        tmp_out = tmp_dir / f"ot_{eid}.json"
        if not tmp_out.exists():
            cmd = [
                "uv",
                "run",
                str(OPENTARGETS_SCRIPT),
                "--output",
                str(tmp_out),
                "--limit",
                "10000",
                "custom-query",
                query,
                "--variables",
                json.dumps({"efoId": eid}),
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0 or not tmp_out.exists():
                continue
        try:
            payload = json.loads(tmp_out.read_text(encoding="utf-8"))
        except Exception:
            continue
        dis = payload.get("disease")
        if not dis:
            continue
        queried_ok += 1
        rows = (dis.get("drugAndClinicalCandidates") or {}).get("rows") or []
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
            lookup[normalize_intervention_drug_name(dname)] = entry
            lookup[dname] = entry
            for syn in drug.get("synonyms") or []:
                if isinstance(syn, dict) and syn.get("label"):
                    sn = normalize_intervention_drug_name(syn["label"])
                    if len(sn) >= 4 and sn not in PLACEBO_COMPARATOR_STOPWORDS and sn not in lookup:
                        lookup[sn] = entry

    out_json_path.write_text(json.dumps(lookup, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[OpenTargets] Queried {queried_ok}/{len(all_efo_ids)} ontologies -> {len(lookup)} normalized drug/synonym -> gene mappings.")
    return lookup


def harvest_non_cns_clinical_trials_macro_atlas(
    raw_jsonl_path: pathlib.Path,
    drug_target_lookup: dict[str, dict[str, Any]],
    max_pages_per_axis: int = 8,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Harvest >= 30,000 real verbatim non-CNS Phase 2-4 drug/biologic clinical trials via API v2."""
    raw_jsonl_path.parent.mkdir(parents=True, exist_ok=True)
    seen_ncts: set[str] = set()
    records: list[dict[str, Any]] = []

    prisma_counts = {
        "api_raw_records_fetched": 0,
        "duplicate_across_axes_removed": 0,
        "cns_neuro_firewall_excluded": 0,
        "final_unique_non_cns_trials": 0,
        "target_resolved_trials": 0,
        "per_axis_counts": {},
    }

    fields = (
        "NCTId,BriefTitle,OfficialTitle,OverallStatus,WhyStopped,HasResults,"
        "Phase,EnrollmentCount,StartDateStruct,ConditionsModule,ArmsInterventionsModule,SponsorCollaboratorsModule"
    )

    with raw_jsonl_path.open("w", encoding="utf-8") as fout:
        for axis in NON_CNS_12_AXES:
            axis_id = axis["axis_id"]
            cond_q = axis["ct_cond"]
            adv_filter = (
                "AREA[StudyType]INTERVENTIONAL AND "
                "(AREA[InterventionType]DRUG OR AREA[InterventionType]BIOLOGICAL) AND "
                "(AREA[Phase]PHASE2 OR AREA[Phase]PHASE3 OR AREA[Phase]PHASE4) AND "
                "(AREA[OverallStatus]COMPLETED OR AREA[OverallStatus]TERMINATED OR AREA[OverallStatus]WITHDRAWN)"
            )
            page_token: str | None = None
            axis_kept = 0
            axis_total_reported = 0

            for page_idx in range(max_pages_per_axis):
                params: dict[str, str] = {
                    "query.cond": cond_q,
                    "filter.advanced": adv_filter,
                    "fields": fields,
                    "pageSize": "1000",
                }
                if page_idx == 0:
                    params["countTotal"] = "true"
                if page_token:
                    params["pageToken"] = page_token

                url = CT_API_V2_URL + "?" + urllib.parse.urlencode(params)
                req = urllib.request.Request(url, headers={"User-Agent": "OmniNiche/2.0"})
                payload: dict[str, Any] | None = None
                for attempt in range(3):
                    try:
                        with urllib.request.urlopen(req, timeout=45) as resp:
                            payload = json.loads(resp.read().decode("utf-8"))
                        break
                    except Exception:
                        time.sleep(1.5 * (attempt + 1))
                if not payload:
                    break

                if page_idx == 0:
                    axis_total_reported = int(payload.get("totalCount") or 0)

                studies = payload.get("studies") or []
                prisma_counts["api_raw_records_fetched"] += len(studies)

                for st in studies:
                    raw_str = json.dumps(st, sort_keys=True, ensure_ascii=False)
                    sha256_hash = hashlib.sha256(raw_str.encode("utf-8")).hexdigest()

                    p = st.get("protocolSection") or {}
                    ident = p.get("identificationModule") or {}
                    nct_id = ident.get("nctId", "")
                    if not nct_id:
                        continue
                    if nct_id in seen_ncts:
                        prisma_counts["duplicate_across_axes_removed"] += 1
                        continue

                    title = ident.get("briefTitle") or ident.get("officialTitle") or ""
                    conds = ((p.get("conditionsModule") or {}).get("conditions")) or []
                    if is_cns_excluded_trial(conds, title):
                        prisma_counts["cns_neuro_firewall_excluded"] += 1
                        continue

                    seen_ncts.add(nct_id)
                    status_mod = p.get("statusModule") or {}
                    overall_status = status_mod.get("overallStatus", "")
                    why_stopped = status_mod.get("whyStopped", "")
                    start_date = ((status_mod.get("startDateStruct") or {}).get("date")) or ""
                    design_mod = p.get("designModule") or {}
                    phases = design_mod.get("phases") or []
                    enrollment = ((design_mod.get("enrollmentInfo") or {}).get("count")) or 0
                    has_results = bool(st.get("hasResults", False))
                    sponsor = (((p.get("sponsorCollaboratorsModule") or {}).get("leadSponsor") or {}).get("class")) or "UNKNOWN"

                    arms_mod = p.get("armsInterventionsModule") or {}
                    interventions_raw = [
                        iv.get("name", "")
                        for iv in (arms_mod.get("interventions") or [])
                        if iv.get("type") in ("DRUG", "BIOLOGICAL") and iv.get("name")
                    ]

                    stop_cat = classify_trial_stop_deconvolution(overall_status, why_stopped)
                    resolved = resolve_intervention_targets(interventions_raw, drug_target_lookup)

                    max_phase_num = 2
                    if "PHASE4" in phases:
                        max_phase_num = 4
                    elif "PHASE3" in phases:
                        max_phase_num = 3

                    rec = {
                        "nct_id": nct_id,
                        "disease_axis": axis_id,
                        "brief_title": title,
                        "overall_status": overall_status,
                        "why_stopped": why_stopped,
                        "stop_deconvolution": stop_cat,
                        "phases": ";".join(phases),
                        "max_phase_num": max_phase_num,
                        "has_results": has_results,
                        "enrollment": int(enrollment or 0),
                        "start_year": start_date[:4] if len(start_date) >= 4 else "",
                        "sponsor_class": sponsor,
                        "conditions": "; ".join(conds[:6]),
                        "interventions_raw": "; ".join(interventions_raw[:6]),
                        "matched_drugs": ";".join(resolved["matched_drugs"]),
                        "target_genes": ";".join(resolved["target_genes"]),
                        "primary_target_gene": resolved["target_genes"][0] if resolved["target_genes"] else "",
                        "drug_modalities": ";".join(resolved["drug_modalities"]),
                        "mechanisms_of_action": "; ".join(resolved["mechanisms_of_action"]),
                        "n_targets": len(resolved["target_genes"]),
                        "sha256_lock": sha256_hash,
                    }
                    records.append(rec)
                    axis_kept += 1
                    fout.write(json.dumps({"nct_id": nct_id, "sha256": sha256_hash, "axis": axis_id, "study": st}, ensure_ascii=False) + "\n")

                page_token = payload.get("nextPageToken")
                if not page_token:
                    break
                time.sleep(0.15)

            prisma_counts["per_axis_counts"][axis_id] = {
                "api_total_available": axis_total_reported,
                "unique_non_cns_kept": axis_kept,
            }
            print(f"[CT.gov v2] {axis_id:25s} | API total={axis_total_reported:6d} | Unique non-CNS kept={axis_kept:5d}")

    df = pd.DataFrame(records)
    prisma_counts["final_unique_non_cns_trials"] = int(len(df))
    prisma_counts["target_resolved_trials"] = int((df["n_targets"] > 0).sum()) if len(df) else 0
    return df, prisma_counts
