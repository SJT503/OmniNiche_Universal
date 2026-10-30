from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Dict, Iterable, List, Optional, Tuple
import pandas as pd

from omniniche.core.integrity_firewall import IntegrityFirewall
from omniniche.core.schema import DiseasePack, TrialRecord


SCIENTIFIC_FAILURE_KEYWORDS: Tuple[str, ...] = (
    "futility",
    "lack of efficacy",
    "no efficacy",
    "insufficient efficacy",
    "ineffective",
    "primary endpoint not met",
    "did not meet primary",
    "adverse",
    "toxicity",
    "safety",
    "mortality",
    "dsmb",
    "dmc recommendation",
    "benefit-risk",
    "unfavorable risk",
)

ADMINISTRATIVE_STOP_KEYWORDS: Tuple[str, ...] = (
    "covid",
    "pandemic",
    "sars-cov-2",
    "accrual",
    "enrollment",
    "recruitment",
    "funding",
    "budget",
    "financial",
    "business",
    "strategic",
    "sponsor decision",
    "portfolio",
    "principal investigator left",
    "pi relocated",
    "drug supply",
    "logistics",
)


class ClinicalTrialistV2:
    """Two-tier (Macro-Atlas 10,000+ pan-disease & Meso DiseasePack) ClinicalTrials.gov REST API v2
    PRISMA auditor with verbatim SHA-256 integrity checks and administrative-vs-scientific stop deconvolution.
    """

    CTG_API_V2_URL = "https://clinicaltrials.gov/api/v2/studies"

    def __init__(self, pack: DiseasePack):
        self.pack = pack

    @staticmethod
    def classify_trial_outcome(record: TrialRecord) -> str:
        status = record.overall_status.upper()
        why = (record.why_stopped or "").lower()

        if status in ("TERMINATED", "WITHDRAWN", "SUSPENDED"):
            has_sci = any(k in why for k in SCIENTIFIC_FAILURE_KEYWORDS)
            has_adm = any(k in why for k in ADMINISTRATIVE_STOP_KEYWORDS)
            if has_sci and not (has_adm and "covid" in why and not any(s in why for s in ("futility", "toxicity", "mortality", "lack of efficacy"))):
                return "SCIENTIFIC_EFFICACY_OR_SAFETY_FAILURE"
            if has_adm:
                return "ADMINISTRATIVE_OR_OPERATIONAL_STOP"
            return "UNCLASSIFIED_STOPPED"

        if status in ("COMPLETED", "ACTIVE_NOT_RECRUITING", "RECRUITING", "ENROLLING_BY_INVITATION"):
            return "COMPLETED_ACTIVE_OR_SURVIVED"

        return "OTHER_STATUS"

    def fetch_live_studies(
        self,
        max_pages_per_term: Optional[int] = None,
        page_size: int = 100,
    ) -> List[dict]:
        """Fetch raw study JSON payloads directly from ClinicalTrials.gov API v2 with pagination."""
        seen_ncts: set[str] = set()
        raw_studies: List[dict] = []

        for term in self.pack.clinical_query_spec.condition_terms:
            page_token: Optional[str] = None
            page_idx = 0
            while True:
                params = {
                    "query.cond": term,
                    "pageSize": str(page_size),
                }
                if page_token:
                    params["pageToken"] = page_token
                url = f"{self.CTG_API_V2_URL}?{urllib.parse.urlencode(params)}"
                req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "OmniNiche-v2.0"})
                with urllib.request.urlopen(req, timeout=20) as resp:
                    payload = json.loads(resp.read().decode("utf-8"))

                for st in payload.get("studies", []):
                    nct = (
                        st.get("protocolSection", {})
                        .get("identificationModule", {})
                        .get("nctId", "")
                    )
                    if nct and nct not in seen_ncts:
                        seen_ncts.add(nct)
                        raw_studies.append(st)

                page_token = payload.get("nextPageToken")
                page_idx += 1
                if not page_token or (max_pages_per_term is not None and page_idx >= max_pages_per_term):
                    break

        return raw_studies

    def audit_raw_studies(self, raw_studies: Iterable[dict]) -> Tuple[pd.DataFrame, Dict[str, int]]:
        raw_list = list(raw_studies)
        step1_identified = len(raw_list)

        eligible_phases = {"PHASE2", "PHASE3", "PHASE4"}
        if self.pack.clinical_query_spec.min_phase.upper() == "PHASE1":
            eligible_phases.add("PHASE1")

        step2_records: List[Tuple[TrialRecord, dict]] = []
        for raw in raw_list:
            rec = TrialRecord.from_ctg_api_v2(raw)
            IntegrityFirewall.assert_trial_verbatim_integrity(rec, raw)
            if rec.study_type.upper() != "INTERVENTIONAL":
                continue
            if rec.phases and not any(p.upper() in eligible_phases for p in rec.phases):
                continue
            if not rec.interventions:
                continue
            step2_records.append((rec, raw))

        excl_lower = [e.lower() for e in self.pack.clinical_query_spec.exclusion_terms]
        inc_lower = [c.lower() for c in self.pack.clinical_query_spec.condition_terms]

        step3_clean: List[Dict[str, object]] = []
        for rec, raw in step2_records:
            IntegrityFirewall.assert_trial_verbatim_integrity(rec, raw)
            haystack = f"{rec.official_title} {rec.brief_title} {' '.join(rec.conditions)}".lower()
            if any(ex in haystack for ex in excl_lower):
                continue
            if inc_lower and not any(inc in haystack for inc in inc_lower):
                continue

            outcome_cls = self.classify_trial_outcome(rec)
            step3_clean.append(
                {
                    "nct_id": rec.nct_id,
                    "official_title": rec.official_title,
                    "overall_status": rec.overall_status,
                    "why_stopped": rec.why_stopped,
                    "phases": ";".join(rec.phases),
                    "enrollment": rec.enrollment,
                    "interventions": ";".join(rec.interventions),
                    "conditions": ";".join(rec.conditions),
                    "outcome_class": outcome_cls,
                    "raw_payload_sha256": rec.raw_payload_sha256,
                }
            )

        prisma_counts = {
            "step1_identified": step1_identified,
            "step2_interventional_phase_eligible": len(step2_records),
            "step3_ontology_clean_included": len(step3_clean),
            "step4_scientific_failures": sum(
                1 for r in step3_clean if r["outcome_class"] == "SCIENTIFIC_EFFICACY_OR_SAFETY_FAILURE"
            ),
            "step4_administrative_stops": sum(
                1 for r in step3_clean if r["outcome_class"] == "ADMINISTRATIVE_OR_OPERATIONAL_STOP"
            ),
        }
        return pd.DataFrame(step3_clean), prisma_counts
