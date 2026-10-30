from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional
import yaml

from omniniche.core.integrity_firewall import IntegrityFirewall
from omniniche.core.schema import (
    BarrierDeliverySpec,
    ClinicalQuerySpec,
    DiseasePack,
    OrganPuritySpec,
    WorkspaceMode,
)


ORGAN_AMBIENT_CONTAMINATION_KB: Dict[str, List[str]] = {
    "kidney": ["UMOD", "SLC12A1", "SLC34A1", "MIOX", "ALB"],
    "liver": ["ALB", "APOA1", "APOA2", "HP", "TTR", "CYP2E1"],
    "lung": ["SFTPC", "SFTPA1", "SFTPB", "SCGB1A1", "HBA1", "HBB"],
    "blood": ["HBA1", "HBB", "HBA2", "PPBP", "PF4"],
    "immune": ["HBA1", "HBB", "PPBP", "PF4", "IGKC"],
    "gut": ["FABP1", "MUC2", "TFF3", "PIGR", "IGHA1"],
    "tumor": ["KRT8", "KRT18", "KRT19", "EPCAM"],
    "heart": ["MYH6", "MYL2", "TNNT2", "MB", "HBA1"],
    "brain": ["TTR", "FOLR1", "PRLR", "CCDC103", "CLIC6"],
    "default": ["HBA1", "HBB", "PPBP", "ALB", "MALAT1"],
}


class AutoPackAgent:
    """Zero-config onboarding agent that drafts and schema-validates a complete DiseasePack YAML
    for any user-requested disease in <30 seconds (optionally enriching synonyms via EMBL-EBI OLS).
    """

    @staticmethod
    def _slugify(name: str) -> str:
        slug = re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")
        return slug or "custom_disease"

    def _infer_ambient_genes(self, organ_hint: str, disease_name: str) -> List[str]:
        combined = f"{organ_hint} {disease_name}".lower()
        genes: List[str] = []
        for key, gset in ORGAN_AMBIENT_CONTAMINATION_KB.items():
            if key != "default" and key in combined:
                for g in gset:
                    if g not in genes:
                        genes.append(g)
        if len(genes) < 3:
            for g in ORGAN_AMBIENT_CONTAMINATION_KB["default"]:
                if g not in genes:
                    genes.append(g)
        return genes

    def _fetch_ols_synonyms(self, disease_name: str) -> List[str]:
        q = urllib.parse.quote(disease_name)
        url = f"https://www.ebi.ac.uk/ols4/api/search?q={q}&ontology=mondo&rows=3"
        try:
            with urllib.request.urlopen(url, timeout=6) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            docs = payload.get("response", {}).get("docs", [])
            terms = [disease_name]
            for d in docs:
                lbl = d.get("label")
                if lbl and lbl not in terms:
                    terms.append(lbl)
            return terms[:5]
        except Exception:
            return [disease_name]

    def generate_disease_pack(
        self,
        disease_name: str,
        organ_hint: str = "Systemic_Tissue",
        mode: WorkspaceMode = WorkspaceMode.LAB_PRIVATE,
        output_yaml_path: Optional[Path | str] = None,
        use_network: bool = True,
    ) -> DiseasePack:
        pack_id = self._slugify(disease_name)
        cond_terms = self._fetch_ols_synonyms(disease_name) if use_network else [disease_name]
        ambient_genes = self._infer_ambient_genes(organ_hint, disease_name)

        is_brain = "brain" in f"{organ_hint} {disease_name}".lower() or "neuro" in f"{organ_hint} {disease_name}".lower()
        barrier_type = (
            "biphasic_blood_brain_barrier"
            if is_brain
            else "microvascular_and_stromal_remodeling_barrier"
        )
        mab_curve = [0.15, 0.58, 0.35, 0.14] if is_brain else [0.70, 0.76, 0.52, 0.30]

        pack = DiseasePack(
            pack_id=pack_id,
            display_name=f"{disease_name} ({mode.value})",
            archetype="auto_onboarded_multicellular_remodeling",
            organ_system=organ_hint,
            time_scale="days_to_months",
            clinical_query_spec=ClinicalQuerySpec(
                condition_terms=cond_terms,
                exclusion_terms=[] if is_brain else ["Traumatic Brain Injury"],
                stages=[
                    "t0_homeostatic_baseline",
                    "t1_acute_initiation",
                    "t2_subacute_amplification",
                    "t3_chronic_remodeling_lockin",
                ],
                min_phase="PHASE2",
            ),
            organ_purity_spec=OrganPuritySpec(
                ambient_contamination_genes=ambient_genes,
            ),
            barrier_delivery_spec=BarrierDeliverySpec(
                barrier_type=barrier_type,
                stage_permeability={
                    "small_molecule": [0.92, 0.90, 0.84, 0.76],
                    "mAb": mab_curve,
                    "ADC": [0.62, 0.68, 0.50, 0.34],
                    "ASO_siRNA": [0.55, 0.60, 0.46, 0.30],
                },
            ),
            cell_ontology_map={
                "Initiator_Immune": ["Monocyte", "Macrophage", "Dendritic Cell"],
                "Amplifier_Stromal": ["Fibroblast", "Myofibroblast", "Pericyte"],
                "Barrier_Vascular": ["Endothelial", "Capillary Endothelial"],
                "Effector_Parenchyma": ["Epithelial", "Parenchymal", "Lymphocyte"],
            },
        )

        if output_yaml_path is not None:
            out_p = Path(output_yaml_path)
            IntegrityFirewall.assert_path_allowed(out_p, mode=mode, is_write=True)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            header = (
                "# PRIVATE_LAB_ONLY: Generated in lab_private workspace mode\n"
                if mode == WorkspaceMode.LAB_PRIVATE
                else "# PUBLIC_BENCHMARK_PACK\n"
            )
            yaml_body = yaml.safe_dump(pack.model_dump(exclude_none=True), sort_keys=False, allow_unicode=True)
            out_p.write_text(header + yaml_body, encoding="utf-8")

        return pack
