from __future__ import annotations

import hashlib
import json
from enum import Enum
from pathlib import Path
from typing import Dict, List, Literal, Optional
from pydantic import BaseModel, Field, field_validator


class WorkspaceMode(str, Enum):
    PUBLIC_BENCHMARK = "public_benchmark"
    LAB_PRIVATE = "lab_private"


class ClinicalQuerySpec(BaseModel):
    condition_terms: List[str] = Field(..., min_length=1)
    exclusion_terms: List[str] = Field(default_factory=list)
    stages: List[str] = Field(..., min_length=2)
    min_phase: str = "PHASE2"


class OrganPuritySpec(BaseModel):
    ambient_contamination_genes: List[str] = Field(default_factory=list)
    dissociation_stress_genes: List[str] = Field(
        default_factory=lambda: ["FOS", "JUN", "HSPA1A", "HSPA1B", "DNAJB1", "EGR1"]
    )
    vital_organ_toxicity_genes: List[str] = Field(
        default_factory=lambda: ["KCNH2", "MYH6", "SCN5A", "CYP3A4"]
    )


class BarrierDeliverySpec(BaseModel):
    barrier_type: str
    stage_permeability: Dict[str, List[float]]

    @field_validator("stage_permeability")
    @classmethod
    def _validate_probabilities(cls, v: Dict[str, List[float]]) -> Dict[str, List[float]]:
        for mod, curve in v.items():
            if not curve:
                raise ValueError(f"Empty permeability curve for modality {mod}")
            for val in curve:
                if not (0.0 <= val <= 1.0):
                    raise ValueError(f"Permeability {val} for {mod} must lie in [0, 1]")
        return v


class CrossOrganTransferMap(BaseModel):
    source_organ: str
    target_organ: str
    compartment_alignment: Dict[str, str]
    gene_ortholog_or_functional_swap: Dict[str, str] = Field(default_factory=dict)


class BenchmarkPositiveNegativeAnchors(BaseModel):
    positive_approved_or_phase3_wins: List[str] = Field(default_factory=list)
    negative_scientific_failures: List[str] = Field(default_factory=list)


class DiseasePack(BaseModel):
    pack_id: str
    display_name: str
    archetype: str
    organ_system: str
    time_scale: str
    clinical_query_spec: ClinicalQuerySpec
    organ_purity_spec: OrganPuritySpec
    barrier_delivery_spec: BarrierDeliverySpec
    cell_ontology_map: Dict[str, List[str]]
    cross_organ_transfer: Optional[CrossOrganTransferMap] = None
    clinical_anchors: Optional[BenchmarkPositiveNegativeAnchors] = None


class TrialRecord(BaseModel):
    nct_id: str
    official_title: str
    brief_title: str = ""
    overall_status: str
    why_stopped: str = ""
    study_type: str = "INTERVENTIONAL"
    phases: List[str] = Field(default_factory=list)
    enrollment: int = 0
    interventions: List[str] = Field(default_factory=list)
    conditions: List[str] = Field(default_factory=list)
    raw_payload_sha256: str

    @classmethod
    def from_ctg_api_v2(cls, raw_study: dict) -> "TrialRecord":
        proto = raw_study.get("protocolSection", {})
        ident = proto.get("identificationModule", {})
        status = proto.get("statusModule", {})
        design = proto.get("designModule", {})
        arms = proto.get("armsInterventionsModule", {})
        conds = proto.get("conditionsModule", {})

        nct_id = ident.get("nctId", "")
        brief = ident.get("briefTitle", "")
        official = ident.get("officialTitle") or brief

        interventions = [
            f"{item.get('type', 'OTHER')}:{item.get('name', '')}"
            for item in arms.get("interventions", [])
            if item.get("name")
        ]
        enrollment = int((design.get("enrollmentInfo") or {}).get("count") or 0)
        canonical_bytes = json.dumps(raw_study, sort_keys=True, ensure_ascii=False).encode("utf-8")
        digest = hashlib.sha256(canonical_bytes).hexdigest()

        return cls(
            nct_id=nct_id,
            official_title=official,
            brief_title=brief,
            overall_status=status.get("overallStatus", "UNKNOWN"),
            why_stopped=status.get("whyStopped", "") or "",
            study_type=design.get("studyType", "UNKNOWN"),
            phases=design.get("phases", []) or [],
            enrollment=enrollment,
            interventions=interventions,
            conditions=conds.get("conditions", []) or [],
            raw_payload_sha256=digest,
        )


class MatrixProvenance(BaseModel):
    accession: str
    file_path: str
    source_matrix_sha256: str
    organism: str
    tissue: str
    condition_summary: str
    n_donors: int = Field(..., ge=1)
    n_cells_or_samples: int = Field(..., ge=1)
    derivation_method: str = "raw_matrix_direct_computation"

    @classmethod
    def from_file(
        cls,
        accession: str,
        file_path: Path | str,
        organism: str,
        tissue: str,
        condition_summary: str,
        n_donors: int,
        n_cells_or_samples: int,
        derivation_method: str = "raw_matrix_direct_computation",
    ) -> "MatrixProvenance":
        p = Path(file_path)
        h = hashlib.sha256()
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return cls(
            accession=accession,
            file_path=str(p),
            source_matrix_sha256=h.hexdigest(),
            organism=organism,
            tissue=tissue,
            condition_summary=condition_summary,
            n_donors=n_donors,
            n_cells_or_samples=n_cells_or_samples,
            derivation_method=derivation_method,
        )


class CompartmentRule(BaseModel):
    name: str
    baseline_production: float = Field(0.02, ge=0.0, le=1.0)
    decay_rate: float = Field(0.25, gt=0.0, le=5.0)
    noise_sigma: float = Field(0.04, ge=0.0, le=0.5)
    marker_genes: List[str] = Field(default_factory=list)


class DirectedEdgeRule(BaseModel):
    source: str
    target: str
    interaction_type: Literal["activation", "inhibition"] = "activation"
    weight_v: float = Field(..., gt=0.0, le=5.0)
    half_max_k: float = Field(0.40, gt=0.01, le=2.0)
    hill_n: float = Field(2.0, ge=1.0, le=6.0)
    ligand_gene: str = ""
    receptor_gene: str = ""
    effector_tf: str = ""
    partial_corr_gain: float = 0.0
    spatial_proximity_weight: float = 1.0


class NicheGrammar(BaseModel):
    grammar_id: str
    source_pack_id: str
    organ_context: str
    compartments: List[CompartmentRule] = Field(..., min_length=2)
    edges: List[DirectedEdgeRule] = Field(default_factory=list)
    remission_compartment: str
    remission_threshold: float = Field(0.30, gt=0.0, lt=1.0)
    bic_score: Optional[float] = None
