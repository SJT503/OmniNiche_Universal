"""TDD tests for Stage II Task 1: Macro-Atlas Non-CNS ClinicalTrials.gov v2 & OpenTargets Drug-Target Engine."""
from __future__ import annotations

import pytest

from omniniche.stage2.macro_clinical_census import (
    classify_trial_stop_deconvolution,
    is_cns_excluded_trial,
    normalize_intervention_drug_name,
    resolve_intervention_targets,
)


def test_cns_exclusion_firewall_strictly_blocks_neurological_and_allows_non_cns() -> None:
    assert is_cns_excluded_trial(["Alzheimer Disease", "Dementia"], "Anti-amyloid mAb Phase 3") is True
    assert is_cns_excluded_trial(["Glioblastoma Multiforme"], "PD-1 blockade in GBM") is True
    assert is_cns_excluded_trial(["Traumatic Brain Injury"], "Neuroimmune modulation after TBI") is True
    assert is_cns_excluded_trial(["Schizophrenia"], "Dopamine receptor modulator") is True

    assert is_cns_excluded_trial(["Non-Small Cell Lung Cancer"], "Pembrolizumab plus chemotherapy") is False
    assert is_cns_excluded_trial(["Non-alcoholic Steatohepatitis", "Liver Fibrosis"], "Resmetirom Phase 3") is False
    assert is_cns_excluded_trial(["Idiopathic Pulmonary Fibrosis"], "Pamrevlumab in IPF") is False


def test_stop_reason_deconvolution_separates_scientific_vs_administrative() -> None:
    assert (
        classify_trial_stop_deconvolution("TERMINATED", "Terminated due to lack of efficacy at interim futility analysis")
        == "Scientific_Efficacy_Futility"
    )
    assert (
        classify_trial_stop_deconvolution("TERMINATED", "Unacceptable hepatotoxicity and serious adverse events")
        == "Scientific_Safety_Toxicity"
    )
    assert (
        classify_trial_stop_deconvolution("WITHDRAWN", "Sponsor strategic business decision and slow patient accrual")
        == "Administrative_Commercial_Enrollment"
    )
    assert classify_trial_stop_deconvolution("COMPLETED", "") == "Completed"


def test_drug_normalization_and_target_resolution() -> None:
    assert normalize_intervention_drug_name("Pembrolizumab 200 mg IV every 3 weeks") == "pembrolizumab"
    assert normalize_intervention_drug_name("Nintedanib esylate 150mg capsule") == "nintedanib"
    assert normalize_intervention_drug_name("Drug: Selonsertib (GS-4997) 18 mg") == "selonsertib"

    lookup = {
        "pembrolizumab": {"targets": ["PDCD1"], "moa": "Programmed cell death protein 1 inhibitor", "max_phase": 4},
        "selonsertib": {"targets": ["MAP3K5"], "moa": "Apoptosis signal-regulating kinase 1 inhibitor", "max_phase": 3},
    }
    res = resolve_intervention_targets(
        ["Pembrolizumab 200 mg IV", "Placebo oral tablet", "Standard of care"],
        lookup,
    )
    assert res["matched_drugs"] == ["pembrolizumab"]
    assert res["target_genes"] == ["PDCD1"]
