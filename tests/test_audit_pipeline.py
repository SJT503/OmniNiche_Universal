import numpy as np
import pandas as pd
import pytest

from omniniche.disease_packs.loader import load_disease_pack
from omniniche.audit.clinical_trialist_v2 import ClinicalTrialistV2
from omniniche.audit.chronotau_tracker import ChronoTauTracker
from omniniche.audit.sixgate_tournament import SixGateTournamentReviewer


def _make_raw_study(nct_id: str, title: str, status: str, why: str, conds: list[str], phase: str = "PHASE3", study_type: str = "INTERVENTIONAL", drug: str = "TestMab"):
    return {
        "protocolSection": {
            "identificationModule": {"nctId": nct_id, "officialTitle": title, "briefTitle": title},
            "statusModule": {"overallStatus": status, "whyStopped": why},
            "designModule": {"studyType": study_type, "phases": [phase], "enrollmentInfo": {"count": 300}},
            "armsInterventionsModule": {"interventions": [{"type": "DRUG", "name": drug}]},
            "conditionsModule": {"conditions": conds},
        }
    }


def test_prisma_funnel_and_stop_reason_deconvolution():
    pack = load_disease_pack("sepsis_ards")
    raw_studies = [
        _make_raw_study("NCT00000001", "Eritoran Phase 3 in Severe Sepsis", "TERMINATED", "Lack of efficacy at interim futility analysis", ["Severe Sepsis"]),
        _make_raw_study("NCT00000002", "Anti-IL6 Trial in ARDS Paused by Pandemic", "TERMINATED", "Stopped early due to COVID-19 pandemic accrual halt", ["Acute Respiratory Distress Syndrome"]),
        _make_raw_study("NCT00000003", "Completed Hydrocortisone Regimen in Septic Shock", "COMPLETED", "", ["Septic Shock"]),
        _make_raw_study("NCT00000004", "Observational Biomarker Registry in Sepsis", "COMPLETED", "", ["Sepsis"], study_type="OBSERVATIONAL"),
        _make_raw_study("NCT00000005", "Off-target Traumatic Brain Injury Trial", "COMPLETED", "", ["Sepsis", "Traumatic Brain Injury"]),
    ]

    agent = ClinicalTrialistV2(pack)
    df, prisma_counts = agent.audit_raw_studies(raw_studies)

    assert prisma_counts["step1_identified"] == 5
    assert prisma_counts["step2_interventional_phase_eligible"] == 4
    assert prisma_counts["step3_ontology_clean_included"] == 3
    assert len(df) == 3

    outcomes = dict(zip(df["nct_id"], df["outcome_class"]))
    assert outcomes["NCT00000001"] == "SCIENTIFIC_EFFICACY_OR_SAFETY_FAILURE"
    assert outcomes["NCT00000002"] == "ADMINISTRATIVE_OR_OPERATIONAL_STOP"
    assert outcomes["NCT00000003"] == "COMPLETED_ACTIVE_OR_SURVIVED"


def test_chronotau_detects_delocalizing_decoys_and_enforces_donor_pseudobulk():
    rng = np.random.default_rng(42)
    compartments = ["Injured_Epithelium", "Profibrotic_Macrophage", "Myofibroblast_ECM", "Microvascular_Endothelium"]
    stages = ["t0", "t1", "t2", "t3"]

    # Build stage x compartment mean expression matrices for 4 genes:
    # 1. POSTN: persistently specific to Myofibroblast_ECM across all stages (persistent elite target)
    # 2. DECOY1: high tau at t0 (specific to Myofibroblast), but spreads pan-cellularly at t1-t3 (delocalizing decoy!)
    # 3. UMOD: ambient kidney contamination gene
    # 4. KCNH2: vital cardiac toxicity gene
    stage_comp_profiles = {
        "t0": pd.DataFrame(
            {
                "POSTN": [0.1, 0.1, 4.5, 0.1],
                "DECOY1": [0.1, 0.1, 4.2, 0.1],
                "UMOD": [6.0, 2.5, 2.2, 2.0],
                "KCNH2": [1.5, 1.4, 1.6, 1.5],
            },
            index=compartments,
        ),
        "t1": pd.DataFrame(
            {
                "POSTN": [0.2, 0.3, 6.2, 0.2],
                "DECOY1": [3.8, 4.1, 4.3, 3.9],
                "UMOD": [5.8, 2.9, 2.6, 2.4],
                "KCNH2": [1.5, 1.5, 1.5, 1.5],
            },
            index=compartments,
        ),
        "t2": pd.DataFrame(
            {
                "POSTN": [0.2, 0.4, 7.8, 0.2],
                "DECOY1": [4.0, 4.2, 4.4, 4.1],
                "UMOD": [5.5, 3.1, 2.8, 2.5],
                "KCNH2": [1.6, 1.5, 1.5, 1.4],
            },
            index=compartments,
        ),
        "t3": pd.DataFrame(
            {
                "POSTN": [0.2, 0.4, 8.5, 0.2],
                "DECOY1": [4.1, 4.0, 4.2, 4.0],
                "UMOD": [5.2, 3.0, 2.9, 2.7],
                "KCNH2": [1.5, 1.6, 1.5, 1.5],
            },
            index=compartments,
        ),
    }

    tracker = ChronoTauTracker()
    traj_df = tracker.compute_dynamic_tau_trajectories(stage_comp_profiles)

    postn_row = traj_df.loc["POSTN"]
    decoy_row = traj_df.loc["DECOY1"]
    assert postn_row["tau_t0"] > 0.85
    assert postn_row["min_disease_tau"] > 0.85
    assert postn_row["delta_tau_drift"] > -0.08

    assert decoy_row["tau_t0"] > 0.85
    assert decoy_row["min_disease_tau"] < 0.20
    assert decoy_row["delta_tau_drift"] < -0.65

    # Bimodality coefficient check on bimodal on/off single-cell vector
    bimodal_vec = np.concatenate([np.zeros(200), rng.normal(6.0, 0.3, size=200)])
    bc = tracker.compute_sarle_bimodality(bimodal_vec)
    assert bc > 0.555

    # Donor-level pseudobulk: must reject n_donors < 2 per group (anti-pseudoreplication gate!)
    single_donor_cells = pd.DataFrame(
        {
            "donor_id": ["D1"] * 50 + ["D2"] * 50,
            "group": ["control"] * 50 + ["disease"] * 50,
            "POSTN": np.concatenate([rng.poisson(1, 50), rng.poisson(8, 50)]),
        }
    )
    with pytest.raises(ValueError, match="pseudoreplication"):
        tracker.compute_donor_pseudobulk_de(single_donor_cells, gene_cols=["POSTN"])

    # Valid multi-donor pseudobulk (4 control donors vs 4 disease donors) passes cleanly
    rows = []
    for d_idx in range(4):
        for _ in range(30):
            rows.append({"donor_id": f"C{d_idx}", "group": "control", "POSTN": rng.normal(1.0, 0.2), "DECOY1": rng.normal(1.0, 0.2), "UMOD": rng.normal(4.0, 0.3), "KCNH2": rng.normal(1.5, 0.2)})
    for d_idx in range(4):
        for _ in range(30):
            rows.append({"donor_id": f"P{d_idx}", "group": "disease", "POSTN": rng.normal(6.5, 0.3), "DECOY1": rng.normal(3.8, 0.3), "UMOD": rng.normal(4.1, 0.3), "KCNH2": rng.normal(1.5, 0.2)})
    multi_donor_df = pd.DataFrame(rows)
    pb_df = tracker.compute_donor_pseudobulk_de(multi_donor_df, gene_cols=["POSTN", "DECOY1", "UMOD", "KCNH2"])
    assert pb_df.loc["POSTN", "log2fc"] > 2.0
    assert pb_df.loc["POSTN", "pvalue"] < 0.01

    # Six-Gate Reviewer + BTL/Elo Tournament
    pack = load_disease_pack("fibrosis_ckm_ipf")
    candidate_table = traj_df.join(pb_df)
    candidate_table["bimodality_bc"] = [0.68, 0.42, 0.48, 0.35]
    candidate_table["preferred_modality"] = ["small_molecule", "mAb", "mAb", "small_molecule"]
    candidate_table["target_stage_progress"] = [0.75, 0.75, 0.75, 0.75]
    candidate_table["external_cohort_concordance"] = [0.89, 0.30, 0.25, 0.15]

    reviewer = SixGateTournamentReviewer(pack)
    ranked_df, ablation_df = reviewer.audit_and_rank(candidate_table)

    assert ranked_df.index[0] == "POSTN"
    assert bool(ranked_df.loc["POSTN", "all_six_gates_pass"]) is True
    assert bool(ranked_df.loc["UMOD", "G1_organ_purity"]) is False
    assert bool(ranked_df.loc["DECOY1", "G3_dynamic_tau_retention"]) is False
    assert bool(ranked_df.loc["KCNH2", "G5_vital_organ_safety"]) is False
    assert ranked_df.loc["POSTN", "elo_score"] > ranked_df.loc["DECOY1", "elo_score"]
    assert len(ablation_df) == 6
