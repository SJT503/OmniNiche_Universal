from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from omniniche.disease_packs.loader import load_disease_pack, transfer_grammar_cross_organ
from omniniche.audit.clinical_trialist_v2 import ClinicalTrialistV2
from omniniche.benchmark.calibrated_arena import CalibratedBenchmarkArena
from omniniche.viz.nature_methods_plotter import NatureMethodsPlotter
from omniniche.cli import main as cli_main
from tests.test_langevin_and_identifiability import _build_bistable_3comp_grammar


def _build_realistic_pan_disease_intervention_benchmark() -> pd.DataFrame:
    """Construct a structured benchmark table of 24 interventions spanning 4 translational regimes:
    - True wins: high dynamic tau, high GRN eta, accessible barrier B(m,t), dosed before lock-in OR dual-unlocking combo
    - Delocalizing decoys: high static healthy tau_t0 (fools B1!), high static LR (fools B4!), but low disease tau
    - Barrier-excluded late mAbs: high GRN eta in vitro (fools B3!), but dense stroma B(mAb, late)=0.22 blocks delivery
    - Post-lockin monotherapies: good target, but dosed late (t=48h > t50) as single agent without feedback unlocking
    """
    rows = [
        # True clinical wins (y_true = 1)
        {"trial_id": "W1", "tau_t0": 0.88, "min_disease_tau": 0.86, "grn_eta": 0.84, "lr_centrality": 0.78, "barrier_access": 0.90, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 1},
        {"trial_id": "W2", "tau_t0": 0.82, "min_disease_tau": 0.80, "grn_eta": 0.79, "lr_centrality": 0.72, "barrier_access": 0.88, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 1},
        {"trial_id": "W3", "tau_t0": 0.76, "min_disease_tau": 0.78, "grn_eta": 0.82, "lr_centrality": 0.68, "barrier_access": 0.84, "dosed_before_t50": 0, "is_unlock_combo": 1, "y_clinical_success": 1},
        {"trial_id": "W4", "tau_t0": 0.90, "min_disease_tau": 0.85, "grn_eta": 0.88, "lr_centrality": 0.75, "barrier_access": 0.78, "dosed_before_t50": 0, "is_unlock_combo": 1, "y_clinical_success": 1},
        {"trial_id": "W5", "tau_t0": 0.79, "min_disease_tau": 0.76, "grn_eta": 0.75, "lr_centrality": 0.70, "barrier_access": 0.86, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 1},
        {"trial_id": "W6", "tau_t0": 0.84, "min_disease_tau": 0.81, "grn_eta": 0.80, "lr_centrality": 0.74, "barrier_access": 0.82, "dosed_before_t50": 0, "is_unlock_combo": 1, "y_clinical_success": 1},
        {"trial_id": "W7", "tau_t0": 0.86, "min_disease_tau": 0.83, "grn_eta": 0.77, "lr_centrality": 0.65, "barrier_access": 0.91, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 1},
        {"trial_id": "W8", "tau_t0": 0.75, "min_disease_tau": 0.74, "grn_eta": 0.85, "lr_centrality": 0.69, "barrier_access": 0.80, "dosed_before_t50": 0, "is_unlock_combo": 1, "y_clinical_success": 1},
        # Failure mode A: Delocalizing decoys (fools B1 static tau_t0 & B4 static LR!)
        {"trial_id": "F_DELOC1", "tau_t0": 0.92, "min_disease_tau": 0.22, "grn_eta": 0.55, "lr_centrality": 0.88, "barrier_access": 0.85, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_DELOC2", "tau_t0": 0.89, "min_disease_tau": 0.25, "grn_eta": 0.60, "lr_centrality": 0.85, "barrier_access": 0.82, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_DELOC3", "tau_t0": 0.91, "min_disease_tau": 0.19, "grn_eta": 0.52, "lr_centrality": 0.84, "barrier_access": 0.88, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_DELOC4", "tau_t0": 0.87, "min_disease_tau": 0.28, "grn_eta": 0.58, "lr_centrality": 0.82, "barrier_access": 0.80, "dosed_before_t50": 1, "is_unlock_combo": 0, "y_clinical_success": 0},
        # Failure mode B: Barrier-excluded late mAbs (fools B3 in-vitro GRN eta!)
        {"trial_id": "F_BAR1", "tau_t0": 0.85, "min_disease_tau": 0.80, "grn_eta": 0.91, "lr_centrality": 0.81, "barrier_access": 0.22, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_BAR2", "tau_t0": 0.83, "min_disease_tau": 0.78, "grn_eta": 0.89, "lr_centrality": 0.79, "barrier_access": 0.24, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_BAR3", "tau_t0": 0.81, "min_disease_tau": 0.77, "grn_eta": 0.87, "lr_centrality": 0.77, "barrier_access": 0.20, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_BAR4", "tau_t0": 0.84, "min_disease_tau": 0.79, "grn_eta": 0.90, "lr_centrality": 0.80, "barrier_access": 0.23, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        # Failure mode C: Post-lockin monotherapies dosed after t50 cliff (fools static B1, B3, B4!)
        {"trial_id": "F_LATE1", "tau_t0": 0.84, "min_disease_tau": 0.80, "grn_eta": 0.82, "lr_centrality": 0.83, "barrier_access": 0.84, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_LATE2", "tau_t0": 0.82, "min_disease_tau": 0.79, "grn_eta": 0.80, "lr_centrality": 0.80, "barrier_access": 0.82, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_LATE3", "tau_t0": 0.86, "min_disease_tau": 0.82, "grn_eta": 0.84, "lr_centrality": 0.85, "barrier_access": 0.85, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
        {"trial_id": "F_LATE4", "tau_t0": 0.80, "min_disease_tau": 0.76, "grn_eta": 0.78, "lr_centrality": 0.78, "barrier_access": 0.81, "dosed_before_t50": 0, "is_unlock_combo": 0, "y_clinical_success": 0},
    ]
    return pd.DataFrame(rows)


def test_five_baseline_calibrated_arena_and_cross_organ_zero_shot_evaluation():
    bench_df = _build_realistic_pan_disease_intervention_benchmark()
    arena = CalibratedBenchmarkArena()
    summary_df = arena.evaluate_five_baselines(bench_df)

    expected_models = {
        "B1_Static_Healthy_Tau",
        "B2_Uninformative_Linear_Mean",
        "B3_SingleCell_GRN_Only",
        "B4_Static_LR_Centrality",
        "OmniNiche_v2_Full",
    }
    assert set(summary_df["model_name"]) == expected_models

    omni_row = summary_df.set_index("model_name").loc["OmniNiche_v2_Full"]
    for b_name in ("B1_Static_Healthy_Tau", "B2_Uninformative_Linear_Mean", "B3_SingleCell_GRN_Only", "B4_Static_LR_Centrality"):
        b_row = summary_df.set_index("model_name").loc[b_name]
        assert omni_row["auroc"] > b_row["auroc"]
        assert omni_row["brier_score"] < b_row["brier_score"]
    assert omni_row["calibrated_delta_auroc_vs_null"] > 0.30

    # Cross-organ zero-shot transfer evaluation: Kidney AKI-to-CKD grammar -> Lung IPF interventions
    fib_pack = load_disease_pack("fibrosis_ckm_ipf")
    kidney_grammar = _build_bistable_3comp_grammar()
    ipf_grammar = transfer_grammar_cross_organ(kidney_grammar, fib_pack.cross_organ_transfer)
    transfer_res = arena.evaluate_cross_organ_fibrosis_transfer(ipf_grammar, fib_pack)

    assert transfer_res["transferred_organ"] == "Lung_Alveolar_Interstitium"
    # Multi-pathway oral small-molecule antifibrotic combo outperforms late dense-ECM excluded mAb monotherapy
    assert transfer_res["dual_small_molecule_remission"] > transfer_res["late_anti_ecm_mab_monotherapy_remission"] + 0.40


def test_nature_methods_plotter_exports_clean_vector_and_raster_dashboard(tmp_path: Path):
    bench_df = _build_realistic_pan_disease_intervention_benchmark()
    arena = CalibratedBenchmarkArena()
    summary_df = arena.evaluate_five_baselines(bench_df)

    plotter = NatureMethodsPlotter()
    out_prefix = tmp_path / "Fig_OmniNiche_v2_Benchmark_Dashboard"
    paths = plotter.render_four_panel_dashboard(
        arena_summary_df=summary_df,
        grammar=_build_bistable_3comp_grammar(),
        pack=load_disease_pack("fibrosis_ckm_ipf"),
        output_prefix=out_prefix,
    )
    for fmt in ("pdf", "svg", "png"):
        p = Path(paths[fmt])
        assert p.exists()
        assert p.stat().st_size > 15_000


def test_cli_commands_and_live_clinicaltrials_gov_api_v2_smoke(tmp_path: Path):
    assert cli_main(["list-packs"]) == 0

    custom_yaml = tmp_path / "pulmonary_arterial_hypertension.yaml"
    rc = cli_main(
        [
            "onboard",
            "--disease",
            "Pulmonary Arterial Hypertension",
            "--organ",
            "Lung_Heart",
            "--output",
            str(custom_yaml),
            "--offline",
        ]
    )
    assert rc == 0
    assert custom_yaml.exists()

    # Live ClinicalTrials.gov REST API v2 smoke check on cardiometabolic_mash pack (10 real live studies)
    mash_pack = load_disease_pack("cardiometabolic_mash")
    trialist = ClinicalTrialistV2(mash_pack)
    live_raw = trialist.fetch_live_studies(max_pages_per_term=1, page_size=8)
    assert len(live_raw) >= 5
    audited_df, prisma = trialist.audit_raw_studies(live_raw)
    assert prisma["step1_identified"] == len(live_raw)
    assert "raw_payload_sha256" in audited_df.columns or len(audited_df) == 0
