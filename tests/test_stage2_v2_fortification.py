"""TDD tests for Stage II v2.0 Statistical, Causal & Biophysical Fortification (Nature Medicine / Science standard)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from omniniche.stage2.v2_statistical_fortification import (
    adjudicate_stoppage_clause_aware,
    bootstrap_classification_metrics,
    build_conserved_fibrosis_grammar_from_organ,
    compute_decision_curve_and_enrichment,
    compute_visium_radial_distance_profile,
    evaluate_six_step_dnacc_ablation,
    exact_permutation_auroc_pvalue,
    fast_delong_paired_test,
)


def test_fast_delong_paired_test_exact_properties() -> None:
    rng = np.random.default_rng(42)
    y = np.array([0] * 120 + [1] * 80, dtype=int)
    p1 = y * 0.55 + rng.uniform(0.0, 0.45, size=len(y))
    p2 = rng.uniform(0.0, 1.0, size=len(y))

    same = fast_delong_paired_test(y, p1, p1)
    assert np.isclose(same["delta_auroc"], 0.0, atol=1e-12)
    assert np.isclose(same["p_value"], 1.0, atol=1e-6)

    diff = fast_delong_paired_test(y, p1, p2)
    assert np.isclose(diff["auroc_1"], roc_auc_score(y, p1), atol=1e-6)
    assert np.isclose(diff["auroc_2"], roc_auc_score(y, p2), atol=1e-6)
    assert diff["delta_auroc"] > 0.20
    assert diff["z_stat"] > 4.0
    assert diff["p_value"] < 1e-4


def test_bootstrap_classification_metrics_stratified_ci_bounds() -> None:
    # Small n=8 cohort (4 pos, 4 neg) must never fail with single-class resample
    y = np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=int)
    p = np.array([0.88, 0.79, 0.65, 0.52, 0.48, 0.25, 0.14, 0.09], dtype=float)
    res = bootstrap_classification_metrics(y, p, n_boot=500, seed=42)
    assert res["auroc_ci_low"] <= res["auroc"] <= res["auroc_ci_high"]
    assert res["auprc_ci_low"] <= res["auprc"] <= res["auprc_ci_high"]
    assert res["brier_ci_low"] <= res["brier_score"] <= res["brier_ci_high"]
    assert res["auroc"] > 0.90


def test_exact_permutation_auroc_pvalue_small_n() -> None:
    y = np.array([1] * 6 + [0] * 8, dtype=int)
    p_good = np.array([0.95, 0.91, 0.88, 0.84, 0.79, 0.75, 0.35, 0.30, 0.25, 0.22, 0.18, 0.15, 0.11, 0.08])
    p_noise = np.array([0.41, 0.52, 0.39, 0.61, 0.44, 0.50, 0.53, 0.48, 0.59, 0.42, 0.51, 0.47, 0.45, 0.55])

    r_good = exact_permutation_auroc_pvalue(y, p_good, y_baseline=p_noise, n_perm=5000, seed=42)
    r_noise = exact_permutation_auroc_pvalue(y, p_noise, n_perm=2000, seed=42)
    assert r_good["p_perm_vs_null"] < 0.005
    assert r_good["p_perm_vs_baseline"] < 0.01
    assert r_noise["p_perm_vs_null"] > 0.25


def test_negation_aware_stoppage_auditor_eliminates_boilerplate_safety_false_positives() -> None:
    assert (
        adjudicate_stoppage_clause_aware(
            "TERMINATED",
            "Sponsor strategic portfolio prioritization; this decision is not related to patient safety or efficacy.",
        )
        == "Administrative_Commercial_Enrollment"
    )
    assert (
        adjudicate_stoppage_clause_aware(
            "TERMINATED",
            "Stopped after interim analysis met futility criteria. No new safety signals were observed.",
        )
        == "Scientific_Efficacy_Futility"
    )
    assert (
        adjudicate_stoppage_clause_aware(
            "TERMINATED",
            "FDA clinical hold following unexpected Grade 4 hepatotoxicity.",
        )
        == "Scientific_Safety_Toxicity"
    )


def test_six_step_orthogonal_dnacc_ablation_and_unseen_target_groupkfold() -> None:
    rng = np.random.default_rng(123)
    n = 300
    genes = [f"GENE_{i:02d}" for i in rng.integers(0, 40, size=n)]
    tau0 = rng.uniform(0.3, 0.95, size=n)
    vb = 1.0 - tau0 + rng.normal(0, 0.05, size=n)
    tau_e = np.clip(tau0 + rng.normal(-0.05, 0.1, size=n), 0.05, 0.99)
    dtau = rng.normal(-0.08, 0.22, size=n)
    dret = np.clip(0.5 * tau_e + 0.5 * (0.5 + dtau), 0.05, 0.99)
    csw = (dtau < -0.18).astype(int)
    st_sh = np.clip(0.25 + 0.35 * csw + rng.normal(0, 0.05, size=n), 0.05, 0.90)
    bar = rng.uniform(0.35, 0.98, size=n)

    logit = -2.0 + 0.3 * tau0 + 2.2 * dtau + 1.8 * dret - 0.9 * csw + 2.4 * bar
    y = (rng.uniform(0, 1, size=n) < (1.0 / (1.0 + np.exp(-logit)))).astype(int)

    df = pd.DataFrame(
        {
            "nct_id": [f"NCT{i:08d}" for i in range(n)],
            "primary_target_gene": genes,
            "n_targets": rng.choice([1, 2], size=n, p=[0.7, 0.3]),
            "tau0_healthy": tau0,
            "vital_bystander_hazard": vb,
            "tau_early": tau_e,
            "delta_tau": dtau,
            "dynamic_retention_score": dret,
            "compartment_switch": csw,
            "stromal_ectopic_share_late": st_sh,
            "barrier_compatibility": bar,
            "clinical_success": y,
        }
    )
    abl_df = evaluate_six_step_dnacc_ablation(df, cohort_name="Synthetic_Test", n_boot=200, seed=42)
    m1_strat = float(abl_df.loc[(abl_df["step_id"] == "M1_Static_Healthy_Tau0") & (abl_df["cv_scheme"] == "Stratified_5Fold"), "auroc"].iloc[0])
    m6_strat = float(abl_df.loc[(abl_df["step_id"] == "M6_Full_DNAcc") & (abl_df["cv_scheme"] == "Stratified_5Fold"), "auroc"].iloc[0])
    m6_group = float(abl_df.loc[(abl_df["step_id"] == "M6_Full_DNAcc") & (abl_df["cv_scheme"] == "UnseenTarget_Group5Fold"), "auroc"].iloc[0])
    assert m6_strat > m1_strat + 0.10
    assert m6_group > 0.65


def test_decision_curve_analysis_and_trial_enrichment_bounds() -> None:
    rng = np.random.default_rng(99)
    y = np.array([1] * 180 + [0] * 60, dtype=int)
    p_dnacc = np.clip(0.55 + 0.30 * y + rng.normal(0, 0.12, size=len(y)), 0.05, 0.98)
    p_tau0 = np.clip(0.72 + rng.normal(0, 0.05, size=len(y)), 0.05, 0.98)

    dca_df, enrich_df = compute_decision_curve_and_enrichment(y, p_dnacc, p_tau0, median_enrollment=250)
    assert len(dca_df) >= 15
    # In actionable decision thresholds around prevalence (0.75), DNAcc net benefit exceeds Static tau0 & Treat-All
    sub = dca_df[(dca_df["threshold_pt"] >= 0.65) & (dca_df["threshold_pt"] <= 0.85)]
    assert float(sub["nb_dnacc"].max()) > float(sub["nb_treat_all"].max())
    assert float(enrich_df["futile_trials_prevented_pct"].max()) > 20.0
    assert float(enrich_df["patients_spared_futile_exposure"].max()) > 1000


def test_visium_continuous_radial_distance_profile_captures_peritumoral_mycaf_ring() -> None:
    xs, ys = np.meshgrid(np.arange(20, dtype=float), np.arange(20, dtype=float))
    coords = np.column_stack([xs.ravel(), ys.ravel()])
    r = np.hypot(coords[:, 0] - 10.0, coords[:, 1] - 10.0)

    scores = pd.DataFrame(
        {
            "Tumor": np.exp(-0.5 * (r / 3.2) ** 2),
            "myCAF": np.exp(-0.5 * ((r - 5.2) / 1.4) ** 2),
            "TAM": 0.9 * np.exp(-0.5 * ((r - 5.2) / 1.5) ** 2),
            "CD8_T": 1.0 / (1.0 + np.exp(-(r - 8.0) * 1.6)),
        }
    )
    prof = compute_visium_radial_distance_profile(coords, scores, max_pitch=6)
    ring_caf = float(prof.loc[prof["pitch_bin"].isin([1, 2]), "myCAF_z"].mean())
    outer_caf = float(prof.loc[prof["pitch_bin"] >= 5, "myCAF_z"].mean())
    core_cd8 = float(prof.loc[prof["pitch_bin"] == 0, "CD8_T_z"].mean())
    outer_cd8 = float(prof.loc[prof["pitch_bin"] >= 5, "CD8_T_z"].mean())
    assert ring_caf > outer_caf + 0.5
    assert outer_cd8 > core_cd8 + 0.5


def test_bidirectional_cross_organ_fibrosis_grammar_calibration() -> None:
    from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver

    lung_donors = pd.DataFrame(
        {
            "Parenchyma_Injury": [0.16, 0.22, 0.50, 0.55, 0.80, 0.87],
            "SPP1_SAMac": [0.14, 0.18, 0.44, 0.48, 0.82, 0.89],
            "Myofibroblast": [0.12, 0.15, 0.42, 0.47, 0.84, 0.92],
        }
    )
    g_lung = build_conserved_fibrosis_grammar_from_organ(lung_donors, organ_label="Human Lung IPF")
    solver = LangevinNicheSolver(g_lung, t_end=96.0, dt=0.5)
    mono = [InterventionRegimen(target_compartment="Parenchyma_Injury", perturbation_eta=0.82, modality="small_molecule")]
    dual = [
        InterventionRegimen(target_compartment="Parenchyma_Injury", perturbation_eta=0.82, modality="small_molecule"),
        InterventionRegimen(target_compartment="Myofibroblast", perturbation_eta=0.80, modality="small_molecule"),
    ]
    r_f2 = solver.simulate_ensemble(s0=np.array([0.51, 0.43, 0.41]), regimens=mono, n_patients=80, seed=1)
    r_f4_mono = solver.simulate_ensemble(s0=np.array([0.84, 0.86, 0.88]), regimens=mono, n_patients=80, seed=2)
    r_f4_dual = solver.simulate_ensemble(s0=np.array([0.84, 0.86, 0.88]), regimens=dual, n_patients=80, seed=3)
    assert float(r_f2["remission_rate"]) > 0.70
    assert float(r_f4_mono["remission_rate"]) < 0.25
    assert float(r_f4_dual["remission_rate"]) > 0.75


def test_profile_likelihood_and_2d_bliss_no_spurious_inactive_synergy() -> None:
    from omniniche.stage2.v2_statistical_fortification import (
        compute_2d_fibrosis_bliss_matrix_and_landscape,
        compute_continuous_profile_likelihood_curves,
    )

    donors = pd.DataFrame(
        {
            "Parenchyma_Injury": [0.16, 0.28, 0.50, 0.46, 0.80, 0.74, 0.35, 0.62],
            "SPP1_SAMac": [0.14, 0.20, 0.44, 0.55, 0.82, 0.69, 0.48, 0.54],
            "Myofibroblast": [0.22, 0.15, 0.42, 0.49, 0.76, 0.88, 0.31, 0.67],
        }
    )
    g = build_conserved_fibrosis_grammar_from_organ(donors, organ_label="Human Liver")
    prof_df = compute_continuous_profile_likelihood_curves(g, n_grid=11)
    assert len(prof_df) == 44
    assert bool(prof_df["identifiable_at_95pct"].all())

    bliss_df, npz_d = compute_2d_fibrosis_bliss_matrix_and_landscape(g)
    assert len(bliss_df) == 36
    assert npz_d["U_2d"].shape == (60, 60)
    # Orthogonal active pair achieves hysteresis escape, while inactive/blocked arms do not spuriously synergize
    active_pair = bliss_df[
        (bliss_df["arm_row"].str.startswith("THRB")) & (bliss_df["arm_col"].str.startswith("TGFBR1"))
    ]["bliss_synergy_excess"].iloc[0]
    blocked_pair = bliss_df[
        (bliss_df["arm_row"].str.startswith("TREM2")) & (bliss_df["arm_col"].str.startswith("CCR2"))
    ]["bliss_synergy_excess"].iloc[0]
    assert float(active_pair) > 0.85
    assert float(blocked_pair) < 0.05


def test_target_cluster_robust_se_and_within_modality_dnacc() -> None:
    from omniniche.stage2.v2_statistical_fortification import (
        evaluate_within_modality_dnacc_stratification,
        fit_multivariable_fixed_effects_logistic,
    )

    rng = np.random.default_rng(7)
    n = 240
    targets = rng.choice([f"TGT_{i}" for i in range(24)], size=n)
    mods = rng.choice(["Small molecule", "Antibody", "Protein"], size=n, p=[0.6, 0.3, 0.1])
    tau0 = rng.uniform(0.3, 0.9, size=n)
    dnacc = rng.uniform(0.15, 0.88, size=n)
    y = (dnacc + rng.normal(0.0, 0.18, size=n) > 0.50).astype(int)
    df = pd.DataFrame(
        {
            "clinical_success": y,
            "tau0_healthy": tau0,
            "drug_modalities": mods,
            "disease_axis": rng.choice(["Solid_Oncology", "Hepatic_MASH_Fibrosis", "Pulmonary_IPF_ARDS"], size=n),
            "start_year": rng.integers(2008, 2025, size=n),
            "sponsor_class": rng.choice(["INDUSTRY", "OTHER"], size=n),
            "n_targets": rng.choice([1, 2], size=n),
            "primary_target_gene": targets,
        }
    )
    fe = fit_multivariable_fixed_effects_logistic(df, dnacc)
    mv_row = fe[fe["model_specification"] == "Multivariable_FixedEffects_DNAcc"].iloc[0]
    assert int(mv_row["n_target_clusters"]) <= 24
    assert float(mv_row["cluster_ci95_low"]) > 1.0
    assert float(mv_row["cluster_p_value"]) < 0.01

    strat = evaluate_within_modality_dnacc_stratification(df, dnacc, tau0, n_boot=100, seed=7)
    sm_row = strat[strat["modality_stratum"] == "Small_Molecule"].iloc[0]
    assert float(sm_row["dnacc_auroc"]) > float(sm_row["static_tau0_auroc"]) + 0.10
