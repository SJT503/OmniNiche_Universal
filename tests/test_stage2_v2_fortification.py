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

    # Realistic cross-compartment human donor variation (r_sm ~ 0.44, matching primary human organs)
    lung_donors = pd.DataFrame(
        {
            "Parenchyma_Injury": [0.16, 0.34, 0.50, 0.42, 0.80, 0.68, 0.25, 0.62],
            "SPP1_SAMac": [0.14, 0.52, 0.30, 0.66, 0.82, 0.45, 0.40, 0.74],
            "Myofibroblast": [0.12, 0.20, 0.58, 0.36, 0.78, 0.88, 0.52, 0.48],
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
            "Parenchyma_Injury": [0.16, 0.34, 0.50, 0.42, 0.80, 0.68, 0.25, 0.62],
            "SPP1_SAMac": [0.14, 0.52, 0.30, 0.66, 0.82, 0.45, 0.40, 0.74],
            "Myofibroblast": [0.12, 0.20, 0.58, 0.36, 0.78, 0.88, 0.52, 0.48],
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


def test_t1_noncircular_240_manual_clinical_stoppage_benchmark() -> None:
    import pathlib
    from omniniche.stage2.v2_statistical_fortification import audit_stoppage_classifier_benchmark

    root = pathlib.Path(r"F:\OmniNiche_Universal")
    ledger = root / "data_vault" / "raw_clinical_trials" / "stoppage_240_manual_clinical_adjudication_ledger.csv"
    master = root / "data_vault" / "raw_clinical_trials" / "non_cns_phase234_trials_master.csv"
    assert ledger.exists()
    df_m = pd.read_csv(master, usecols=["nct_id", "disease_axis", "overall_status", "why_stopped"], low_memory=False)
    val_df, cm_df, sample_df = audit_stoppage_classifier_benchmark(df_m, manual_ledger_path=ledger)
    assert len(sample_df) == 240
    assert cm_df.shape == (4, 4)

    cl_col = "auditor_name" if "auditor_name" in val_df.columns else "classifier"
    cl = val_df[val_df[cl_col] == "Clause_Aware_Auditor_v2_2"].iloc[0]
    nv = val_df[val_df[cl_col] == "Naive_Keyword_Baseline"].iloc[0]
    assert float(cl["overall_accuracy"]) >= 0.94
    assert float(cl["cohen_kappa"]) >= 0.92
    assert float(cl["overall_accuracy"]) > float(nv["overall_accuracy"]) + 0.20


def test_t2_intrahepatic_glp1r_absent_and_objective_compartment_weights() -> None:
    from omniniche.engine.grn_perturb import RidgeGRNPerturber
    from omniniche.stage2.fibrosis_mash_deconvolution import (
        FIBROSIS_COMPARTMENT_GENES,
        LANDMARK_FIBROSIS_TRIALS,
        compute_unsigned_grn_leverage,
        derive_objective_compartment_weights,
    )

    # Jara et al. Nat. Med. 2025: absent intrahepatic GLP1R in organ-resident macrophages
    assert "GLP1R" not in FIBROSIS_COMPARTMENT_GENES["SPP1_SAMac"]
    assert len(LANDMARK_FIBROSIS_TRIALS) == 24
    assert all("compartments_hit" not in t for t in LANDMARK_FIBROSIS_TRIALS)

    w_glp = derive_objective_compartment_weights(["GLP1R"])
    assert w_glp["Parenchyma_Injury"] == 1.0
    assert w_glp.get("SPP1_SAMac", 0.0) == 0.0

    # Protective nuclear agonist THRB (negatively correlated with COL1A1) still yields positive unsigned leverage
    expr = pd.DataFrame(
        {
            "THRB": [4.2, 3.8, 2.5, 1.4, 0.9],
            "GDF15": [1.1, 1.5, 2.8, 3.9, 4.5],
            "COL1A1": [0.8, 1.2, 2.6, 4.0, 4.9],
        }
    )
    pert = RidgeGRNPerturber(ridge_lambda=0.5, propagation_steps=2).fit(
        expr, prior_edges=[("THRB", "GDF15"), ("GDF15", "COL1A1")], pathogenic_program_genes=["GDF15", "COL1A1"]
    )
    assert compute_unsigned_grn_leverage(pert, ["THRB"]) > 0.10


def test_t4_opentargets_human_genetics_orthogonality_and_joint_lift() -> None:
    from omniniche.stage2.v2_statistical_fortification import evaluate_human_genetics_vs_dnacc_orthogonality

    rng = np.random.default_rng(2026)
    n = 320
    genes = [f"G_{i:02d}" for i in rng.integers(0, 45, size=n)]
    g_assoc = np.clip(rng.beta(0.7, 1.8, size=n), 0.0, 1.0)
    g_som = np.zeros(n)
    tau0 = rng.uniform(0.3, 0.9, size=n)
    dtau = rng.normal(-0.05, 0.25, size=n)
    bar = rng.uniform(0.35, 0.96, size=n)
    logit = -1.6 + 2.1 * g_assoc + 2.6 * dtau + 2.4 * bar
    y = (rng.uniform(0, 1, size=n) < (1.0 / (1.0 + np.exp(-logit)))).astype(int)

    df = pd.DataFrame(
        {
            "primary_target_gene": genes,
            "ot_genetic_assoc": g_assoc,
            "ot_somatic_mut": g_som,
            "ot_human_genetic_max": g_assoc,
            "tau0_healthy": tau0,
            "vital_bystander_hazard": 1.0 - tau0,
            "tau_early": tau0,
            "delta_tau": dtau,
            "dynamic_retention_score": np.clip(0.5 + dtau, 0.05, 0.95),
            "compartment_switch": (dtau < -0.15).astype(int),
            "stromal_ectopic_share_late": np.clip(0.3 - 0.5 * dtau, 0.05, 0.90),
            "barrier_compatibility": bar,
            "clinical_success": y,
        }
    )
    res = evaluate_human_genetics_vs_dnacc_orthogonality(df, "Synthetic_Genetics", n_boot=100, seed=42)
    st = res[res["evaluation_stratum"] == "Stratified_5Fold"].set_index("model_id")
    assert float(st.loc["G5_Joint_Human_Genetics_Plus_DNAcc", "auroc"]) > float(st.loc["G1_OpenTargets_Human_Genetics_Only", "auroc"]) + 0.04
    sup = res[res["evaluation_stratum"].str.startswith("Within_Genetically_Supported")].iloc[0]
    assert float(sup["delta_auroc_vs_tau0"]) > 0.08


def test_t5_empirical_ecm_density_and_singletarget_group5fold() -> None:
    import pathlib
    from omniniche.stage2.pan_disease_drift_law import compute_empirical_axis_ecm_density_map

    vault_pb = pathlib.Path(r"F:\OmniNiche_Universal\data_vault\processed_pseudobulk")
    ecm_map = compute_empirical_axis_ecm_density_map(vault_pb)
    assert len(ecm_map) >= 6
    assert ecm_map["Pulmonary_IPF_ARDS"] > ecm_map["Infectious_Sepsis"]
    assert all(0.05 <= v <= 0.40 for v in ecm_map.values())


def test_t7_lopo_reference_scaling_zero_test_patient_leakage() -> None:
    import sys

    if r"F:\OmniNiche_Universal\scripts" not in sys.path:
        sys.path.insert(0, r"F:\OmniNiche_Universal\scripts")
    from run_stage2_task4_solid_tumor_icb import _lopo_robust_01, _lopo_zscore

    arr1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 100.0])
    arr2 = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 500.0])
    # Patient index 5's LOPO reference median/MAD depends ONLY on indices 0..4!
    z1 = _lopo_zscore(arr1)
    z2 = _lopo_zscore(arr2)
    # Scaling factor (median=3.0, MAD*1.4826=1.4826) for patient 5 is strictly identical regardless of patient 5's own value
    assert np.isclose((arr1[5] - 3.0) / z1[5], (arr2[5] - 3.0) / z2[5], rtol=1e-10)
    p1 = _lopo_robust_01(arr1, 0.14, 0.88)
    assert np.all((p1 >= 0.14) & (p1 <= 0.88))


def test_v23_continuous_2d_synergy_surface_and_isobologram() -> None:
    import pathlib

    p = pathlib.Path(r"F:\OmniNiche_Universal\outputs\stage2_tables\Table_S4h_Continuous_2D_Dose_Synergy_Surface.csv")
    assert p.exists()
    df = pd.read_csv(p)
    assert len(df) == 169
    mono_mask = (df["eta_parenchyma"] == 0.0) | (df["eta_myofibroblast"] == 0.0)
    assert np.allclose(df.loc[mono_mask, "bliss_synergy_excess"].to_numpy(), 0.0, atol=1e-8)
    assert float(df["bliss_synergy_excess"].max()) > 0.95
    assert int((df["bliss_synergy_excess"] > 0.01).sum()) >= 80


def test_v23_antipd1_monotherapy_74patient_survival_and_roc() -> None:
    import pathlib

    t_dir = pathlib.Path(r"F:\OmniNiche_Universal\outputs\stage2_tables")
    df_surv = pd.read_csv(t_dir / "Table_S5l_AntiPD1_Monotherapy_Survival_74Patients.csv")
    assert len(df_surv) == 74
    assert int(df_surv["event_dead"].sum()) == 43
    assert int(df_surv["stromal_barrier_high"].sum()) == 37
    assert int((df_surv["stromal_barrier_high"] == 0).sum()) == 37

    df_sum = pd.read_csv(t_dir / "Table_S5m_AntiPD1_Survival_Summary.csv").set_index("stratification")
    row = df_sum.loc["Stromal_CAF_Exclusion_Barrier_High_vs_Low"]
    assert float(row["hazard_ratio"]) > 1.80
    assert float(row["hr_ci95_low"]) > 1.00
    assert float(row["logrank_p"]) < 0.05

    df_roc = pd.read_csv(t_dir / "Table_S5n_PreTreatment_107Patients_ROC_Curves.csv")
    assert len(df_roc) == 505
    for _, grp in df_roc.groupby("method"):
        assert len(grp) == 101


def test_v23_pure_vector_fig1_svg_editable_text_and_zero_forbidden_jargon() -> None:
    import pathlib

    fig_dir = pathlib.Path(r"F:\OmniNiche_Universal\outputs\stage2_figures")
    fig1_svg = fig_dir / "Fig1_OmniNiche_Architecture_and_Framework.svg"
    assert fig1_svg.exists()
    svg_text = fig1_svg.read_text(encoding="utf-8")
    assert svg_text.count("<text") >= 30
    assert svg_text.count("<image") <= 1

    svg_files = sorted(fig_dir.glob("*.svg"))
    assert len(svg_files) >= 9
    forbidden_tokens = ("non" + "-cns", "firewall" + "-2", "non" + "-neurological")
    for svg_path in svg_files:
        low = svg_path.read_text(encoding="utf-8").lower()
        for tok in forbidden_tokens:
            assert tok not in low, f"Forbidden jargon '{tok}' found in {svg_path.name}"


def test_v24_peer_review_closure_cohort_stratified_os_and_fig_parity() -> None:
    import pathlib

    t_dir = pathlib.Path(r"F:\OmniNiche_Universal\outputs\stage2_tables")
    fig_dir = pathlib.Path(r"F:\OmniNiche_Universal\outputs\stage2_figures")

    # 1. T11: Cohort-stratified Cox HR & Mantel-Cox log-rank P on n = 74 monotherapy OS
    df_sum = pd.read_csv(t_dir / "Table_S5m_AntiPD1_Survival_Summary.csv").set_index("stratification")
    r_bar = df_sum.loc["Stromal_CAF_Exclusion_Barrier_High_vs_Low"]
    assert float(r_bar["stratified_hazard_ratio"]) > 1.80
    assert float(r_bar["stratified_hr_ci95_low"]) > 1.00
    assert float(r_bar["stratified_logrank_p"]) < 0.05

    # 2. T13a: Fig. 1d Card 2 dynamic CCR2 6.19 -> 1.47 (-76%) parity (no stale CCR2 7.99)
    fig1_svg = (fig_dir / "Fig1_OmniNiche_Architecture_and_Framework.svg").read_text(encoding="utf-8")
    assert "CCR2 6.19" in fig1_svg and "1.47 (-76%)" in fig1_svg
    assert "CCR2 7.99" not in fig1_svg

    # 3. T15a: Fig. 2a exact modality counts (18.1k, 7.3k, 4.4k)
    fig2_svg = (fig_dir / "Fig2_MacroAtlas_41396_Systemic_Organ_ClinicalTrials.svg").read_text(encoding="utf-8")
    assert "18.1k" in fig2_svg and "7.3k" in fig2_svg and "4.4k" in fig2_svg

    # 4. T12 & T13b: Fig. 4a 64-Donor Closure-Free Triad & Fig. 5g 3-arm rescue (+ SPP1/TREM2)
    fig4_svg = (fig_dir / "Fig4_PanOrgan_Fibrosis_ZeroShot_and_Hysteresis.svg").read_text(encoding="utf-8")
    assert "64-Donor Closure-Free Triad" in fig4_svg
    fig5_svg = (fig_dir / "Fig5_SolidTumor_Visium_Exclusion_and_PreICB_Responder.svg").read_text(encoding="utf-8")
    assert "74.6%" in fig5_svg and "TREM2" in fig5_svg


def test_v26_harmonized_ecm_m6_vs_m3_and_manuscript_sync() -> None:
    import json
    import pathlib

    root = pathlib.Path(r"F:\OmniNiche_Universal")
    t_dir = root / "outputs" / "stage2_tables"
    fig_dir = root / "outputs" / "stage2_figures"
    ms_dir = root / "manuscript"

    # 1. T1: Harmonized log2(CPM+1) 10-gene ECM density scale across all 12 axes
    s3m = pd.read_csv(t_dir / "Table_S3m_Empirical_Axis_ECM_Density.csv").set_index("disease_axis")
    assert abs(float(s3m.loc["Solid_Oncology", "empirical_ecm_density"]) - 0.3500) < 1e-4
    assert float(s3m.loc["Pulmonary_IPF_ARDS", "empirical_ecm_density"]) > 0.30
    assert abs(float(s3m.loc["Infectious_Sepsis", "empirical_ecm_density"]) - 0.0700) < 1e-4

    # 2. T2: Pairwise M6 vs M3 & sequential DeLong tests in Table_S3f
    s3f = pd.read_csv(t_dir / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv")
    assert {"delta_auroc_vs_m3", "delong_p_vs_m3", "delta_auroc_vs_prev", "delong_p_vs_prev"}.issubset(s3f.columns)
    mat_m6 = s3f[
        (s3f["cohort"] == "Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials")
        & (s3f["cv_scheme"] == "Stratified_5Fold")
        & (s3f["step_id"] == "M6_Full_DNAcc")
    ].iloc[0]
    assert float(mat_m6["delta_auroc_vs_m3"]) > 0.01
    assert float(mat_m6["delong_p_vs_m3"]) < 1e-10

    # 3. T4/T5/T10: Figure SVG parity for Youden threshold, CD8:myCAF 6th series, and 20,880 QC spots
    fig4_svg = (fig_dir / "Fig4_PanOrgan_Fibrosis_ZeroShot_and_Hysteresis.svg").read_text(encoding="utf-8")
    fig5_svg = (fig_dir / "Fig5_SolidTumor_Visium_Exclusion_and_PreICB_Responder.svg").read_text(encoding="utf-8")
    assert "Youden" in fig4_svg
    assert "CD8:myCAF" in fig5_svg
    assert "20,880" in fig5_svg

    # 4. Manuscript EN/ZH 1:1 line-parity, 61 verified references (PMID: 36934740), and zero forbidden jargon
    en_lines = (ms_dir / "OmniNiche_Flagship_Manuscript_EN.md").read_text(encoding="utf-8").splitlines()
    zh_lines = (ms_dir / "OmniNiche_Flagship_Manuscript_ZH_Preview.md").read_text(encoding="utf-8").splitlines()
    assert len(en_lines) == len(zh_lines) and len(en_lines) >= 335
    refs = json.loads((t_dir / "verified_pubmed_references.json").read_text(encoding="utf-8"))
    assert len(refs) >= 61 and any(str(r.get("pmid")) == "36934740" for r in refs)
    assert "PMID: 36934740" in "\n".join(en_lines) and "PMID: 36934740" in "\n".join(zh_lines)
    forbidden_tokens = ("non" + "-cns", "firewall" + "-2", "non" + "-neurological")
    for doc_txt in ("\n".join(en_lines).lower(), "\n".join(zh_lines).lower()):
        for tok in forbidden_tokens:
            assert tok not in doc_txt


