import numpy as np
import pandas as pd
import pytest

from omniniche.core.schema import (
    CompartmentRule,
    DirectedEdgeRule,
    NicheGrammar,
)
from omniniche.disease_packs.loader import load_disease_pack
from omniniche.engine.langevin_solver import (
    InterventionRegimen,
    LangevinNicheSolver,
)


def _build_bistable_3comp_grammar() -> NicheGrammar:
    # Realistic 72-hour cellular turnover rates (h^-1): half-life ~ 6-10 hours, lock-in tipping window ~ 12-22 hours
    return NicheGrammar(
        grammar_id="bistable_fibrosis_test",
        source_pack_id="fibrosis_ckm_ipf",
        organ_context="Kidney_Cortex",
        compartments=[
            CompartmentRule(name="Injured_Epithelium", baseline_production=0.0015, decay_rate=0.055, noise_sigma=0.018),
            CompartmentRule(name="Profibrotic_Macrophage", baseline_production=0.0015, decay_rate=0.050, noise_sigma=0.018),
            CompartmentRule(name="Myofibroblast_ECM", baseline_production=0.0010, decay_rate=0.042, noise_sigma=0.016),
        ],
        edges=[
            DirectedEdgeRule(
                source="Injured_Epithelium",
                target="Profibrotic_Macrophage",
                interaction_type="activation",
                weight_v=0.115,
                half_max_k=0.38,
                hill_n=3.0,
            ),
            DirectedEdgeRule(
                source="Profibrotic_Macrophage",
                target="Myofibroblast_ECM",
                interaction_type="activation",
                weight_v=0.115,
                half_max_k=0.40,
                hill_n=3.0,
            ),
            DirectedEdgeRule(
                source="Myofibroblast_ECM",
                target="Injured_Epithelium",
                interaction_type="activation",
                weight_v=0.125,
                half_max_k=0.40,
                hill_n=3.0,
            ),
            DirectedEdgeRule(
                source="Myofibroblast_ECM",
                target="Myofibroblast_ECM",
                interaction_type="activation",
                weight_v=0.088,
                half_max_k=0.52,
                hill_n=4.0,
            ),
        ],
        remission_compartment="Myofibroblast_ECM",
        remission_threshold=0.30,
    )


def test_bistable_langevin_sde_monotherapy_cliff_t50_and_combination_synergy():
    pack = load_disease_pack("fibrosis_ckm_ipf")
    grammar = _build_bistable_3comp_grammar()
    solver = LangevinNicheSolver(grammar=grammar, pack=pack, t_end=96.0, dt=0.5)

    # Acute injury trigger initial state: high primary epithelial injury, quiescent downstream macrophage & ECM at t=0
    s0 = np.array([0.82, 0.08, 0.04])

    untreated = solver.simulate_ensemble(s0=s0, regimens=[], n_patients=200, seed=42)
    assert np.all(untreated["trajectories"] >= 0.0) and np.all(untreated["trajectories"] <= 1.0)
    assert untreated["remission_rate"] < 0.15

    # Early monotherapy targeting Injured_Epithelium at t=1h intercepts before ECM self-lock crosses separatrix
    early_mono = [
        InterventionRegimen(
            target_compartment="Injured_Epithelium",
            perturbation_eta=0.90,
            modality="small_molecule",
            start_time=1.0,
        )
    ]
    res_early = solver.simulate_ensemble(s0=s0, regimens=early_mono, n_patients=200, seed=42)
    assert res_early["remission_rate"] > 0.75

    # Late monotherapy targeting only upstream Injured_Epithelium at t=40h fails due to Myofibroblast_ECM self-lock
    late_mono_a = [
        InterventionRegimen(
            target_compartment="Injured_Epithelium",
            perturbation_eta=0.90,
            modality="small_molecule",
            start_time=40.0,
        )
    ]
    late_mono_b = [
        InterventionRegimen(
            target_compartment="Myofibroblast_ECM",
            perturbation_eta=0.55,
            modality="small_molecule",
            start_time=40.0,
        )
    ]
    res_late_a = solver.simulate_ensemble(s0=s0, regimens=late_mono_a, n_patients=200, seed=42)
    res_late_b = solver.simulate_ensemble(s0=s0, regimens=late_mono_b, n_patients=200, seed=42)
    assert res_late_a["remission_rate"] < 0.25

    # Compute monotherapy window cliff t50
    window_res = solver.compute_monotherapy_window_t50(
        s0=s0,
        target_compartment="Injured_Epithelium",
        perturbation_eta=0.90,
        modality="small_molecule",
        onset_grid=np.linspace(0.0, 48.0, 13),
        n_patients=120,
    )
    assert 3.0 < window_res["t50_window_cliff"] < 38.0

    # Late dual-compartment rescue combination breaks both cross-cell feedback and ECM self-lock (Bliss synergy!)
    synergy = solver.evaluate_bliss_synergy(
        s0=s0,
        regimen_a=late_mono_a[0],
        regimen_b=late_mono_b[0],
        n_patients=200,
        seed=42,
    )
    assert synergy["combo_remission_rate"] > 0.65
    assert synergy["bliss_excess"] > 0.20


def test_waddington_quasi_potential_landscape_and_escape_barrier_drop():
    pack = load_disease_pack("fibrosis_ckm_ipf")
    grammar = _build_bistable_3comp_grammar()
    solver = LangevinNicheSolver(grammar=grammar, pack=pack, t_end=60.0, dt=0.5)

    landscape_untreated = solver.compute_quasi_potential_landscape(
        regimens=[],
        comp_x="Injured_Epithelium",
        comp_y="Myofibroblast_ECM",
        n_grid=25,
        n_trajectories=160,
        seed=7,
    )
    rescue_regimens = [
        InterventionRegimen(target_compartment="Injured_Epithelium", perturbation_eta=0.88, modality="small_molecule", start_time=0.0),
        InterventionRegimen(target_compartment="Myofibroblast_ECM", perturbation_eta=0.80, modality="small_molecule", start_time=0.0),
    ]
    landscape_rescue = solver.compute_quasi_potential_landscape(
        regimens=rescue_regimens,
        comp_x="Injured_Epithelium",
        comp_y="Myofibroblast_ECM",
        n_grid=25,
        n_trajectories=160,
        seed=7,
    )

    assert landscape_untreated["U_grid"].shape == (25, 25)
    assert np.all(np.isfinite(landscape_untreated["U_grid"]))
    # Pathological-to-healthy escape barrier drops under dual rescue perturbation
    assert landscape_rescue["pathological_lock_depth"] < landscape_untreated["pathological_lock_depth"]


def test_profile_likelihood_identifiability_and_loco_cv():
    pack = load_disease_pack("fibrosis_ckm_ipf")
    grammar = _build_bistable_3comp_grammar()
    # Add a disconnected 4th bystander compartment with near-zero state so an edge from it is structurally sloppy/unidentifiable
    grammar_with_sloppy = grammar.model_copy(
        update={
            "compartments": grammar.compartments
            + [CompartmentRule(name="Silent_Bystander", baseline_production=0.0, decay_rate=2.0, noise_sigma=0.005)],
            "edges": grammar.edges
            + [
                DirectedEdgeRule(
                    source="Silent_Bystander",
                    target="Myofibroblast_ECM",
                    interaction_type="activation",
                    weight_v=0.50,
                    half_max_k=0.85,
                    hill_n=4.0,
                )
            ],
        }
    )
    solver = LangevinNicheSolver(grammar=grammar_with_sloppy, pack=pack, t_end=48.0, dt=0.5)
    s0 = np.array([0.80, 0.20, 0.12, 0.001])

    # Generate 3 realistic cohort observation tables at t = [0, 12, 24, 48] with biological measurement noise
    rng = np.random.default_rng(99)
    det_ref = solver.simulate_deterministic(s0=s0, regimens=[])
    obs_times = [0.0, 12.0, 24.0, 48.0]
    cohort_rows = []
    for cid in ("Cohort_A", "Cohort_B", "Cohort_C"):
        for t_pt in obs_times:
            idx = int(np.argmin(np.abs(det_ref["time"] - t_pt)))
            noisy_state = np.clip(det_ref["state"][idx] + rng.normal(0, 0.025, size=4), 0.0, 1.0)
            row = {"cohort_id": cid, "time_hr": t_pt}
            for c_idx, comp in enumerate(grammar_with_sloppy.compartments):
                row[comp.name] = float(noisy_state[c_idx])
            cohort_rows.append(row)
    obs_df = pd.DataFrame(cohort_rows)

    loco_res = solver.run_loco_cv(obs_df, s0=s0)
    assert len(loco_res["per_cohort_rmse"]) == 3
    assert loco_res["mean_out_of_sample_rmse"] < 0.15

    pl_df = solver.profile_likelihood_identifiability(obs_df, s0=s0, obs_sigma=0.04)
    core_edge_row = pl_df[
        (pl_df["source"] == "Injured_Epithelium") & (pl_df["target"] == "Profibrotic_Macrophage")
    ].iloc[0]
    sloppy_edge_row = pl_df[
        (pl_df["source"] == "Silent_Bystander") & (pl_df["target"] == "Myofibroblast_ECM")
    ].iloc[0]

    assert bool(core_edge_row["identifiable_95ci"]) is True
    assert bool(sloppy_edge_row["identifiable_95ci"]) is False
