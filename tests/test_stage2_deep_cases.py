"""TDD tests for Stage II Tasks 3 & 4: Visium Spatial Collagen Barrier, Blind ICB Responder Prediction, and Cross-Organ Fibrosis Transfer."""
from __future__ import annotations

import numpy as np
import pandas as pd

from omniniche.stage2.solid_tumor_spatial_icb import (
    compute_visium_spatial_barrier_and_proximity,
    score_baseline_patients_via_langevin_icb,
)


def test_visium_spatial_barrier_detects_mycaf_exclusion_of_cd8_t_cells() -> None:
    rng = np.random.default_rng(7)
    # Create a 20x20 spatial grid (400 spots)
    xs, ys = np.meshgrid(np.arange(20, dtype=float), np.arange(20, dtype=float))
    coords = np.column_stack([xs.ravel(), ys.ravel()])
    dist_center = np.hypot(coords[:, 0] - 10.0, coords[:, 1] - 10.0)

    # Tumor core (r < 5), myCAF + TAM dense ring (4 <= r <= 8), CD8_T excluded to outer periphery (r > 8)
    tumor = np.exp(-0.5 * (dist_center / 3.8) ** 2) + rng.normal(0, 0.03, len(coords))
    mycaf = np.exp(-0.5 * ((dist_center - 6.0) / 1.6) ** 2) + rng.normal(0, 0.03, len(coords))
    tam = 0.85 * mycaf + 0.15 * rng.uniform(0, 1, len(coords))
    cd8 = (1.0 / (1.0 + np.exp(-(dist_center - 8.5) * 1.5))) + rng.normal(0, 0.03, len(coords))

    scores = pd.DataFrame(
        {
            "Tumor": np.clip(tumor, 0, None),
            "TAM": np.clip(tam, 0, None),
            "myCAF": np.clip(mycaf, 0, None),
            "CD8_T": np.clip(cd8, 0, None),
        }
    )
    res = compute_visium_spatial_barrier_and_proximity(coords, scores, k_neighbors=6)
    prox = res["proximity_weights_df"]
    # TAM and myCAF colocalize strongly, whereas Tumor and CD8_T are physically separated by the myCAF ring
    assert prox.loc["TAM", "myCAF"] > prox.loc["Tumor", "CD8_T"] + 0.20
    assert res["stromal_exclusion_index"] > 0.15


def test_score_baseline_patients_via_langevin_icb_separates_low_and_locked_stroma() -> None:
    # 6 Responder-like baseline patients (low myCAF/TAM lock, preserved T effector potential)
    # 6 Non-Responder-like baseline patients (high myCAF + SPP1 TAM lock + deep T exhaustion)
    patients = pd.DataFrame(
        {
            "patient_id": [f"P{i:02d}" for i in range(12)],
            "Tumor": [0.35, 0.40, 0.30, 0.42, 0.38, 0.33, 0.72, 0.78, 0.81, 0.75, 0.84, 0.79],
            "TAM": [0.22, 0.25, 0.18, 0.28, 0.24, 0.20, 0.76, 0.82, 0.85, 0.79, 0.88, 0.80],
            "myCAF": [0.18, 0.22, 0.15, 0.26, 0.20, 0.19, 0.80, 0.85, 0.88, 0.77, 0.90, 0.83],
            "CD8_Exhaustion": [0.24, 0.28, 0.21, 0.30, 0.25, 0.22, 0.74, 0.81, 0.86, 0.78, 0.89, 0.82],
            "responder": [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        }
    ).set_index("patient_id")

    preds = score_baseline_patients_via_langevin_icb(patients)
    assert preds.loc["P00", "simulated_icb_remission_prob"] > preds.loc["P08", "simulated_icb_remission_prob"] + 0.35


def test_conserved_fibrosis_grammar_exhibits_f2_vs_f4_hysteresis_and_combo_rescue() -> None:
    from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver
    from omniniche.stage2.fibrosis_mash_deconvolution import build_conserved_fibrosis_grammar_from_liver

    mock_liver = pd.DataFrame(
        {
            "Parenchyma_Injury": [0.16, 0.34, 0.50, 0.42, 0.80, 0.68, 0.25, 0.62],
            "SPP1_SAMac": [0.14, 0.52, 0.30, 0.66, 0.82, 0.45, 0.40, 0.74],
            "Myofibroblast": [0.12, 0.20, 0.58, 0.36, 0.78, 0.88, 0.52, 0.48],
        }
    )
    grammar = build_conserved_fibrosis_grammar_from_liver(mock_liver)
    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)

    # Upstream parenchymal monotherapy (e.g. Resmetirom in F2/F3 vs Selonsertib in F4 Cirrhosis)
    reg_parenchyma = [InterventionRegimen(target_compartment="Parenchyma_Injury", perturbation_eta=0.82, modality="small_molecule")]
    f2_res = solver.simulate_ensemble(s0=np.array([0.52, 0.43, 0.41]), regimens=reg_parenchyma, n_patients=120, seed=11)
    f4_res = solver.simulate_ensemble(s0=np.array([0.84, 0.86, 0.88]), regimens=reg_parenchyma, n_patients=120, seed=12)

    # Barrier-breaking dual therapy in F4 Cirrhosis (Parenchyma + Myofibroblast small-molecule loop breaker)
    reg_combo = [
        InterventionRegimen(target_compartment="Parenchyma_Injury", perturbation_eta=0.82, modality="small_molecule"),
        InterventionRegimen(target_compartment="Myofibroblast", perturbation_eta=0.80, modality="small_molecule"),
    ]
    f4_combo = solver.simulate_ensemble(s0=np.array([0.84, 0.86, 0.88]), regimens=reg_combo, n_patients=120, seed=13)

    assert float(f2_res["remission_rate"]) > 0.70
    assert float(f4_res["remission_rate"]) < 0.25
    assert float(f4_combo["remission_rate"]) > float(f4_res["remission_rate"]) + 0.50
