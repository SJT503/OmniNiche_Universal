"""TDD tests for Stage II Task 2: HPA Pan-Organ Baseline tau_0 and 6-Axis Dynamic Drift Law."""
from __future__ import annotations

import numpy as np
import pandas as pd

from omniniche.stage2.pan_disease_drift_law import (
    compute_dynamic_drift_metrics,
    compute_yanai_tau_matrix,
    evaluate_translational_progression_models,
)


def test_compute_yanai_tau_matrix_exact_bounds() -> None:
    # Gene A: strictly specific to cell type 0 -> tau = 1.0
    # Gene B: uniform across all 4 cell types -> tau = 0.0
    mat = np.array(
        [
            [10.0, 0.0, 0.0, 0.0],
            [5.0, 5.0, 5.0, 5.0],
        ],
        dtype=float,
    )
    tau = compute_yanai_tau_matrix(mat)
    assert np.isclose(tau[0], 1.0, atol=1e-6)
    assert np.isclose(tau[1], 0.0, atol=1e-6)


def test_compute_dynamic_drift_metrics_detects_late_stage_collapse_and_retention() -> None:
    genes = ["GENE_STABLE_HIGH", "GENE_COLLAPSED_LATE", "GENE_UBIQUITOUS"]
    # 3 compartments: [Parenchyma, Macrophage, Myofibroblast]
    early = pd.DataFrame(
        [[0.1, 8.0, 0.2], [0.1, 9.0, 0.2], [4.0, 4.2, 3.9]],
        index=genes,
        columns=["Parenchyma", "Macrophage", "Myofibroblast"],
    )
    late = pd.DataFrame(
        [[0.1, 8.5, 0.2], [3.5, 3.2, 3.6], [4.1, 4.0, 4.2]],
        index=genes,
        columns=["Parenchyma", "Macrophage", "Myofibroblast"],
    )
    df = compute_dynamic_drift_metrics(early, late)
    assert df.loc["GENE_STABLE_HIGH", "tau_late"] > 0.85
    assert df.loc["GENE_COLLAPSED_LATE", "delta_tau"] < -0.45
    assert df.loc["GENE_STABLE_HIGH", "dynamic_retention_score"] > df.loc["GENE_COLLAPSED_LATE", "dynamic_retention_score"]


def test_evaluate_translational_progression_models_outperforms_static_tau0_under_drift() -> None:
    rng = np.random.default_rng(42)
    n = 240
    tau0 = rng.uniform(0.3, 0.95, size=n)
    delta_tau = rng.normal(-0.1, 0.25, size=n)
    barrier_compat = rng.uniform(0.35, 1.0, size=n)
    # True clinical success depends on high tau0 + positive/stable retention (delta_tau) + barrier compatibility
    logit = -1.8 + 1.4 * tau0 + 2.6 * delta_tau + 2.2 * barrier_compat
    prob = 1.0 / (1.0 + np.exp(-logit))
    y = (rng.uniform(0, 1, size=n) < prob).astype(int)

    trials_df = pd.DataFrame(
        {
            "nct_id": [f"NCT{10000000+i}" for i in range(n)],
            "clinical_success": y,
            "tau0_healthy": tau0,
            "tau_early": np.clip(tau0 + rng.normal(0, 0.05, size=n), 0, 1),
            "delta_tau": delta_tau,
            "barrier_compatibility": barrier_compat,
        }
    )
    res = evaluate_translational_progression_models(trials_df)
    assert res["omniniche_dynamic_auroc"] > res["static_tau0_auroc"] + 0.05
    assert res["omniniche_odds_ratio_per_sd"] > 1.3
