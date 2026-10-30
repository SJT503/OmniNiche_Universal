from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from omniniche.engine.grn_perturb import RidgeGRNPerturber
from omniniche.engine.autogrammar_compiler import DeconfoundedAutoGrammarCompiler


def test_ridge_grn_perturber_propagates_tf_knockout_and_polypharmacology():
    rng = np.random.default_rng(123)
    n_cells = 400
    genes = ["TGFBR2", "SMAD3", "POSTN", "COL1A1", "ACTA2", "HOMEOS_TF", "HOUSE_TF", "ACTB"]

    tgfbr2 = rng.uniform(0.5, 3.0, size=n_cells)
    homeos = rng.uniform(0.5, 3.0, size=n_cells)
    house = rng.uniform(0.5, 3.0, size=n_cells)

    smad3 = 0.85 * tgfbr2 - 0.55 * homeos + rng.normal(0, 0.15, size=n_cells)
    postn = 0.90 * smad3 + rng.normal(0, 0.15, size=n_cells)
    col1a1 = 0.95 * smad3 + rng.normal(0, 0.15, size=n_cells)
    acta2 = 0.80 * smad3 - 0.40 * homeos + rng.normal(0, 0.15, size=n_cells)
    actb = 0.75 * house + rng.normal(0, 0.15, size=n_cells)

    expr_df = pd.DataFrame(
        np.column_stack([tgfbr2, smad3, postn, col1a1, acta2, homeos, house, actb]),
        columns=genes,
    )

    prior_edges = [
        ("TGFBR2", "SMAD3"),
        ("SMAD3", "POSTN"),
        ("SMAD3", "COL1A1"),
        ("SMAD3", "ACTA2"),
        ("HOMEOS_TF", "SMAD3"),
        ("HOMEOS_TF", "ACTA2"),
        ("HOUSE_TF", "ACTB"),
    ]
    pathogenic_signature = ["POSTN", "COL1A1", "ACTA2"]

    perturber = RidgeGRNPerturber(ridge_lambda=1.0, propagation_steps=3)
    perturber.fit(expr_df, prior_edges=prior_edges, pathogenic_program_genes=pathogenic_signature)

    eta_smad3_ko = perturber.predict_knockout_efficacy("SMAD3")
    eta_upstream_ko = perturber.predict_knockout_efficacy("TGFBR2")
    eta_house_ko = perturber.predict_knockout_efficacy("HOUSE_TF")
    eta_homeos_ko = perturber.predict_knockout_efficacy("HOMEOS_TF")

    # Inhibiting SMAD3 or upstream receptor TGFBR2 strongly suppresses pathogenic program
    assert eta_smad3_ko > 0.50
    assert eta_upstream_ko > 0.35
    # Inhibiting housekeeping TF has negligible effect on pathogenic program
    assert abs(eta_house_ko) < 0.15
    # Knocking out protective HOMEOS_TF worsens pathogenic program (negative suppression efficacy)
    assert eta_homeos_ko < -0.30

    # Polypharmacology drug profile (pIC50 >= 6.0 active)
    eta_drug = perturber.predict_drug_efficacy({"TGFBR2": 8.2, "SMAD3": 7.5})
    assert -1.0 <= eta_drug <= 1.0
    assert eta_drug > eta_upstream_ko * 0.9


def test_deconfounded_autogrammar_eliminates_spurious_confounder_and_recovers_feedback_loop(tmp_path: Path):
    rng = np.random.default_rng(2026)
    # Simulate 24 longitudinal samples across 6 donors x 4 ordered time stages (t=0,1,2,3)
    n_donors = 12
    n_stages = 4
    rows = []

    for d in range(n_donors):
        z_prev = rng.uniform(0.1, 0.4)
        a_prev = rng.uniform(0.1, 0.3)
        b_prev = rng.uniform(0.1, 0.3)
        c_prev = rng.uniform(0.1, 0.3)
        for st in range(n_stages):
            # Shared systemic inflammation confounder Z(t)
            z_t = np.clip(0.35 * z_prev + 0.22 * st + rng.normal(0, 0.05), 0.05, 0.98)
            if st == 0:
                a_t, b_t, c_t = a_prev, b_prev, c_prev
            else:
                # True bidirectional positive feedback between Myeloid (A) and Fibroblast (B)
                a_t = np.clip(0.25 * a_prev + 0.62 * b_prev + 0.12 * z_t + rng.normal(0, 0.03), 0.02, 0.98)
                b_t = np.clip(0.25 * b_prev + 0.65 * a_prev + 0.10 * z_t + rng.normal(0, 0.03), 0.02, 0.98)
                # Bystander_Epithelium (C) is ONLY driven by systemic confounder Z(t), zero causal input from A or B!
                c_t = np.clip(0.15 * c_prev + 0.82 * z_t + rng.normal(0, 0.03), 0.02, 0.98)

            rows.append(
                {
                    "donor_id": f"D{d}",
                    "stage_idx": st,
                    "Initiator_Myeloid": a_t,
                    "Amplifier_Fibroblast": b_t,
                    "Confounded_Bystander": c_t,
                    "Systemic_Confounder_Proxy": z_t,
                }
            )
            z_prev, a_prev, b_prev, c_prev = z_t, a_t, b_t, c_t

    longitudinal_df = pd.DataFrame(rows)
    compartment_cols = ["Initiator_Myeloid", "Amplifier_Fibroblast", "Confounded_Bystander"]

    # Spatial proximity matrix: Myeloid & Fibroblast co-localize (1.5), Bystander is physically segregated (-1.0)
    spatial_matrix = pd.DataFrame(
        [
            [0.0, 1.5, -1.0],
            [1.5, 0.0, -1.0],
            [-1.0, -1.0, 0.0],
        ],
        index=compartment_cols,
        columns=compartment_cols,
    )

    compiler = DeconfoundedAutoGrammarCompiler(min_partial_lag_gain=0.22)
    grammar = compiler.compile_from_longitudinal_states(
        longitudinal_df=longitudinal_df,
        compartment_cols=compartment_cols,
        pack_id="autoimmune_ibd_ra",
        organ_context="Gut_Mucosa",
        remission_compartment="Amplifier_Fibroblast",
        spatial_enrichment=spatial_matrix,
        confounder_col="Systemic_Confounder_Proxy",
    )

    edge_pairs = {(e.source, e.target) for e in grammar.edges}
    assert ("Initiator_Myeloid", "Amplifier_Fibroblast") in edge_pairs
    assert ("Amplifier_Fibroblast", "Initiator_Myeloid") in edge_pairs
    # Spurious confounder-driven edge to Confounded_Bystander must be eliminated!
    assert ("Initiator_Myeloid", "Confounded_Bystander") not in edge_pairs
    assert ("Amplifier_Fibroblast", "Confounded_Bystander") not in edge_pairs
    assert grammar.bic_score is not None and np.isfinite(grammar.bic_score)

    yaml_out = tmp_path / "compiled_ibd_grammar.yaml"
    compiler.export_yaml(grammar, yaml_out)
    assert yaml_out.exists()
    reloaded = compiler.load_yaml(yaml_out)
    assert len(reloaded.edges) == len(grammar.edges)
