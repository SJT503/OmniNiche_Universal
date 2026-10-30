"""Execute Stage II Task 3 (Deep Case #1 / Fig. 4):
MASH & Pan-Organ Fibrosis Phase III Failure Deconvolution + Zero-Shot Cross-Organ Transfer.

100% real primary human data:
- Organ 1 (Source Grammar): Pure Human Liver Healthy vs NASH vs F4 Cirrhosis + NAFLD Blood
  (GSE136103 Ramachandran et al., Nature 2019 + GSE212837 Wang et al., Sci. Transl. Med. 2023)
- Organ 2 (Zero-Shot Transfer #1): Human Lung Control vs ILD vs IPF
  (GSE135893 Habermann et al., Sci. Adv. 2020; 114,396 annotated single cells)
- Organ 3 (Zero-Shot Transfer #2): Human Kidney Control vs Glomerular CKD & Diabetic Kidney Disease
  (GSE271188; 17 human kidney snRNA-seq biopsies including Control, MCD, MN, IgAN, DKD)
- 14 Verbatim SHA-256-Locked Phase IIb/III Landmark Antifibrotic Clinical Trials (ClinicalTrials.gov API v2)
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.engine.autogrammar_compiler import DeconfoundedAutoGrammarCompiler
from omniniche.engine.grn_perturb import RidgeGRNPerturber
from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver
from omniniche.stage2.fibrosis_mash_deconvolution import (
    FIBROSIS_COMPARTMENT_GENES,
    LANDMARK_FIBROSIS_TRIALS,
    build_conserved_fibrosis_grammar_from_liver,
    evaluate_14_landmark_fibrosis_trials,
    extract_human_kidney_ckd_gse271188,
    extract_human_lung_ipf_gse135893,
    extract_pure_human_liver_gse136103_and_gse212837,
    fetch_verbatim_landmark_ncts,
)
from omniniche.stage2.pan_disease_drift_law import compute_dynamic_drift_metrics

# Directed biological receptor/nuclear-receptor -> pathogenic effector prior edges for RidgeGRNPerturber
FIBROSIS_PRIOR_EDGES: list[tuple[str, str]] = [
    # Metabolic & nuclear receptors regulating parenchymal stress and myofibroblast activation
    ("THRB", "GDF15"), ("THRB", "CCL2"), ("THRB", "TGFB2"), ("THRB", "SOX9"),
    ("GLP1R", "GDF15"), ("GLP1R", "CCL2"), ("GLP1R", "IL1B"), ("GLP1R", "SPP1"),
    ("PPARA", "GDF15"), ("PPARA", "CCL2"),
    ("PPARD", "CCL2"), ("PPARD", "IL1B"),
    ("PPARG", "ACTA2"), ("PPARG", "COL1A1"), ("PPARG", "COL1A2"), ("PPARG", "POSTN"), ("PPARG", "TGFB1"),
    ("NR3C2", "SPP1"), ("NR3C2", "TGFB1"), ("NR3C2", "COL1A1"), ("NR3C2", "FN1"), ("NR3C2", "TIMP1"),
    # Receptor tyrosine kinases & TGF-beta axis driving myofibroblast ECM deposition
    ("FGFR1", "COL1A1"), ("FGFR1", "ACTA2"), ("FGFR1", "POSTN"),
    ("PDGFRA", "COL1A1"), ("PDGFRA", "COL3A1"), ("PDGFRA", "ACTA2"), ("PDGFRA", "CTHRC1"),
    ("PDGFRB", "COL1A1"), ("PDGFRB", "COL1A2"), ("PDGFRB", "POSTN"), ("PDGFRB", "FN1"),
    ("KDR", "POSTN"), ("KDR", "VCAM1"),
    ("TGFB1", "TGFBR1"), ("TGFB1", "COL1A1"), ("TGFB1", "COL1A2"), ("TGFB1", "COL3A1"), ("TGFB1", "POSTN"), ("TGFB1", "ACTA2"),
    ("TGFBR1", "COL1A1"), ("TGFBR1", "COL1A2"), ("TGFBR1", "POSTN"), ("TGFBR1", "CTHRC1"),
    # SAMac myeloid axis (TREM2/SPP1 vs upstream monocyte chemokine CCR2/CCR5 & apoptosis kinase MAP3K5)
    ("TREM2", "SPP1"), ("TREM2", "GPNMB"), ("TREM2", "CD9"), ("TREM2", "TGFB1"),
    ("SPP1", "COL1A1"), ("SPP1", "POSTN"), ("SPP1", "ACTA2"),
    ("MAP3K5", "GDF15"), ("MAP3K5", "ICAM1"),
    ("CCR2", "IL1B"), ("CCR5", "IL1B"),
    # Terminal ECM cross-linking / matricellular proteins (downstream effectors with minimal upstream TF control)
    ("COL1A1", "LOXL2"), ("TGFB1", "CCN2"),
]


def main() -> None:
    force_reextract = "--reextract" in sys.argv
    vault_sc = ROOT / "data_vault" / "raw_single_cell"
    vault_pb = ROOT / "data_vault" / "processed_pseudobulk"
    vault_ct = ROOT / "data_vault" / "raw_clinical_trials"
    out_tables = ROOT / "outputs" / "stage2_tables"
    vault_pb.mkdir(parents=True, exist_ok=True)
    out_tables.mkdir(parents=True, exist_ok=True)

    t2 = pd.read_csv(out_tables / "Table_S2_Target_Resolved_Trials_Deconvolution.csv", usecols=["target_genes"])
    clinical_genes: set[str] = set()
    for val in t2["target_genes"].dropna().astype(str):
        for g in val.split(";"):
            gs = g.strip().upper()
            if gs:
                clinical_genes.add(gs)
    for gset in FIBROSIS_COMPARTMENT_GENES.values():
        clinical_genes |= gset
    for tr in LANDMARK_FIBROSIS_TRIALS:
        clinical_genes |= {g.upper() for g in tr["targets"]}

    print(f">>> Tracking {len(clinical_genes)} human target & fibrotic niche genes across 3 human organs...")
    tau0_df = pd.read_csv(vault_pb / "hpa_pan_organ_tau0_table.csv", index_col=0)

    print("\n>>> [1/6] Fetching & SHA-256-locking 14 verbatim Phase IIb/III antifibrotic trials (ClinicalTrials.gov API v2)...")
    verified_ncts = fetch_verbatim_landmark_ncts(vault_ct / "landmark_14_fibrosis_ncts_sha256.json")
    print(f"    Verified {len(verified_ncts)} landmark NCT trials with SHA-256 provenance.")

    # Organ 1: Pure Human Liver
    liv_e_p = vault_pb / "axis1_liver_fibrosis_early_comp.csv"
    liv_l_p = vault_pb / "axis1_liver_fibrosis_late_comp.csv"
    liv_d_p = out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Liver.csv"
    ccr2_p = out_tables / "Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv"
    if not force_reextract and all(p.exists() for p in (liv_e_p, liv_l_p, liv_d_p, ccr2_p)):
        print("\n>>> [2/6] Loading verified Pure Human Liver pseudobulk & donor profiles from data_vault cache...")
        liv = {
            "early_comp_df": pd.read_csv(liv_e_p, index_col=0),
            "late_comp_df": pd.read_csv(liv_l_p, index_col=0),
            "donors_df": pd.read_csv(liv_d_p),
            "ccr2_paradox_df": pd.read_csv(ccr2_p),
        }
    else:
        print("\n>>> [2/6] Streaming Pure Human Liver (GSE136103 Nature 2019 + GSE212837 Sci. Transl. Med. 2023)...")
        liv = extract_pure_human_liver_gse136103_and_gse212837(vault_sc, clinical_genes)
        liv["early_comp_df"].to_csv(liv_e_p, encoding="utf-8")
        liv["late_comp_df"].to_csv(liv_l_p, encoding="utf-8")
        liv["donors_df"].to_csv(liv_d_p, index=False, encoding="utf-8")
        liv["ccr2_paradox_df"].to_csv(ccr2_p, index=False, encoding="utf-8")
    liver_drift = compute_dynamic_drift_metrics(liv["early_comp_df"], liv["late_comp_df"])
    liver_drift.to_csv(vault_pb / "axis1_liver_fibrosis_drift.csv", encoding="utf-8")

    # Organ 2: Human Lung IPF
    lng_e_p = vault_pb / "axis2_lung_ipf_early_comp.csv"
    lng_l_p = vault_pb / "axis2_lung_ipf_late_comp.csv"
    lng_d_p = out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Lung.csv"
    if not force_reextract and all(p.exists() for p in (lng_e_p, lng_l_p, lng_d_p)):
        print("\n>>> [3/6] Loading verified Human Lung IPF (GSE135893; 114,396 annotated cells) from data_vault cache...")
        lng = {
            "early_comp_df": pd.read_csv(lng_e_p, index_col=0),
            "late_comp_df": pd.read_csv(lng_l_p, index_col=0),
            "donors_df": pd.read_csv(lng_d_p),
        }
    else:
        print("\n>>> [3/6] Streaming Human Lung IPF (GSE135893 Science Advances 2020; 114,396 annotated cells)...")
        lng = extract_human_lung_ipf_gse135893(vault_sc, clinical_genes)
        lng["early_comp_df"].to_csv(lng_e_p, encoding="utf-8")
        lng["late_comp_df"].to_csv(lng_l_p, encoding="utf-8")
        lng["donors_df"].to_csv(lng_d_p, index=False, encoding="utf-8")
    lung_drift = compute_dynamic_drift_metrics(lng["early_comp_df"], lng["late_comp_df"])
    lung_drift.to_csv(vault_pb / "axis2_lung_ipf_drift.csv", encoding="utf-8")

    # Organ 3: Human Kidney CKD/DKD
    kid_e_p = vault_pb / "axis3_kidney_ckd_early_comp.csv"
    kid_l_p = vault_pb / "axis3_kidney_ckd_late_comp.csv"
    kid_d_p = out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Kidney.csv"
    if not force_reextract and all(p.exists() for p in (kid_e_p, kid_l_p, kid_d_p)):
        print("\n>>> [4/6] Loading verified Human Kidney CKD/DKD (GSE271188; 17 biopsies) from data_vault cache...")
        kid = {
            "early_comp_df": pd.read_csv(kid_e_p, index_col=0),
            "late_comp_df": pd.read_csv(kid_l_p, index_col=0),
            "donors_df": pd.read_csv(kid_d_p),
        }
    else:
        print("\n>>> [4/6] Streaming Human Kidney CKD & Diabetic Kidney Disease (GSE271188; 17 biopsies)...")
        kid = extract_human_kidney_ckd_gse271188(vault_sc, clinical_genes)
        kid["early_comp_df"].to_csv(kid_e_p, encoding="utf-8")
        kid["late_comp_df"].to_csv(kid_l_p, encoding="utf-8")
        kid["donors_df"].to_csv(kid_d_p, index=False, encoding="utf-8")
    kidney_drift = compute_dynamic_drift_metrics(kid["early_comp_df"], kid["late_comp_df"])
    kidney_drift.to_csv(vault_pb / "axis3_kidney_ckd_drift.csv", encoding="utf-8")

    all_donors = pd.concat([liv["donors_df"], lng["donors_df"], kid["donors_df"]], ignore_index=True)
    all_donors.to_csv(out_tables / "Table_S4a_PanOrgan_Fibrosis_All_Human_Donors.csv", index=False, encoding="utf-8")

    print("\n>>> [5/6] Compiling Human Liver Bistable Fibrotic Grammar + Directed RidgeGRNPerturber + Zero-Shot Arena...")
    grammar = build_conserved_fibrosis_grammar_from_liver(liv["donors_df"])
    DeconfoundedAutoGrammarCompiler.export_yaml(
        grammar, ROOT / "outputs" / "compiled_grammars" / "conserved_pan_organ_fibrosis_v2.yaml"
    )

    grn_genes = sorted(
        (
            FIBROSIS_COMPARTMENT_GENES["Parenchyma_Injury"]
            | FIBROSIS_COMPARTMENT_GENES["SPP1_SAMac"]
            | FIBROSIS_COMPARTMENT_GENES["Myofibroblast"]
        )
        & set(liv["late_comp_df"].index)
    )
    # Fit RidgeGRNPerturber STRICTLY on Human Liver (Source Organ) so Lung IPF & Kidney CKD remain 100% ZERO-SHOT!
    expr_samples = pd.concat(
        [
            liv["early_comp_df"].loc[grn_genes].T,
            liv["late_comp_df"].loc[grn_genes].T,
            0.5 * (liv["early_comp_df"].loc[grn_genes].T + liv["late_comp_df"].loc[grn_genes].T),
        ],
        ignore_index=True,
    )
    pathogenic_core = sorted(
        {"COL1A1", "COL1A2", "COL3A1", "POSTN", "CTHRC1", "ACTA2", "SPP1", "TREM2", "TGFB1", "GDF15", "CCL2"}
        & set(grn_genes)
    )
    perturber = RidgeGRNPerturber(ridge_lambda=0.6, propagation_steps=3).fit(
        expr_samples,
        prior_edges=FIBROSIS_PRIOR_EDGES,
        pathogenic_program_genes=pathogenic_core,
    )

    df_trials, df_arena = evaluate_14_landmark_fibrosis_trials(
        grammar=grammar,
        tau0_df=tau0_df,
        liver_drift_df=liver_drift,
        lung_drift_df=lung_drift,
        kidney_drift_df=kidney_drift,
        verified_ncts=verified_ncts,
        grn_perturber=perturber,
        ccr2_paradox_df=liv["ccr2_paradox_df"],
    )
    df_trials.to_csv(out_tables / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv", index=False, encoding="utf-8")
    df_arena.to_csv(out_tables / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv", index=False, encoding="utf-8")
    print(
        df_trials[
            [
                "nct_id", "trial_acronym", "drug", "organ", "stage_enrolled",
                "clinical_win", "tau0_healthy", "organ_delta_tau", "grn_causal_eta",
                "scar_barrier_access", "effective_perturbation_eta", "omniniche_remission_prob",
            ]
        ].to_string(index=False)
    )
    print("\n--- 5-Baseline Calibrated Arena (14 Phase IIb/III Antifibrotic Trials) ---")
    print(df_arena.to_string(index=False))

    print("\n>>> [6/6] Running Profile-Likelihood Identifiability, t50 Window Cliff & Locked-F4 Bliss Synergy Scan...")
    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)

    # Profile Likelihood on Human Liver GSE136103 Healthy -> F4 Cirrhosis progression
    comps = ["Parenchyma_Injury", "SPP1_SAMac", "Myofibroblast"]
    g136 = liv["donors_df"][liv["donors_df"]["organ"] == "Liver_GSE136103"].copy()
    for c in comps:
        v = g136[c].to_numpy(dtype=float)
        g136[c] = 0.14 + 0.68 * (v - np.min(v)) / (np.max(v) - np.min(v) + 1e-8)
    s_means = g136.groupby("stage")[comps].mean()
    h_vec = s_means.loc["Control_Healthy"].to_numpy(dtype=float)
    c_vec = s_means.loc["F4_Cirrhosis"].to_numpy(dtype=float)
    mid_vec = 0.5 * (h_vec + c_vec)
    obs_df = pd.DataFrame(
        [
            {"cohort_id": "GSE136103_Human_Liver", "time_hr": 4.0, **dict(zip(comps, h_vec))},
            {"cohort_id": "GSE136103_Human_Liver", "time_hr": 38.0, **dict(zip(comps, mid_vec))},
            {"cohort_id": "GSE136103_Human_Liver", "time_hr": 88.0, **dict(zip(comps, c_vec))},
        ]
    )
    prof_df = solver.profile_likelihood_identifiability(obs_df=obs_df, s0=np.array([0.48, 0.42, 0.40]), obs_sigma=0.04)
    prof_df.to_csv(out_tables / "Table_S4d_Fibrosis_Profile_Likelihood_Identifiability.csv", index=False, encoding="utf-8")

    # Monotherapy window cliff t50 across delayed onset (simulating natural progression from F1 -> F2 -> F3 -> F4)
    t50_res = solver.compute_monotherapy_window_t50(
        s0=np.array([0.42, 0.22, 0.18], dtype=float),
        target_compartment="Parenchyma_Injury",
        perturbation_eta=0.84,
        modality="small_molecule",
        onset_grid=np.linspace(0.0, 48.0, 13),
        n_patients=220,
    )
    pd.DataFrame(
        {
            "onset_time_hr": t50_res["onset_times"],
            "parenchymal_monotherapy_remission_rate": t50_res["remission_rates"],
            "t50_window_cliff_hr": t50_res["t50_window_cliff"],
        }
    ).to_csv(out_tables / "Fig4_Monotherapy_Window_Cliff_Curve.csv", index=False, encoding="utf-8")

    # Bliss dual-target synergy matrix in deeply locked F4 Cirrhosis (s0 = [0.85, 0.88, 0.90])
    # Evaluate realistic clinical partial-occupancy regimens where monotherapies fail to cross the F4 saddle
    s0_locked = np.array([0.85, 0.88, 0.90], dtype=float)
    candidate_arms = [
        ("THRB_Resmetirom_Parenchyma_SM", InterventionRegimen("Parenchyma_Injury", 0.72, "small_molecule")),
        ("GLP1R_Semaglutide_Metabolic_Peptide", InterventionRegimen("Parenchyma_Injury", 0.68, "peptide")),
        ("TGFBR1_PDGFRB_Myofibro_LowDose_SM", InterventionRegimen("Myofibroblast", 0.145, "small_molecule")),
        ("TREM2_SPP1_SAMac_Modulator_SM", InterventionRegimen("SPP1_SAMac", 0.140, "small_molecule")),
        ("LOXL2_Simtuzumab_StericBlocked_mAb", InterventionRegimen("Myofibroblast", 0.004, "antibody")),
        ("CCR2_Cenicriviroc_TargetLost_SM", InterventionRegimen("SPP1_SAMac", 0.005, "small_molecule")),
    ]
    syn_rows = []
    for i in range(len(candidate_arms)):
        for j in range(i + 1, len(candidate_arms)):
            na, ra = candidate_arms[i]
            nb, rb = candidate_arms[j]
            if ra.target_compartment == rb.target_compartment:
                continue
            syn = solver.evaluate_bliss_synergy(s0=s0_locked, regimen_a=ra, regimen_b=rb, n_patients=220, seed=700 + i * 10 + j)
            syn_rows.append(
                {
                    "arm_a": na,
                    "arm_b": nb,
                    "comp_a": ra.target_compartment,
                    "comp_b": rb.target_compartment,
                    **{k: round(float(v), 4) for k, v in syn.items()},
                }
            )
    syn_df = pd.DataFrame(syn_rows).sort_values(["bliss_excess", "combo_remission_rate"], ascending=False)
    syn_df.to_csv(out_tables / "Table_S4e_Fibrosis_LockedStage_Bliss_Synergy.csv", index=False, encoding="utf-8")

    qp_untr = solver.compute_quasi_potential_landscape([], "SPP1_SAMac", "Myofibroblast", n_grid=32, n_trajectories=240, seed=77)
    qp_combo = solver.compute_quasi_potential_landscape(
        [
            InterventionRegimen("Parenchyma_Injury", 0.82, "small_molecule"),
            InterventionRegimen("Myofibroblast", 0.45, "small_molecule"),
        ],
        "SPP1_SAMac",
        "Myofibroblast",
        n_grid=32,
        n_trajectories=240,
        seed=78,
    )
    np.savez_compressed(
        out_tables / "Fig4_Waddington_QuasiPotential_Fibrosis.npz",
        grid_centers=qp_untr["grid_centers"],
        U_untreated=qp_untr["U_grid"],
        U_combo_rescue=qp_combo["U_grid"],
        lock_depth_untreated=qp_untr["pathological_lock_depth"],
        lock_depth_combo=qp_combo["pathological_lock_depth"],
    )
    print(f"    Profile Likelihood Identifiable Edges: {int(prof_df['identifiable_95ci'].sum())}/{len(prof_df)}")
    print(f"    t50 Window Cliff = {t50_res['t50_window_cliff']:.2f}h | Top Locked-F4 Dual Bliss Synergy:")
    print(syn_df.head(6).to_string(index=False))
    print("\n=== STAGE II TASK 3 (PAN-ORGAN FIBROSIS & MASH DECONVOLUTION) COMPLETE ===")


if __name__ == "__main__":
    main()
