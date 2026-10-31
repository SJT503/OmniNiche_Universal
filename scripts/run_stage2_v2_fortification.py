"""Execute Stage II v2.2 Complete Statistical, Causal, Genetic & Biophysical Fortification (T1-T10).

Runs deterministically (< 2.5 GB peak RAM) on SHA-256-locked empirical caches in F:/OmniNiche_Universal,
exporting all fortified CSV/NPZ tables required for Nature Medicine / Science v2.2 figures and manuscript:
- T1: Non-circular 240-trial manual clinical stoppage adjudication benchmark & unified clause-aware stoppage labels.
- T2: De-heuristicized data-driven fibrosis Langevin grammar (>25% empirical coupling span, objective single-cell
      compartment weights, absent intrahepatic GLP1R) & expanded 24-trial Phase IIb/III RCT benchmark (+ StageBlind/StagePrior).
- T3: Exact Fisher's exact test ORs, exact Wilson score 95% binomial CIs, and 100% code-to-Methods mathematical parity.
- T4: Open Targets Platform v4 human genetic support orthogonality benchmark & ChEMBL/OpenTargets regulatory approval Regime R6.
- T5: Empirical single-cell/HPA ECM density in B_acc, single-target monotherapy Group5Fold CV, and 6-vs-12 axis audit.
- T6/T7: Response-label-free Axis-4 tumor drift, exact umi_stride=1 streaming, 4-cohort pre-treatment anti-PD-1 benchmark
      (N=83 pure monotherapy + N=107 all pre-treatment biopsies) with LOPO scaling, published biomarkers (18-gene GEP,
      CAF/TIDE Exclusion, CD8:myCAF, PD-L1), and 10x Visium P_ij spatial ablation.
- T8/T9: Strictly peripheral solid-parenchymal interstitial 9-axis Regime R7 & dual Scientific/Maturation Table S3b blocks.
"""
from __future__ import annotations

import json
import pathlib
import sys
import tracemalloc
from typing import Any

import numpy as np
import pandas as pd

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.engine.grn_perturb import RidgeGRNPerturber
from omniniche.stage2.fibrosis_mash_deconvolution import (
    FIBROSIS_COMPARTMENT_GENES,
    build_conserved_fibrosis_grammar_from_liver,
    evaluate_14_landmark_fibrosis_trials,
)
from omniniche.stage2.pan_disease_drift_law import (
    compute_empirical_axis_ecm_density_map,
    compute_modality_barrier_compatibility,
)
from omniniche.stage2.solid_tumor_spatial_icb import (
    build_default_solid_tumor_icb_grammar,
    score_baseline_patients_via_langevin_icb,
)
from omniniche.stage2.v2_statistical_fortification import (
    _oof_logistic_predictions,
    _wilson_ci95,
    audit_stoppage_classifier_benchmark,
    bootstrap_classification_metrics,
    compute_2d_fibrosis_bliss_matrix_and_landscape,
    compute_2d_tumor_waddington_surface_shift,
    compute_calibration_deciles,
    compute_continuous_fibrosis_hysteresis_curve,
    compute_continuous_profile_likelihood_curves,
    compute_decision_curve_and_enrichment,
    compute_macro_longitudinal_and_modality_odds,
    evaluate_bidirectional_3x3_fibrosis_transfer,
    evaluate_human_genetics_vs_dnacc_orthogonality,
    evaluate_six_step_dnacc_ablation,
    evaluate_unreported_stop_sensitivity_regimes,
    evaluate_within_modality_dnacc_stratification,
    exact_permutation_auroc_pvalue,
    extract_12section_visium_radial_and_contrasting_archetypes,
    fast_delong_paired_test,
    fit_multivariable_fixed_effects_logistic,
)
from scripts.run_stage2_task2_drift_law import unify_clause_aware_stoppage_tables
from scripts.run_stage2_task3_fibrosis_mash import FIBROSIS_PRIOR_EDGES
from scripts.run_stage2_task4_solid_tumor_icb import (
    ALL_PANEL_GENES,
    PANEL_GENES,
    _lopo_robust_01,
    fast_stream_tsv_gz_gene_subset,
)


def _assemble_full_target_resolved_drift_df(
    master_df: pd.DataFrame,
    tau0_df: pd.DataFrame,
    atlas_6ax: pd.DataFrame,
    emp_ecm_map: dict[str, float],
    ot_payload: dict[str, Any],
) -> pd.DataFrame:
    t2 = master_df[master_df["n_targets"] > 0].copy()
    axis_drift_tables = {ax: sub.set_index("gene_symbol") for ax, sub in atlas_6ax.groupby("disease_axis")}
    num_cols = [
        "tau_early",
        "tau_late",
        "delta_tau",
        "compartment_switch",
        "primary_comp_delta_expr",
        "stromal_ectopic_share_late",
        "dynamic_retention_score",
    ]
    consensus_drift = atlas_6ax.groupby("gene_symbol")[num_cols].mean()
    consensus_drift["compartment_switch"] = (consensus_drift["compartment_switch"] >= 0.5).astype(int)

    axis_alias = {
        "Hepatic_MASH_Fibrosis": axis_drift_tables["Hepatic_MASH_Fibrosis"],
        "Pulmonary_IPF_ARDS": axis_drift_tables["Pulmonary_IPF_ARDS"],
        "Renal_CKD_AKI": axis_drift_tables["Renal_CKD_AKI"],
        "Solid_Oncology": axis_drift_tables["Solid_Oncology"],
        "Autoimmune_IBD_Rheum": axis_drift_tables["Autoimmune_IBD_Rheum"],
        "Dermatologic_Immune": axis_drift_tables["Autoimmune_IBD_Rheum"],
        "Infectious_Sepsis": axis_drift_tables["Infectious_Sepsis"],
    }

    axis_ot_scores: dict[str, dict[str, dict[str, float]]] = ot_payload.get("axis_gene_scores", {})
    pan_ot_scores: dict[str, dict[str, float]] = ot_payload.get("pan_gene_scores", {})
    axis_app_targets: dict[str, set[str]] = {
        k: {str(x).upper() for x in v} for k, v in ot_payload.get("axis_approved_targets", {}).items()
    }
    axis_app_drugs: dict[str, set[str]] = {
        k: {str(x).strip().lower() for x in v} for k, v in ot_payload.get("axis_approved_drugs", {}).items()
    }

    tau0_dict = tau0_df["tau0_healthy"].to_dict()
    vb_dict = tau0_df["vital_bystander_hazard"].to_dict()
    records: list[dict[str, Any]] = []
    for row in t2.itertuples(index=False):
        ax = str(row.disease_axis)
        dtbl = axis_alias.get(ax, consensus_drift)
        tgenes = [g.strip().upper() for g in str(row.target_genes).split(";") if g.strip()]
        if not tgenes:
            continue
        t0_list = [tau0_dict[g] for g in tgenes if g in tau0_dict]
        vb_list = [vb_dict[g] for g in tgenes if g in vb_dict]
        tau0_v = float(np.mean(t0_list)) if t0_list else 0.50
        vb_v = float(np.mean(vb_list)) if vb_list else 0.35

        matched_g = [g for g in tgenes if g in dtbl.index]
        if matched_g:
            sub_g = dtbl.loc[matched_g]
            te_v = float(sub_g["tau_early"].mean())
            tl_v = float(sub_g["tau_late"].mean())
            dt_v = float(sub_g["delta_tau"].mean())
            sw_v = int(float(sub_g["compartment_switch"].mean()) >= 0.5)
            st_v = float(sub_g["stromal_ectopic_share_late"].mean())
            dr_v = float(sub_g["dynamic_retention_score"].mean())
        else:
            te_v, tl_v, dt_v, sw_v, st_v, dr_v = tau0_v, tau0_v, 0.0, 0, 0.25, 0.50 * tau0_v

        bar_v = compute_modality_barrier_compatibility(
            str(row.drug_modalities),
            ax,
            st_v,
            sw_v,
            axis_ecm_density_map=emp_ecm_map,
        )

        # Open Targets Platform v4 Human Genetic Support & ChEMBL/OpenTargets Axis Regulatory Approval
        ax_ot = axis_ot_scores.get(ax, {})
        ga_vals, sm_vals, hm_vals = [], [], []
        for g in tgenes:
            grec = ax_ot.get(g) or pan_ot_scores.get(g)
            if grec:
                ga_vals.append(float(grec.get("genetic_association", 0.0)))
                sm_vals.append(float(grec.get("somatic_mutation", 0.0)))
                hm_vals.append(float(grec.get("human_genetic_max", 0.0)))
        ot_ga = float(max(ga_vals)) if ga_vals else 0.0
        ot_sm = float(max(sm_vals)) if sm_vals else 0.0
        ot_hm = float(max(hm_vals)) if hm_vals else 0.0

        mdrugs = [d.strip().lower() for d in str(getattr(row, "matched_drugs", "") or "").split(";") if d.strip()]
        app_tgts = axis_app_targets.get(ax, set())
        app_drgs = axis_app_drugs.get(ax, set())
        is_app = int(any(d in app_drgs for d in mdrugs) or any(g in app_tgts for g in tgenes))

        records.append(
            {
                "nct_id": row.nct_id,
                "disease_axis": ax,
                "is_core_6axis": ax in axis_drift_tables,
                "overall_status": row.overall_status,
                "why_stopped": row.why_stopped,
                "stop_deconvolution": row.stop_deconvolution,
                "naive_keyword_stop": getattr(row, "naive_keyword_stop", row.stop_deconvolution),
                "max_phase_num": int(row.max_phase_num),
                "has_results": row.has_results,
                "enrollment": row.enrollment,
                "start_year": row.start_year,
                "sponsor_class": row.sponsor_class,
                "target_genes": str(row.target_genes),
                "primary_target_gene": row.primary_target_gene,
                "matched_drugs": str(getattr(row, "matched_drugs", "") or ""),
                "drug_modalities": row.drug_modalities,
                "n_targets": int(row.n_targets),
                "tau0_healthy": round(tau0_v, 5),
                "vital_bystander_hazard": round(vb_v, 5),
                "tau_early": round(te_v, 5),
                "tau_late": round(tl_v, 5),
                "delta_tau": round(dt_v, 5),
                "compartment_switch": sw_v,
                "stromal_ectopic_share_late": round(st_v, 5),
                "dynamic_retention_score": round(dr_v, 5),
                "barrier_compatibility": round(bar_v, 5),
                "ot_genetic_assoc": round(ot_ga, 5),
                "ot_somatic_mut": round(ot_sm, 5),
                "ot_human_genetic_max": round(ot_hm, 5),
                "is_approved_in_axis": is_app,
            }
        )
    return pd.DataFrame(records)


def _compute_gse123813_facs_deconfounding_sensitivity(
    sc_dir: pathlib.Path,
    grammar: Any,
) -> pd.DataFrame:
    """Compare naive CD3+-diluted pooling vs non-CD3+ FACS-deconfounded TME evaluation on GSE123813 BCC (n=8, exact umi_stride=1)."""
    cdir = sc_dir / "GSE123813_BCC_SCC_AntiPD1"
    clin = pd.read_excel(cdir / "Yost_2019_NatMed_SuppTables.xlsx", sheet_name="SuppTable1", skiprows=3).dropna(subset=["Patient"])
    clin["Patient"] = clin["Patient"].astype(str).str.strip()
    bcc_clin = clin[clin["Tumor Type"] == "BCC"].set_index("Patient")

    bcc_meta = pd.read_csv(cdir / "GSE123813_bcc_all_metadata.txt.gz", sep="\t")
    bcc_meta["cell.id"] = bcc_meta["cell.id"].astype(str)
    barcodes, gene_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE123813_bcc_scRNA_counts.txt.gz", ALL_PANEL_GENES, umi_stride=1
    )
    bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    bcc_meta = bcc_meta[bcc_meta["cell.id"].isin(bc_to_idx)].copy()
    bcc_meta["col_idx"] = [bc_to_idx[c] for c in bcc_meta["cell.id"]]
    pre_meta = bcc_meta[bcc_meta["treatment"] == "pre"].copy()

    def _mp(idxs: np.ndarray, genes: set[str]) -> float:
        present = [gene_rows[g][idxs] for g in genes if g in gene_rows]
        if not present or len(idxs) == 0:
            return 0.0
        return float(np.mean([np.mean(v) for v in present]))

    rows = []
    for pat in sorted(bcc_clin.index):
        p_pre = pre_meta[pre_meta["patient"] == pat]
        non_t = p_pre[p_pre["sort"] != "CD45+ CD3+"]
        if len(non_t) < 15:
            continue
        is_resp = 1 if str(bcc_clin.loc[pat, "Response"]).strip().lower().startswith("yes") else 0
        cd3_sort_frac = float((p_pre["sort"] == "CD45+ CD3+").mean())

        naive_caf_frac = float(p_pre["cluster"].isin(["CAFs", "Myofibroblasts"]).sum()) / max(len(p_pre), 1)
        deconf_caf_frac = float(non_t["cluster"].isin(["CAFs", "Myofibroblasts"]).sum()) / max(len(non_t), 1)
        naive_tam_frac = float(p_pre["cluster"].isin(["Macrophages", "DCs"]).sum()) / max(len(p_pre), 1)
        deconf_tam_frac = float(non_t["cluster"].isin(["Macrophages", "DCs"]).sum()) / max(len(non_t), 1)

        all_idx = p_pre["col_idx"].to_numpy(dtype=int)
        nt_idx = non_t["col_idx"].to_numpy(dtype=int)
        cd8_pre = p_pre[p_pre["cluster"].str.contains("CD8", na=False)]["col_idx"].to_numpy(dtype=int)
        n_t_all = float(p_pre["cluster"].isin(["CD4_T_cells", "Tregs", "CD8_act_T_cells", "CD8_ex_T_cells", "CD8_mem_T_cells", "Tcell_prolif"]).sum())
        cd8_reac = float(p_pre["cluster"].isin(["CD8_act_T_cells", "CD8_ex_T_cells", "Tcell_prolif"]).sum()) / max(n_t_all, 1.0)
        cd8_ex = float((p_pre["cluster"] == "CD8_ex_T_cells").sum()) / max(float(len(cd8_pre)), 1.0)
        eff_prog = _mp(cd8_pre if len(cd8_pre) >= 3 else all_idx, PANEL_GENES["CD8_Effector"]) * (0.45 + 2.5 * cd8_reac)
        exh_prog = _mp(cd8_pre if len(cd8_pre) >= 3 else all_idx, PANEL_GENES["CD8_Exhaustion"]) + 0.8 * cd8_ex

        rows.append(
            {
                "patient_id": pat,
                "responder": is_resp,
                "n_pre_total": int(len(p_pre)),
                "n_pre_non_cd3": int(len(non_t)),
                "cd3_facs_sort_share": round(cd3_sort_frac, 4),
                "naive_pooled_caf_frac": round(naive_caf_frac, 4),
                "deconfounded_nonT_caf_frac": round(deconf_caf_frac, 4),
                "caf_dilution_fold": round(deconf_caf_frac / max(naive_caf_frac, 1e-5), 2),
                "naive_pooled_tam_frac": round(naive_tam_frac, 4),
                "deconfounded_nonT_tam_frac": round(deconf_tam_frac, 4),
                "naive_tum": _mp(all_idx, PANEL_GENES["Tumor"]) + 2.8 * (float(p_pre["cluster"].isin(["Tumor_1", "Tumor_2"]).sum()) / len(p_pre)),
                "naive_tam": _mp(all_idx, PANEL_GENES["TAM"]) + 2.2 * naive_tam_frac,
                "naive_caf": _mp(all_idx, PANEL_GENES["myCAF"]) + 2.8 * naive_caf_frac,
                "deconf_tum": _mp(nt_idx, PANEL_GENES["Tumor"]) + 2.8 * (float(non_t["cluster"].isin(["Tumor_1", "Tumor_2"]).sum()) / len(non_t)),
                "deconf_tam": _mp(nt_idx, PANEL_GENES["TAM"]) + 2.2 * deconf_tam_frac,
                "deconf_caf": _mp(nt_idx, PANEL_GENES["myCAF"]) + 2.8 * deconf_caf_frac,
                "raw_eff": eff_prog,
                "raw_exh_net": exh_prog - 0.45 * eff_prog,
            }
        )
    pdf = pd.DataFrame(rows).set_index("patient_id")
    for pref in ("naive", "deconf"):
        tmp = pd.DataFrame(index=pdf.index)
        tmp["Tumor"] = _lopo_robust_01(pdf[f"{pref}_tum"].to_numpy(), 0.16, 0.86)
        tmp["TAM"] = _lopo_robust_01(pdf[f"{pref}_tam"].to_numpy(), 0.14, 0.88)
        tmp["myCAF"] = _lopo_robust_01(pdf[f"{pref}_caf"].to_numpy(), 0.14, 0.88)
        tmp["CD8_Exhaustion"] = _lopo_robust_01(pdf["raw_exh_net"].to_numpy(), 0.16, 0.88)
        tmp["CD8_Effector_Reserve"] = _lopo_robust_01(pdf["raw_eff"].to_numpy(), 0.16, 0.90)
        sim = score_baseline_patients_via_langevin_icb(tmp, grammar=grammar, n_ensemble=250, seed=101, compute_combos=False)
        pdf[f"{pref}_icb_prob"] = sim["simulated_icb_remission_prob"]

    keep = [
        "responder",
        "n_pre_total",
        "n_pre_non_cd3",
        "cd3_facs_sort_share",
        "naive_pooled_caf_frac",
        "deconfounded_nonT_caf_frac",
        "caf_dilution_fold",
        "naive_pooled_tam_frac",
        "deconfounded_nonT_tam_frac",
        "naive_icb_prob",
        "deconf_icb_prob",
    ]
    return pdf[keep].reset_index()


def main() -> None:
    tracemalloc.start()
    vault_ct = ROOT / "data_vault" / "raw_clinical_trials"
    vault_pb = ROOT / "data_vault" / "processed_pseudobulk"
    vault_sc = ROOT / "data_vault" / "raw_single_cell"
    out_tables = ROOT / "outputs" / "stage2_tables"
    out_tables.mkdir(parents=True, exist_ok=True)

    print(">>> [Step 1/7] Unifying Clause-Aware Stoppage Deconvolution & Non-Circular 240-Trial Manual Adjudication...", flush=True)
    unify_clause_aware_stoppage_tables(vault_ct, out_tables)

    master_df = pd.read_csv(vault_ct / "non_cns_phase234_trials_master.csv", low_memory=False)
    tau0_df = pd.read_csv(vault_pb / "hpa_pan_organ_tau0_table.csv", index_col=0)
    atlas_6ax = pd.read_csv(out_tables / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv")

    s1_p = out_tables / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv"
    s1_df = pd.read_csv(s1_p)
    print(f"    Unified Table S1 ({len(s1_df)} peripheral organ axes; total trials = {len(master_df):,}).", flush=True)

    audit_metrics_df, audit_cm_df, audit_sample_df = audit_stoppage_classifier_benchmark(master_df, n_per_stratum=60, seed=42)
    audit_metrics_df.to_csv(out_tables / "Table_S2a_Stoppage_Auditor_Validation.csv", index=False, encoding="utf-8")
    audit_cm_df.to_csv(out_tables / "Table_S2a_Stoppage_Confusion_Matrix.csv", encoding="utf-8")
    audit_sample_df.to_csv(out_tables / "Table_S2a_Stoppage_Audit_240_Sample.csv", index=False, encoding="utf-8")
    print("    Non-Circular Manual Stoppage Adjudication Benchmark (n=240):")
    print(audit_metrics_df.to_string(index=False), flush=True)

    trend_df, odds_df, enr_dist_df = compute_macro_longitudinal_and_modality_odds(master_df)
    trend_df.to_csv(out_tables / "Table_S2b_Longitudinal_2005_2025_Trends.csv", index=False, encoding="utf-8")
    odds_df.to_csv(out_tables / "Table_S2c_Modality_Termination_Odds.csv", index=False, encoding="utf-8")
    enr_dist_df.to_csv(out_tables / "Table_S2d_Phase_Enrollment_Distributions.csv", index=False, encoding="utf-8")

    print("\n>>> [Step 2/7] Empirical ECM Density, Open Targets Human Genetics vs DNAcc, 6-Step Ablation & GroupKFold...", flush=True)
    emp_ecm_map = compute_empirical_axis_ecm_density_map(
        vault_pb,
        vault_sc / "hpa_pan_organ_normal" / "rna_single_cell_type.tsv",
    )
    ecm_rows = [{"disease_axis": k, "empirical_ecm_density": round(v, 5)} for k, v in sorted(emp_ecm_map.items())]
    pd.DataFrame(ecm_rows).to_csv(out_tables / "Table_S3m_Empirical_Axis_ECM_Density.csv", index=False, encoding="utf-8")

    ot_payload = json.loads((vault_ct / "opentargets_target_genetics_and_approval.json").read_text(encoding="utf-8"))
    full_drift_df = _assemble_full_target_resolved_drift_df(master_df, tau0_df, atlas_6ax, emp_ecm_map, ot_payload)

    sci_mask = ((full_drift_df["stop_deconvolution"] == "Completed") & (full_drift_df["max_phase_num"] >= 3)) | (
        full_drift_df["stop_deconvolution"].isin(["Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"])
    )
    sci_df = full_drift_df[sci_mask].copy()
    sci_df["clinical_success"] = ((sci_df["stop_deconvolution"] == "Completed") & (sci_df["max_phase_num"] >= 3)).astype(int)
    sci_6ax = sci_df[sci_df["is_core_6axis"]].copy()

    prog_mask = full_drift_df["stop_deconvolution"].isin(
        ["Completed", "Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"]
    )
    prog_df = full_drift_df[prog_mask].copy()
    prog_df["clinical_success"] = ((prog_df["stop_deconvolution"] == "Completed") & (prog_df["max_phase_num"] >= 3)).astype(int)

    abl_parts = []
    gen_parts = []
    summary_s3a_rows = []
    dnacc_oof_sci: np.ndarray | None = None
    tau0_oof_sci: np.ndarray | None = None
    dnacc_oof_prog: np.ndarray | None = None
    tau0_oof_prog: np.ndarray | None = None

    for c_name, c_df in [
        ("Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes", sci_df),
        ("Deconfounded_Scientific_Ph34_vs_SciFail_Core6SingleCellAxes", sci_6ax),
        ("Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials", prog_df),
    ]:
        abl_c = evaluate_six_step_dnacc_ablation(c_df, cohort_name=c_name, n_boot=2000, seed=42)
        abl_parts.append(abl_c)

        gen_c = evaluate_human_genetics_vs_dnacc_orthogonality(c_df, cohort_name=c_name, n_boot=2000, seed=42)
        gen_parts.append(gen_c)

        y_c = c_df["clinical_success"].to_numpy(dtype=int)
        X1_c = c_df[["tau0_healthy"]].to_numpy(dtype=float)
        X6_c = np.column_stack(
            [
                c_df["tau0_healthy"].to_numpy(dtype=float),
                c_df["vital_bystander_hazard"].to_numpy(dtype=float),
                c_df["tau_early"].to_numpy(dtype=float),
                c_df["delta_tau"].to_numpy(dtype=float),
                c_df["dynamic_retention_score"].to_numpy(dtype=float),
                c_df["compartment_switch"].to_numpy(dtype=float),
                c_df["stromal_ectopic_share_late"].to_numpy(dtype=float),
                c_df["barrier_compatibility"].to_numpy(dtype=float),
                c_df["tau_early"].to_numpy(dtype=float) * c_df["barrier_compatibility"].to_numpy(dtype=float),
            ]
        )
        p1_c = _oof_logistic_predictions(X1_c, y_c, seed=42)
        p6_c = _oof_logistic_predictions(X6_c, y_c, seed=42)

        if c_name == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes":
            dnacc_oof_sci, tau0_oof_sci = p6_c, p1_c
            c_out = c_df.copy()
            c_out["oof_dnacc_prob"] = np.round(p6_c, 5)
            c_out["oof_omniniche_prob"] = np.round(p6_c, 5)
            c_out["oof_static_tau0_prob"] = np.round(p1_c, 5)
            c_out.to_csv(out_tables / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv", index=False, encoding="utf-8")
            mod_strat_df = evaluate_within_modality_dnacc_stratification(c_df, p6_c, p1_c, cohort_name=c_name, n_boot=2000, seed=42)
            mod_strat_df.to_csv(out_tables / "Table_S3k_WithinModality_DNAcc_Stratified_Audit.csv", index=False, encoding="utf-8")
        elif c_name == "Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials":
            dnacc_oof_prog, tau0_oof_prog = p6_c, p1_c

        m1_row = abl_c[(abl_c["cv_scheme"] == "Stratified_5Fold") & (abl_c["step_id"] == "M1_Static_Healthy_Tau0")].iloc[0]
        m6_row = abl_c[(abl_c["cv_scheme"] == "Stratified_5Fold") & (abl_c["step_id"] == "M6_Full_DNAcc")].iloc[0]
        m6_grp = abl_c[(abl_c["cv_scheme"] == "UnseenTarget_Group5Fold") & (abl_c["step_id"] == "M6_Full_DNAcc")].iloc[0]
        st_grp_df = abl_c[abl_c["cv_scheme"] == "Subgroup_SingleTarget_UnseenTarget_Group5Fold"]
        st_grp = st_grp_df.iloc[0] if not st_grp_df.empty else m6_grp

        g1_row = gen_c[(gen_c["evaluation_stratum"] == "Stratified_5Fold") & (gen_c["model_id"] == "G1_OpenTargets_Human_Genetics_Only")].iloc[0]
        g5_row = gen_c[(gen_c["evaluation_stratum"] == "Stratified_5Fold") & (gen_c["model_id"] == "G5_Joint_Human_Genetics_Plus_DNAcc")].iloc[0]

        fe_c = fit_multivariable_fixed_effects_logistic(c_df, p6_c)
        or_u_dnacc = fe_c[fe_c["model_specification"] == "Unadjusted_Univariate_DNAcc"].iloc[0]
        or_u_tau0 = fe_c[fe_c["model_specification"] == "Unadjusted_Univariate_Static_Tau0"].iloc[0]
        or_adj_dnacc = fe_c[fe_c["model_specification"] == "Multivariable_FixedEffects_DNAcc"].iloc[0]

        q_o = pd.qcut(pd.Series(p6_c).rank(method="first"), q=4, labels=["Q1", "Q2", "Q3", "Q4"])
        q_s = pd.qcut(c_df["tau0_healthy"].rank(method="first"), q=4, labels=["Q1", "Q2", "Q3", "Q4"])
        rr_o = float(y_c[q_o == "Q4"].mean() / max(y_c[q_o == "Q1"].mean(), 1e-6))
        rr_s = float(y_c[q_s == "Q4"].mean() / max(y_c[q_s == "Q1"].mean(), 1e-6))

        summary_s3a_rows.append(
            {
                "benchmark_cohort": c_name,
                "n_trials": int(len(y_c)),
                "n_positive": int(np.sum(y_c)),
                "human_genetics_auroc": g1_row["auroc"],
                "human_genetics_ci_low": g1_row["auroc_ci_low"],
                "human_genetics_ci_high": g1_row["auroc_ci_high"],
                "static_tau0_auroc": m1_row["auroc"],
                "static_tau0_ci_low": m1_row["auroc_ci_low"],
                "static_tau0_ci_high": m1_row["auroc_ci_high"],
                "dnacc_auroc": m6_row["auroc"],
                "dnacc_ci_low": m6_row["auroc_ci_low"],
                "dnacc_ci_high": m6_row["auroc_ci_high"],
                "joint_genetics_dnacc_auroc": g5_row["auroc"],
                "joint_genetics_dnacc_ci_low": g5_row["auroc_ci_low"],
                "joint_genetics_dnacc_ci_high": g5_row["auroc_ci_high"],
                "unseen_target_groupkfold_auroc": m6_grp["auroc"],
                "unseen_target_ci_low": m6_grp["auroc_ci_low"],
                "unseen_target_ci_high": m6_grp["auroc_ci_high"],
                "single_target_groupkfold_auroc": st_grp["auroc"],
                "single_target_groupkfold_ci_low": st_grp["auroc_ci_low"],
                "single_target_groupkfold_ci_high": st_grp["auroc_ci_high"],
                "delta_auroc_vs_tau0": m6_row["delta_auroc_vs_m1"],
                "delong_p_value": m6_row["delong_p_vs_m1"],
                "delong_log10_p": m6_row["delong_log10p_vs_m1"],
                "static_tau0_auprc": m1_row["auprc"],
                "dnacc_auprc": m6_row["auprc"],
                "static_tau0_brier": m1_row["brier_score"],
                "dnacc_brier": m6_row["brier_score"],
                "static_tau0_or_per_sd": or_u_tau0["odds_ratio_per_sd"],
                "static_tau0_or_ci_low": or_u_tau0["ci95_low"],
                "static_tau0_or_ci_high": or_u_tau0["ci95_high"],
                "static_tau0_or_p": or_u_tau0["p_value"],
                "dnacc_or_per_sd": or_u_dnacc["odds_ratio_per_sd"],
                "dnacc_or_ci_low": or_u_dnacc["ci95_low"],
                "dnacc_or_ci_high": or_u_dnacc["ci95_high"],
                "dnacc_or_p": or_u_dnacc["p_value"],
                "dnacc_adjusted_or_per_sd": or_adj_dnacc["odds_ratio_per_sd"],
                "dnacc_adjusted_or_ci_low": or_adj_dnacc["ci95_low"],
                "dnacc_adjusted_or_ci_high": or_adj_dnacc["ci95_high"],
                "dnacc_adjusted_or_p": or_adj_dnacc["p_value"],
                "q4_vs_q1_ratio_static_tau0": round(rr_s, 4),
                "q4_vs_q1_ratio_dnacc": round(rr_o, 4),
            }
        )

    assert dnacc_oof_sci is not None and tau0_oof_sci is not None
    assert dnacc_oof_prog is not None and tau0_oof_prog is not None

    # Exact Wilson score 95% binomial CIs on quartile maturation & scientific success rates (T3)
    y_prog = prog_df["clinical_success"].to_numpy(dtype=int)
    y_sci = sci_df["clinical_success"].to_numpy(dtype=int)
    q_omni = pd.qcut(pd.Series(dnacc_oof_prog).rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
    q_stat = pd.qcut(prog_df["tau0_healthy"].rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
    q_omni_sci = pd.qcut(pd.Series(dnacc_oof_sci).rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
    q_stat_sci = pd.qcut(sci_df["tau0_healthy"].rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
    q_rows = []
    for ql in ["Q1_Low", "Q2", "Q3", "Q4_High"]:
        ms = q_stat == ql
        ns = int(np.sum(ms))
        ks = int(np.sum(y_prog[ms]))
        slo, shi = _wilson_ci95(ks, ns)

        mo = q_omni == ql
        no = int(np.sum(mo))
        ko = int(np.sum(y_prog[mo]))
        olo, ohi = _wilson_ci95(ko, no)

        mss = q_stat_sci == ql
        nss = int(np.sum(mss))
        kss = int(np.sum(y_sci[mss]))
        sslo, sshi = _wilson_ci95(kss, nss)

        mos = q_omni_sci == ql
        nos = int(np.sum(mos))
        kos = int(np.sum(y_sci[mos]))
        oslo, oshi = _wilson_ci95(kos, nos)

        q_rows.append(
            {
                "quartile": ql,
                "n_trials_static": ns,
                "static_tau0_maturation_pct": round(100.0 * ks / max(ns, 1), 3),
                "static_ci95_low": round(100.0 * slo, 3),
                "static_ci95_high": round(100.0 * shi, 3),
                "n_trials_dnacc": no,
                "dnacc_maturation_pct": round(100.0 * ko / max(no, 1), 3),
                "dnacc_ci95_low": round(100.0 * olo, 3),
                "dnacc_ci95_high": round(100.0 * ohi, 3),
                "n_sci_static": nss,
                "sci_static_tau0_success_pct": round(100.0 * kss / max(nss, 1), 3),
                "sci_static_ci95_low": round(100.0 * sslo, 3),
                "sci_static_ci95_high": round(100.0 * sshi, 3),
                "n_sci_dnacc": nos,
                "sci_dnacc_success_pct": round(100.0 * kos / max(nos, 1), 3),
                "sci_dnacc_ci95_low": round(100.0 * oslo, 3),
                "sci_dnacc_ci95_high": round(100.0 * oshi, 3),
            }
        )
    pd.DataFrame(q_rows).to_csv(out_tables / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv", index=False, encoding="utf-8")

    abl_all_df = pd.concat(abl_parts, ignore_index=True)
    abl_all_df.to_csv(out_tables / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv", index=False, encoding="utf-8")

    gen_all_df = pd.concat(gen_parts, ignore_index=True)
    gen_all_df.to_csv(out_tables / "Table_S3l_Genetics_vs_DNAcc_Orthogonality.csv", index=False, encoding="utf-8")
    print("    Open Targets Human Genetics vs Static tau_0 vs DNAcc (Primary Scientific Cohort):")
    print(
        gen_all_df[gen_all_df["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"][
            ["evaluation_stratum", "model_id", "n_trials", "auroc", "delta_auroc_vs_genetics", "delong_log10p_vs_genetics"]
        ].to_string(index=False),
        flush=True,
    )

    s3a_df = pd.DataFrame(summary_s3a_rows)
    s3a_df.to_csv(out_tables / "Table_S3a_PanDisease_Dynamic_Drift_Model_Comparison.csv", index=False, encoding="utf-8")
    print(
        s3a_df[
            [
                "benchmark_cohort",
                "n_trials",
                "human_genetics_auroc",
                "static_tau0_auroc",
                "dnacc_auroc",
                "joint_genetics_dnacc_auroc",
                "unseen_target_groupkfold_auroc",
                "single_target_groupkfold_auroc",
                "delong_log10_p",
            ]
        ].to_string(index=False),
        flush=True,
    )

    # Per-Axis Performance across BOTH Definitive_Maturation (prog_df) AND Deconfounded_Scientific (sci_df) row blocks (T9)
    per_ax_rows = []
    for block_name, block_df in [
        ("Definitive_Maturation", prog_df),
        ("Deconfounded_Scientific", sci_df),
    ]:
        for ax_name in [
            "Hepatic_MASH_Fibrosis",
            "Pulmonary_IPF_ARDS",
            "Renal_CKD_AKI",
            "Solid_Oncology",
            "Autoimmune_IBD_Rheum",
            "Infectious_Sepsis",
        ]:
            sub = block_df[block_df["disease_axis"] == ax_name].copy()
            ya = sub["clinical_success"].to_numpy(dtype=int)
            if len(ya) < 25 or len(np.unique(ya)) < 2:
                continue
            X1a = sub[["tau0_healthy"]].to_numpy(dtype=float)
            X6a = np.column_stack(
                [
                    sub["tau0_healthy"].to_numpy(dtype=float),
                    sub["vital_bystander_hazard"].to_numpy(dtype=float),
                    sub["tau_early"].to_numpy(dtype=float),
                    sub["delta_tau"].to_numpy(dtype=float),
                    sub["dynamic_retention_score"].to_numpy(dtype=float),
                    sub["compartment_switch"].to_numpy(dtype=float),
                    sub["stromal_ectopic_share_late"].to_numpy(dtype=float),
                    sub["barrier_compatibility"].to_numpy(dtype=float),
                    sub["tau_early"].to_numpy(dtype=float) * sub["barrier_compatibility"].to_numpy(dtype=float),
                ]
            )
            p1a = _oof_logistic_predictions(X1a, ya, seed=42)
            p6a = _oof_logistic_predictions(X6a, ya, seed=42)
            b1a = bootstrap_classification_metrics(ya, p1a, n_boot=2000, seed=42)
            b6a = bootstrap_classification_metrics(ya, p6a, n_boot=2000, seed=43)
            dla = fast_delong_paired_test(ya, p6a, p1a)
            fea = fit_multivariable_fixed_effects_logistic(sub, p6a)
            ora_tau = fea[fea["model_specification"] == "Unadjusted_Univariate_Static_Tau0"].iloc[0]
            ora_dna = fea[fea["model_specification"] == "Unadjusted_Univariate_DNAcc"].iloc[0]
            per_ax_rows.append(
                {
                    "cohort_block": block_name,
                    "disease_axis": ax_name,
                    "n_definitive_trials": int(len(ya)),
                    "n_phase34_completed": int(np.sum(ya)),
                    "static_tau0_auroc": round(b1a["auroc"], 4),
                    "static_tau0_ci_low": round(b1a["auroc_ci_low"], 4),
                    "static_tau0_ci_high": round(b1a["auroc_ci_high"], 4),
                    "dnacc_auroc": round(b6a["auroc"], 4),
                    "dnacc_ci_low": round(b6a["auroc_ci_low"], 4),
                    "dnacc_ci_high": round(b6a["auroc_ci_high"], 4),
                    "delta_auroc": round(dla["delta_auroc"], 4),
                    "delong_p_value": dla["p_value"],
                    "static_tau0_or_per_sd": ora_tau["odds_ratio_per_sd"],
                    "static_tau0_or_ci_low": ora_tau["ci95_low"],
                    "static_tau0_or_ci_high": ora_tau["ci95_high"],
                    "dnacc_or_per_sd": ora_dna["odds_ratio_per_sd"],
                    "dnacc_or_ci_low": ora_dna["ci95_low"],
                    "dnacc_or_ci_high": ora_dna["ci95_high"],
                    "dnacc_or_p_value": ora_dna["p_value"],
                }
            )
    pd.DataFrame(per_ax_rows).to_csv(out_tables / "Table_S3b_PerAxis_Dynamic_Drift_Performance.csv", index=False, encoding="utf-8")

    print("\n>>> [Step 3/7] Running Multivariable Fixed-Effects (+Human Genetics), 7-Regime Sensitivity (R1-R7), Calibration & DCA...", flush=True)
    fe_sci = fit_multivariable_fixed_effects_logistic(sci_df, dnacc_oof_sci)
    fe_sci["cohort"] = "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"
    fe_prog = fit_multivariable_fixed_effects_logistic(prog_df, dnacc_oof_prog)
    fe_prog["cohort"] = "Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials"
    fe_all = pd.concat([fe_sci, fe_prog], ignore_index=True)
    fe_all.to_csv(out_tables / "Table_S3g_Multivariable_FixedEffects_Control.csv", index=False, encoding="utf-8")

    sens_df = evaluate_unreported_stop_sensitivity_regimes(full_drift_df, n_boot=2000, seed=42)
    sens_df.to_csv(out_tables / "Table_S3e_Unreported_Stop_Sensitivity_Bounds.csv", index=False, encoding="utf-8")
    print("    7-Regime Sensitivity & Boundary Audit (R1-R7):")
    print(sens_df[["regime_id", "n_trials", "n_positive", "static_tau0_auroc", "dnacc_auroc", "delta_auroc", "delong_log10_p"]].to_string(index=False), flush=True)

    cal_prog = compute_calibration_deciles(prog_df["clinical_success"].to_numpy(dtype=int), dnacc_oof_prog, tau0_oof_prog, n_bins=10)
    cal_prog["cohort"] = f"Phase2_to_Phase34_Maturation_n{len(prog_df)}"
    cal_sci = compute_calibration_deciles(sci_df["clinical_success"].to_numpy(dtype=int), dnacc_oof_sci, tau0_oof_sci, n_bins=10)
    cal_sci["cohort"] = f"Deconfounded_Scientific_Outcomes_n{len(sci_df)}"
    pd.concat([cal_prog, cal_sci], ignore_index=True).to_csv(out_tables / "Table_S3h_Calibration_Deciles.csv", index=False, encoding="utf-8")

    ph3_enr = pd.to_numeric(master_df[master_df["max_phase_num"] == 3]["enrollment"], errors="coerce").dropna()
    med_ph3_n = int(round(float(np.median(ph3_enr[ph3_enr > 0])))) if len(ph3_enr) > 0 else 250
    dca_df, enrich_df = compute_decision_curve_and_enrichment(
        prog_df["clinical_success"].to_numpy(dtype=int),
        dnacc_oof_prog,
        tau0_oof_prog,
        median_enrollment=med_ph3_n,
    )
    dca_df.to_csv(out_tables / "Table_S3i_Decision_Curve_Analysis.csv", index=False, encoding="utf-8")
    enrich_df["median_phase3_enrollment_used"] = med_ph3_n
    enrich_df.to_csv(out_tables / "Table_S3j_Simulated_Trial_Enrichment.csv", index=False, encoding="utf-8")

    print("\n>>> [Step 4/7] Running Fibrosis 24-Trial Benchmark (+StageBlind/StagePrior), Data-Driven 3x3 Transfer & Profile Likelihood...", flush=True)
    verified_ncts = json.loads((vault_ct / "landmark_24_fibrosis_ncts_sha256.json").read_text(encoding="utf-8"))
    organ_data = {
        "Liver_MASH": {
            "donors": pd.read_csv(out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Liver.csv"),
            "early": pd.read_csv(vault_pb / "axis1_liver_fibrosis_early_comp.csv", index_col=0),
            "late": pd.read_csv(vault_pb / "axis1_liver_fibrosis_late_comp.csv", index_col=0),
            "drift": pd.read_csv(vault_pb / "axis1_liver_fibrosis_drift.csv", index_col=0),
        },
        "Lung_IPF": {
            "donors": pd.read_csv(out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Lung.csv"),
            "early": pd.read_csv(vault_pb / "axis2_lung_ipf_early_comp.csv", index_col=0),
            "late": pd.read_csv(vault_pb / "axis2_lung_ipf_late_comp.csv", index_col=0),
            "drift": pd.read_csv(vault_pb / "axis2_lung_ipf_drift.csv", index_col=0),
        },
        "Kidney_CKD": {
            "donors": pd.read_csv(out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Kidney.csv"),
            "early": pd.read_csv(vault_pb / "axis3_kidney_ckd_early_comp.csv", index_col=0),
            "late": pd.read_csv(vault_pb / "axis3_kidney_ckd_late_comp.csv", index_col=0),
            "drift": pd.read_csv(vault_pb / "axis3_kidney_ckd_drift.csv", index_col=0),
        },
    }
    organ_late_map = {k: v["late"] for k, v in organ_data.items()}
    ccr2_df = pd.read_csv(out_tables / "Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv")

    grammar_liv = build_conserved_fibrosis_grammar_from_liver(organ_data["Liver_MASH"]["donors"])
    grn_genes = sorted(
        (
            FIBROSIS_COMPARTMENT_GENES["Parenchyma_Injury"]
            | FIBROSIS_COMPARTMENT_GENES["SPP1_SAMac"]
            | FIBROSIS_COMPARTMENT_GENES["Myofibroblast"]
        )
        & set(organ_data["Liver_MASH"]["late"].index)
    )
    e_liv = organ_data["Liver_MASH"]["early"].loc[grn_genes].T
    l_liv = organ_data["Liver_MASH"]["late"].loc[grn_genes].T
    expr_liv = pd.concat([e_liv, l_liv, 0.5 * (e_liv + l_liv)], ignore_index=True)
    pathogenic_core = sorted(
        {"COL1A1", "COL1A2", "COL3A1", "POSTN", "CTHRC1", "ACTA2", "SPP1", "TREM2", "TGFB1", "GDF15", "CCL2"} & set(grn_genes)
    )
    perturber_liv = RidgeGRNPerturber(ridge_lambda=0.6, propagation_steps=3).fit(
        expr_liv, prior_edges=FIBROSIS_PRIOR_EDGES, pathogenic_program_genes=pathogenic_core
    )

    df_fib_tr, _ = evaluate_14_landmark_fibrosis_trials(
        grammar=grammar_liv,
        tau0_df=tau0_df,
        liver_drift_df=organ_data["Liver_MASH"]["drift"],
        lung_drift_df=organ_data["Lung_IPF"]["drift"],
        kidney_drift_df=organ_data["Kidney_CKD"]["drift"],
        verified_ncts=verified_ncts,
        grn_perturber=perturber_liv,
        ccr2_paradox_df=ccr2_df,
        organ_late_comp_map=organ_late_map,
    )
    df_fib_tr.to_csv(out_tables / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv", index=False, encoding="utf-8")

    fib_methods = [
        ("OmniNiche_v2_ZeroShot_Langevin", "omniniche_remission_prob", True),
        ("OmniNiche_StageBlind_FixedS0", "omniniche_stageblind_prob", True),
        ("B0_Enrolled_Stage_Prior_Only", "b0_enrolled_stage_only", True),
        ("B1_Static_Healthy_Tau0", "b1_static_tau0", True),
        ("B2_Polypharmacology_Target_Count", "b2_linear_target_count", False),
        ("B3_SingleCell_GRN_Retention_Only", "b3_single_cell_grn_only", True),
        ("B4_Static_LR_Compartment_Degree", "b4_static_lr_degree", False),
        ("B5_Uninformative_Modality_Prior", "b5_modality_prior_only", True),
    ]
    arena_v2_rows = []
    for bench_label, sub_f in [
        ("PanOrgan_Fibrosis_24_Phase2b3_Trials", df_fib_tr),
        ("PanOrgan_Fibrosis_Core14_Phase2b3_Trials", df_fib_tr[df_fib_tr["is_core14"]].copy()),
    ]:
        yf = sub_f["clinical_win"].to_numpy(dtype=int)
        p_tau_f = sub_f["b1_static_tau0"].to_numpy(dtype=float)
        p_stg_f = sub_f["b0_enrolled_stage_only"].to_numpy(dtype=float)
        for mname, col, is_p in fib_methods:
            raw_f = sub_f[col].to_numpy(dtype=float)
            prob_f = (
                np.clip(raw_f, 0.01, 0.99)
                if is_p
                else np.clip((raw_f - raw_f.min()) / (raw_f.max() - raw_f.min() + 1e-8), 0.05, 0.95)
            )
            bf = bootstrap_classification_metrics(yf, prob_f, n_boot=2000, seed=42)
            pf_tau = exact_permutation_auroc_pvalue(yf, prob_f, y_baseline=p_tau_f, n_perm=10000, seed=42)
            pf_stg = exact_permutation_auroc_pvalue(yf, prob_f, y_baseline=p_stg_f, n_perm=10000, seed=43)
            arena_v2_rows.append(
                {
                    "benchmark": bench_label,
                    "method": mname,
                    "n_trials": int(len(yf)),
                    "n_winners": int(np.sum(yf)),
                    "auroc": round(bf["auroc"], 4),
                    "auroc_ci_low": round(bf["auroc_ci_low"], 4),
                    "auroc_ci_high": round(bf["auroc_ci_high"], 4),
                    "auprc": round(bf["auprc"], 4),
                    "auprc_ci_low": round(bf["auprc_ci_low"], 4),
                    "auprc_ci_high": round(bf["auprc_ci_high"], 4),
                    "brier_score": round(bf["brier_score"], 4),
                    "brier_ci_low": round(bf["brier_ci_low"], 4),
                    "brier_ci_high": round(bf["brier_ci_high"], 4),
                    "exact_perm_p_vs_null": round(pf_tau["p_perm_vs_null"], 5),
                    "exact_perm_p_vs_tau0": round(pf_tau["p_perm_vs_baseline"], 5),
                    "exact_perm_p_vs_stage_prior": round(pf_stg["p_perm_vs_baseline"], 5),
                }
            )
    s4b_df = pd.DataFrame(arena_v2_rows)
    s4b_df.to_csv(out_tables / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv", index=False, encoding="utf-8")
    print("    Fibrosis 24-Trial & Core-14 Arena (with 95% CIs & 10,000-perm P):")
    print(
        s4b_df[s4b_df["benchmark"] == "PanOrgan_Fibrosis_24_Phase2b3_Trials"][
            ["method", "auroc", "auroc_ci_low", "auroc_ci_high", "brier_score", "exact_perm_p_vs_null", "exact_perm_p_vs_tau0"]
        ].to_string(index=False),
        flush=True,
    )

    bi_3x3_df = evaluate_bidirectional_3x3_fibrosis_transfer(
        organ_data=organ_data,
        tau0_df=tau0_df,
        verified_ncts=verified_ncts,
        prior_edges=FIBROSIS_PRIOR_EDGES,
        ccr2_paradox_df=ccr2_df,
    )
    bi_3x3_df.to_csv(out_tables / "Table_S4f_Bidirectional_3x3_CrossOrgan_Transfer.csv", index=False, encoding="utf-8")
    print("    Data-Driven Bidirectional 3x3 & Leave-One-Organ-Out Cross-Organ Transfer:")
    print(bi_3x3_df.to_string(index=False), flush=True)

    print("\n>>> [Step 5/7] Computing 25-Point Continuous Profile Likelihood, Hysteresis Curve & 6x6 Bliss Matrix...", flush=True)
    prof_curves_df = compute_continuous_profile_likelihood_curves(grammar_liv, n_grid=25)
    prof_curves_df.to_csv(out_tables / "Table_S4g_Profile_Likelihood_Curves.csv", index=False, encoding="utf-8")

    hyst_df = compute_continuous_fibrosis_hysteresis_curve(grammar_liv, n_points=25)
    hyst_df.to_csv(out_tables / "Fig4_Continuous_Hysteresis_Bifurcation_Curve.csv", index=False, encoding="utf-8")

    bliss_2d_df, fib_quiver_npz = compute_2d_fibrosis_bliss_matrix_and_landscape(grammar_liv)
    bliss_2d_df.to_csv(out_tables / "Table_S4e_Fibrosis_2D_Bliss_Matrix.csv", index=False, encoding="utf-8")
    np.savez_compressed(out_tables / "Fig4_Waddington_2D_Quiver_Fibrosis.npz", **fib_quiver_npz)

    print("\n>>> [Step 6/7] Fortifying 4-Cohort Solid-Tumor Anti-PD-1 Benchmarks (N=83 Pure Mono & N=107 All Pre-Treatment)...", flush=True)
    df_bcc = pd.read_csv(out_tables / "Table_S5c_GSE123813_PreAntiPD1_Patient_Predictions.csv")
    df_m91 = pd.read_csv(out_tables / "Table_S5j_GSE91061_Melanoma_PreNivo_Mono_Predictions.csv")
    df_m78 = pd.read_csv(out_tables / "Table_S5k_GSE78220_Melanoma_PrePembro_Mono_Predictions.csv")
    df_ns_bk = pd.read_csv(out_tables / "Table_S5d_GSE207422_NSCLC_PreBiopsy_Bulk_Predictions.csv")
    df_pool107 = pd.read_csv(out_tables / "Table_S5i_Pooled_PreTreatment_107Patients_Master.csv")
    df_pool83 = df_pool107[df_pool107["regimen_class"] == "AntiPD1_Monotherapy"].copy()
    df_ns_sc = pd.read_csv(out_tables / "Table_S5e_GSE207422_NSCLC_scRNA_Predictions.csv")

    icb_specs = [
        ("GSE123813_BCC_PreTreatment_AntiPD1_Mono_scRNA", df_bcc),
        ("GSE91061_Melanoma_PreTreatment_Nivolumab_Mono", df_m91),
        ("GSE78220_Melanoma_PreTreatment_Pembrolizumab_Mono", df_m78),
        ("Pooled_Pure_AntiPD1_Monotherapy_PreTreatment_n83", df_pool83),
        ("GSE207422_NSCLC_PreTreatment_ChemoAntiPD1_Biopsy", df_ns_bk),
        ("Pooled_All_Strictly_PreTreatment_Biopsies_n107", df_pool107),
        ("GSE207422_NSCLC_Resected_TME_scRNA", df_ns_sc),
    ]
    method_cols = [
        ("OmniNiche_v2_Spatial_Langevin", "simulated_icb_remission_prob"),
        ("OmniNiche_NoVisium_UniformPij", "no_visium_icb_remission_prob"),
        ("B1_Static_PDL1_CD274", "b1_static_pdl1_prob"),
        ("B2_Static_Healthy_Tau0", "b2_static_tau0_prob"),
        ("B3_SingleComp_TCell_GRN", "b3_single_tcell_grn_prob"),
        ("B4_Static_LR_Centrality", "b4_static_lr_centrality_prob"),
        ("B5_Uninformative_Linear_Mean", "b5_linear_mean_prob"),
        ("B6_Published_18Gene_Tcell_GEP", "b6_published_gep18_prob"),
        ("B7_Published_CAF_TIDE_Exclusion", "b7_published_caf_exclusion_prob"),
        ("B8_Algebraic_CD8_to_myCAF_Ratio", "b8_algebraic_cd8_mycaf_ratio_prob"),
    ]

    s5_rows = []
    for cname, cdf in icb_specs:
        yi = cdf["responder"].to_numpy(dtype=int)
        norm_cols = {
            col: np.clip(cdf[col].to_numpy(dtype=float), 1e-4, 1.0 - 1e-4)
            for _, col in method_cols
        }
        pdl1_ref = norm_cols["b1_static_pdl1_prob"]
        gep_ref = norm_cols["b6_published_gep18_prob"]

        for m_idx, (mlabel, col) in enumerate(method_cols):
            pi = norm_cols[col]
            bi = bootstrap_classification_metrics(yi, pi, n_boot=2000, seed=100 + m_idx)
            pm_pdl1 = exact_permutation_auroc_pvalue(yi, pi, y_baseline=pdl1_ref, n_perm=10000, seed=200 + m_idx)
            pm_gep = exact_permutation_auroc_pvalue(yi, pi, y_baseline=gep_ref, n_perm=10000, seed=300 + m_idx)
            s5_rows.append(
                {
                    "cohort": cname,
                    "n_patients": int(len(yi)),
                    "n_responders": int(np.sum(yi)),
                    "method": mlabel,
                    "auroc": round(bi["auroc"], 4),
                    "auroc_ci_low": round(bi["auroc_ci_low"], 4),
                    "auroc_ci_high": round(bi["auroc_ci_high"], 4),
                    "auprc": round(bi["auprc"], 4),
                    "auprc_ci_low": round(bi["auprc_ci_low"], 4),
                    "auprc_ci_high": round(bi["auprc_ci_high"], 4),
                    "brier_score": round(bi["brier_score"], 4),
                    "brier_ci_low": round(bi["brier_ci_low"], 4),
                    "brier_ci_high": round(bi["brier_ci_high"], 4),
                    "exact_perm_p_vs_null": round(pm_pdl1["p_perm_vs_null"], 5),
                    "exact_perm_p_vs_pdl1": round(pm_pdl1["p_perm_vs_baseline"], 5),
                    "exact_perm_p_vs_gep18": round(pm_gep["p_perm_vs_baseline"], 5),
                }
            )
    s5_df = pd.DataFrame(s5_rows)
    s5_df.to_csv(out_tables / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv", index=False, encoding="utf-8")
    print("    Solid-Tumor Anti-PD-1 Multi-Cohort Benchmarks (OmniNiche vs NoVisium vs Published GEP/PD-L1):")
    print(
        s5_df[
            s5_df["method"].isin(
                [
                    "OmniNiche_v2_Spatial_Langevin",
                    "OmniNiche_NoVisium_UniformPij",
                    "B6_Published_18Gene_Tcell_GEP",
                    "B1_Static_PDL1_CD274",
                ]
            )
        ][["cohort", "n_patients", "method", "auroc", "auroc_ci_low", "auroc_ci_high", "auprc", "brier_score", "exact_perm_p_vs_null"]].to_string(index=False),
        flush=True,
    )

    print("\n>>> [Step 7/7] Streaming 12 Human Visium Sections (21,796 Spots) Radial Profiles + Exact-UMI FACS Sensitivity + 2D Tumor Waddington...", flush=True)
    radial_12_df, contrasting_vis_df = extract_12section_visium_radial_and_contrasting_archetypes(vault_sc)
    radial_12_df.to_csv(out_tables / "Table_S5g_Visium_12Section_Radial_Profiles.csv", index=False, encoding="utf-8")
    contrasting_vis_df.to_csv(out_tables / "Fig5_Contrasting_Visium_Sections.csv", index=False, encoding="utf-8")

    rad_summary = (
        radial_12_df.groupby(["pitch_bin", "approx_distance_um"])[["Tumor_z", "myCAF_z", "TAM_z", "CD8_T_z"]]
        .agg(["mean", "sem"])
    )
    rad_summary.columns = [f"{c}_{st}" for c, st in rad_summary.columns]
    rad_summary = rad_summary.reset_index()
    rad_summary.to_csv(out_tables / "Table_S5g_Visium_Radial_Consensus_Summary.csv", index=False, encoding="utf-8")

    prox_4 = pd.read_csv(out_tables / "Table_S5b_Visium_Consensus_Proximity_Matrix.csv", index_col=0)
    prox_4_mapped = prox_4.rename(index={"CD8_T": "CD8_Exhaustion"}, columns={"CD8_T": "CD8_Exhaustion"})
    icb_grammar = build_default_solid_tumor_icb_grammar(prox_4_mapped)
    facs_sens_df = _compute_gse123813_facs_deconfounding_sensitivity(vault_sc, icb_grammar)
    facs_sens_df.to_csv(out_tables / "Table_S5h_FACS_Deconfounding_Sensitivity.csv", index=False, encoding="utf-8")

    tum_quiver_npz = compute_2d_tumor_waddington_surface_shift(prox_4_mapped)
    np.savez_compressed(out_tables / "Fig5_Waddington_2D_Quiver_SolidTumor.npz", **tum_quiver_npz)

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(
        f"\n=== STAGE II v2.2 STATISTICAL, GENETIC & BIOPHYSICAL FORTIFICATION COMPLETE (Peak Python Memory = {peak_bytes / 1e6:.1f} MB << 7,000 MB) ===",
        flush=True,
    )


if __name__ == "__main__":
    main()
