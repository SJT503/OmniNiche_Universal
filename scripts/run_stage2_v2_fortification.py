"""Execute Stage II v2.0 Complete Statistical, Causal & Biophysical Fortification (Phase 1 & Phase 2: G01-G16).

Runs deterministically (< 2 GB peak RAM) on SHA-256-locked empirical caches in F:/OmniNiche_Universal,
exporting all fortified CSV/NPZ tables required for Nature Medicine / Science v2.0 figures and manuscript.
"""
from __future__ import annotations

import json
import math
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
from omniniche.stage2.pan_disease_drift_law import compute_modality_barrier_compatibility
from omniniche.stage2.solid_tumor_spatial_icb import (
    build_default_solid_tumor_icb_grammar,
    score_baseline_patients_via_langevin_icb,
)
from omniniche.stage2.v2_statistical_fortification import (
    _oof_logistic_predictions,
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
    evaluate_six_step_dnacc_ablation,
    evaluate_unreported_stop_sensitivity_regimes,
    evaluate_within_modality_dnacc_stratification,
    exact_permutation_auroc_pvalue,
    extract_12section_visium_radial_and_contrasting_archetypes,
    fast_delong_paired_test,
    fit_multivariable_fixed_effects_logistic,
)
from scripts.run_stage2_task3_fibrosis_mash import FIBROSIS_PRIOR_EDGES
from scripts.run_stage2_task4_solid_tumor_icb import ALL_PANEL_GENES, PANEL_GENES, _robust_01, fast_stream_tsv_gz_gene_subset


def _assemble_full_target_resolved_drift_df(
    master_df: pd.DataFrame,
    tau0_df: pd.DataFrame,
    atlas_6ax: pd.DataFrame,
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

        bar_v = compute_modality_barrier_compatibility(str(row.drug_modalities), ax, st_v, sw_v)
        records.append(
            {
                "nct_id": row.nct_id,
                "disease_axis": ax,
                "is_core_6axis": ax in axis_drift_tables,
                "overall_status": row.overall_status,
                "why_stopped": row.why_stopped,
                "stop_deconvolution": row.stop_deconvolution,
                "max_phase_num": int(row.max_phase_num),
                "has_results": row.has_results,
                "enrollment": row.enrollment,
                "start_year": row.start_year,
                "sponsor_class": row.sponsor_class,
                "primary_target_gene": row.primary_target_gene,
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
            }
        )
    return pd.DataFrame(records)


def _compute_gse123813_facs_deconfounding_sensitivity(
    sc_dir: pathlib.Path,
    grammar: Any,
) -> pd.DataFrame:
    """Compare naive CD3+-diluted pooling vs non-CD3+ FACS-deconfounded TME evaluation on GSE123813 BCC (n=8)."""
    cdir = sc_dir / "GSE123813_BCC_SCC_AntiPD1"
    clin = pd.read_excel(cdir / "Yost_2019_NatMed_SuppTables.xlsx", sheet_name="SuppTable1", skiprows=3).dropna(subset=["Patient"])
    clin["Patient"] = clin["Patient"].astype(str).str.strip()
    bcc_clin = clin[clin["Tumor Type"] == "BCC"].set_index("Patient")

    bcc_meta = pd.read_csv(cdir / "GSE123813_bcc_all_metadata.txt.gz", sep="\t")
    bcc_meta["cell.id"] = bcc_meta["cell.id"].astype(str)
    barcodes, gene_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE123813_bcc_scRNA_counts.txt.gz", ALL_PANEL_GENES, umi_stride=20
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

        # Naive all-cell CAF fraction vs FACS-deconfounded non-T CAF fraction
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
        tmp["Tumor"] = _robust_01(pdf[f"{pref}_tum"].to_numpy(), 0.16, 0.86)
        tmp["TAM"] = _robust_01(pdf[f"{pref}_tam"].to_numpy(), 0.14, 0.88)
        tmp["myCAF"] = _robust_01(pdf[f"{pref}_caf"].to_numpy(), 0.14, 0.88)
        tmp["CD8_Exhaustion"] = _robust_01(pdf["raw_exh_net"].to_numpy(), 0.16, 0.88)
        tmp["CD8_Effector_Reserve"] = _robust_01(pdf["raw_eff"].to_numpy(), 0.16, 0.90)
        sim = score_baseline_patients_via_langevin_icb(tmp, grammar=grammar, n_ensemble=250, seed=101)
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
    resume_step5 = "--resume-from-step5" in sys.argv
    resume_step3g = "--resume-from-step3g" in sys.argv
    tracemalloc.start()
    vault_ct = ROOT / "data_vault" / "raw_clinical_trials"
    vault_pb = ROOT / "data_vault" / "processed_pseudobulk"
    vault_sc = ROOT / "data_vault" / "raw_single_cell"
    out_tables = ROOT / "outputs" / "stage2_tables"
    out_tables.mkdir(parents=True, exist_ok=True)

    if resume_step5:
        print(">>> [Resume Step 5-7] Loading Liver fibrosis grammar & completing Steps 5-7...", flush=True)
        liv_donors = pd.read_csv(out_tables / "Table_S4a_PanOrgan_Fibrosis_Donor_States_Liver.csv")
        grammar_liv = build_conserved_fibrosis_grammar_from_liver(liv_donors)
        _run_steps_5_to_7(grammar_liv, vault_sc, out_tables)
        return

    master_df = pd.read_csv(vault_ct / "non_cns_phase234_trials_master.csv", low_memory=False)
    tau0_df = pd.read_csv(vault_pb / "hpa_pan_organ_tau0_table.csv", index_col=0)
    atlas_6ax = pd.read_csv(out_tables / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv")

    if not resume_step3g:
        print(">>> [Phase 1 & 2 | Step 1/7] De-scarring Table S1 & Auditing Macro-Clinical Stoppage Deconvolution...", flush=True)
        old_s1 = out_tables / "Table_S1_NonCNS_ClinicalTrials_PRISMA_Summary.csv"
        new_s1 = out_tables / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv"
        if old_s1.exists():
            s1_df = pd.read_csv(old_s1).rename(columns={"unique_non_cns_trials": "unique_systemic_organ_trials"})
            s1_df.to_csv(new_s1, index=False, encoding="utf-8")
            old_s1.unlink()
        elif new_s1.exists():
            s1_df = pd.read_csv(new_s1)
        print(f"    Saved clean Table S1 ({len(s1_df)} systemic/parenchymal organ axes; total trials = {len(master_df):,}).", flush=True)

        audit_metrics_df, audit_cm_df, audit_sample_df = audit_stoppage_classifier_benchmark(master_df, n_per_stratum=60, seed=42)
        audit_metrics_df.to_csv(out_tables / "Table_S2a_Stoppage_Auditor_Validation.csv", index=False, encoding="utf-8")
        audit_cm_df.to_csv(out_tables / "Table_S2a_Stoppage_Confusion_Matrix.csv", encoding="utf-8")
        audit_sample_df.to_csv(out_tables / "Table_S2a_Stoppage_Audit_240_Sample.csv", index=False, encoding="utf-8")
        print("    Stoppage Validation (n=240):")
        print(audit_metrics_df.to_string(index=False), flush=True)

        trend_df, odds_df, enr_dist_df = compute_macro_longitudinal_and_modality_odds(master_df)
        trend_df.to_csv(out_tables / "Table_S2b_Longitudinal_2005_2025_Trends.csv", index=False, encoding="utf-8")
        odds_df.to_csv(out_tables / "Table_S2c_Modality_Termination_Odds.csv", index=False, encoding="utf-8")
        enr_dist_df.to_csv(out_tables / "Table_S2d_Phase_Enrollment_Distributions.csv", index=False, encoding="utf-8")

        print("\n>>> [Step 2/7] Running 6-Step Orthogonal DNAcc Ablation, Unseen-Target GroupKFold & DeLong Tests...", flush=True)
        full_drift_df = _assemble_full_target_resolved_drift_df(master_df, tau0_df, atlas_6ax)
    else:
        print(">>> [Resume Mode] Recomputing fixed-effects Table S3g/S3a & completing Steps 4-7...", flush=True)
        full_drift_df = _assemble_full_target_resolved_drift_df(master_df, tau0_df, atlas_6ax)

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

    if resume_step3g:
        s3a_df = pd.read_csv(out_tables / "Table_S3a_PanDisease_Dynamic_Drift_Model_Comparison.csv")
        fe_parts = []
        for idx_c, (c_name, c_df) in enumerate([
            ("Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes", sci_df),
            ("Deconfounded_Scientific_Ph34_vs_SciFail_Core6SingleCellAxes", sci_6ax),
            ("Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials", prog_df),
        ]):
            y_c = c_df["clinical_success"].to_numpy(dtype=int)
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
            p6_c = _oof_logistic_predictions(X6_c, y_c, seed=42)
            fe_c = fit_multivariable_fixed_effects_logistic(c_df, p6_c)
            fe_c["cohort"] = c_name
            if c_name != "Deconfounded_Scientific_Ph34_vs_SciFail_Core6SingleCellAxes":
                fe_parts.append(fe_c)
            or_adj_dnacc = fe_c[fe_c["model_specification"] == "Multivariable_FixedEffects_DNAcc"].iloc[0]
            mask_row = s3a_df["benchmark_cohort"] == c_name
            s3a_df.loc[mask_row, "dnacc_adjusted_or_per_sd"] = or_adj_dnacc["odds_ratio_per_sd"]
            s3a_df.loc[mask_row, "dnacc_adjusted_or_ci_low"] = or_adj_dnacc["ci95_low"]
            s3a_df.loc[mask_row, "dnacc_adjusted_or_ci_high"] = or_adj_dnacc["ci95_high"]
            s3a_df.loc[mask_row, "dnacc_adjusted_or_p"] = or_adj_dnacc["p_value"]
        s3a_df.to_csv(out_tables / "Table_S3a_PanDisease_Dynamic_Drift_Model_Comparison.csv", index=False, encoding="utf-8")
        fe_all = pd.concat(fe_parts, ignore_index=True)
        fe_all.to_csv(out_tables / "Table_S3g_Multivariable_FixedEffects_Control.csv", index=False, encoding="utf-8")
        print("    Updated Table S3a & Table S3g (5 pre-outcome baseline confounders, collider-free):", flush=True)
        print(fe_all.to_string(index=False), flush=True)
    else:
        abl_parts = []
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

            # Extract M1 and M6 OOF probabilities under Stratified 5-fold CV for calibration, DCA, and Table S3a
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
                q_omni = pd.qcut(pd.Series(p6_c).rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
                q_stat = pd.qcut(c_df["tau0_healthy"].rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
                q_rows = []
                for ql in ["Q1_Low", "Q2", "Q3", "Q4_High"]:
                    ns = int(np.sum(q_stat == ql))
                    ps = float(y_c[q_stat == ql].mean())
                    ses = math.sqrt(max(ps * (1.0 - ps) / max(ns, 1), 1e-9))
                    no = int(np.sum(q_omni == ql))
                    po = float(y_c[q_omni == ql].mean())
                    seo = math.sqrt(max(po * (1.0 - po) / max(no, 1), 1e-9))
                    q_rows.append(
                        {
                            "quartile": ql,
                            "n_trials_static": ns,
                            "static_tau0_maturation_pct": round(100.0 * ps, 3),
                            "static_ci95_low": round(100.0 * max(0.0, ps - 1.96 * ses), 3),
                            "static_ci95_high": round(100.0 * min(1.0, ps + 1.96 * ses), 3),
                            "n_trials_dnacc": no,
                            "dnacc_maturation_pct": round(100.0 * po, 3),
                            "dnacc_ci95_low": round(100.0 * max(0.0, po - 1.96 * seo), 3),
                            "dnacc_ci95_high": round(100.0 * min(1.0, po + 1.96 * seo), 3),
                        }
                    )
                pd.DataFrame(q_rows).to_csv(out_tables / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv", index=False, encoding="utf-8")

            m1_row = abl_c[(abl_c["cv_scheme"] == "Stratified_5Fold") & (abl_c["step_id"] == "M1_Static_Healthy_Tau0")].iloc[0]
            m6_row = abl_c[(abl_c["cv_scheme"] == "Stratified_5Fold") & (abl_c["step_id"] == "M6_Full_DNAcc")].iloc[0]
            m6_grp = abl_c[(abl_c["cv_scheme"] == "UnseenTarget_Group5Fold") & (abl_c["step_id"] == "M6_Full_DNAcc")].iloc[0]
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
                    "static_tau0_auroc": m1_row["auroc"],
                    "static_tau0_ci_low": m1_row["auroc_ci_low"],
                    "static_tau0_ci_high": m1_row["auroc_ci_high"],
                    "dnacc_auroc": m6_row["auroc"],
                    "dnacc_ci_low": m6_row["auroc_ci_low"],
                    "dnacc_ci_high": m6_row["auroc_ci_high"],
                    "unseen_target_groupkfold_auroc": m6_grp["auroc"],
                    "unseen_target_ci_low": m6_grp["auroc_ci_low"],
                    "unseen_target_ci_high": m6_grp["auroc_ci_high"],
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

        abl_all_df = pd.concat(abl_parts, ignore_index=True)
        abl_all_df.to_csv(out_tables / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv", index=False, encoding="utf-8")
        s3a_df = pd.DataFrame(summary_s3a_rows)
        s3a_df.to_csv(out_tables / "Table_S3a_PanDisease_Dynamic_Drift_Model_Comparison.csv", index=False, encoding="utf-8")
        print(s3a_df[["benchmark_cohort", "n_trials", "static_tau0_auroc", "dnacc_auroc", "unseen_target_groupkfold_auroc", "delong_log10_p", "dnacc_adjusted_or_per_sd"]].to_string(index=False))

        # Per-Axis Performance with 2,000-bootstrap CIs, OR 95% CIs, and DeLong P-values
        per_ax_rows = []
        for ax_name in [
            "Hepatic_MASH_Fibrosis",
            "Pulmonary_IPF_ARDS",
            "Renal_CKD_AKI",
            "Solid_Oncology",
            "Autoimmune_IBD_Rheum",
            "Infectious_Sepsis",
        ]:
            sub = prog_df[prog_df["disease_axis"] == ax_name].copy()
            ya = sub["clinical_success"].to_numpy(dtype=int)
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

        print("\n>>> [Step 3/7] Running Multivariable Fixed-Effects, 5-Regime Sensitivity Bounds, Calibration & DCA...")
        assert dnacc_oof_sci is not None and tau0_oof_sci is not None
        assert dnacc_oof_prog is not None and tau0_oof_prog is not None
        fe_sci = fit_multivariable_fixed_effects_logistic(sci_df, dnacc_oof_sci)
        fe_sci["cohort"] = "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"
        fe_prog = fit_multivariable_fixed_effects_logistic(prog_df, dnacc_oof_prog)
        fe_prog["cohort"] = "Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials"
        fe_all = pd.concat([fe_sci, fe_prog], ignore_index=True)
        fe_all.to_csv(out_tables / "Table_S3g_Multivariable_FixedEffects_Control.csv", index=False, encoding="utf-8")

        sens_df = evaluate_unreported_stop_sensitivity_regimes(full_drift_df, n_boot=2000, seed=42)
        sens_df.to_csv(out_tables / "Table_S3e_Unreported_Stop_Sensitivity_Bounds.csv", index=False, encoding="utf-8")
        print("    5-Regime Unreported/Negation Sensitivity Bounds:")
        print(sens_df[["regime_id", "n_trials", "static_tau0_auroc", "dnacc_auroc", "delta_auroc", "delong_log10_p"]].to_string(index=False))

        # Calibration, DCA & Trial Portfolio Enrichment on Definitive Phase II->III/IV Maturation (n=23,967) and Primary Scientific (n=11,843)
        cal_prog = compute_calibration_deciles(prog_df["clinical_success"].to_numpy(dtype=int), dnacc_oof_prog, tau0_oof_prog, n_bins=10)
        cal_prog["cohort"] = "Phase2_to_Phase34_Maturation_n23967"
        cal_sci = compute_calibration_deciles(sci_df["clinical_success"].to_numpy(dtype=int), dnacc_oof_sci, tau0_oof_sci, n_bins=10)
        cal_sci["cohort"] = "Deconfounded_Scientific_Outcomes_n11843"
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
        print("    Trial Portfolio Enrichment (median Phase III enrollment =", med_ph3_n, "patients):")
        print(enrich_df.to_string(index=False))

    print("\n>>> [Step 4/7] Running Fibrosis 14-Trial Benchmark (Semaglutide=Peptide), 3x3 Cross-Organ Transfer & Profile Likelihood...")
    verified_ncts = json.loads((vault_ct / "landmark_14_fibrosis_ncts_sha256.json").read_text(encoding="utf-8"))
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
    pathogenic_core = sorted({"COL1A1", "COL1A2", "COL3A1", "POSTN", "CTHRC1", "ACTA2", "SPP1", "TREM2", "TGFB1", "GDF15", "CCL2"} & set(grn_genes))
    perturber_liv = RidgeGRNPerturber(ridge_lambda=0.6, propagation_steps=3).fit(
        expr_liv, prior_edges=FIBROSIS_PRIOR_EDGES, pathogenic_program_genes=pathogenic_core
    )

    df_fib_tr, df_fib_arena = evaluate_14_landmark_fibrosis_trials(
        grammar=grammar_liv,
        tau0_df=tau0_df,
        liver_drift_df=organ_data["Liver_MASH"]["drift"],
        lung_drift_df=organ_data["Lung_IPF"]["drift"],
        kidney_drift_df=organ_data["Kidney_CKD"]["drift"],
        verified_ncts=verified_ncts,
        grn_perturber=perturber_liv,
        ccr2_paradox_df=ccr2_df,
    )
    df_fib_tr.to_csv(out_tables / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv", index=False, encoding="utf-8")

    # Add 2,000-bootstrap 95% CIs + 10,000-permutation exact P-values to Table S4b
    yf = df_fib_tr["clinical_win"].to_numpy(dtype=int)
    p_tau_f = df_fib_tr["b1_static_tau0"].to_numpy(dtype=float)
    arena_v2_rows = []
    for mname, col, is_p in [
        ("OmniNiche_v2_ZeroShot_Langevin", "omniniche_remission_prob", True),
        ("B1_Static_Healthy_Tau0", "b1_static_tau0", True),
        ("B2_Polypharmacology_Target_Count", "b2_linear_target_count", False),
        ("B3_SingleCell_GRN_Retention_Only", "b3_single_cell_grn_only", True),
        ("B4_Static_LR_Compartment_Degree", "b4_static_lr_degree", False),
        ("B5_Uninformative_Modality_Prior", "b5_modality_prior_only", True),
    ]:
        raw_f = df_fib_tr[col].to_numpy(dtype=float)
        prob_f = np.clip(raw_f, 0.01, 0.99) if is_p else np.clip((raw_f - raw_f.min()) / (raw_f.max() - raw_f.min() + 1e-8), 0.05, 0.95)
        bf = bootstrap_classification_metrics(yf, prob_f, n_boot=2000, seed=42)
        pf = exact_permutation_auroc_pvalue(yf, prob_f, y_baseline=p_tau_f, n_perm=10000, seed=42)
        arena_v2_rows.append(
            {
                "benchmark": "PanOrgan_Fibrosis_14_Phase2b3_Trials",
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
                "exact_perm_p_vs_null": round(pf["p_perm_vs_null"], 5),
                "exact_perm_p_vs_tau0": round(pf["p_perm_vs_baseline"], 5),
            }
        )
    s4b_df = pd.DataFrame(arena_v2_rows)
    s4b_df.to_csv(out_tables / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv", index=False, encoding="utf-8")
    print("    Fibrosis 14-Trial Arena (with 95% CIs & 10,000-perm P):")
    print(s4b_df[["method", "auroc", "auroc_ci_low", "auroc_ci_high", "brier_score", "exact_perm_p_vs_null"]].to_string(index=False))

    bi_3x3_df = evaluate_bidirectional_3x3_fibrosis_transfer(
        organ_data=organ_data,
        tau0_df=tau0_df,
        verified_ncts=verified_ncts,
        prior_edges=FIBROSIS_PRIOR_EDGES,
        ccr2_paradox_df=ccr2_df,
    )
    bi_3x3_df.to_csv(out_tables / "Table_S4f_Bidirectional_3x3_CrossOrgan_Transfer.csv", index=False, encoding="utf-8")
    print("    Bidirectional 3x3 & Leave-One-Organ-Out Cross-Organ Transfer:")
    print(bi_3x3_df.to_string(index=False))

    _run_steps_5_to_7(grammar_liv, vault_sc, out_tables)


def _run_steps_5_to_7(grammar_liv: Any, vault_sc: pathlib.Path, out_tables: pathlib.Path) -> None:
    print("\n>>> [Step 5/7] Computing 25-Point Continuous Profile Likelihood, Hysteresis Curve & 6x6 Bliss Matrix...", flush=True)
    prof_curves_df = compute_continuous_profile_likelihood_curves(grammar_liv, n_grid=25)
    prof_curves_df.to_csv(out_tables / "Table_S4g_Profile_Likelihood_Curves.csv", index=False, encoding="utf-8")

    hyst_df = compute_continuous_fibrosis_hysteresis_curve(grammar_liv, n_points=25)
    hyst_df.to_csv(out_tables / "Fig4_Continuous_Hysteresis_Bifurcation_Curve.csv", index=False, encoding="utf-8")

    bliss_2d_df, fib_quiver_npz = compute_2d_fibrosis_bliss_matrix_and_landscape(grammar_liv)
    bliss_2d_df.to_csv(out_tables / "Table_S4e_Fibrosis_2D_Bliss_Matrix.csv", index=False, encoding="utf-8")
    np.savez_compressed(out_tables / "Fig4_Waddington_2D_Quiver_Fibrosis.npz", **fib_quiver_npz)

    print("\n>>> [Step 6/7] Fortifying Solid-Tumor ICB Benchmarks (CIs + Exact Permutations) & Pooled n=32 Waterfall...")
    df_bcc = pd.read_csv(out_tables / "Table_S5c_GSE123813_PreAntiPD1_Patient_Predictions.csv")
    df_ns_bk = pd.read_csv(out_tables / "Table_S5d_GSE207422_NSCLC_PreBiopsy_Bulk_Predictions.csv")
    df_ns_sc = pd.read_csv(out_tables / "Table_S5e_GSE207422_NSCLC_scRNA_Predictions.csv")
    df_pool32 = pd.concat(
        [
            df_bcc.assign(cancer_cohort="BCC_scRNA_GSE123813"),
            df_ns_bk.assign(cancer_cohort="NSCLC_PreBiopsy_GSE207422"),
        ],
        ignore_index=True,
    )
    df_pool32.to_csv(out_tables / "Table_S5i_Pooled32_PreTreatment_Waterfall_and_Rescue.csv", index=False, encoding="utf-8")

    icb_specs = [
        ("GSE123813_BCC_PreTreatment_FullTME_scRNA", df_bcc),
        ("GSE207422_NSCLC_PreTreatment_Biopsy_Bulk", df_ns_bk),
        ("GSE207422_NSCLC_TME_Validation_scRNA", df_ns_sc),
        ("Pooled_Strictly_PreTreatment_Biopsies_n32", df_pool32),
    ]
    method_cols = [
        ("OmniNiche_v2_Spatial_Langevin", "simulated_icb_remission_prob", True),
        ("B1_Static_PDL1_CD274", "b1_static_pdl1", False),
        ("B2_Static_Healthy_Tau0", "b2_static_tau0", False),
        ("B3_SingleComp_TCell_GRN", "b3_single_tcell_grn", False),
        ("B4_Static_LR_Centrality", "b4_static_lr_centrality", False),
        ("B5_Uninformative_Linear_Mean", "b5_linear_mean", False),
    ]

    s5_rows = []
    for cname, cdf in icb_specs:
        yi = cdf["responder"].to_numpy(dtype=int)
        norm_cols: dict[str, np.ndarray] = {}
        for mlabel, col, is_p in method_cols:
            v = cdf[col].to_numpy(dtype=float)
            if is_p:
                norm_cols[col] = np.clip(v, 0.02, 0.98)
            else:
                mn, mx = float(np.min(v)), float(np.max(v))
                norm_cols[col] = 0.08 + 0.84 * (v - mn) / (mx - mn + 1e-8)
        base_ref = norm_cols["b1_static_pdl1"]

        for m_idx, (mlabel, col, _) in enumerate(method_cols):
            pi = norm_cols[col]
            bi = bootstrap_classification_metrics(yi, pi, n_boot=2000, seed=100 + m_idx)
            pm = exact_permutation_auroc_pvalue(yi, pi, y_baseline=base_ref, n_perm=10000, seed=200 + m_idx)
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
                    "exact_perm_p_vs_null": round(pm["p_perm_vs_null"], 5),
                    "exact_perm_p_vs_pdl1": round(pm["p_perm_vs_baseline"], 5),
                }
            )
    s5_df = pd.DataFrame(s5_rows)
    s5_df.to_csv(out_tables / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv", index=False, encoding="utf-8")
    print("    Solid-Tumor ICB Benchmarks (OmniNiche rows):")
    print(s5_df[s5_df["method"] == "OmniNiche_v2_Spatial_Langevin"][["cohort", "n_patients", "auroc", "auroc_ci_low", "auroc_ci_high", "brier_score", "exact_perm_p_vs_null"]].to_string(index=False))

    print("\n>>> [Step 7/7] Streaming 12 Human Visium Sections (21,796 Spots) Radial Profiles + FACS Sensitivity + 2D Tumor Waddington...")
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
    print("    12-Section Visium Radial Consensus Profile (0 to 600 um):")
    print(rad_summary.round(4).to_string(index=False))

    prox_4 = pd.read_csv(out_tables / "Table_S5b_Visium_Consensus_Proximity_Matrix.csv", index_col=0)
    prox_4_mapped = prox_4.rename(index={"CD8_T": "CD8_Exhaustion"}, columns={"CD8_T": "CD8_Exhaustion"})
    icb_grammar = build_default_solid_tumor_icb_grammar(prox_4_mapped)
    facs_sens_df = _compute_gse123813_facs_deconfounding_sensitivity(vault_sc, icb_grammar)
    facs_sens_df.to_csv(out_tables / "Table_S5h_FACS_Deconfounding_Sensitivity.csv", index=False, encoding="utf-8")

    tum_quiver_npz = compute_2d_tumor_waddington_surface_shift(prox_4_mapped)
    np.savez_compressed(out_tables / "Fig5_Waddington_2D_Quiver_SolidTumor.npz", **tum_quiver_npz)

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"\n=== PHASE 1 & PHASE 2 STATISTICAL/BIOPHYSICAL FORTIFICATION COMPLETE (Peak Python Memory = {peak_bytes / 1e6:.1f} MB << 7,000 MB) ===")


if __name__ == "__main__":
    main()
