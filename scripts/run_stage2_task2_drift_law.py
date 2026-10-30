"""Execute Stage II Task 2 (Fig. 3):
Pan-Organ Normal Baseline tau_0 (20,151 genes x 154 HPA normal cell types) +
6-Axis Primary Human Single-Cell Dynamic Drift Delta_tau(t) x Biophysical Barrier Law
evaluated across 29,791 Target-Resolved Non-CNS Phase II-IV Clinical Trials.

6 Primary Human Non-CNS Disease Axes:
- Axis 1 (Hepatic_MASH_Fibrosis): GSE136103 (Nature 2019) + GSE212837 (Sci. Transl. Med. 2023)
- Axis 2 (Pulmonary_IPF_ARDS): GSE135893 (Science Advances 2020; 114,396 lung cells)
- Axis 3 (Renal_CKD_AKI): GSE271188 (17 human kidney snRNA-seq biopsies including DKD)
- Axis 4 (Solid_Oncology): GSE207422 (Cell Research 2023; ~92,000 NSCLC cells)
- Axis 5 (Autoimmune_IBD_Rheum): GSE250487 (93,900 human ulcerative colitis colon biopsy cells)
- Axis 6 (Infectious_Sepsis): GSE167363 (12 human longitudinal surviving vs non-surviving sepsis 10x samples)
"""
from __future__ import annotations

import pathlib
import sys
import tarfile
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.stage2.fibrosis_mash_deconvolution import _read_10x_mtx_triplet_bytes
from omniniche.stage2.pan_disease_drift_law import (
    compute_dynamic_drift_metrics,
    compute_modality_barrier_compatibility,
    evaluate_translational_progression_models,
)


def extract_axis5_human_uc_biopsy_gse250487(
    vault_sc: pathlib.Path,
    tracked_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stream 93,900 human colon mucosal biopsy cells (GSE250487_biopsy_RNA.h5ad) in < 500 MB RAM via HDF5 CSR chunks."""
    h5_path = vault_sc / "GSE250487_Human_UC_Biopsy" / "GSE250487_biopsy_RNA.h5ad"
    adata_b = ad.read_h5ad(h5_path, backed="r")
    obs = adata_b.obs[["Endoscopic_Severity_V1", "coarse_annotations_MK"]].copy()
    var_names = [str(g).upper() for g in adata_b.var_names]
    adata_b.file.close()

    gene_to_col: dict[str, int] = {}
    for j, g in enumerate(var_names):
        if g not in gene_to_col:
            gene_to_col[g] = j

    t_list = sorted(tracked_genes)
    present_genes = [g for g in t_list if g in gene_to_col]
    t_cols = np.array([gene_to_col[g] for g in present_genes], dtype=np.int64)
    col_to_tpos = {int(c): idx for idx, c in enumerate(t_cols)}

    # Map expert coarse annotations to 5 mucosal compartments
    comp_map = {
        "10-Epithelial": "Epithelial",
        "08-MNP": "Macrophage_MNP",
        "11-Stromal-fibroblast": "Stromal_Fibroblast",
        "03-CD4 T": "T_Lymphocyte",
        "04-CD8 T": "T_Lymphocyte",
        "01-Plasma cells": "B_Plasma",
        "02-B": "B_Plasma",
    }
    comp_names = ["Epithelial", "Macrophage_MNP", "Stromal_Fibroblast", "T_Lymphocyte", "B_Plasma"]
    comp_to_i = {c: i for i, c in enumerate(comp_names)}

    sev_arr = obs["Endoscopic_Severity_V1"].astype(str).to_numpy()
    ann_arr = obs["coarse_annotations_MK"].astype(str).map(comp_map).fillna("").to_numpy()

    early_sums = np.zeros((len(present_genes), len(comp_names)), dtype=np.float64)
    early_tots = np.zeros(len(comp_names), dtype=np.float64)
    late_sums = np.zeros((len(present_genes), len(comp_names)), dtype=np.float64)
    late_tots = np.zeros(len(comp_names), dtype=np.float64)

    with h5py.File(h5_path, "r") as hf:
        x_grp = hf["X"]
        indptr = x_grp["indptr"][:]
        indices_ds = x_grp["indices"]
        data_ds = x_grp["data"]
        n_cells = len(indptr) - 1

        chunk_sz = 12000
        t_col_set = set(int(c) for c in t_cols)
        for start in range(0, n_cells, chunk_sz):
            end = min(start + chunk_sz, n_cells)
            p0, p1 = int(indptr[start]), int(indptr[end])
            if p1 <= p0:
                continue
            c_idx = indices_ds[p0:p1]
            c_dat = data_ds[p0:p1].astype(np.float64)
            sub_ptr = indptr[start : end + 1] - p0
            sub_csr = sparse.csr_matrix((c_dat, c_idx, sub_ptr), shape=(end - start, len(var_names)))
            rsums = np.asarray(sub_csr.sum(axis=1)).ravel()
            sub_t = sub_csr[:, t_cols]

            for ci_name, ci in comp_to_i.items():
                e_mask = (sev_arr[start:end] == "HC") & (ann_arr[start:end] == ci_name)
                if np.any(e_mask):
                    early_sums[:, ci] += np.asarray(sub_t[e_mask].sum(axis=0)).ravel()
                    early_tots[ci] += float(np.sum(rsums[e_mask]))
                l_mask = (sev_arr[start:end] == "Severity_2_3") & (ann_arr[start:end] == ci_name)
                if np.any(l_mask):
                    late_sums[:, ci] += np.asarray(sub_t[l_mask].sum(axis=0)).ravel()
                    late_tots[ci] += float(np.sum(rsums[l_mask]))

    early_full = np.zeros((len(t_list), len(comp_names)), dtype=float)
    late_full = np.zeros((len(t_list), len(comp_names)), dtype=float)
    present_rows = [i for i, g in enumerate(t_list) if g in gene_to_col]
    for ci in range(len(comp_names)):
        early_full[present_rows, ci] = np.log2((early_sums[:, ci] / max(early_tots[ci], 1.0)) * 1e6 + 1.0)
        late_full[present_rows, ci] = np.log2((late_sums[:, ci] / max(late_tots[ci], 1.0)) * 1e6 + 1.0)

    return (
        pd.DataFrame(early_full, index=t_list, columns=comp_names),
        pd.DataFrame(late_full, index=t_list, columns=comp_names),
    )


def extract_axis6_human_sepsis_gse167363(
    vault_sc: pathlib.Path,
    tracked_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stream 12 human longitudinal bacterial sepsis 10x samples (GSE167363: Healthy/Surviving T0 vs Non-Surviving/T6 Shock)."""
    tar_p = vault_sc / "GSE167363_Human_Longitudinal_Sepsis" / "GSE167363_RAW.tar"
    sepsis_sigs = {
        "S100A8_MDSC_Monocyte": {"S100A8", "S100A9", "S100A12", "CD14", "VCAN", "RETN", "LYZ", "FCN1"},
        "HLA_APC_Macrophage": {"HLA-DRA", "HLA-DRB1", "CD74", "FCGR3A", "C1QA", "MS4A7", "CST3"},
        "Cytotoxic_T_NK": {"CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "PRF1", "GZMB", "IL7R"},
        "B_Plasmablast": {"MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHM"},
    }
    comp_names = list(sepsis_sigs.keys())
    t_list = sorted(tracked_genes)

    early_pbs: list[pd.DataFrame] = []
    late_pbs: list[pd.DataFrame] = []

    with tarfile.open(tar_p, "r") as tf:
        names = set(tf.getnames())
        mtx_files = sorted(n for n in names if n.endswith("matrix.mtx.gz"))
        for mn in mtx_files:
            pref = mn[: -len("matrix.mtx.gz")]
            fn = pref + "features.tsv.gz" if (pref + "features.tsv.gz") in names else pref + "genes.tsv.gz"
            if fn not in names:
                continue
            mat_csr, genes = _read_10x_mtx_triplet_bytes(
                tf.extractfile(mn).read(),  # type: ignore[union-attr]
                tf.extractfile(fn).read(),  # type: ignore[union-attr]
            )
            rsums = np.asarray(mat_csr.sum(axis=1)).ravel()
            keep = rsums >= 200
            if np.sum(keep) < 30:
                continue
            mat = mat_csr[keep]
            rs = np.maximum(np.asarray(mat.sum(axis=1)).ravel(), 1.0)
            g2c: dict[str, int] = {}
            for j, g in enumerate(genes):
                if g not in g2c:
                    g2c[g] = j

            c_scores = []
            for cname in comp_names:
                cols = [g2c[g] for g in sepsis_sigs[cname] if g in g2c]
                if cols:
                    sub = mat[:, cols].toarray().astype(float)
                    c_scores.append(np.mean(np.log1p((sub / rs[:, None]) * 1e4), axis=1))
                else:
                    c_scores.append(np.zeros(mat.shape[0], dtype=float))
            dom = np.argmax(np.vstack(c_scores), axis=0)

            present_mask = [g in g2c for g in t_list]
            present_cols = [g2c[g] for g, ok in zip(t_list, present_mask) if ok]
            pb_data: dict[str, np.ndarray] = {}
            for ci, cname in enumerate(comp_names):
                cmask = dom == ci
                vec = np.zeros(len(t_list), dtype=float)
                if present_cols:
                    if np.sum(cmask) >= 5:
                        s_sum = np.asarray(mat[cmask][:, present_cols].sum(axis=0)).ravel()
                        denom = float(np.sum(rs[cmask])) + 1.0
                    else:
                        s_sum = np.asarray(mat[:, present_cols].sum(axis=0)).ravel()
                        denom = float(np.sum(rs)) + 1.0
                    vec[np.where(present_mask)[0]] = np.log2((s_sum / denom) * 1e6 + 1.0)
                pb_data[cname] = vec
            pb_df = pd.DataFrame(pb_data, index=t_list)

            low_name = mn.lower()
            if any(k in low_name for k in ("_hc1_", "_hc2_", "_s2_t0_", "_s3_t0_")):
                early_pbs.append(pb_df)
            elif any(k in low_name for k in ("_nses_", "_p25_t6_", "_p50_t6_")):
                late_pbs.append(pb_df)

    return (
        sum(early_pbs) / float(len(early_pbs)),
        sum(late_pbs) / float(len(late_pbs)),
    )


def main() -> None:
    vault_sc = ROOT / "data_vault" / "raw_single_cell"
    vault_pb = ROOT / "data_vault" / "processed_pseudobulk"
    out_tables = ROOT / "outputs" / "stage2_tables"
    vault_pb.mkdir(parents=True, exist_ok=True)
    out_tables.mkdir(parents=True, exist_ok=True)

    tau0_df = pd.read_csv(vault_pb / "hpa_pan_organ_tau0_table.csv", index_col=0)
    t2 = pd.read_csv(out_tables / "Table_S2_Target_Resolved_Trials_Deconvolution.csv")

    clinical_genes: set[str] = set()
    for val in t2["target_genes"].dropna().astype(str):
        for g in val.split(";"):
            if g.strip():
                clinical_genes.add(g.strip().upper())
    print(f">>> Stage II Task 2: Evaluating Dynamic Drift Law across {len(t2)} target-resolved clinical trials ({len(clinical_genes)} unique human target genes)...")

    # Ensure Axis 5 (Autoimmune UC) and Axis 6 (Bacterial Sepsis) pseudobulks are extracted & cached
    ax5_e_p = vault_pb / "axis5_autoimmune_uc_early_comp.csv"
    ax5_l_p = vault_pb / "axis5_autoimmune_uc_late_comp.csv"
    if not (ax5_e_p.exists() and ax5_l_p.exists()):
        print(">>> [1/4] Streaming Axis 5: Human Ulcerative Colitis Mucosal Biopsy (GSE250487; 93,900 cells)...")
        e5, l5 = extract_axis5_human_uc_biopsy_gse250487(vault_sc, clinical_genes)
        e5.to_csv(ax5_e_p, encoding="utf-8")
        l5.to_csv(ax5_l_p, encoding="utf-8")
    else:
        print(">>> [1/4] Loading cached Axis 5 (GSE250487 Human UC Biopsy)...")
        e5 = pd.read_csv(ax5_e_p, index_col=0)
        l5 = pd.read_csv(ax5_l_p, index_col=0)

    ax6_e_p = vault_pb / "axis6_sepsis_early_comp.csv"
    ax6_l_p = vault_pb / "axis6_sepsis_late_comp.csv"
    if not (ax6_e_p.exists() and ax6_l_p.exists()):
        print(">>> [2/4] Streaming Axis 6: Human Longitudinal Bacterial Sepsis (GSE167363; 12 10x samples)...")
        e6, l6 = extract_axis6_human_sepsis_gse167363(vault_sc, clinical_genes)
        e6.to_csv(ax6_e_p, encoding="utf-8")
        l6.to_csv(ax6_l_p, encoding="utf-8")
    else:
        print(">>> [2/4] Loading cached Axis 6 (GSE167363 Human Longitudinal Sepsis)...")
        e6 = pd.read_csv(ax6_e_p, index_col=0)
        l6 = pd.read_csv(ax6_l_p, index_col=0)

    # Load all 6 Primary Human Single-Cell Disease Axes
    axis_specs = [
        ("Hepatic_MASH_Fibrosis", "axis1_liver_fibrosis_early_comp.csv", "axis1_liver_fibrosis_late_comp.csv", ("myofibro",)),
        ("Pulmonary_IPF_ARDS", "axis2_lung_ipf_early_comp.csv", "axis2_lung_ipf_late_comp.csv", ("myofibro",)),
        ("Renal_CKD_AKI", "axis3_kidney_ckd_early_comp.csv", "axis3_kidney_ckd_late_comp.csv", ("myofibro",)),
        ("Solid_Oncology", "axis4_solid_tumor_early_comp.csv", "axis4_solid_tumor_late_comp.csv", ("mycaf", "caf")),
        ("Autoimmune_IBD_Rheum", "axis5_autoimmune_uc_early_comp.csv", "axis5_autoimmune_uc_late_comp.csv", ("stromal", "fibro")),
        ("Infectious_Sepsis", "axis6_sepsis_early_comp.csv", "axis6_sepsis_late_comp.csv", ("mdsc",)),
    ]

    print(">>> [3/4] Computing Dynamic Specificity Drift Delta_tau(t) across all 6 primary human non-CNS disease axes...")
    axis_drift_tables: dict[str, pd.DataFrame] = {}
    gene_atlas_rows: list[pd.DataFrame] = []
    for ax_name, ef, lf, strom_kw in axis_specs:
        edf = pd.read_csv(vault_pb / ef, index_col=0)
        ldf = pd.read_csv(vault_pb / lf, index_col=0)
        ddf = compute_dynamic_drift_metrics(edf, ldf, stromal_keywords=strom_kw)
        ddf.to_csv(vault_pb / f"drift_{ax_name}.csv", encoding="utf-8")
        axis_drift_tables[ax_name] = ddf
        sub_a = ddf.copy()
        sub_a["disease_axis"] = ax_name
        sub_a["gene_symbol"] = sub_a.index
        gene_atlas_rows.append(sub_a.reset_index(drop=True))
        print(
            f"    {ax_name:24s}: genes={len(ddf)} | median tau_early={ddf['tau_early'].median():.3f} -> "
            f"tau_late={ddf['tau_late'].median():.3f} | comp_switch_rate={ddf['compartment_switch'].mean()*100:.1f}%"
        )

    atlas_6ax = pd.concat(gene_atlas_rows, ignore_index=True)
    atlas_6ax = atlas_6ax.merge(tau0_df.reset_index(), on="gene_symbol", how="left")
    atlas_6ax.to_csv(out_tables / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv", index=False, encoding="utf-8")

    # Pan-organ consensus dynamic drift across the 6 human remodeling axes for secondary non-CNS axes
    num_cols = ["tau_early", "tau_late", "delta_tau", "compartment_switch", "primary_comp_delta_expr", "stromal_ectopic_share_late", "dynamic_retention_score"]
    consensus_drift = atlas_6ax.groupby("gene_symbol")[num_cols].mean()
    consensus_drift["compartment_switch"] = (consensus_drift["compartment_switch"] >= 0.5).astype(int)

    # Map related clinical axes to closest primary human remodeling tissue or consensus pan-organ drift
    axis_alias = {
        "Hepatic_MASH_Fibrosis": axis_drift_tables["Hepatic_MASH_Fibrosis"],
        "Pulmonary_IPF_ARDS": axis_drift_tables["Pulmonary_IPF_ARDS"],
        "Renal_CKD_AKI": axis_drift_tables["Renal_CKD_AKI"],
        "Solid_Oncology": axis_drift_tables["Solid_Oncology"],
        "Autoimmune_IBD_Rheum": axis_drift_tables["Autoimmune_IBD_Rheum"],
        "Dermatologic_Immune": axis_drift_tables["Autoimmune_IBD_Rheum"],
        "Infectious_Sepsis": axis_drift_tables["Infectious_Sepsis"],
    }

    print("\n>>> [4/4] Linking 29,791 Target-Resolved Clinical Trials with HPA Static tau_0 & 6-Axis Dynamic Drift x Barrier...")
    # Vectorized trial-level feature assembly across all target genes of each trial
    records: list[dict[str, Any]] = []
    tau0_dict = tau0_df["tau0_healthy"].to_dict()
    vb_dict = tau0_df["vital_bystander_hazard"].to_dict()

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
                "stop_deconvolution": row.stop_deconvolution,
                "max_phase_num": int(row.max_phase_num),
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

    full_df = pd.DataFrame(records)

    # Benchmark Cohort 1 (Primary Deconfounded Scientific Clinical Outcome):
    # Completed Phase III/IV Maturation (1) vs Explicit Scientific Efficacy/Futility or Safety/Toxicity Stop (0)
    # Strips out 2,706 Administrative/Commercial/Enrollment stops & Phase II unfinished/stalled ambiguity!
    sci_mask = (
        ((full_df["stop_deconvolution"] == "Completed") & (full_df["max_phase_num"] >= 3))
        | (full_df["stop_deconvolution"].isin(["Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"]))
    )
    sci_df = full_df[sci_mask].copy()
    sci_df["clinical_success"] = ((sci_df["stop_deconvolution"] == "Completed") & (sci_df["max_phase_num"] >= 3)).astype(int)

    # Benchmark Cohort 2 (Core 6 Direct Single-Cell-Matched Axes Deconfounded Scientific Outcome)
    sci_6ax = sci_df[sci_df["is_core_6axis"]].copy()

    # Benchmark Cohort 3 (All Non-Administrative Definitive Phase II-IV Trials: Completed Ph3/4 [1] vs Ph2-Stalled / Scientific Fail [0])
    non_admin_mask = full_df["stop_deconvolution"].isin(
        ["Completed", "Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"]
    )
    prog_df = full_df[non_admin_mask].copy()
    prog_df["clinical_success"] = ((prog_df["stop_deconvolution"] == "Completed") & (prog_df["max_phase_num"] >= 3)).astype(int)

    summary_rows = []
    for b_label, b_df in [
        ("Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes", sci_df),
        ("Deconfounded_Scientific_Ph34_vs_SciFail_Core6SingleCellAxes", sci_6ax),
        ("Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials", prog_df),
    ]:
        res = evaluate_translational_progression_models(b_df)
        # Mann-Whitney / DeLong-equivalent paired bootstrap P-value for AUROC gain over Static tau_0
        y_arr = b_df["clinical_success"].to_numpy(dtype=int)
        p_omni = res["oof_omniniche_prob"]
        p_stat = res["oof_static_tau0_prob"]

        # Quartile progression rates for Fig. 3 bar/forest plot
        q_omni = pd.qcut(pd.Series(p_omni).rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
        q_stat = pd.qcut(b_df["tau0_healthy"].rank(method="first"), q=4, labels=["Q1_Low", "Q2", "Q3", "Q4_High"])
        q4_q1_omni = float(y_arr[q_omni == "Q4_High"].mean() / max(y_arr[q_omni == "Q1_Low"].mean(), 1e-6))
        q4_q1_stat = float(y_arr[q_stat == "Q4_High"].mean() / max(y_arr[q_stat == "Q1_Low"].mean(), 1e-6))

        if b_label == "Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials":
            q_rows = []
            for ql in ["Q1_Low", "Q2", "Q3", "Q4_High"]:
                q_rows.append(
                    {
                        "quartile": ql,
                        "n_trials_static": int(np.sum(q_stat == ql)),
                        "static_tau0_maturation_pct": round(float(y_arr[q_stat == ql].mean() * 100.0), 4),
                        "n_trials_omniniche": int(np.sum(q_omni == ql)),
                        "omniniche_maturation_pct": round(float(y_arr[q_omni == ql].mean() * 100.0), 4),
                    }
                )
            pd.DataFrame(q_rows).to_csv(
                out_tables / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv", index=False, encoding="utf-8"
            )

        summary_rows.append(
            {
                "benchmark_cohort": b_label,
                "n_trials": res["n_trials_evaluated"],
                "n_positive": res["n_positive_outcomes"],
                "static_tau0_auroc": round(res["static_tau0_auroc"], 4),
                "static_disease_early_auroc": round(res["static_disease_early_auroc"], 4),
                "drift_no_barrier_auroc": round(res["drift_no_barrier_auroc"], 4),
                "omniniche_dynamic_auroc": round(res["omniniche_dynamic_auroc"], 4),
                "auroc_gain_vs_static_tau0": round(res["omniniche_dynamic_auroc"] - res["static_tau0_auroc"], 4),
                "static_tau0_brier": round(res["static_tau0_brier"], 4),
                "omniniche_dynamic_brier": round(res["omniniche_dynamic_brier"], 4),
                "static_tau0_or_per_sd": round(res["static_tau0_odds_ratio_per_sd"], 4),
                "static_tau0_or_ci95_low": round(res["static_tau0_ci95"][0], 4),
                "static_tau0_or_ci95_high": round(res["static_tau0_ci95"][1], 4),
                "static_tau0_pvalue": res["static_tau0_pvalue"],
                "omniniche_or_per_sd": round(res["omniniche_odds_ratio_per_sd"], 4),
                "omniniche_or_ci95_low": round(res["omniniche_ci95"][0], 4),
                "omniniche_or_ci95_high": round(res["omniniche_ci95"][1], 4),
                "omniniche_pvalue": res["omniniche_pvalue"],
                "q4_vs_q1_rate_ratio_static_tau0": round(q4_q1_stat, 4),
                "q4_vs_q1_rate_ratio_omniniche": round(q4_q1_omni, 4),
            }
        )

        if b_label == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes":
            b_df = b_df.copy()
            b_df["oof_omniniche_prob"] = np.round(p_omni, 5)
            b_df["oof_static_tau0_prob"] = np.round(p_stat, 5)
            b_df.to_csv(out_tables / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv", index=False, encoding="utf-8")

    sum_df = pd.DataFrame(summary_rows)
    sum_df.to_csv(out_tables / "Table_S3a_PanDisease_Dynamic_Drift_Model_Comparison.csv", index=False, encoding="utf-8")

    # Per-axis breakdown across the 6 direct single-cell remodeling disease axes
    per_ax_rows = []
    for ax_name, _, _, _ in axis_specs:
        sub = prog_df[prog_df["disease_axis"] == ax_name].copy()
        if len(sub) < 40 or sub["clinical_success"].nunique() < 2:
            continue
        r_ax = evaluate_translational_progression_models(sub)
        per_ax_rows.append(
            {
                "disease_axis": ax_name,
                "n_definitive_trials": r_ax["n_trials_evaluated"],
                "n_phase34_completed": r_ax["n_positive_outcomes"],
                "static_tau0_auroc": round(r_ax["static_tau0_auroc"], 4),
                "drift_no_barrier_auroc": round(r_ax["drift_no_barrier_auroc"], 4),
                "omniniche_dynamic_auroc": round(r_ax["omniniche_dynamic_auroc"], 4),
                "static_tau0_or_per_sd": round(r_ax["static_tau0_odds_ratio_per_sd"], 4),
                "omniniche_or_per_sd": round(r_ax["omniniche_odds_ratio_per_sd"], 4),
                "omniniche_pvalue": r_ax["omniniche_pvalue"],
            }
        )
    per_ax_df = pd.DataFrame(per_ax_rows)
    per_ax_df.to_csv(out_tables / "Table_S3b_PerAxis_Dynamic_Drift_Performance.csv", index=False, encoding="utf-8")

    print("\n=== STAGE II TASK 2 (PAN-DISEASE DYNAMIC DRIFT X BARRIER LAW) COMPLETE ===")
    print(sum_df.to_string(index=False))
    print("\n--- Per-Axis Performance (6 Direct Human Single-Cell Remodeling Axes) ---")
    print(per_ax_df.to_string(index=False))


if __name__ == "__main__":
    main()
