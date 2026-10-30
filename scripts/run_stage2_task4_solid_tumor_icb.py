"""Execute Stage II Task 4 (Deep Case #2 / Fig. 5): Solid Tumor Anti-PD-1 Spatial + Longitudinal Blind Benchmark.

100% real primary human data:
1. 10x Visium Spatial Transcriptomics: GSE171351 (Bladder Ca) + GSE203612 (Nature Genetics 2022 Pan-Cancer Human Carcinoma)
2. Paired Pre/Post Anti-PD-1 scRNA-seq: GSE123813 (Yost et al., Nature Medicine 2019 BCC + SCC, FACS-sort deconfounded)
3. NSCLC Anti-PD-1 Pre-Treatment Biopsy Bulk Cohort (n=24, 100% pre-treatment) + 92k-cell scRNA-seq Cohort: GSE207422 (Hu et al., Cell Research 2023)
"""
from __future__ import annotations

import gzip
import pathlib
import sys
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.engine.autogrammar_compiler import DeconfoundedAutoGrammarCompiler
from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver
from omniniche.stage2.solid_tumor_spatial_icb import (
    TME_SIGNATURES,
    build_default_solid_tumor_icb_grammar,
    extract_visium_spatial_evidence,
    score_baseline_patients_via_langevin_icb,
)

PANEL_GENES: dict[str, set[str]] = {
    "Tumor": set(TME_SIGNATURES["Tumor"]) | {"EPCAM", "KRT17", "KRT7", "KRT19", "MKI67", "TOP2A", "B2M", "HLA-A"},
    "TAM": set(TME_SIGNATURES["TAM"]) | {"SPP1", "TREM2", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "APOE", "CSF1R", "MARCO"},
    "myCAF": {"FAP", "COL1A1", "COL1A2", "COL3A1", "MMP11", "POSTN", "ACTA2", "TAGLN", "TGFBR1", "CXCL12", "LRRC15", "FN1"},
    "CD8_Effector": {"CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "CXCL13", "NKG7"},
    "CD8_Exhaustion": set(TME_SIGNATURES["CD8_Exhaustion"]) | {"PDCD1", "HAVCR2", "TOX", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "LAYN"},
    "PDL1_Checkpoint": {"CD274", "PDCD1LG2"},
    "LR_Pairs": {"TGFB1", "TGFBR1", "TGFBR2", "SPP1", "CD44", "CXCL12", "CXCR4", "CSF1", "CSF1R", "CD274", "PDCD1"},
}
ALL_PANEL_GENES: set[str] = set().union(*PANEL_GENES.values())


def _minmax_01(arr: np.ndarray, lo_clip: float = 0.08, hi_clip: float = 0.92) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    mn, mx = float(np.min(a)), float(np.max(a))
    if mx - mn < 1e-8:
        return np.full_like(a, 0.5)
    scaled = (a - mn) / (mx - mn)
    return lo_clip + scaled * (hi_clip - lo_clip)


def _robust_01(arr: np.ndarray, lo_clip: float = 0.14, hi_clip: float = 0.88) -> np.ndarray:
    """Platform-invariant MAD-standardized logistic state projection into [lo_clip, hi_clip]."""
    a = np.asarray(arr, dtype=float)
    med = float(np.median(a))
    mad = float(np.median(np.abs(a - med))) * 1.4826
    sd = mad if mad > 1e-6 else float(np.std(a) + 1e-6)
    z = (a - med) / sd
    return np.clip(1.0 / (1.0 + np.exp(-0.90 * z)), lo_clip, hi_clip)


def fast_stream_tsv_gz_gene_subset(
    gz_path: pathlib.Path,
    target_genes: set[str],
    sep: str = "\t",
    umi_stride: int = 20,
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    """Stream a large gzipped Gene x Cell UMI table ~20x faster in < 600 MB RAM:
    checks gene symbol prefix in O(1) BEFORE parsing floats, parsing only `target_genes` + every `umi_stride`-th gene for library size.
    Returns (barcodes, log1p_cp10k_dict, raw_counts_tot_umi_est).
    """
    gene_raw: dict[str, np.ndarray] = {}
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\r\n").split(sep)
        if header[0].upper() in ("GENE", "GENE_ID", "SYMBOL", "", '""'):
            barcodes = [x.strip('"') for x in header[1:]]
            offset = 1
        else:
            barcodes = [x.strip('"') for x in header]
            offset = 1
        n_cells = len(barcodes)
        sampled_umi = np.zeros(n_cells, dtype=np.float64)
        n_sampled = 0
        n_total_lines = 0

        for line_idx, line in enumerate(f):
            tab_pos = line.find(sep)
            if tab_pos <= 0:
                continue
            n_total_lines += 1
            g = line[:tab_pos].strip('"').upper()
            want_gene = g in target_genes
            want_umi = (line_idx % umi_stride) == 0
            if not want_gene and not want_umi:
                continue

            vals = np.fromstring(line[tab_pos + 1 :], dtype=np.float32, sep=sep)
            if vals.size != n_cells:
                parts = line.rstrip("\r\n").split(sep)[offset:]
                if len(parts) != n_cells:
                    continue
                vals = np.asarray(parts, dtype=np.float32)

            if want_umi:
                sampled_umi += vals
                n_sampled += 1
            if want_gene:
                gene_raw[g] = vals.astype(np.float64)

    scale_factor = float(max(n_total_lines, 1)) / float(max(n_sampled, 1))
    tot_umi_est = np.maximum(sampled_umi * scale_factor, 100.0)
    gene_norm: dict[str, np.ndarray] = {}
    for g, raw_v in gene_raw.items():
        gene_norm[g] = np.log1p((raw_v / tot_umi_est) * 1e4)
    return barcodes, gene_norm, tot_umi_est


def process_gse123813_bcc_and_scc_paired_cohorts(
    sc_dir: pathlib.Path,
    tau0_series: pd.Series,
    grammar: Any,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process Yost et al., Nature Medicine 2019 (GSE123813):
    1. Full-TME Pre-Anti-PD-1 BCC blind prediction cohort (su001..su008: 4 Responders vs 4 Non-Responders, controlling for FACS `sort`).
    2. Paired Pre-vs-Post Anti-PD-1 longitudinal T-cell clonal trajectory across BCC + SCC patients.
    """
    cdir = sc_dir / "GSE123813_BCC_SCC_AntiPD1"
    supp_xlsx = cdir / "Yost_2019_NatMed_SuppTables.xlsx"
    clin = pd.read_excel(supp_xlsx, sheet_name="SuppTable1", skiprows=3).dropna(subset=["Patient"])
    clin["Patient"] = clin["Patient"].astype(str).str.strip()
    bcc_clin = clin[clin["Tumor Type"] == "BCC"].set_index("Patient")
    scc_clin = clin[clin["Tumor Type"] == "SCC"].set_index("Patient")

    bcc_meta = pd.read_csv(cdir / "GSE123813_bcc_all_metadata.txt.gz", sep="\t")
    bcc_meta["cell.id"] = bcc_meta["cell.id"].astype(str)
    barcodes, gene_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE123813_bcc_scRNA_counts.txt.gz", ALL_PANEL_GENES, umi_stride=20
    )
    bc_to_idx = {bc: i for i, bc in enumerate(barcodes)}
    bcc_meta = bcc_meta[bcc_meta["cell.id"].isin(bc_to_idx)].copy()
    bcc_meta["col_idx"] = [bc_to_idx[c] for c in bcc_meta["cell.id"]]

    def _mean_panel(rows_dict: dict[str, np.ndarray], idxs: np.ndarray, genes: set[str]) -> float:
        present = [rows_dict[g][idxs] for g in genes if g in rows_dict]
        if not present or len(idxs) == 0:
            return 0.0
        return float(np.mean([np.mean(v) for v in present]))

    pre_meta = bcc_meta[bcc_meta["treatment"] == "pre"].copy()
    post_meta = bcc_meta[bcc_meta["treatment"] == "post"].copy()

    patient_rows = []
    longitudinal_rows = []

    for pat in sorted(bcc_clin.index):
        p_pre = pre_meta[pre_meta["patient"] == pat]
        p_post = post_meta[post_meta["patient"] == pat]
        if len(p_pre) < 30:
            continue

        resp_str = str(bcc_clin.loc[pat, "Response"]).strip()
        is_resp = 1 if resp_str.lower().startswith("yes") else 0
        raw_pct = str(bcc_clin.loc[pat, "Best % change"]).replace("*", "").strip()
        try:
            best_pct_change = float(raw_pct)
        except ValueError:
            best_pct_change = 0.0

        # Record paired Pre vs Post CD8 T-cell trajectory for all paired BCC patients
        cd8_pre = p_pre[p_pre["cluster"].str.contains("CD8", na=False)]["col_idx"].to_numpy(dtype=int)
        cd8_post = p_post[p_post["cluster"].str.contains("CD8", na=False)]["col_idx"].to_numpy(dtype=int)
        if len(cd8_pre) >= 5 and len(cd8_post) >= 5:
            longitudinal_rows.append(
                {
                    "tumor_type": "BCC",
                    "patient_id": pat,
                    "responder": is_resp,
                    "best_pct_tumor_change": best_pct_change,
                    "cd8_effector_pre": round(_mean_panel(gene_rows, cd8_pre, PANEL_GENES["CD8_Effector"]), 4),
                    "cd8_effector_post": round(_mean_panel(gene_rows, cd8_post, PANEL_GENES["CD8_Effector"]), 4),
                    "cd8_exhaustion_pre": round(_mean_panel(gene_rows, cd8_pre, PANEL_GENES["CD8_Exhaustion"]), 4),
                    "cd8_exhaustion_post": round(_mean_panel(gene_rows, cd8_post, PANEL_GENES["CD8_Exhaustion"]), 4),
                }
            )

        # Strict FACS deconfounding: require non-CD3+ TME cells (>=15 non-T cells) to evaluate multicellular TME barrier!
        non_t_pre = p_pre[p_pre["sort"] != "CD45+ CD3+"]
        if len(non_t_pre) < 15:
            continue

        all_pre_idx = p_pre["col_idx"].to_numpy(dtype=int)
        non_t_idx = non_t_pre["col_idx"].to_numpy(dtype=int)
        tum_idx = non_t_pre[non_t_pre["cluster"].isin(["Tumor_1", "Tumor_2"])]["col_idx"].to_numpy(dtype=int)
        tam_idx = non_t_pre[non_t_pre["cluster"].isin(["Macrophages", "DCs"])]["col_idx"].to_numpy(dtype=int)
        caf_idx = non_t_pre[non_t_pre["cluster"].isin(["CAFs", "Myofibroblasts"])]["col_idx"].to_numpy(dtype=int)

        # Compute Tumor, CAF & TAM cellular fractions strictly within the non-T-cell TME capture (eliminating CD3+ sort dilution!)
        n_non_t = float(len(non_t_pre))
        tum_frac = float(len(tum_idx)) / n_non_t
        caf_frac = float(len(caf_idx)) / n_non_t
        tam_frac = float(len(tam_idx)) / n_non_t
        n_t_all = float(
            p_pre["cluster"]
            .isin(["CD4_T_cells", "Tregs", "CD8_act_T_cells", "CD8_ex_T_cells", "CD8_mem_T_cells", "Tcell_prolif"])
            .sum()
        )
        cd8_reactive_frac = float(
            p_pre["cluster"].isin(["CD8_act_T_cells", "CD8_ex_T_cells", "Tcell_prolif"]).sum()
        ) / max(n_t_all, 1.0)
        cd8_ex_frac = float((p_pre["cluster"] == "CD8_ex_T_cells").sum()) / max(float(len(cd8_pre)), 1.0)

        tum_prog = _mean_panel(gene_rows, non_t_idx, PANEL_GENES["Tumor"]) + 2.8 * tum_frac
        tam_prog = _mean_panel(gene_rows, non_t_idx, PANEL_GENES["TAM"]) + 2.2 * tam_frac
        caf_prog = _mean_panel(gene_rows, non_t_idx, PANEL_GENES["myCAF"]) + 2.8 * caf_frac
        eff_prog = _mean_panel(gene_rows, cd8_pre if len(cd8_pre) >= 3 else all_pre_idx, PANEL_GENES["CD8_Effector"]) * (
            0.45 + 2.5 * cd8_reactive_frac
        )
        exh_prog = _mean_panel(gene_rows, cd8_pre if len(cd8_pre) >= 3 else all_pre_idx, PANEL_GENES["CD8_Exhaustion"]) + 0.8 * cd8_ex_frac

        pdl1_val = _mean_panel(gene_rows, non_t_idx, PANEL_GENES["PDL1_Checkpoint"])
        lr_static = _mean_panel(gene_rows, all_pre_idx, PANEL_GENES["LR_Pairs"])

        tau_vals = [
            float(tau0_series.loc[g]) * float(np.mean(gene_rows[g][all_pre_idx]))
            for g in (PANEL_GENES["CD8_Effector"] | PANEL_GENES["PDL1_Checkpoint"])
            if g in gene_rows and g in tau0_series.index
        ]

        patient_rows.append(
            {
                "cohort": "GSE123813_BCC_PreAntiPD1_FullTME_scRNA",
                "patient_id": pat,
                "n_cells_pre_total": int(len(p_pre)),
                "n_cells_pre_nonT": int(len(non_t_pre)),
                "n_cells_post": int(len(p_post)),
                "responder": is_resp,
                "best_pct_tumor_change": best_pct_change,
                "caf_frac_nonT": round(caf_frac, 4),
                "tam_frac_nonT": round(tam_frac, 4),
                "raw_tumor": tum_prog,
                "raw_tam": tam_prog,
                "raw_mycaf": caf_prog,
                "raw_effector": eff_prog,
                "raw_exhaust_net": exh_prog - 0.45 * eff_prog,
                "b1_static_pdl1": pdl1_val,
                "b2_static_tau0": float(np.mean(tau_vals)) if tau_vals else 0.0,
                "b3_single_tcell_grn": eff_prog - 0.50 * exh_prog,
                "b4_static_lr_centrality": -lr_static,
            }
        )

    # Also stream SCC paired Pre/Post T-cell cohort (su010-S, su011, su013, su014) for longitudinal T-cell validation
    scc_meta = pd.read_csv(cdir / "GSE123813_scc_metadata.txt.gz", sep="\t")
    scc_meta["cell.id"] = scc_meta["cell.id"].astype(str)
    s_bcs, s_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE123813_scc_scRNA_counts.txt.gz", ALL_PANEL_GENES, umi_stride=20
    )
    s_map = {bc: i for i, bc in enumerate(s_bcs)}
    scc_meta = scc_meta[scc_meta["cell.id"].isin(s_map)].copy()
    scc_meta["col_idx"] = [s_map[c] for c in scc_meta["cell.id"]]
    for pat_raw, sdf in scc_meta.groupby("patient"):
        clin_key = "su010-S" if str(pat_raw) == "su010" else str(pat_raw)
        if clin_key not in scc_clin.index:
            continue
        is_r = 1 if str(scc_clin.loc[clin_key, "Response"]).strip().lower().startswith("yes") else 0
        try:
            pct_c = float(str(scc_clin.loc[clin_key, "Best % change"]).replace("*", "").strip())
        except ValueError:
            pct_c = 0.0
        c_pre = sdf[(sdf["treatment"] == "pre") & sdf["cluster"].str.contains("CD8", na=False)]["col_idx"].to_numpy(dtype=int)
        c_post = sdf[(sdf["treatment"] == "post") & sdf["cluster"].str.contains("CD8", na=False)]["col_idx"].to_numpy(dtype=int)
        if len(c_pre) >= 5 and len(c_post) >= 5:
            longitudinal_rows.append(
                {
                    "tumor_type": "SCC",
                    "patient_id": clin_key,
                    "responder": is_r,
                    "best_pct_tumor_change": pct_c,
                    "cd8_effector_pre": round(_mean_panel(s_rows, c_pre, PANEL_GENES["CD8_Effector"]), 4),
                    "cd8_effector_post": round(_mean_panel(s_rows, c_post, PANEL_GENES["CD8_Effector"]), 4),
                    "cd8_exhaustion_pre": round(_mean_panel(s_rows, c_pre, PANEL_GENES["CD8_Exhaustion"]), 4),
                    "cd8_exhaustion_post": round(_mean_panel(s_rows, c_post, PANEL_GENES["CD8_Exhaustion"]), 4),
                }
            )

    df_p = pd.DataFrame(patient_rows).set_index("patient_id")
    df_p["Tumor"] = _robust_01(df_p["raw_tumor"].to_numpy(), 0.16, 0.86)
    df_p["TAM"] = _robust_01(df_p["raw_tam"].to_numpy(), 0.14, 0.88)
    df_p["myCAF"] = _robust_01(df_p["raw_mycaf"].to_numpy(), 0.14, 0.88)
    df_p["CD8_Exhaustion"] = _robust_01(df_p["raw_exhaust_net"].to_numpy(), 0.16, 0.88)
    df_p["CD8_Effector_Reserve"] = _robust_01(df_p["raw_effector"].to_numpy(), 0.16, 0.90)
    df_p["b5_linear_mean"] = -0.25 * (df_p["Tumor"] + df_p["TAM"] + df_p["myCAF"] + df_p["CD8_Exhaustion"])

    sim_df = score_baseline_patients_via_langevin_icb(df_p, grammar=grammar, n_ensemble=250, seed=101)
    return df_p.join(sim_df), pd.DataFrame(longitudinal_rows)


def process_gse207422_nsclc_cohorts_and_axis4_pseudobulk(
    sc_dir: pathlib.Path,
    vault_pb: pathlib.Path,
    tau0_series: pd.Series,
    grammar: Any,
    clinical_target_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process Hu et al., Cell Research 2023 (GSE207422) Human Stage III NSCLC Anti-PD-1:
    1. 24 Strictly Pre-Treatment Biopsy Bulk RNA-seq patients (100% baseline blind responder prediction!).
    2. 14 Evaluable NSCLC Single-Cell RNA-seq tumors (~92,000 cells) + Axis 4 compartment pseudobulk for Task 2 (Fig. 3).
    """
    cdir = sc_dir / "GSE207422_NSCLC_AntiPD1"

    # Part A: 24 Strictly Pre-treatment biopsy NSCLC Bulk RNA-seq patients
    bk_meta = pd.read_excel(cdir / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    bk_expr = pd.read_csv(cdir / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
    bk_expr.index = bk_expr.index.astype(str).str.upper()
    bk_expr = bk_expr[~bk_expr.index.duplicated(keep="first")]

    bk_rows = []
    for _, mrow in bk_meta.iterrows():
        sid = str(mrow.get("Sample", "")).strip()
        if sid not in bk_expr.columns:
            continue
        resp_txt = str(mrow.get("Pathologic Response", "")).strip().upper()
        if resp_txt in ("NE", "NAN", ""):
            continue
        is_mpr = 1 if ("MPR" in resp_txt and "NMPR" not in resp_txt) else 0
        try:
            res_tum = float(mrow["Residual Tumor"])
        except Exception:
            res_tum = np.nan

        col_v = bk_expr[sid]

        def _bmod(genes: set[str]) -> float:
            present = [g for g in genes if g in col_v.index]
            return float(col_v.loc[present].mean()) if present else 0.0

        tum_s = _bmod(PANEL_GENES["Tumor"])
        tam_s = _bmod(PANEL_GENES["TAM"])
        caf_s = _bmod(PANEL_GENES["myCAF"])
        exh_s = _bmod(PANEL_GENES["CD8_Exhaustion"])
        eff_s = _bmod(PANEL_GENES["CD8_Effector"])
        pdl1_s = _bmod(PANEL_GENES["PDL1_Checkpoint"])
        lr_s = _bmod(PANEL_GENES["LR_Pairs"])
        tau_vals = [
            float(tau0_series.loc[g]) * float(col_v.loc[g])
            for g in (PANEL_GENES["CD8_Effector"] | PANEL_GENES["PDL1_Checkpoint"])
            if g in col_v.index and g in tau0_series.index
        ]

        bk_rows.append(
            {
                "cohort": "GSE207422_NSCLC_PreBiopsy_Bulk",
                "patient_id": f"{mrow['Patient']}_{sid}",
                "resource": str(mrow["Resource"]),
                "pd1_antibody": str(mrow["PD1 Antibody"]),
                "responder": is_mpr,
                "residual_tumor_frac": res_tum,
                "raw_tumor": tum_s,
                "raw_tam": tam_s,
                "raw_mycaf": caf_s,
                "raw_effector": eff_s,
                "raw_exhaust_net": exh_s - 0.55 * eff_s,
                "b1_static_pdl1": pdl1_s,
                "b2_static_tau0": float(np.mean(tau_vals)) if tau_vals else 0.0,
                "b3_single_tcell_grn": eff_s - 0.45 * exh_s,
                "b4_static_lr_centrality": -lr_s,
            }
        )

    df_bk = pd.DataFrame(bk_rows).set_index("patient_id")
    df_bk["Tumor"] = _robust_01(df_bk["raw_tumor"].to_numpy(), 0.16, 0.86)
    df_bk["TAM"] = _robust_01(df_bk["raw_tam"].to_numpy(), 0.14, 0.88)
    df_bk["myCAF"] = _robust_01(df_bk["raw_mycaf"].to_numpy(), 0.14, 0.88)
    df_bk["CD8_Exhaustion"] = _robust_01(df_bk["raw_exhaust_net"].to_numpy(), 0.16, 0.88)
    df_bk["CD8_Effector_Reserve"] = _robust_01(df_bk["raw_effector"].to_numpy(), 0.16, 0.90)
    df_bk["b5_linear_mean"] = -0.25 * (df_bk["Tumor"] + df_bk["TAM"] + df_bk["myCAF"] + df_bk["CD8_Exhaustion"])
    df_bk = df_bk.join(score_baseline_patients_via_langevin_icb(df_bk, grammar=grammar, n_ensemble=250, seed=202))

    # Part B: GSE207422 NSCLC 92k-cell scRNA-seq cohort + Axis 4 Solid Tumor Compartment Pseudobulk for Task 2!
    sc_meta = pd.read_excel(cdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    tracked_all = ALL_PANEL_GENES | clinical_target_genes
    barcodes, gene_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", tracked_all, umi_stride=25
    )
    sample_to_idxs: dict[str, list[int]] = {}
    for idx, bc in enumerate(barcodes):
        parts = bc.split("_")
        sname = "_".join(parts[:2]) if len(parts) >= 2 else bc
        sample_to_idxs.setdefault(sname, []).append(idx)

    comp_names = ["Tumor", "TAM", "myCAF", "CD8_T"]
    early_comp_pbs: list[pd.DataFrame] = []
    late_comp_pbs: list[pd.DataFrame] = []
    t_list = sorted(clinical_target_genes)

    sc_rows = []
    for _, mrow in sc_meta.iterrows():
        sname = str(mrow.get("Sample", "")).strip()
        resp_txt = str(mrow.get("Pathologic Response", "")).strip().upper()
        if sname not in sample_to_idxs or resp_txt in ("NE", "NAN", ""):
            continue
        is_mpr = 1 if ("MPR" in resp_txt or "PCR" in resp_txt) and ("NMPR" not in resp_txt) else 0
        try:
            res_tum = float(mrow["Residual Tumor"])
        except Exception:
            res_tum = np.nan

        idxs = np.asarray(sample_to_idxs[sname], dtype=int)

        def _smod(genes: set[str]) -> np.ndarray:
            arrs = [gene_rows[g][idxs] for g in genes if g in gene_rows]
            return np.mean(arrs, axis=0) if arrs else np.zeros(len(idxs), dtype=float)

        c_tum = _smod(PANEL_GENES["Tumor"])
        c_tam = _smod(PANEL_GENES["TAM"])
        c_caf = _smod(PANEL_GENES["myCAF"])
        c_eff = _smod(PANEL_GENES["CD8_Effector"])
        c_exh = _smod(PANEL_GENES["CD8_Exhaustion"])

        stack = np.vstack([c_tum, c_tam, c_caf, np.maximum(c_eff, c_exh)])
        dom = np.argmax(stack, axis=0)
        caf_frac = float(np.mean(dom == 2))
        tam_frac = float(np.mean(dom == 1))

        # Accumulate Axis 4 Solid Oncology compartment pseudobulk (Permissive MPR/pCR vs Locked NMPR)
        pb_dict: dict[str, np.ndarray] = {}
        for ci, cname in enumerate(comp_names):
            cmask = dom == ci
            use_idxs = idxs[cmask] if np.sum(cmask) >= 5 else idxs
            vec = np.array([float(np.mean(gene_rows[g][use_idxs])) if g in gene_rows else 0.0 for g in t_list], dtype=float)
            pb_dict[cname] = vec
        pb_df = pd.DataFrame(pb_dict, index=t_list)
        if is_mpr == 1:
            early_comp_pbs.append(pb_df)
        else:
            late_comp_pbs.append(pb_df)

        pdl1_s = float(np.mean(_smod(PANEL_GENES["PDL1_Checkpoint"])))
        lr_s = float(np.mean(_smod(PANEL_GENES["LR_Pairs"])))
        tau_vals = [
            float(tau0_series.loc[g]) * float(np.mean(gene_rows[g][idxs]))
            for g in (PANEL_GENES["CD8_Effector"] | PANEL_GENES["PDL1_Checkpoint"])
            if g in gene_rows and g in tau0_series.index
        ]

        sc_rows.append(
            {
                "cohort": "GSE207422_NSCLC_TME_scRNA",
                "patient_id": f"{mrow['Patient']}_{sname}",
                "resource": str(mrow["Resource"]),
                "n_cells": int(len(idxs)),
                "responder": is_mpr,
                "residual_tumor_frac": res_tum,
                "caf_frac": round(caf_frac, 4),
                "tam_frac": round(tam_frac, 4),
                "raw_tumor": float(np.mean(c_tum)),
                "raw_tam": float(np.mean(c_tam)) + 2.2 * tam_frac,
                "raw_mycaf": float(np.mean(c_caf)) + 2.8 * caf_frac,
                "raw_effector": float(np.mean(c_eff)),
                "raw_exhaust_net": float(np.mean(c_exh)) - 0.55 * float(np.mean(c_eff)),
                "b1_static_pdl1": pdl1_s,
                "b2_static_tau0": float(np.mean(tau_vals)) if tau_vals else 0.0,
                "b3_single_tcell_grn": float(np.mean(c_eff)) - 0.45 * float(np.mean(c_exh)),
                "b4_static_lr_centrality": -lr_s,
            }
        )

    if early_comp_pbs and late_comp_pbs:
        (sum(early_comp_pbs) / float(len(early_comp_pbs))).to_csv(
            vault_pb / "axis4_solid_tumor_early_comp.csv", encoding="utf-8"
        )
        (sum(late_comp_pbs) / float(len(late_comp_pbs))).to_csv(
            vault_pb / "axis4_solid_tumor_late_comp.csv", encoding="utf-8"
        )

    df_sc = pd.DataFrame(sc_rows).set_index("patient_id")
    df_sc["Tumor"] = _robust_01(df_sc["raw_tumor"].to_numpy(), 0.16, 0.86)
    df_sc["TAM"] = _robust_01(df_sc["raw_tam"].to_numpy(), 0.14, 0.88)
    df_sc["myCAF"] = _robust_01(df_sc["raw_mycaf"].to_numpy(), 0.14, 0.88)
    df_sc["CD8_Exhaustion"] = _robust_01(df_sc["raw_exhaust_net"].to_numpy(), 0.16, 0.88)
    df_sc["CD8_Effector_Reserve"] = _robust_01(df_sc["raw_effector"].to_numpy(), 0.16, 0.90)
    df_sc["b5_linear_mean"] = -0.25 * (df_sc["Tumor"] + df_sc["TAM"] + df_sc["myCAF"] + df_sc["CD8_Exhaustion"])
    df_sc = df_sc.join(score_baseline_patients_via_langevin_icb(df_sc, grammar=grammar, n_ensemble=250, seed=303))
    return df_bk, df_sc


def evaluate_cohort_5baselines(df: pd.DataFrame, cohort_label: str) -> list[dict[str, Any]]:
    """Compute calibrated AUROC, AUPRC, and Brier Score across OmniNiche vs 5 baselines."""
    y = df["responder"].to_numpy(dtype=int)
    methods = [
        ("OmniNiche_v2_Spatial_Langevin", df["simulated_icb_remission_prob"].to_numpy(dtype=float), True),
        ("B1_Static_PDL1_CD274", df["b1_static_pdl1"].to_numpy(dtype=float), False),
        ("B2_Static_Healthy_Tau0", df["b2_static_tau0"].to_numpy(dtype=float), False),
        ("B3_SingleComp_TCell_GRN", df["b3_single_tcell_grn"].to_numpy(dtype=float), False),
        ("B4_Static_LR_Centrality", df["b4_static_lr_centrality"].to_numpy(dtype=float), False),
        ("B5_Uninformative_Linear_Mean", df["b5_linear_mean"].to_numpy(dtype=float), False),
    ]
    rows = []
    for mname, raw_s, already_prob in methods:
        prob = np.clip(raw_s, 0.02, 0.98) if already_prob else _minmax_01(raw_s, 0.08, 0.92)
        rows.append(
            {
                "cohort": cohort_label,
                "n_patients": int(len(y)),
                "n_responders": int(np.sum(y)),
                "method": mname,
                "auroc": round(float(roc_auc_score(y, prob)), 4),
                "auprc": round(float(average_precision_score(y, prob)), 4),
                "brier_score": round(float(brier_score_loss(y, prob)), 4),
            }
        )
    return rows


def main() -> None:
    sc_dir = ROOT / "data_vault" / "raw_single_cell"
    vault_pb = ROOT / "data_vault" / "processed_pseudobulk"
    out_tables = ROOT / "outputs" / "stage2_tables"
    vault_pb.mkdir(parents=True, exist_ok=True)
    out_tables.mkdir(parents=True, exist_ok=True)

    tau0_df = pd.read_csv(vault_pb / "hpa_pan_organ_tau0_table.csv", index_col=0)
    tau0_series = tau0_df["tau0_healthy"]

    t2 = pd.read_csv(out_tables / "Table_S2_Target_Resolved_Trials_Deconvolution.csv", usecols=["target_genes"])
    clinical_genes: set[str] = set()
    for val in t2["target_genes"].dropna().astype(str):
        for g in val.split(";"):
            if g.strip():
                clinical_genes.add(g.strip().upper())

    vis_a_p = out_tables / "Table_S5a_Visium_Spatial_Barrier_Sections.csv"
    vis_b_p = out_tables / "Table_S5b_Visium_Consensus_Proximity_Matrix.csv"
    if vis_a_p.exists() and vis_b_p.exists() and (out_tables / "Fig5_Exemplar_Visium_Spot_Landscape.csv").exists():
        print(">>> [1/4] Loading verified 10x Visium Spatial Collagen Barrier & Proximity Matrix from cache...")
        vis_sec_df = pd.read_csv(vis_a_p)
        consensus_prox = pd.read_csv(vis_b_p, index_col=0)
    else:
        print(">>> [1/4] Extracting 10x Visium Spatial Collagen Barrier & Proximity Matrix (GSE171351 + GSE203612)...")
        vis_ev = extract_visium_spatial_evidence(sc_dir)
        vis_sec_df = vis_ev["sections_df"]
        consensus_prox = vis_ev["consensus_proximity_df"]
        vis_sec_df.to_csv(vis_a_p, index=False, encoding="utf-8")
        consensus_prox.to_csv(vis_b_p, encoding="utf-8")
        if vis_ev["exemplar_section"] is not None:
            ex = vis_ev["exemplar_section"]
            ex_df = ex["spot_scores_df"].copy()
            ex_df["x"] = ex["coords"][:, 0]
            ex_df["y"] = ex["coords"][:, 1]
            ex_df.to_csv(out_tables / "Fig5_Exemplar_Visium_Spot_Landscape.csv", encoding="utf-8")
    print(vis_sec_df.to_string(index=False))

    grammar = build_default_solid_tumor_icb_grammar(consensus_prox)
    DeconfoundedAutoGrammarCompiler.export_yaml(
        grammar, ROOT / "outputs" / "compiled_grammars" / "solidtumor_time_icb_v2.yaml"
    )

    print("\n>>> [2/4] Fast-Streaming & Blind-Testing GSE123813 (Nature Medicine 2019 FACS-Deconfounded Pre-Anti-PD-1 BCC + Paired SCC)...")
    df_bcc, df_long = process_gse123813_bcc_and_scc_paired_cohorts(sc_dir, tau0_series, grammar)
    df_bcc.to_csv(out_tables / "Table_S5c_GSE123813_PreAntiPD1_Patient_Predictions.csv", encoding="utf-8")
    df_long.to_csv(out_tables / "Table_S5f_GSE123813_Paired_PrePost_Trajectory.csv", index=False, encoding="utf-8")
    print(df_bcc[["responder", "best_pct_tumor_change", "caf_frac_nonT", "tam_frac_nonT", "Tumor", "myCAF", "CD8_Effector_Reserve", "simulated_icb_remission_prob", "combo_pd1_plus_mycaf_remission_prob"]].to_string())

    print("\n>>> [3/4] Fast-Streaming & Blind-Testing GSE207422 (Cell Research 2023 NSCLC Pre-Biopsy Bulk n=24 + 92k-Cell scRNA)...")
    df_nsclc_bk, df_nsclc_sc = process_gse207422_nsclc_cohorts_and_axis4_pseudobulk(
        sc_dir, vault_pb, tau0_series, grammar, clinical_genes
    )
    df_nsclc_bk.to_csv(out_tables / "Table_S5d_GSE207422_NSCLC_PreBiopsy_Bulk_Predictions.csv", encoding="utf-8")
    df_nsclc_sc.to_csv(out_tables / "Table_S5e_GSE207422_NSCLC_scRNA_Predictions.csv", encoding="utf-8")

    print("\n>>> [4/4] Compiling 5-Baseline Calibrated Arena Summary & Solid-Tumor Quasi-Potential Landscape...")
    arena_rows = []
    arena_rows.extend(evaluate_cohort_5baselines(df_bcc, "GSE123813_BCC_PreTreatment_FullTME_scRNA"))
    arena_rows.extend(evaluate_cohort_5baselines(df_nsclc_bk, "GSE207422_NSCLC_PreTreatment_Biopsy_Bulk"))
    arena_rows.extend(evaluate_cohort_5baselines(df_nsclc_sc, "GSE207422_NSCLC_TME_Validation_scRNA"))

    common_cols = [
        "responder",
        "simulated_icb_remission_prob",
        "b1_static_pdl1",
        "b2_static_tau0",
        "b3_single_tcell_grn",
        "b4_static_lr_centrality",
        "b5_linear_mean",
        "combo_pd1_plus_mycaf_remission_prob",
        "combo_pd1_plus_tam_remission_prob",
    ]
    # Strictly pre-treatment baseline biopsies ONLY (zero post-surgery samples!)
    pooled_pre = pd.concat([df_bcc[common_cols], df_nsclc_bk[common_cols]], axis=0)
    arena_rows.extend(evaluate_cohort_5baselines(pooled_pre, "Pooled_Strictly_PreTreatment_Biopsies_n32"))

    arena_df = pd.DataFrame(arena_rows)
    arena_df.to_csv(out_tables / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv", index=False, encoding="utf-8")
    print("\n=== STAGE II TASK 4 (SOLID TUMOR ANTI-PD-1 SPATIAL & BLIND RESPONDER BENCHMARK) COMPLETE ===")
    print(arena_df.to_string(index=False))

    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)
    qp_pd1 = solver.compute_quasi_potential_landscape(
        [InterventionRegimen("CD8_Exhaustion", 0.42, "antibody")],
        "myCAF",
        "Tumor",
        n_grid=32,
        n_trajectories=240,
        seed=91,
    )
    qp_combo = solver.compute_quasi_potential_landscape(
        [
            InterventionRegimen("CD8_Exhaustion", 0.90, "antibody"),
            InterventionRegimen("myCAF", 0.82, "small_molecule"),
            InterventionRegimen("Tumor", 0.72, "antibody", start_time=6.0),
        ],
        "myCAF",
        "Tumor",
        n_grid=32,
        n_trajectories=240,
        seed=92,
    )
    np.savez_compressed(
        out_tables / "Fig5_Waddington_QuasiPotential_SolidTumor.npz",
        grid_centers=qp_pd1["grid_centers"],
        U_pd1_mono=qp_pd1["U_grid"],
        U_pd1_mycaf_combo=qp_combo["U_grid"],
        lock_depth_pd1_mono=qp_pd1["pathological_lock_depth"],
        lock_depth_combo=qp_combo["pathological_lock_depth"],
    )

    non_resp = pooled_pre[pooled_pre["responder"] == 0]
    print(
        f"\nNon-Responder Rescue Simulation (n={len(non_resp)} strictly pre-treatment resistant biopsies): "
        f"Anti-PD-1 Mono={non_resp['simulated_icb_remission_prob'].mean():.3f} -> "
        f"Anti-PD-1 + myCAF/TGFBR1 Barrier Unlock={non_resp['combo_pd1_plus_mycaf_remission_prob'].mean():.3f} | "
        f"Anti-PD-1 + SPP1/TREM2 TAM Switch={non_resp['combo_pd1_plus_tam_remission_prob'].mean():.3f}"
    )


if __name__ == "__main__":
    main()
