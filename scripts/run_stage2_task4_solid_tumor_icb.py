"""Execute Stage II Task 4 (Pillar III / Fig. 5): Multi-Cohort Solid Tumor Anti-PD-1 Spatial + Pre-Treatment Blind Benchmark.

100% real primary human data:
1. 10x Visium Spatial Transcriptomics (12 sections, 21,796 spots): GSE171351 (Bladder Ca) + GSE203612 (Pan-Cancer Human Carcinoma)
2. Pure Anti-PD-1 Monotherapy Pre-Treatment Single-Cell Cohort: GSE123813 (Yost et al., Nature Medicine 2019 BCC, n=8 FACS-deconfounded + Paired Pre/Post BCC/SCC)
3. Pure Anti-PD-1 Monotherapy Pre-Treatment Melanoma Cohort #1: GSE91061 (Riaz et al., Cell 2017 Nivolumab Monotherapy, n=49 strictly Pre biopsies)
4. Pure Anti-PD-1 Monotherapy Pre-Treatment Melanoma Cohort #2: GSE78220 (Hugo et al., Cell 2016 Pembrolizumab Monotherapy, n=26 strictly Pre unique patients)
5. Neoadjuvant Anti-PD-1 + Platinum Chemo-IO Pre-Treatment NSCLC Biopsy Cohort (n=24) + Resected NSCLC TME scRNA-seq Cohort (n=14, ~92k cells): GSE207422 (Hu et al., Cell Research 2023)

Key Stage II v2.2 Methodological Safeguards:
- Exact full-transcriptome streaming cell UMI normalization (`umi_stride = 1`, zero gene-sampling approximation, < 700 MB RAM).
- Response-label-free Axis-4 Solid Oncology compartment pseudobulk (`early_comp` vs `late_comp` stratified strictly by intrinsic desmoplastic stromal remodeling burden, 100% blind to clinical response).
- Leave-One-Patient-Out (LOPO) reference scaling (`_lopo_robust_01`, zero test-patient moment leakage).
- Head-to-head benchmarking against published immunotherapy signatures (`18-gene T-cell-inflamed GEP`, Ayers et al. JCI 2017; `CAF/TIDE Exclusion`, Jiang et al. Nat. Med. 2018; static `PD-L1`, algebraic `CD8:myCAF` ratio, and `OmniNiche_NoVisium_UniformPij` spatial ablation).
"""
from __future__ import annotations

import gzip
import json
import pathlib
import sys
from typing import Any

import numpy as np
import openpyxl
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

# Published 18-gene T-cell-inflamed Gene Expression Profile (Ayers et al., J. Clin. Invest. 2017; PMID 28650338)
PUBLISHED_GEP_18_GENES: set[str] = {
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1",
    "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT",
}

# Published TIDE Cytotoxic T-Lymphocyte vs CAF/Myeloid Exclusion signatures (Jiang et al., Nature Medicine 2018; PMID 30127393)
PUBLISHED_TIDE_CTL_GENES: set[str] = {"CD8A", "CD8B", "GZMA", "GZMB", "PRF1"}
PUBLISHED_CAF_EXCLUSION_GENES: set[str] = {
    "FAP", "COL1A1", "COL1A2", "COL3A1", "POSTN", "ACTA2", "TAGLN", "TGFB1", "MMP11", "CXCL12", "SPP1", "S100A8", "S100A9",
}

MELANOMA_TUMOR_GENES: set[str] = {
    "SOX10", "MITF", "MLANA", "PMEL", "TYR", "DCT", "SOX9", "MKI67", "TOP2A", "VEGFA", "TGFB1", "B2M",
}
# Purity- & leukocyte-deconfounded tumor-intrinsic oncogenic/angiogenic resistance vs MHC-I antigen presentation for bulk biopsies
BULK_TUMOR_PROLIF_RESIST_GENES: set[str] = {
    "VEGFA", "SOX9", "MYC",
}
BULK_MHC1_PRESENTATION_GENES: set[str] = {
    "B2M", "HLA-A", "HLA-E", "PSMB10",
}
# SPP1/TREM2 profibrotic-immunosuppressive TAM polarity vs CXCL9 M1 polarity in bulk mixtures (Bill et al., Science 2023)
BULK_PROFIBROTIC_TAM_GENES: set[str] = {
    "SPP1", "TREM2", "CD163", "MRC1", "MARCO", "APOE", "CSF1R", "TGFB1", "IL10",
}
BULK_M1_POLARITY_GENES: set[str] = {
    "CXCL9", "STAT1",
}
# Terminal exhaustion TFs/ectonucleotidases distinct from pan-activation checkpoints in bulk mixtures
BULK_TERMINAL_EXHAUSTION_GENES: set[str] = {
    "HAVCR2", "TOX", "ENTPD1", "LAYN",
}

PANEL_GENES: dict[str, set[str]] = {
    "Tumor": set(TME_SIGNATURES["Tumor"]) | {"EPCAM", "KRT17", "KRT7", "KRT19", "MKI67", "TOP2A", "MYC", "B2M", "HLA-A"},
    "TAM": set(TME_SIGNATURES["TAM"]) | {"SPP1", "TREM2", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "APOE", "CSF1R", "MARCO"},
    "myCAF": {"FAP", "COL1A1", "COL1A2", "COL3A1", "MMP11", "POSTN", "ACTA2", "TAGLN", "TGFBR1", "CXCL12", "LRRC15", "FN1"},
    "CD8_Effector": {"CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "CXCL13", "NKG7"},
    "CD8_Exhaustion": set(TME_SIGNATURES["CD8_Exhaustion"]) | {"PDCD1", "HAVCR2", "TOX", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "LAYN"},
    "PDL1_Checkpoint": {"CD274", "PDCD1LG2"},
    "LR_Pairs": {"TGFB1", "TGFBR1", "TGFBR2", "SPP1", "CD44", "CXCL12", "CXCR4", "CSF1", "CSF1R", "CD274", "PDCD1"},
    "GEP_18": PUBLISHED_GEP_18_GENES,
    "TIDE_CTL": PUBLISHED_TIDE_CTL_GENES,
    "CAF_Exclusion": PUBLISHED_CAF_EXCLUSION_GENES,
    "Melanoma_Tumor": MELANOMA_TUMOR_GENES,
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


def _lopo_zscore(arr: np.ndarray) -> np.ndarray:
    """Leave-One-Patient-Out (LOPO) robust MAD z-score (patient i scaled strictly against the N-1 other patients)."""
    a = np.asarray(arr, dtype=float)
    n = len(a)
    if n <= 2:
        med = float(np.median(a))
        mad = float(np.median(np.abs(a - med))) * 1.4826
        sd = mad if mad > 1e-6 else float(np.std(a) + 1e-6)
        return (a - med) / sd
    z = np.empty(n, dtype=float)
    for i in range(n):
        ref = np.delete(a, i)
        med = float(np.median(ref))
        mad = float(np.median(np.abs(ref - med))) * 1.4826
        sd = mad if mad > 1e-6 else float(np.std(ref) + 1e-6)
        z[i] = (a[i] - med) / sd
    return z


def _lopo_robust_01(arr: np.ndarray, lo_clip: float = 0.14, hi_clip: float = 0.88) -> np.ndarray:
    """Leave-One-Patient-Out (LOPO) MAD-standardized logistic state projection into [lo_clip, hi_clip]."""
    z = _lopo_zscore(arr)
    return np.clip(1.0 / (1.0 + np.exp(-0.90 * z)), lo_clip, hi_clip)


def _attach_lopo_states_and_simulations(
    df: pd.DataFrame,
    grammar_visium: Any,
    grammar_no_visium: Any,
    seed: int,
) -> pd.DataFrame:
    """Project raw compartment scores via LOPO scaling, run Visium vs No-Visium Langevin SDE, and standardize baselines."""
    out = df.copy()
    out["Tumor"] = _lopo_robust_01(out["raw_tumor"].to_numpy(), 0.16, 0.86)
    out["TAM"] = _lopo_robust_01(out["raw_tam"].to_numpy(), 0.14, 0.88)
    out["myCAF"] = _lopo_robust_01(out["raw_mycaf"].to_numpy(), 0.14, 0.88)
    out["CD8_Exhaustion"] = _lopo_robust_01(out["raw_exhaust_net"].to_numpy(), 0.16, 0.88)
    out["CD8_Effector_Reserve"] = _lopo_robust_01(out["raw_effector"].to_numpy(), 0.16, 0.90)
    out["b5_linear_mean"] = -0.25 * (out["Tumor"] + out["TAM"] + out["myCAF"] + out["CD8_Exhaustion"])

    # Published CAF/TIDE Exclusion score (LOPO-standardized CTL infiltration minus CAF/myeloid exclusion)
    out["b7_published_caf_exclusion"] = _lopo_zscore(out["raw_tide_ctl"].to_numpy()) - 0.85 * _lopo_zscore(
        out["raw_caf_exclusion"].to_numpy()
    )
    # Simple algebraic CD8:myCAF ratio on identical LOPO state variables (tests value of nonlinear Langevin + Visium dynamics)
    out["b8_algebraic_cd8_mycaf_ratio"] = out["CD8_Effector_Reserve"] / np.maximum(out["myCAF"], 0.05)

    # Platform-invariant LOPO probability projections [0.08, 0.92] for fair cross-platform pooled calibration & AUROC
    for col in [
        "b1_static_pdl1",
        "b2_static_tau0",
        "b3_single_tcell_grn",
        "b4_static_lr_centrality",
        "b5_linear_mean",
        "b6_published_gep18",
        "b7_published_caf_exclusion",
        "b8_algebraic_cd8_mycaf_ratio",
    ]:
        out[f"{col}_prob"] = _lopo_robust_01(out[col].to_numpy(), 0.08, 0.92)

    sim_vis = score_baseline_patients_via_langevin_icb(out, grammar=grammar_visium, n_ensemble=250, seed=seed, compute_combos=True)
    sim_novis = score_baseline_patients_via_langevin_icb(
        out, grammar=grammar_no_visium, n_ensemble=250, seed=seed + 5000, compute_combos=False
    )
    out = out.join(sim_vis)
    out["raw_icb_remission_prob"] = out["simulated_icb_remission_prob"]
    out["raw_no_visium_remission_prob"] = sim_novis["simulated_icb_remission_prob"]
    # Apply identical platform-invariant LOPO [0.08, 0.92] calibration across scRNA / FPKM / TPM cohorts as baselines b1..b8
    out["simulated_icb_remission_prob"] = np.round(_lopo_robust_01(out["raw_icb_remission_prob"].to_numpy(), 0.08, 0.92), 4)
    out["no_visium_icb_remission_prob"] = np.round(_lopo_robust_01(out["raw_no_visium_remission_prob"].to_numpy(), 0.08, 0.92), 4)
    return out


def fast_stream_tsv_gz_gene_subset(
    gz_path: pathlib.Path,
    target_genes: set[str],
    sep: str = "\t",
    umi_stride: int = 1,
) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    """Stream a large gzipped Gene x Cell UMI table in < 650 MB RAM with exact full-transcriptome cell UMI totals (`umi_stride=1`)."""
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

    scale_factor = 1.0 if umi_stride == 1 else float(max(n_total_lines, 1)) / float(max(n_sampled, 1))
    tot_umi_exact = np.maximum(sampled_umi * scale_factor, 100.0)
    gene_norm: dict[str, np.ndarray] = {}
    for g, raw_v in gene_raw.items():
        gene_norm[g] = np.log1p((raw_v / tot_umi_exact) * 1e4)
    return barcodes, gene_norm, tot_umi_exact


def process_gse123813_bcc_and_scc_paired_cohorts(
    sc_dir: pathlib.Path,
    vault_pb: pathlib.Path,
    tau0_series: pd.Series,
    grammar_visium: Any,
    grammar_no_visium: Any,
    force_reextract: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process Yost et al., Nature Medicine 2019 (GSE123813):
    1. Full-TME Pre-Anti-PD-1 Monotherapy BCC blind prediction cohort (su001..su008: 4 Responders vs 4 Non-Responders, FACS-deconfounded).
    2. Paired Pre-vs-Post Anti-PD-1 longitudinal T-cell clonal trajectory across BCC + SCC patients.
    """
    cache_raw = vault_pb / "gse123813_bcc_raw_features_exact_umi1.csv"
    cache_long = vault_pb / "gse123813_paired_prepost_exact_umi1.csv"
    if not force_reextract and cache_raw.exists() and cache_long.exists():
        df_p = pd.read_csv(cache_raw).set_index("patient_id")
        df_p = _attach_lopo_states_and_simulations(df_p, grammar_visium, grammar_no_visium, seed=101)
        return df_p, pd.read_csv(cache_long)

    cdir = sc_dir / "GSE123813_BCC_SCC_AntiPD1"
    supp_xlsx = cdir / "Yost_2019_NatMed_SuppTables.xlsx"
    clin = pd.read_excel(supp_xlsx, sheet_name="SuppTable1", skiprows=3).dropna(subset=["Patient"])
    clin["Patient"] = clin["Patient"].astype(str).str.strip()
    bcc_clin = clin[clin["Tumor Type"] == "BCC"].set_index("Patient")
    scc_clin = clin[clin["Tumor Type"] == "SCC"].set_index("Patient")

    bcc_meta = pd.read_csv(cdir / "GSE123813_bcc_all_metadata.txt.gz", sep="\t")
    bcc_meta["cell.id"] = bcc_meta["cell.id"].astype(str)
    barcodes, gene_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE123813_bcc_scRNA_counts.txt.gz", ALL_PANEL_GENES, umi_stride=1
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

        non_t_pre = p_pre[p_pre["sort"] != "CD45+ CD3+"]
        if len(non_t_pre) < 15:
            continue

        all_pre_idx = p_pre["col_idx"].to_numpy(dtype=int)
        non_t_idx = non_t_pre["col_idx"].to_numpy(dtype=int)
        tum_idx = non_t_pre[non_t_pre["cluster"].isin(["Tumor_1", "Tumor_2"])]["col_idx"].to_numpy(dtype=int)
        tam_idx = non_t_pre[non_t_pre["cluster"].isin(["Macrophages", "DCs"])]["col_idx"].to_numpy(dtype=int)
        caf_idx = non_t_pre[non_t_pre["cluster"].isin(["CAFs", "Myofibroblasts"])]["col_idx"].to_numpy(dtype=int)

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
        gep18_val = _mean_panel(gene_rows, all_pre_idx, PUBLISHED_GEP_18_GENES)
        tide_ctl = _mean_panel(gene_rows, cd8_pre if len(cd8_pre) >= 3 else all_pre_idx, PUBLISHED_TIDE_CTL_GENES)
        caf_excl = _mean_panel(gene_rows, non_t_idx, PUBLISHED_CAF_EXCLUSION_GENES) + 2.0 * caf_frac

        tau_vals = [
            float(tau0_series.loc[g]) * float(np.mean(gene_rows[g][all_pre_idx]))
            for g in (PANEL_GENES["CD8_Effector"] | PANEL_GENES["PDL1_Checkpoint"])
            if g in gene_rows and g in tau0_series.index
        ]

        patient_rows.append(
            {
                "cohort": "GSE123813_BCC_PreTreatment_AntiPD1_Mono_scRNA",
                "regimen_class": "AntiPD1_Monotherapy",
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
                "raw_tide_ctl": tide_ctl,
                "raw_caf_exclusion": caf_excl,
                "b1_static_pdl1": pdl1_val,
                "b2_static_tau0": float(np.mean(tau_vals)) if tau_vals else 0.0,
                "b3_single_tcell_grn": eff_prog - 0.50 * exh_prog,
                "b4_static_lr_centrality": -lr_static,
                "b6_published_gep18": gep18_val,
            }
        )

    scc_meta = pd.read_csv(cdir / "GSE123813_scc_metadata.txt.gz", sep="\t")
    scc_meta["cell.id"] = scc_meta["cell.id"].astype(str)
    s_bcs, s_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE123813_scc_scRNA_counts.txt.gz", ALL_PANEL_GENES, umi_stride=1
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

    df_p_raw = pd.DataFrame(patient_rows)
    df_long = pd.DataFrame(longitudinal_rows)
    df_p_raw.to_csv(cache_raw, index=False, encoding="utf-8")
    df_long.to_csv(cache_long, index=False, encoding="utf-8")
    df_p = df_p_raw.set_index("patient_id")
    df_p = _attach_lopo_states_and_simulations(df_p, grammar_visium, grammar_no_visium, seed=101)
    return df_p, df_long


def _score_bulk_logexpr_cohort(
    expr_df: pd.DataFrame,
    sample_records: list[dict[str, Any]],
    tau0_series: pd.Series,
    tumor_gene_set: set[str],
    grammar_visium: Any,
    grammar_no_visium: Any,
    seed: int,
) -> pd.DataFrame:
    """Compute baseline TME state variables & published biomarkers on log-scale bulk RNA-seq cohorts using LOPO scaling.

    Uses purity-deconfounded proliferative/immune-evasive tumor drive (`BULK_TUMOR_PROLIF_RESIST_GENES`), M2/M1 polarity-deconfounded
    profibrotic TAM drive, non-tumor admixture-deconfounded stromal myCAF (`caf_s - 0.35 * eff_s`), and per-T-cell terminal
    exhaustion (`0.55 * term_exh_s - 0.85 * eff_s`, cancelling total T-cell abundance in unsorted bulk mixtures).
    """
    rows = []
    for rec in sample_records:
        sid = rec["sample_col"]
        col_v = expr_df[sid]

        def _bmod(genes: set[str]) -> float:
            present = [g for g in genes if g in col_v.index]
            return float(col_v.loc[present].mean()) if present else 0.0

        tum_s = _bmod(tumor_gene_set) - 0.55 * _bmod(BULK_MHC1_PRESENTATION_GENES)
        tam_s = _bmod(BULK_PROFIBROTIC_TAM_GENES) - 0.45 * _bmod(BULK_M1_POLARITY_GENES)
        caf_s = _bmod(PANEL_GENES["myCAF"])
        exh_s = _bmod(PANEL_GENES["CD8_Exhaustion"])
        term_exh_s = _bmod(BULK_TERMINAL_EXHAUSTION_GENES)
        eff_s = _bmod(PANEL_GENES["CD8_Effector"])
        pdl1_s = _bmod(PANEL_GENES["PDL1_Checkpoint"])
        lr_s = _bmod(PANEL_GENES["LR_Pairs"])
        gep18_s = _bmod(PUBLISHED_GEP_18_GENES)
        ctl_s = _bmod(PUBLISHED_TIDE_CTL_GENES)
        excl_s = _bmod(PUBLISHED_CAF_EXCLUSION_GENES)

        tau_vals = [
            float(tau0_series.loc[g]) * float(col_v.loc[g])
            for g in (PANEL_GENES["CD8_Effector"] | PANEL_GENES["PDL1_Checkpoint"])
            if g in col_v.index and g in tau0_series.index
        ]

        rows.append(
            {
                **{k: v for k, v in rec.items() if k != "sample_col"},
                "raw_tumor": tum_s,
                "raw_tam": tam_s,
                "raw_mycaf": caf_s - 0.35 * eff_s,
                "raw_effector": eff_s,
                "raw_exhaust_net": 0.55 * term_exh_s - 0.85 * eff_s,
                "raw_tide_ctl": ctl_s,
                "raw_caf_exclusion": excl_s,
                "b1_static_pdl1": pdl1_s,
                "b2_static_tau0": float(np.mean(tau_vals)) if tau_vals else 0.0,
                "b3_single_tcell_grn": eff_s - 0.45 * exh_s,
                "b4_static_lr_centrality": -lr_s,
                "b6_published_gep18": gep18_s,
            }
        )
    df_out = pd.DataFrame(rows).set_index("patient_id")
    return _attach_lopo_states_and_simulations(df_out, grammar_visium, grammar_no_visium, seed=seed)


def process_gse91061_melanoma_nivolumab_monotherapy(
    sc_dir: pathlib.Path,
    tau0_series: pd.Series,
    grammar_visium: Any,
    grammar_no_visium: Any,
) -> pd.DataFrame:
    """Process Riaz et al., Cell 2017 (GSE91061): strictly pre-treatment Nivolumab monotherapy melanoma biopsies (n=49)."""
    cdir = sc_dir / "GSE91061_Melanoma_Nivolumab_Monotherapy"
    entrez_map: dict[str, str] = json.loads((cdir / "ncbi_entrez_to_hgnc.json").read_text(encoding="utf-8"))

    with gzip.open(cdir / "GSE91061_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore") as f:
        meta_lines = [l.rstrip("\r\n") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
    titles = [x.strip('"') for x in meta_lines[0].split("\t")[1:]]
    chars = [[x.strip('"') for x in l.split("\t")[1:]] for l in meta_lines[1:]]

    raw_fpkm = pd.read_csv(cdir / "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz", index_col=0)
    raw_fpkm.index = [entrez_map.get(str(idx).strip(), "") for idx in raw_fpkm.index]
    raw_fpkm = raw_fpkm[raw_fpkm.index != ""]
    raw_fpkm = raw_fpkm[~raw_fpkm.index.duplicated(keep="first")]
    log_expr = np.log2(np.maximum(raw_fpkm.astype(float), 0.0) + 1.0)

    records = []
    seen_pts: set[str] = set()
    for idx_s, s_title in enumerate(titles):
        if s_title not in log_expr.columns:
            continue
        row_c = [c[idx_s] for c in chars if idx_s < len(c)]
        visit = next((x.split(":", 1)[1].strip() for x in row_c if "visit" in x.lower()), "")
        resp = next((x.split(":", 1)[1].strip().upper() for x in row_c if "response" in x.lower()), "")
        if visit.lower() != "pre" or resp not in ("PRCR", "CR", "PR", "SD", "PD"):
            continue
        pt_id = s_title.split("_")[0]
        if pt_id in seen_pts:
            continue
        seen_pts.add(pt_id)
        records.append(
            {
                "cohort": "GSE91061_Melanoma_PreTreatment_Nivolumab_Mono",
                "regimen_class": "AntiPD1_Monotherapy",
                "patient_id": f"GSE91061_{pt_id}",
                "sample_col": s_title,
                "pd1_antibody": "Nivolumab",
                "recist_response": resp,
                "responder": 1 if resp in ("PRCR", "CR", "PR") else 0,
            }
        )
    return _score_bulk_logexpr_cohort(
        log_expr,
        records,
        tau0_series,
        BULK_TUMOR_PROLIF_RESIST_GENES,
        grammar_visium,
        grammar_no_visium,
        seed=404,
    )


def process_gse78220_melanoma_pembrolizumab_monotherapy(
    sc_dir: pathlib.Path,
    tau0_series: pd.Series,
    grammar_visium: Any,
    grammar_no_visium: Any,
) -> pd.DataFrame:
    """Process Hugo et al., Cell 2016 (GSE78220): strictly pre-treatment Pembrolizumab monotherapy melanoma patients (n=26 unique)."""
    cdir = sc_dir / "GSE78220_Melanoma_Pembrolizumab_Monotherapy"
    with gzip.open(cdir / "GSE78220_series_matrix.txt.gz", "rt", encoding="utf-8", errors="ignore") as f:
        meta_lines = [l.rstrip("\r\n") for l in f if l.startswith("!Sample_title") or l.startswith("!Sample_characteristics_ch1")]
    titles = [x.strip('"') for x in meta_lines[0].split("\t")[1:]]
    chars = [[x.strip('"') for x in l.split("\t")[1:]] for l in meta_lines[1:]]

    wb = openpyxl.load_workbook(cdir / "GSE78220_PatientFPKM.xlsx", read_only=True, data_only=True)
    ws = wb.active
    rows_iter = ws.iter_rows(values_only=True)
    header = list(next(rows_iter))
    col_names = [str(c).strip() for c in header[1:] if c is not None]
    n_cols = len(col_names)
    genes = []
    mat_rows = []
    for r in rows_iter:
        if r[0] is None:
            continue
        genes.append(str(r[0]).strip().upper())
        mat_rows.append([float(x) if x is not None else 0.0 for x in r[1 : 1 + n_cols]])
    wb.close()

    fpkm_df = pd.DataFrame(mat_rows, index=genes, columns=col_names)
    fpkm_df = fpkm_df[~fpkm_df.index.duplicated(keep="first")]
    log_raw = np.log2(np.maximum(fpkm_df.astype(float), 0.0) + 1.0)

    # Map sample title prefix (e.g. 'Pt1', 'Pt27A', 'Pt27B') -> xlsx baseline column, strictly excluding on-treatment ('Pt16')
    title_to_col = {c.split(".")[0]: c for c in log_raw.columns}
    patient_cols: dict[str, list[str]] = {}
    patient_meta: dict[str, dict[str, Any]] = {}

    for idx_s, s_title in enumerate(titles):
        row_c = [c[idx_s] for c in chars if idx_s < len(c) and c[idx_s]]
        btime = next((x.split(":", 1)[1].strip().lower() for x in row_c if "biopsy time" in x.lower()), "")
        if btime != "pre-treatment":
            continue
        xcol = title_to_col.get(s_title)
        if not xcol or "ontx" in xcol.lower():
            continue
        pid_clean = next((x.split(":", 1)[1].strip() for x in row_c if "patient id" in x.lower()), s_title)
        resp_str = next((x.split(":", 1)[1].strip() for x in row_c if "anti-pd-1 response" in x.lower()), "")
        if not resp_str:
            continue
        is_resp = 1 if ("complete" in resp_str.lower() or "partial" in resp_str.lower()) else 0
        patient_cols.setdefault(pid_clean, []).append(xcol)
        patient_meta[pid_clean] = {
            "recist_response": resp_str,
            "responder": is_resp,
        }

    # De-duplicate synchronous pre-treatment biopsies from the same patient (Pt27A + Pt27B -> Pt27 mean log2FPKM)
    merged_cols = {}
    records = []
    for pid_clean in sorted(patient_cols.keys()):
        cols = patient_cols[pid_clean]
        merged_cols[pid_clean] = log_raw[cols].mean(axis=1)
        records.append(
            {
                "cohort": "GSE78220_Melanoma_PreTreatment_Pembrolizumab_Mono",
                "regimen_class": "AntiPD1_Monotherapy",
                "patient_id": f"GSE78220_{pid_clean}",
                "sample_col": pid_clean,
                "pd1_antibody": "Pembrolizumab",
                "n_merged_lesions": len(cols),
                "recist_response": patient_meta[pid_clean]["recist_response"],
                "responder": patient_meta[pid_clean]["responder"],
            }
        )
    merged_expr = pd.DataFrame(merged_cols, index=log_raw.index)
    return _score_bulk_logexpr_cohort(
        merged_expr,
        records,
        tau0_series,
        BULK_TUMOR_PROLIF_RESIST_GENES,
        grammar_visium,
        grammar_no_visium,
        seed=505,
    )


def process_gse207422_nsclc_cohorts_and_axis4_pseudobulk(
    sc_dir: pathlib.Path,
    vault_pb: pathlib.Path,
    tau0_series: pd.Series,
    grammar_visium: Any,
    grammar_no_visium: Any,
    clinical_target_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Process Hu et al., Cell Research 2023 (GSE207422) Human Stage III NSCLC Neoadjuvant Anti-PD-1 + Platinum Chemo-IO:
    1. 24 Strictly Pre-Treatment Biopsy Bulk RNA-seq patients.
    2. 14 Evaluable Resected NSCLC TME Single-Cell RNA-seq tumors (~92,000 cells).
    3. Response-Label-Free Axis 4 Solid Oncology compartment pseudobulk (stratified strictly by intrinsic desmoplastic stromal remodeling!).
    """
    cdir = sc_dir / "GSE207422_NSCLC_AntiPD1"

    # Part A: 24 Strictly Pre-treatment biopsy NSCLC Bulk RNA-seq patients (anti-PD-1 + platinum-doublet chemo-IO)
    bk_meta = pd.read_excel(cdir / "GSE207422_NSCLC_bulk_RNAseq_metadata.xlsx")
    bk_expr = pd.read_csv(cdir / "GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz", sep="\t", index_col=0)
    bk_expr.index = bk_expr.index.astype(str).str.upper()
    bk_expr = bk_expr[~bk_expr.index.duplicated(keep="first")]

    bk_records = []
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
        bk_records.append(
            {
                "cohort": "GSE207422_NSCLC_PreTreatment_ChemoAntiPD1_Biopsy",
                "regimen_class": "AntiPD1_Plus_Platinum_ChemoIO",
                "patient_id": f"{mrow['Patient']}_{sid}",
                "sample_col": sid,
                "resource": str(mrow["Resource"]),
                "pd1_antibody": f"{mrow['PD1 Antibody']} + Platinum-Doublet",
                "responder": is_mpr,
                "residual_tumor_frac": res_tum,
            }
        )
    df_bk = _score_bulk_logexpr_cohort(
        bk_expr,
        bk_records,
        tau0_series,
        BULK_TUMOR_PROLIF_RESIST_GENES,
        grammar_visium,
        grammar_no_visium,
        seed=202,
    )

    # Part B: GSE207422 NSCLC 92k-cell scRNA-seq cohort (exact umi_stride=1) + Response-Label-Free Axis 4 Pseudobulk!
    cache_ns_sc = vault_pb / "gse207422_nsclc_sc_raw_features_exact_umi1.csv"
    ax4_e_p = vault_pb / "axis4_solid_tumor_early_comp.csv"
    ax4_l_p = vault_pb / "axis4_solid_tumor_late_comp.csv"
    if "--reextract" not in sys.argv and cache_ns_sc.exists() and ax4_e_p.exists() and ax4_l_p.exists():
        df_sc = pd.read_csv(cache_ns_sc).set_index("patient_id")
        df_sc = _attach_lopo_states_and_simulations(df_sc, grammar_visium, grammar_no_visium, seed=303)
        return df_bk, df_sc

    sc_meta = pd.read_excel(cdir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    tracked_all = ALL_PANEL_GENES | clinical_target_genes
    barcodes, gene_rows, _ = fast_stream_tsv_gz_gene_subset(
        cdir / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz", tracked_all, umi_stride=1
    )
    sample_to_idxs: dict[str, list[int]] = {}
    for idx, bc in enumerate(barcodes):
        parts = bc.split("_")
        sname = "_".join(parts[:2]) if len(parts) >= 2 else bc
        sample_to_idxs.setdefault(sname, []).append(idx)

    comp_names = ["Tumor", "TAM", "myCAF", "CD8_T"]
    t_list = sorted(clinical_target_genes)

    # First pass over ALL scRNA-seq tumor samples in sample_to_idxs (100% BLIND to Pathologic Response!) to build Axis 4 Early vs Late
    blind_tumor_pbs: list[tuple[float, pd.DataFrame]] = []
    for sname in sorted(sample_to_idxs.keys()):
        idxs_all = np.asarray(sample_to_idxs[sname], dtype=int)
        if len(idxs_all) < 50:
            continue

        def _smod_blind(genes: set[str]) -> np.ndarray:
            arrs = [gene_rows[g][idxs_all] for g in genes if g in gene_rows]
            return np.mean(arrs, axis=0) if arrs else np.zeros(len(idxs_all), dtype=float)

        bt_tum = _smod_blind(PANEL_GENES["Tumor"])
        bt_tam = _smod_blind(PANEL_GENES["TAM"])
        bt_caf = _smod_blind(PANEL_GENES["myCAF"])
        bt_t = np.maximum(_smod_blind(PANEL_GENES["CD8_Effector"]), _smod_blind(PANEL_GENES["CD8_Exhaustion"]))
        dom_b = np.argmax(np.vstack([bt_tum, bt_tam, bt_caf, bt_t]), axis=0)
        desmoplastic_burden = float(np.mean(bt_caf)) + 2.8 * float(np.mean(dom_b == 2))

        pb_dict: dict[str, np.ndarray] = {}
        for ci, cname in enumerate(comp_names):
            cmask = dom_b == ci
            use_idxs = idxs_all[cmask] if np.sum(cmask) >= 5 else idxs_all
            vec = np.array([float(np.mean(gene_rows[g][use_idxs])) if g in gene_rows else 0.0 for g in t_list], dtype=float)
            pb_dict[cname] = vec
        blind_tumor_pbs.append((desmoplastic_burden, pd.DataFrame(pb_dict, index=t_list)))

    if len(blind_tumor_pbs) >= 4:
        med_desmo = float(np.median([x[0] for x in blind_tumor_pbs]))
        early_pbs = [pb for score_d, pb in blind_tumor_pbs if score_d <= med_desmo]
        late_pbs = [pb for score_d, pb in blind_tumor_pbs if score_d > med_desmo]
        if early_pbs and late_pbs:
            (sum(early_pbs) / float(len(early_pbs))).to_csv(
                vault_pb / "axis4_solid_tumor_early_comp.csv", encoding="utf-8"
            )
            (sum(late_pbs) / float(len(late_pbs))).to_csv(
                vault_pb / "axis4_solid_tumor_late_comp.csv", encoding="utf-8"
            )

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

        pdl1_s = float(np.mean(_smod(PANEL_GENES["PDL1_Checkpoint"])))
        lr_s = float(np.mean(_smod(PANEL_GENES["LR_Pairs"])))
        gep18_s = float(np.mean(_smod(PUBLISHED_GEP_18_GENES)))
        ctl_s = float(np.mean(_smod(PUBLISHED_TIDE_CTL_GENES)))
        excl_s = float(np.mean(_smod(PUBLISHED_CAF_EXCLUSION_GENES))) + 2.0 * caf_frac
        tau_vals = [
            float(tau0_series.loc[g]) * float(np.mean(gene_rows[g][idxs]))
            for g in (PANEL_GENES["CD8_Effector"] | PANEL_GENES["PDL1_Checkpoint"])
            if g in gene_rows and g in tau0_series.index
        ]

        sc_rows.append(
            {
                "cohort": "GSE207422_NSCLC_Resected_TME_scRNA",
                "regimen_class": "AntiPD1_Plus_Platinum_ChemoIO_ResectedTME",
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
                "raw_tide_ctl": ctl_s,
                "raw_caf_exclusion": excl_s,
                "b1_static_pdl1": pdl1_s,
                "b2_static_tau0": float(np.mean(tau_vals)) if tau_vals else 0.0,
                "b3_single_tcell_grn": float(np.mean(c_eff)) - 0.45 * float(np.mean(c_exh)),
                "b4_static_lr_centrality": -lr_s,
                "b6_published_gep18": gep18_s,
            }
        )

    df_sc_raw = pd.DataFrame(sc_rows)
    df_sc_raw.to_csv(cache_ns_sc, index=False, encoding="utf-8")
    df_sc = df_sc_raw.set_index("patient_id")
    df_sc = _attach_lopo_states_and_simulations(df_sc, grammar_visium, grammar_no_visium, seed=303)
    return df_bk, df_sc


def evaluate_cohort_methods(df: pd.DataFrame, cohort_label: str) -> list[dict[str, Any]]:
    """Compute calibrated AUROC, AUPRC, and Brier Score across OmniNiche, Visium spatial ablation, and 8 clinical/biological baselines."""
    y = df["responder"].to_numpy(dtype=int)
    methods = [
        ("OmniNiche_v2_Spatial_Langevin", df["simulated_icb_remission_prob"].to_numpy(dtype=float)),
        ("OmniNiche_NoVisium_UniformPij", df["no_visium_icb_remission_prob"].to_numpy(dtype=float)),
        ("B1_Static_PDL1_CD274", df["b1_static_pdl1_prob"].to_numpy(dtype=float)),
        ("B2_Static_Healthy_Tau0", df["b2_static_tau0_prob"].to_numpy(dtype=float)),
        ("B3_SingleComp_TCell_GRN", df["b3_single_tcell_grn_prob"].to_numpy(dtype=float)),
        ("B4_Static_LR_Centrality", df["b4_static_lr_centrality_prob"].to_numpy(dtype=float)),
        ("B5_Uninformative_Linear_Mean", df["b5_linear_mean_prob"].to_numpy(dtype=float)),
        ("B6_Published_18Gene_Tcell_GEP", df["b6_published_gep18_prob"].to_numpy(dtype=float)),
        ("B7_Published_CAF_TIDE_Exclusion", df["b7_published_caf_exclusion_prob"].to_numpy(dtype=float)),
        ("B8_Algebraic_CD8_to_myCAF_Ratio", df["b8_algebraic_cd8_mycaf_ratio_prob"].to_numpy(dtype=float)),
    ]
    rows = []
    for mname, prob_s in methods:
        prob = np.clip(prob_s, 1e-4, 1.0 - 1e-4)
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
        print(">>> [1/6] Loading verified 10x Visium Spatial Collagen Barrier & Proximity Matrix from cache...")
        vis_sec_df = pd.read_csv(vis_a_p)
        consensus_prox = pd.read_csv(vis_b_p, index_col=0)
    else:
        print(">>> [1/6] Extracting 10x Visium Spatial Collagen Barrier & Proximity Matrix (GSE171351 + GSE203612)...")
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

    grammar_visium = build_default_solid_tumor_icb_grammar(consensus_prox)
    uniform_prox = pd.DataFrame(
        np.ones((4, 4), dtype=float),
        index=["Tumor", "TAM", "myCAF", "CD8_Exhaustion"],
        columns=["Tumor", "TAM", "myCAF", "CD8_Exhaustion"],
    )
    grammar_no_visium = build_default_solid_tumor_icb_grammar(uniform_prox)
    DeconfoundedAutoGrammarCompiler.export_yaml(
        grammar_visium, ROOT / "outputs" / "compiled_grammars" / "solidtumor_time_icb_v2.yaml"
    )

    print("\n>>> [2/6] Exact-UMI Streaming & LOPO Blind-Testing GSE123813 (Nat. Med. 2019 Anti-PD-1 Monotherapy BCC n=8 + Paired SCC)...")
    df_bcc, df_long = process_gse123813_bcc_and_scc_paired_cohorts(
        sc_dir, vault_pb, tau0_series, grammar_visium, grammar_no_visium, force_reextract=("--reextract" in sys.argv)
    )
    df_bcc.to_csv(out_tables / "Table_S5c_GSE123813_PreAntiPD1_Patient_Predictions.csv", encoding="utf-8")
    df_long.to_csv(out_tables / "Table_S5f_GSE123813_Paired_PrePost_Trajectory.csv", index=False, encoding="utf-8")

    print("\n>>> [3/6] LOPO Blind-Testing GSE91061 (Riaz et al., Cell 2017 Pre-Nivolumab Monotherapy Melanoma n=49)...")
    df_91061 = process_gse91061_melanoma_nivolumab_monotherapy(sc_dir, tau0_series, grammar_visium, grammar_no_visium)
    df_91061.to_csv(out_tables / "Table_S5j_GSE91061_Melanoma_PreNivo_Mono_Predictions.csv", encoding="utf-8")

    print("\n>>> [4/6] LOPO Blind-Testing GSE78220 (Hugo et al., Cell 2016 Pre-Pembrolizumab Monotherapy Melanoma n=26 unique)...")
    df_78220 = process_gse78220_melanoma_pembrolizumab_monotherapy(sc_dir, tau0_series, grammar_visium, grammar_no_visium)
    df_78220.to_csv(out_tables / "Table_S5k_GSE78220_Melanoma_PrePembro_Mono_Predictions.csv", encoding="utf-8")

    print("\n>>> [5/6] Exact-UMI Streaming & LOPO Blind-Testing GSE207422 (Cell Res. 2023 NSCLC Chemo-Anti-PD-1 Pre-Biopsy n=24 + Resected TME scRNA n=14)...")
    df_nsclc_bk, df_nsclc_sc = process_gse207422_nsclc_cohorts_and_axis4_pseudobulk(
        sc_dir, vault_pb, tau0_series, grammar_visium, grammar_no_visium, clinical_genes
    )
    df_nsclc_bk.to_csv(out_tables / "Table_S5d_GSE207422_NSCLC_PreBiopsy_Bulk_Predictions.csv", encoding="utf-8")
    df_nsclc_sc.to_csv(out_tables / "Table_S5e_GSE207422_NSCLC_scRNA_Predictions.csv", encoding="utf-8")

    print("\n>>> [6/6] Compiling 10-Method Calibrated Multi-Cohort Arena Summary & Solid-Tumor Quasi-Potential Landscape...")
    common_cols = [
        "cohort",
        "regimen_class",
        "responder",
        "Tumor",
        "TAM",
        "myCAF",
        "CD8_Exhaustion",
        "CD8_Effector_Reserve",
        "raw_icb_remission_prob",
        "simulated_icb_remission_prob",
        "no_visium_icb_remission_prob",
        "b1_static_pdl1_prob",
        "b2_static_tau0_prob",
        "b3_single_tcell_grn_prob",
        "b4_static_lr_centrality_prob",
        "b5_linear_mean_prob",
        "b6_published_gep18_prob",
        "b7_published_caf_exclusion_prob",
        "b8_algebraic_cd8_mycaf_ratio_prob",
        "combo_pd1_plus_mycaf_remission_prob",
        "combo_pd1_plus_tam_remission_prob",
    ]
    pooled_mono_pre = pd.concat([df_bcc[common_cols], df_91061[common_cols], df_78220[common_cols]], axis=0)
    pooled_all_pre = pd.concat([pooled_mono_pre, df_nsclc_bk[common_cols]], axis=0)
    pooled_all_pre.to_csv(out_tables / "Table_S5i_Pooled_PreTreatment_107Patients_Master.csv", encoding="utf-8")

    arena_rows = []
    arena_rows.extend(evaluate_cohort_methods(df_bcc, "GSE123813_BCC_PreTreatment_AntiPD1_Mono_scRNA"))
    arena_rows.extend(evaluate_cohort_methods(df_91061, "GSE91061_Melanoma_PreTreatment_Nivolumab_Mono"))
    arena_rows.extend(evaluate_cohort_methods(df_78220, "GSE78220_Melanoma_PreTreatment_Pembrolizumab_Mono"))
    arena_rows.extend(evaluate_cohort_methods(pooled_mono_pre, f"Pooled_Pure_AntiPD1_Monotherapy_PreTreatment_n{len(pooled_mono_pre)}"))
    arena_rows.extend(evaluate_cohort_methods(df_nsclc_bk, "GSE207422_NSCLC_PreTreatment_ChemoAntiPD1_Biopsy"))
    arena_rows.extend(evaluate_cohort_methods(pooled_all_pre, f"Pooled_All_Strictly_PreTreatment_Biopsies_n{len(pooled_all_pre)}"))
    arena_rows.extend(evaluate_cohort_methods(df_nsclc_sc, "GSE207422_NSCLC_Resected_TME_scRNA"))

    arena_df = pd.DataFrame(arena_rows)
    arena_df.to_csv(out_tables / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv", index=False, encoding="utf-8")
    print("\n=== STAGE II v2.2 TASK 4 (MULTI-COHORT SOLID TUMOR ANTI-PD-1 SPATIAL & BLIND RESPONDER BENCHMARK) COMPLETE ===")
    print(arena_df.to_string(index=False))

    solver = LangevinNicheSolver(grammar_visium, t_end=96.0, dt=0.5)
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

    non_resp_mono = pooled_mono_pre[pooled_mono_pre["responder"] == 0]
    non_resp_all = pooled_all_pre[pooled_all_pre["responder"] == 0]
    print(
        f"\nPure Anti-PD-1 Monotherapy Non-Responder Rescue (n={len(non_resp_mono)} resistant pre-treatment patients): "
        f"Anti-PD-1 Mono={non_resp_mono['simulated_icb_remission_prob'].mean():.3f} -> "
        f"Anti-PD-1 + myCAF/TGFBR1 Unlock={non_resp_mono['combo_pd1_plus_mycaf_remission_prob'].mean():.3f} | "
        f"Anti-PD-1 + SPP1/TREM2 TAM Switch={non_resp_mono['combo_pd1_plus_tam_remission_prob'].mean():.3f}"
    )
    print(
        f"All Pre-Treatment Non-Responder Rescue (n={len(non_resp_all)} resistant pre-treatment biopsies): "
        f"Anti-PD-1 Mono={non_resp_all['simulated_icb_remission_prob'].mean():.3f} -> "
        f"Anti-PD-1 + myCAF/TGFBR1 Unlock={non_resp_all['combo_pd1_plus_mycaf_remission_prob'].mean():.3f} | "
        f"Anti-PD-1 + SPP1/TREM2 TAM Switch={non_resp_all['combo_pd1_plus_tam_remission_prob'].mean():.3f}"
    )


if __name__ == "__main__":
    main()
