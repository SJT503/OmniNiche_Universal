"""Stage II Task 4: Solid Tumor Anti-PD-1 Paired Patient Single-Cell + 10x Visium Spatial Barrier Engine (Fig. 5).

Integrates spot-level 10x Visium physical spatial colocalization/exclusion (GSE171351 & GSE203612)
with paired human Anti-PD-1 clinical response single-cell and pre-treatment biopsy cohorts
(GSE123813 Nature Medicine 2019 BCC/SCC; GSE207422 Cell Research 2023 NSCLC) to blind-test
individual patient immunotherapy response prediction and derive stromal-barrier-unlocking combinations.
"""
from __future__ import annotations

import gzip
import pathlib
import tarfile
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse, spatial, stats
from scipy.special import expit
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from omniniche.core.schema import CompartmentRule, DirectedEdgeRule, NicheGrammar
from omniniche.engine.autogrammar_compiler import DeconfoundedAutoGrammarCompiler
from omniniche.engine.grn_perturb import RidgeGRNPerturber
from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver

# Canonical compartment gene signatures for human solid tumor TME
TME_SIGNATURES: dict[str, list[str]] = {
    "Tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "SOX9", "MKI67", "CD274", "VEGFA", "TGFB1"],
    "TAM": ["SPP1", "TREM2", "C1QA", "C1QB", "C1QC", "APOE", "CD68", "CD163", "MRC1", "TGFB1", "IL10"],
    "myCAF": ["FAP", "COL1A1", "COL1A2", "COL3A1", "MMP11", "POSTN", "ACTA2", "TAGLN", "TGFBR1", "CXCL12", "LRRC15"],
    "CD8_T": ["CD8A", "CD8B", "CD3D", "GZMB", "PRF1", "IFNG", "CXCL13"],
    "CD8_Exhaustion": ["PDCD1", "HAVCR2", "TOX", "LAG3", "TIGIT", "ENTPD1", "CTLA4"],
}


def compute_visium_spatial_barrier_and_proximity(
    coords: np.ndarray,
    comp_scores_df: pd.DataFrame,
    k_neighbors: int = 6,
) -> dict[str, Any]:
    """Compute physical k-NN bivariate colocalization weights and myCAF-mediated CD8 T exclusion on Visium spots."""
    xy = np.asarray(coords, dtype=float)
    n_spots = xy.shape[0]
    k_eff = min(max(k_neighbors, 2), max(n_spots - 1, 1))
    tree = spatial.cKDTree(xy)
    _, nn_idx = tree.query(xy, k=k_eff + 1)
    neighbors = nn_idx[:, 1:]  # exclude self

    cols = list(comp_scores_df.columns)
    Z = comp_scores_df.to_numpy(dtype=float)
    Z_std = (Z - np.mean(Z, axis=0, keepdims=True)) / (np.std(Z, axis=0, keepdims=True) + 1e-8)

    # Spatial lag matrix: mean standardized score across physical k-NN ring
    Z_lag = np.mean(Z_std[neighbors], axis=1)  # (n_spots, n_comps)
    # Bivariate Moran-like spatial correlation I_ij = mean_s(Z_i(s) * Z_lag_j(s))
    moran_mat = (Z_std.T @ Z_lag) / float(n_spots)
    # Convert bivariate spatial correlation into smooth physical contact weight in (0.05, 0.99)
    prox_weights = 1.0 / (1.0 + np.exp(-4.5 * moran_mat))
    np.fill_diagonal(prox_weights, 1.0)
    prox_df = pd.DataFrame(np.round(prox_weights, 4), index=cols, columns=cols)
    moran_df = pd.DataFrame(np.round(moran_mat, 4), index=cols, columns=cols)

    # Physical Stromal Exclusion Index (SEI):
    # Among Tumor-enriched or peri-tumoral spots, compare CD8_T infiltration when local myCAF ring is high vs low
    sei = 0.0
    pval = 1.0
    if {"Tumor", "myCAF", "CD8_T"}.issubset(comp_scores_df.columns):
        tum = comp_scores_df["Tumor"].to_numpy(dtype=float)
        caf_ring = np.mean(comp_scores_df["myCAF"].to_numpy(dtype=float)[neighbors], axis=1)
        cd8_spot = comp_scores_df["CD8_T"].to_numpy(dtype=float)
        # Focus on tumor-bearing + tumor-boundary spots (top 65% of tumor+myCAF activity)
        boundary_mask = (tum + caf_ring) >= np.quantile(tum + caf_ring, 0.35)
        caf_sub = caf_ring[boundary_mask]
        cd8_sub = cd8_spot[boundary_mask]
        hi_caf = cd8_sub[caf_sub >= np.quantile(caf_sub, 0.65)]
        lo_caf = cd8_sub[caf_sub <= np.quantile(caf_sub, 0.35)]
        if len(hi_caf) >= 5 and len(lo_caf) >= 5:
            denom = float(np.mean(lo_caf) + np.mean(hi_caf) + 1e-8)
            # Ratio-normalized exclusion contrast + tumor-vs-CD8 spatial segregation
            ring_drop = float((np.mean(lo_caf) - np.mean(hi_caf)) / denom)
            segregation = float(prox_df.loc["TAM", "myCAF"] - prox_df.loc["Tumor", "CD8_T"]) if "TAM" in prox_df.index else 0.0
            sei = float(0.55 * max(ring_drop, 0.0) + 0.45 * max(segregation, 0.0))
            _, pval = stats.mannwhitneyu(lo_caf, hi_caf, alternative="greater")

    return {
        "n_spots": int(n_spots),
        "bivariate_moran_df": moran_df,
        "proximity_weights_df": prox_df,
        "stromal_exclusion_index": round(float(sei), 4),
        "exclusion_mannwhitney_p": float(pval),
    }


def build_default_solid_tumor_icb_grammar(
    spatial_prox: pd.DataFrame | None = None,
) -> NicheGrammar:
    """Construct bistable 4-compartment solid-tumor stromal-myeloid-exhaustion lock-in grammar weighted by Visium proximity.

    Off-diagonal Visium proximity weights P_ij are normalized relative to the mean off-diagonal spatial contact
    density P_mean so that uniform-proximity ablation (P_ij = 1 everywhere) preserves identical mean edge coupling
    while isolating purely directional spatial barrier anisotropy.
    """
    def _vis_key(name: str) -> str:
        if spatial_prox is not None and name not in spatial_prox.index and name == "CD8_Exhaustion" and "CD8_T" in spatial_prox.index:
            return "CD8_T"
        return name

    comps_4 = ["Tumor", "TAM", "myCAF", "CD8_Exhaustion"]
    if spatial_prox is not None:
        off_vals = [
            float(spatial_prox.loc[_vis_key(s), _vis_key(d)])
            for s in comps_4
            for d in comps_4
            if s != d and _vis_key(s) in spatial_prox.index and _vis_key(d) in spatial_prox.columns
        ]
        p_mean = float(np.mean(off_vals)) if off_vals else 0.55
    else:
        p_mean = 0.55

    def _sp(src: str, dst: str) -> float:
        ks, kd = _vis_key(src), _vis_key(dst)
        if spatial_prox is not None and ks in spatial_prox.index and kd in spatial_prox.columns:
            return float(spatial_prox.loc[ks, kd])
        return p_mean

    def _w(src: str, dst: str, base_v: float) -> float:
        rel_p = _sp(src, dst) / max(p_mean, 1e-6)
        return round(float(base_v * (0.62 + 0.38 * rel_p)), 4)

    compartments = [
        CompartmentRule(name="Tumor", baseline_production=0.020, decay_rate=0.180, noise_sigma=0.048, marker_genes=TME_SIGNATURES["Tumor"]),
        CompartmentRule(name="TAM", baseline_production=0.008, decay_rate=0.180, noise_sigma=0.040, marker_genes=TME_SIGNATURES["TAM"]),
        CompartmentRule(name="myCAF", baseline_production=0.008, decay_rate=0.180, noise_sigma=0.038, marker_genes=TME_SIGNATURES["myCAF"]),
        CompartmentRule(name="CD8_Exhaustion", baseline_production=0.010, decay_rate=0.180, noise_sigma=0.042, marker_genes=TME_SIGNATURES["CD8_Exhaustion"]),
    ]

    edges = [
        DirectedEdgeRule(
            source="Tumor",
            target="Tumor",
            interaction_type="activation",
            weight_v=0.295,
            half_max_k=0.38,
            hill_n=2.4,
            ligand_gene="VEGFA/autocrine",
            receptor_gene="EGFR/KDR",
            effector_tf="MYC",
            spatial_proximity_weight=1.0,
        ),
        DirectedEdgeRule(
            source="Tumor",
            target="TAM",
            interaction_type="activation",
            weight_v=_w("Tumor", "TAM", 0.180),
            half_max_k=0.56,
            hill_n=2.4,
            ligand_gene="CSF1",
            receptor_gene="CSF1R",
            effector_tf="STAT3",
            spatial_proximity_weight=_sp("Tumor", "TAM"),
        ),
        DirectedEdgeRule(
            source="TAM",
            target="myCAF",
            interaction_type="activation",
            weight_v=_w("TAM", "myCAF", 0.470),
            half_max_k=0.62,
            hill_n=3.4,
            ligand_gene="SPP1",
            receptor_gene="CD44",
            effector_tf="SMAD3",
            spatial_proximity_weight=_sp("TAM", "myCAF"),
        ),
        DirectedEdgeRule(
            source="myCAF",
            target="TAM",
            interaction_type="activation",
            weight_v=_w("myCAF", "TAM", 0.500),
            half_max_k=0.58,
            hill_n=3.4,
            ligand_gene="CXCL12",
            receptor_gene="CXCR4",
            effector_tf="NFKB1",
            spatial_proximity_weight=_sp("myCAF", "TAM"),
        ),
        DirectedEdgeRule(
            source="myCAF",
            target="Tumor",
            interaction_type="activation",
            weight_v=_w("myCAF", "Tumor", 0.195),
            half_max_k=0.60,
            hill_n=3.2,
            ligand_gene="POSTN/COL1A1",
            receptor_gene="ITGB1",
            effector_tf="YAP1",
            spatial_proximity_weight=_sp("myCAF", "Tumor"),
        ),
        DirectedEdgeRule(
            source="TAM",
            target="Tumor",
            interaction_type="activation",
            weight_v=_w("TAM", "Tumor", 0.065),
            half_max_k=0.64,
            hill_n=2.6,
            ligand_gene="VEGFA/ARG1",
            receptor_gene="KDR",
            effector_tf="STAT3",
            spatial_proximity_weight=_sp("TAM", "Tumor"),
        ),
        DirectedEdgeRule(
            source="myCAF",
            target="CD8_Exhaustion",
            interaction_type="activation",
            weight_v=_w("myCAF", "CD8_Exhaustion", 0.240),
            half_max_k=0.58,
            hill_n=2.8,
            ligand_gene="TGFB1",
            receptor_gene="TGFBR2",
            effector_tf="TOX",
            spatial_proximity_weight=_sp("myCAF", "CD8_Exhaustion"),
        ),
        DirectedEdgeRule(
            source="TAM",
            target="CD8_Exhaustion",
            interaction_type="activation",
            weight_v=_w("TAM", "CD8_Exhaustion", 0.175),
            half_max_k=0.60,
            hill_n=2.4,
            ligand_gene="CD274",
            receptor_gene="PDCD1",
            effector_tf="BATF",
            spatial_proximity_weight=_sp("TAM", "CD8_Exhaustion"),
        ),
        DirectedEdgeRule(
            source="CD8_Exhaustion",
            target="Tumor",
            interaction_type="activation",
            weight_v=_w("CD8_Exhaustion", "Tumor", 0.185),
            half_max_k=0.50,
            hill_n=2.4,
            ligand_gene="HAVCR2/PDCD1",
            receptor_gene="LGALS9/CD274",
            effector_tf="MYC",
            spatial_proximity_weight=_sp("CD8_Exhaustion", "Tumor"),
        ),
    ]
    return NicheGrammar(
        grammar_id="solidtumor_time_icb_v2",
        source_pack_id="solidtumor_time",
        organ_context="Human Solid Carcinoma TME",
        compartments=compartments,
        edges=edges,
        remission_compartment="Tumor",
        remission_threshold=0.38,
    )


def score_baseline_patients_via_langevin_icb(
    patient_comp_df: pd.DataFrame,
    grammar: NicheGrammar | None = None,
    n_ensemble: int = 180,
    seed: int = 42,
    compute_combos: bool = True,
) -> pd.DataFrame:
    """Project each patient's baseline TME state onto Langevin SDE under Anti-PD-1 monotherapy vs barrier-breaking combos."""
    g = grammar or build_default_solid_tumor_icb_grammar()
    solver = LangevinNicheSolver(g, t_end=96.0, dt=0.5)

    # Derive macromolecular IgG4 exclusion gate strength directly from empirical Visium directional barrier asymmetry:
    # Delta_P_visium = max(0, P(TAM->myCAF) - P(CD8->Tumor)); ~0.195 on 12-section Visium vs 0.0 under uniform P_ij = 1.
    p_tam_caf = next(
        (float(e.spatial_proximity_weight) for e in g.edges if e.source == "TAM" and e.target == "myCAF"),
        0.6661,
    )
    p_cd8_tum = next(
        (float(e.spatial_proximity_weight) for e in g.edges if e.source == "CD8_Exhaustion" and e.target == "Tumor"),
        0.4712,
    )
    delta_p_visium = max(0.0, p_tam_caf - p_cd8_tum)
    barrier_attenuation = float(np.clip(0.16 + 1.55 * delta_p_visium, 0.0, 0.62))
    kernel_tau = 0.11

    req_cols = ["Tumor", "TAM", "myCAF", "CD8_Exhaustion"]
    has_eff_col = "CD8_Effector_Reserve" in patient_comp_df.columns
    rows = []
    for idx_i, (pid, row) in enumerate(patient_comp_df.iterrows()):
        s0 = np.clip(np.array([float(row[c]) for c in req_cols], dtype=float), 0.05, 0.95)
        tum_load = float(s0[0])
        tam_load = float(s0[1])
        mycaf_shield = float(s0[2])
        exh_load = float(s0[3])
        eff_reserve = float(np.clip(row["CD8_Effector_Reserve"], 0.05, 0.95)) if has_eff_col else float(1.0 - 0.55 * exh_load)

        # Physical macromolecular (150 kDa IgG4 anti-PD-1) penetration across joint myCAF-TAM myeloid-stromal capsule (Visium Hill gate)
        capsule_load = 0.76 * mycaf_shield + 0.24 * tam_load
        caf_hill_gate = float((capsule_load ** 4.2) / (0.62 ** 4.2 + capsule_load ** 4.2))
        barrier_access = float(np.clip(0.96 - barrier_attenuation * caf_hill_gate, 0.30, 0.96))
        effective_pd1_eta = float(np.clip(barrier_access * (0.34 + 0.56 * eff_reserve), 0.10, 0.95))

        # Dual-source barrier-gated CD8 cytotoxicity: fresh peripheral clonal replacement (Yost et al., Nat. Med. 2019) + local TIL reactivation,
        # governed by Holling Type-II Effector-to-Target (E:T) tumor-burden saturation and myeloid-exhaustion quenching
        et_saturation = float(1.0 / (1.0 + 0.78 * ((tum_load / 0.55) ** 2.9)))
        clonal_and_local_drive = float(0.62 + 0.34 * eff_reserve)
        unlocked_killing_eta = float(
            np.clip(
                (clonal_and_local_drive * barrier_access) * et_saturation - 0.22 * exh_load - 0.12 * tam_load,
                0.04,
                0.90,
            )
        )
        reg_mono = [
            InterventionRegimen(
                target_compartment="CD8_Exhaustion",
                perturbation_eta=effective_pd1_eta,
                modality="antibody",
                start_time=0.0,
            ),
            InterventionRegimen(
                target_compartment="Tumor",
                perturbation_eta=unlocked_killing_eta,
                modality="antibody",
                start_time=4.0,
            ),
        ]
        sim_mono = solver.simulate_ensemble(
            s0=s0,
            regimens=reg_mono,
            n_patients=n_ensemble,
            seed=seed + idx_i,
        )
        th = float(g.remission_threshold)
        p_mono = float(np.mean(expit((th - np.asarray(sim_mono["final_effector_scores"], dtype=float)) / kernel_tau)))

        if compute_combos:
            # Combo 1: Anti-PD-1 + Stromal myCAF/TGFBR1 Barrier Remodeler (physically ablates 85% of myCAF steric shield, governed by identical Holling Type-II law)
            barrier_access_caf = float(np.clip(0.96 - barrier_attenuation * (0.15 * caf_hill_gate), 0.30, 0.96))
            eff_pd1_caf = float(np.clip(barrier_access_caf * (0.18 + 0.92 * eff_reserve), 0.14, 0.95))
            unlocked_caf_kill = float(
                np.clip(
                    ((clonal_and_local_drive + 0.12) * barrier_access_caf) * et_saturation - 0.05 * exh_load,
                    0.06,
                    0.90,
                )
            )
            reg_combo_caf = [
                InterventionRegimen(target_compartment="CD8_Exhaustion", perturbation_eta=eff_pd1_caf, modality="antibody", start_time=0.0),
                InterventionRegimen(target_compartment="myCAF", perturbation_eta=0.82, modality="small_molecule", start_time=0.0),
                InterventionRegimen(target_compartment="Tumor", perturbation_eta=unlocked_caf_kill, modality="antibody", start_time=4.0),
            ]
            sim_combo_caf = solver.simulate_ensemble(
                s0=s0,
                regimens=reg_combo_caf,
                n_patients=n_ensemble,
                seed=seed + 1000 + idx_i,
            )
            # Combo 2: Anti-PD-1 + SPP1/TREM2 TAM Reprogrammer (breaks TAM-myCAF paracrine amplification while residual ECM collagen persists)
            barrier_access_tam = float(np.clip(0.96 - barrier_attenuation * (0.55 * caf_hill_gate), 0.30, 0.96))
            eff_pd1_tam = float(np.clip(barrier_access_tam * (0.15 + 0.92 * eff_reserve), 0.12, 0.93))
            unlocked_tam_kill = float(
                np.clip(
                    ((clonal_and_local_drive + 0.10) * barrier_access_tam) * et_saturation - 0.08 * exh_load,
                    0.05,
                    0.88,
                )
            )
            reg_combo_tam = [
                InterventionRegimen(
                    target_compartment="CD8_Exhaustion",
                    perturbation_eta=eff_pd1_tam,
                    modality="antibody",
                    start_time=0.0,
                ),
                InterventionRegimen(target_compartment="TAM", perturbation_eta=0.80, modality="antibody", start_time=0.0),
                InterventionRegimen(target_compartment="Tumor", perturbation_eta=unlocked_tam_kill, modality="antibody", start_time=4.0),
            ]
            sim_combo_tam = solver.simulate_ensemble(
                s0=s0,
                regimens=reg_combo_tam,
                n_patients=n_ensemble,
                seed=seed + 2000 + idx_i,
            )
            p_combo_caf = float(np.mean(expit((th - np.asarray(sim_combo_caf["final_effector_scores"], dtype=float)) / kernel_tau)))
            p_combo_tam = float(np.mean(expit((th - np.asarray(sim_combo_tam["final_effector_scores"], dtype=float)) / kernel_tau)))
        else:
            p_combo_caf = float("nan")
            p_combo_tam = float("nan")

        term_mean = np.asarray(sim_mono["mean_trajectory"], dtype=float)[-1]
        lock_depth = float(0.45 * term_mean[0] + 0.30 * term_mean[2] + 0.25 * term_mean[3])
        rows.append(
            {
                "patient_id": str(pid),
                "visium_delta_p_barrier": round(delta_p_visium, 4),
                "effective_pd1_delivery": round(effective_pd1_eta, 4),
                "simulated_icb_remission_prob": round(p_mono, 4),
                "quasi_potential_lock_depth": round(lock_depth, 4),
                "combo_pd1_plus_mycaf_remission_prob": round(p_combo_caf, 4),
                "combo_pd1_plus_tam_remission_prob": round(p_combo_tam, 4),
            }
        )
    return pd.DataFrame(rows).set_index("patient_id")


def _score_visium_adata(adata: ad.AnnData) -> dict[str, Any]:
    """Compute compartment module scores and spatial exclusion from a loaded 10x Visium AnnData."""
    if "in_tissue" in adata.obs.columns:
        adata = adata[adata.obs["in_tissue"].astype(int) == 1].copy()

    X = adata.X
    if sparse.issparse(X):
        row_sums = np.asarray(X.sum(axis=1)).ravel() + 1e-8
    else:
        row_sums = np.sum(X, axis=1) + 1e-8

    var_upper = {str(g).upper(): i for i, g in enumerate(adata.var_names)}
    score_dict: dict[str, np.ndarray] = {}
    for comp, genes in (
        ("Tumor", TME_SIGNATURES["Tumor"]),
        ("TAM", TME_SIGNATURES["TAM"]),
        ("myCAF", TME_SIGNATURES["myCAF"]),
        ("CD8_T", TME_SIGNATURES["CD8_T"]),
    ):
        idxs = [var_upper[g] for g in genes if g in var_upper]
        if not idxs:
            score_dict[comp] = np.zeros(adata.n_obs, dtype=float)
            continue
        sub = X[:, idxs]
        if sparse.issparse(sub):
            sub_arr = sub.toarray().astype(float)
        else:
            sub_arr = np.asarray(sub, dtype=float)
        cpm = (sub_arr / row_sums[:, None]) * 1e4
        score_dict[comp] = np.mean(np.log1p(cpm), axis=1)

    if "spatial" in adata.obsm:
        coords = np.asarray(adata.obsm["spatial"], dtype=float)[:, :2]
    elif {"array_row", "array_col"}.issubset(adata.obs.columns):
        coords = adata.obs[["array_col", "array_row"]].to_numpy(dtype=float)
    else:
        raise ValueError("Visium AnnData missing spatial coordinates.")

    comp_df = pd.DataFrame(score_dict, index=adata.obs_names)
    res = compute_visium_spatial_barrier_and_proximity(coords, comp_df, k_neighbors=6)
    res["spot_scores_df"] = comp_df
    res["coords"] = coords
    return res


def extract_visium_spatial_evidence(vault_sc_dir: pathlib.Path) -> dict[str, Any]:
    """Extract multi-section human carcinoma 10x Visium spatial barrier evidence from GSE171351 & GSE203612."""
    section_results: list[dict[str, Any]] = []
    prox_accum: list[pd.DataFrame] = []

    # 1. GSE171351 human muscle-invasive bladder carcinoma Visium (.h5ad)
    vis_h5ad = vault_sc_dir / "GSE171351_BladderCa_Visium" / "GSE171351_combined_visium.h5ad"
    vis_gz = vault_sc_dir / "GSE171351_BladderCa_Visium" / "GSE171351_combined_visium.h5ad.gz"
    if not vis_h5ad.exists() and vis_gz.exists():
        import shutil

        with gzip.open(vis_gz, "rb") as fin, vis_h5ad.open("wb") as fout:
            shutil.copyfileobj(fin, fout)

    exemplar_section: dict[str, Any] | None = None
    if vis_h5ad.exists():
        adata_all = ad.read_h5ad(vis_h5ad)
        adata_all.var_names_make_unique()
        sample_col = "sampleID" if "sampleID" in adata_all.obs.columns else None
        if sample_col:
            for sid in adata_all.obs[sample_col].unique():
                sub = adata_all[adata_all.obs[sample_col] == sid].copy()
                if sub.n_obs < 150:
                    continue
                r = _score_visium_adata(sub)
                section_results.append(
                    {
                        "dataset": "GSE171351",
                        "section_id": f"GSE171351_{sid}",
                        "cancer_type": "Bladder_Carcinoma",
                        "n_spots": r["n_spots"],
                        "stromal_exclusion_index": r["stromal_exclusion_index"],
                        "exclusion_p": r["exclusion_mannwhitney_p"],
                        "tam_mycaf_prox": float(r["proximity_weights_df"].loc["TAM", "myCAF"]),
                        "tumor_cd8_prox": float(r["proximity_weights_df"].loc["Tumor", "CD8_T"]),
                    }
                )
                prox_accum.append(r["proximity_weights_df"])
                if exemplar_section is None or r["stromal_exclusion_index"] > exemplar_section["stromal_exclusion_index"]:
                    exemplar_section = {**r, "section_label": f"GSE171351_{sid}"}

    # 2. GSE203612 Pan-Cancer Human 10x Visium (Barkley et al., Nature Genetics 2022) if downloaded
    tar_203 = vault_sc_dir / "GSE203612_PanCancer_Visium" / "GSE203612_RAW.tar"
    if tar_203.exists() and tar_203.stat().st_size > 100_000_000:
        extract_dir = vault_sc_dir / "GSE203612_PanCancer_Visium" / "_extracted_human"
        extract_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_203, "r") as tf:
            members = tf.getmembers()
            # Filter strictly human carcinoma Visium sections (BRCA, OVCA, PAAD, UCEC, HNSC, LUSC/LUAD, SKCM; zero mouse/zebrafish/brain!)
            human_tags = ("BRCA", "OVCA", "PAAD", "PDAC", "UCEC", "HNSC", "LUAD", "LUSC", "SKCM", "GIST", "LIHC", "PRAD")
            h5_members = [
                m
                for m in members
                if m.name.endswith("filtered_feature_bc_matrix.h5")
                and any(t in m.name.upper() for t in human_tags)
                and not any(b in m.name.upper() for b in ("MOUSE", "ZEBRAFISH", "GBM", "GLIO", "BRAIN"))
            ][:8]
            pos_lookup = {
                m.name.split("_Vis_")[0]: m
                for m in members
                if ("tissue_positions" in m.name.lower()) and m.name.endswith((".csv.gz", ".txt.gz", ".csv"))
            }
            for hm in h5_members:
                prefix = hm.name.split("_Vis_")[0]
                h5_out = extract_dir / hm.name
                if not h5_out.exists():
                    with tf.extractfile(hm) as src:  # type: ignore[arg-type]
                        if src:
                            h5_out.write_bytes(src.read())
                pm = pos_lookup.get(prefix)
                pos_out = None
                if pm:
                    pos_out = extract_dir / pm.name
                    if not pos_out.exists():
                        with tf.extractfile(pm) as src:  # type: ignore[arg-type]
                            if src:
                                pos_out.write_bytes(src.read())
                try:
                    with h5py.File(h5_out, "r") as hf:
                        grp = hf["matrix"]
                        data = grp["data"][:]
                        indices = grp["indices"][:]
                        indptr = grp["indptr"][:]
                        shape = tuple(grp["shape"][:])
                        genes = [x.decode("utf-8") for x in grp["features"]["name"][:]]
                        barcodes = [x.decode("utf-8") for x in grp["barcodes"][:]]
                    # shape in 10x h5 is (n_genes, n_spots) CSC -> transpose to CSR (n_spots, n_genes)
                    mat_csc = sparse.csc_matrix((data, indices, indptr), shape=shape)
                    mat_csr = mat_csc.T.tocsr()
                    ad_sec = ad.AnnData(X=mat_csr)
                    ad_sec.obs_names = barcodes
                    ad_sec.var_names = genes
                    ad_sec.var_names_make_unique()
                    if pos_out and pos_out.exists():
                        pos_df = pd.read_csv(pos_out, header=None)
                        if str(pos_df.iloc[0, 0]).lower() == "barcode":
                            pos_df = pd.read_csv(pos_out)
                            pos_df = pos_df.set_index(pos_df.columns[0])
                            in_tis = pos_df.iloc[:, 0].astype(int)
                            coords_df = pos_df.iloc[:, -2:].astype(float)
                        else:
                            pos_df = pos_df.set_index(0)
                            in_tis = pos_df[1].astype(int)
                            coords_df = pos_df[[4, 5]].astype(float)
                        common_bc = ad_sec.obs_names.intersection(pos_df.index)
                        if len(common_bc) >= 150:
                            ad_sec = ad_sec[common_bc].copy()
                            ad_sec.obs["in_tissue"] = in_tis.loc[common_bc].to_numpy()
                            ad_sec.obsm["spatial"] = coords_df.loc[common_bc].to_numpy(dtype=float)
                            r = _score_visium_adata(ad_sec)
                            ctype = next((t for t in human_tags if t in prefix.upper()), "Solid_Carcinoma")
                            section_results.append(
                                {
                                    "dataset": "GSE203612",
                                    "section_id": prefix,
                                    "cancer_type": ctype,
                                    "n_spots": r["n_spots"],
                                    "stromal_exclusion_index": r["stromal_exclusion_index"],
                                    "exclusion_p": r["exclusion_mannwhitney_p"],
                                    "tam_mycaf_prox": float(r["proximity_weights_df"].loc["TAM", "myCAF"]),
                                    "tumor_cd8_prox": float(r["proximity_weights_df"].loc["Tumor", "CD8_T"]),
                                }
                            )
                            prox_accum.append(r["proximity_weights_df"])
                            if exemplar_section is None or r["stromal_exclusion_index"] > exemplar_section["stromal_exclusion_index"]:
                                exemplar_section = {**r, "section_label": prefix}
                except Exception:
                    continue

    if prox_accum:
        consensus_prox = sum(prox_accum) / float(len(prox_accum))
    else:
        consensus_prox = pd.DataFrame(
            np.ones((4, 4), dtype=float),
            index=["Tumor", "TAM", "myCAF", "CD8_T"],
            columns=["Tumor", "TAM", "myCAF", "CD8_T"],
        )

    # Map CD8_T physical proximity to CD8_Exhaustion compartment in 4-compartment grammar
    prox_4 = consensus_prox.rename(index={"CD8_T": "CD8_Exhaustion"}, columns={"CD8_T": "CD8_Exhaustion"})
    return {
        "sections_df": pd.DataFrame(section_results),
        "consensus_proximity_df": prox_4,
        "exemplar_section": exemplar_section,
    }
