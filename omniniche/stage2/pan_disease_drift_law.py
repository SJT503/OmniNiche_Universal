"""Stage II Task 2: Pan-Organ Normal Baseline tau_0 and Multi-Axis Dynamic Drift Statistical Law (Fig. 3).

Quantifies healthy pan-organ single-cell specificity tau_0 across 20,151 human genes x 154 cell types
(Human Protein Atlas) alongside stage-resolved disease single-cell specificity drift Delta_tau(t)
and organ delivery barrier compatibility across 6 primary human non-CNS remodeling disease axes,
testing predictive power across 29,791 target-resolved Phase II-IV clinical trials.
"""
from __future__ import annotations

import math
import pathlib
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

VITAL_BYSTANDER_KEYWORDS = (
    "cardiomyocyte",
    "hepatocyte",
    "proximal tubular",
    "alveolar cells",
    "podocyte",
    "enterocyte",
    "endothelial",
)


def compute_yanai_tau_matrix(expr_matrix: np.ndarray) -> np.ndarray:
    """Compute exact Yanai specificity index tau in [0, 1] for shape (n_genes, n_compartments)."""
    mat = np.maximum(np.asarray(expr_matrix, dtype=float), 0.0)
    if mat.ndim != 2 or mat.shape[1] < 2:
        raise ValueError("expr_matrix must be 2D with >= 2 columns.")
    n_cols = mat.shape[1]
    row_max = np.max(mat, axis=1)
    tau = np.zeros(mat.shape[0], dtype=float)
    nz = row_max > 1e-9
    normed = mat[nz] / row_max[nz, None]
    tau[nz] = np.sum(1.0 - normed, axis=1) / float(n_cols - 1)
    return np.clip(tau, 0.0, 1.0)


def compute_hpa_pan_organ_tau0(hpa_tsv_path: pathlib.Path) -> pd.DataFrame:
    """Compute healthy human pan-organ single-cell tau_0 across 154 normal cell types from HPA."""
    raw = pd.read_csv(hpa_tsv_path, sep="\t")
    val_col = "nTPM" if "nTPM" in raw.columns else raw.columns[-1]
    piv = raw.pivot_table(index="Gene name", columns="Cell type", values=val_col, aggfunc="mean").fillna(0.0)
    log_mat = np.log2(np.maximum(piv.to_numpy(dtype=float), 0.0) + 1.0)
    tau0 = compute_yanai_tau_matrix(log_mat)

    cell_types = list(piv.columns)
    max_idx = np.argmax(log_mat, axis=1)
    max_expr = np.max(log_mat, axis=1)
    mean_expr = np.mean(log_mat, axis=1)

    vital_cols = [
        i for i, ct in enumerate(cell_types) if any(k in ct.lower() for k in VITAL_BYSTANDER_KEYWORDS)
    ]
    if vital_cols:
        vital_mean = np.mean(log_mat[:, vital_cols], axis=1)
    else:
        vital_mean = mean_expr

    out = pd.DataFrame(
        {
            "gene_symbol": piv.index.astype(str),
            "tau0_healthy": np.round(tau0, 5),
            "hpa_primary_cell_type": [cell_types[i] for i in max_idx],
            "hpa_max_log2_ntpm": np.round(max_expr, 4),
            "hpa_mean_log2_ntpm": np.round(mean_expr, 4),
            "vital_bystander_hazard": np.round(vital_mean / (max_expr + 0.5), 4),
        }
    ).set_index("gene_symbol")
    return out


def compute_dynamic_drift_metrics(
    early_comp_df: pd.DataFrame,
    late_comp_df: pd.DataFrame,
    stromal_keywords: tuple[str, ...] = ("fibro", "myofibro", "stellate", "caf", "stroma", "mesenchym"),
) -> pd.DataFrame:
    """Compute cross-stage specificity drift Delta_tau(t), compartment switch, and dynamic retention."""
    common_genes = early_comp_df.index.intersection(late_comp_df.index)
    common_cols = [c for c in early_comp_df.columns if c in late_comp_df.columns]
    if len(common_cols) < 2:
        raise ValueError("Need >= 2 shared compartments between early and late stages.")

    e_mat = np.maximum(early_comp_df.loc[common_genes, common_cols].to_numpy(dtype=float), 0.0)
    l_mat = np.maximum(late_comp_df.loc[common_genes, common_cols].to_numpy(dtype=float), 0.0)

    tau_early = compute_yanai_tau_matrix(e_mat)
    tau_late = compute_yanai_tau_matrix(l_mat)
    delta_tau = tau_late - tau_early

    e_idx = np.argmax(e_mat, axis=1)
    l_idx = np.argmax(l_mat, axis=1)
    comp_switch = (e_idx != l_idx).astype(int)

    # Expression change inside the target's intended primary compartment
    rows_arange = np.arange(len(common_genes))
    e_prim_val = e_mat[rows_arange, e_idx]
    l_in_e_prim = l_mat[rows_arange, e_idx]
    prim_log2fc = l_in_e_prim - e_prim_val

    stromal_cols = [i for i, c in enumerate(common_cols) if any(k in c.lower() for k in stromal_keywords)]
    if stromal_cols:
        stromal_sum = np.sum(l_mat[:, stromal_cols], axis=1)
        total_sum = np.sum(l_mat, axis=1) + 1e-6
        stromal_share_late = stromal_sum / total_sum
    else:
        stromal_share_late = np.zeros(len(common_genes), dtype=float)

    # Dynamic retention rewards high late specificity, non-negative delta_tau, and preserved effector expression
    retention_raw = (
        0.50 * tau_late
        + 0.30 * (0.5 + 0.5 * np.tanh(2.5 * delta_tau))
        + 0.20 * (0.5 + 0.5 * np.tanh(prim_log2fc))
    )
    dynamic_retention = np.clip(retention_raw * (1.0 - 0.25 * comp_switch), 0.0, 1.0)

    return pd.DataFrame(
        {
            "tau_early": np.round(tau_early, 5),
            "tau_late": np.round(tau_late, 5),
            "delta_tau": np.round(delta_tau, 5),
            "primary_comp_early": [common_cols[i] for i in e_idx],
            "primary_comp_late": [common_cols[i] for i in l_idx],
            "compartment_switch": comp_switch,
            "primary_comp_delta_expr": np.round(prim_log2fc, 4),
            "stromal_ectopic_share_late": np.round(stromal_share_late, 4),
            "dynamic_retention_score": np.round(dynamic_retention, 5),
        },
        index=common_genes,
    )


def compute_modality_barrier_compatibility(
    drug_modality: str,
    disease_axis: str,
    stromal_share_late: float,
    compartment_switch: int,
) -> float:
    """Compute biophysical delivery barrier compatibility B(m, axis) in [0.15, 1.0]."""
    mod = (drug_modality or "").lower()
    is_macromolecule = any(k in mod for k in ("antibody", "protein", "cell", "gene", "enzyme"))
    is_oligo = "oligo" in mod

    # High ECM/stromal remodeling axes impose physical steric exclusion & TMDD sink on 150 kDa macromolecules
    dense_fibrotic_axes = {
        "Solid_Oncology": 0.28,
        "Hepatic_MASH_Fibrosis": 0.32,
        "Pulmonary_IPF_ARDS": 0.34,
        "Renal_CKD_AKI": 0.26,
        "Autoimmune_IBD_Rheum": 0.18,
    }
    base_penalty = dense_fibrotic_axes.get(disease_axis, 0.12)
    if is_macromolecule:
        penalty = base_penalty * (0.55 + 0.75 * float(stromal_share_late)) + 0.12 * int(compartment_switch)
    elif is_oligo:
        penalty = (0.10 if disease_axis == "Hepatic_MASH_Fibrosis" else 0.28) + 0.08 * int(compartment_switch)
    else:
        # Small molecule penetrates collagen matrix better, mild penalty only under compartment switch
        penalty = 0.06 * float(stromal_share_late) + 0.10 * int(compartment_switch)
    return float(np.clip(1.0 - penalty, 0.15, 1.0))


def _fit_cv_logistic_probs(X: np.ndarray, y: np.ndarray, seed: int = 42) -> np.ndarray:
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.zeros(len(y), dtype=float)
    for tr, va in skf.split(X, y):
        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X[tr])
        X_va = scaler.transform(X[va])
        clf = LogisticRegression(C=1.0, max_iter=500, solver="lbfgs")
        clf.fit(X_tr, y[tr])
        oof[va] = clf.predict_proba(X_va)[:, 1]
    return oof


def _wald_or_per_sd(x_1d: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    """Compute univariate logistic Odds Ratio per 1-SD increase + 95% Wald CI + two-sided P-value."""
    z = (x_1d - np.mean(x_1d)) / (np.std(x_1d) + 1e-9)
    X_des = np.column_stack([np.ones(len(z)), z])
    clf = LogisticRegression(C=1e4, max_iter=500, solver="lbfgs", fit_intercept=False)
    clf.fit(X_des, y)
    beta = clf.coef_[0]
    p = 1.0 / (1.0 + np.exp(-(X_des @ beta)))
    w = np.maximum(p * (1.0 - p), 1e-6)
    cov = np.linalg.pinv((X_des.T * w) @ X_des)
    se_b1 = math.sqrt(max(float(cov[1, 1]), 1e-12))
    b1 = float(beta[1])
    z_stat = b1 / se_b1
    pval = float(2.0 * stats.norm.sf(abs(z_stat)))
    or_val = math.exp(np.clip(b1, -10.0, 10.0))
    ci_lo = math.exp(np.clip(b1 - 1.96 * se_b1, -10.0, 10.0))
    ci_hi = math.exp(np.clip(b1 + 1.96 * se_b1, -10.0, 10.0))
    return or_val, ci_lo, ci_hi, pval


def evaluate_translational_progression_models(trials_df: pd.DataFrame) -> dict[str, Any]:
    """Compare Static Healthy tau_0 (Science 2026 baseline) vs OmniNiche Dynamic Drift x Barrier Law."""
    y = trials_df["clinical_success"].to_numpy(dtype=int)
    tau0 = trials_df["tau0_healthy"].to_numpy(dtype=float)
    tau_early = trials_df["tau_early"].to_numpy(dtype=float) if "tau_early" in trials_df else tau0
    delta_tau = trials_df["delta_tau"].to_numpy(dtype=float)
    barrier = trials_df["barrier_compatibility"].to_numpy(dtype=float)
    vital_haz = (
        trials_df["vital_bystander_hazard"].to_numpy(dtype=float)
        if "vital_bystander_hazard" in trials_df
        else np.zeros(len(y), dtype=float)
    )

    p_null = np.full(len(y), float(np.mean(y)))
    p_static_tau0 = _fit_cv_logistic_probs(tau0[:, None], y)
    p_static_early = _fit_cv_logistic_probs(np.column_stack([tau0, tau_early]), y)
    p_drift_only = _fit_cv_logistic_probs(np.column_stack([tau0, tau_early, delta_tau]), y)
    X_omni = np.column_stack([tau0, tau_early, delta_tau, barrier, tau_early * barrier, vital_haz])
    p_omniniche = _fit_cv_logistic_probs(X_omni, y)

    or_static, lo_s, hi_s, p_s = _wald_or_per_sd(tau0, y)
    # Composite OmniNiche linear predictor (out-of-fold logit)
    omni_logit = np.log(np.clip(p_omniniche, 1e-5, 1.0 - 1e-5) / (1.0 - np.clip(p_omniniche, 1e-5, 1.0 - 1e-5)))
    or_omni, lo_o, hi_o, p_o = _wald_or_per_sd(omni_logit, y)

    return {
        "n_trials_evaluated": int(len(y)),
        "n_positive_outcomes": int(np.sum(y)),
        "null_auroc": 0.5,
        "null_brier": float(brier_score_loss(y, p_null)),
        "static_tau0_auroc": float(roc_auc_score(y, p_static_tau0)),
        "static_tau0_auprc": float(average_precision_score(y, p_static_tau0)),
        "static_tau0_brier": float(brier_score_loss(y, p_static_tau0)),
        "static_tau0_odds_ratio_per_sd": float(or_static),
        "static_tau0_ci95": [float(lo_s), float(hi_s)],
        "static_tau0_pvalue": float(p_s),
        "static_disease_early_auroc": float(roc_auc_score(y, p_static_early)),
        "drift_no_barrier_auroc": float(roc_auc_score(y, p_drift_only)),
        "omniniche_dynamic_auroc": float(roc_auc_score(y, p_omniniche)),
        "omniniche_dynamic_auprc": float(average_precision_score(y, p_omniniche)),
        "omniniche_dynamic_brier": float(brier_score_loss(y, p_omniniche)),
        "omniniche_odds_ratio_per_sd": float(or_omni),
        "omniniche_ci95": [float(lo_o), float(hi_o)],
        "omniniche_pvalue": float(p_o),
        "oof_omniniche_prob": p_omniniche,
        "oof_static_tau0_prob": p_static_tau0,
    }
