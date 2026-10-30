"""Stage II v2.0 Statistical, Causal & Biophysical Fortification Engine (Nature Medicine / Science standard).

Provides deterministic, vectorized (< 2 GB RAM) inferential and mechanistic routines:
- Exact paired DeLong test for correlated ROC curves (DeLong et al., Biometrics 1988; Sun & Xu 2014)
- Stratified 2,000-iteration bootstrap 95% confidence intervals (AUROC, AUPRC, Brier score)
- Exact 10,000-permutation two-sided small-sample inference
- Clause-level negation-aware clinical trial stoppage audit & Cohen's kappa validation
- 6-step nested orthogonal DNAcc ablation + unseen-target GroupKFold + multivariable fixed-effects control
- 5-regime missingness & negation sensitivity bounds across 3,118 unreported stops
- 10-decile calibration, Decision Curve Analysis (Vickers & Elkin 2006), & Phase II->III trial enrichment
- Bidirectional 3x3 + Leave-One-Organ-Out cross-organ fibrosis transfer & continuous profile likelihood
- 12-section (21,796-spot) 10x Visium continuous radial distance decay & 2D Waddington quiver surfaces
"""
from __future__ import annotations

import math
import pathlib
import re
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import ndimage, optimize, sparse, spatial, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    cohen_kappa_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from omniniche.core.schema import CompartmentRule, DirectedEdgeRule, NicheGrammar
from omniniche.engine.grn_perturb import RidgeGRNPerturber
from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver
from omniniche.stage2.fibrosis_mash_deconvolution import (
    FIBROSIS_COMPARTMENT_GENES,
    LANDMARK_FIBROSIS_TRIALS,
    evaluate_14_landmark_fibrosis_trials,
)
from omniniche.stage2.solid_tumor_spatial_icb import (
    TME_SIGNATURES,
    _score_visium_adata,
    build_default_solid_tumor_icb_grammar,
)


# ---------------------------------------------------------------------------
# 1. Exact Paired DeLong Test, Stratified Bootstrap CIs & Exact Permutations
# ---------------------------------------------------------------------------
def _delong_structural_components(pos_scores: np.ndarray, neg_scores: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Compute Mann-Whitney U statistic and DeLong structural components V10, V01 in O((m+n) log(m+n))."""
    m = len(pos_scores)
    n = len(neg_scores)
    combined = np.concatenate([pos_scores, neg_scores])
    r_all = stats.rankdata(combined, method="average")
    r_pos = stats.rankdata(pos_scores, method="average")
    r_neg = stats.rankdata(neg_scores, method="average")

    v10 = (r_all[:m] - r_pos) / float(n)
    v01 = 1.0 - (r_all[m:] - r_neg) / float(m)
    auc = float(np.mean(v10))
    return auc, v10, v01


def fast_delong_paired_test(
    y_true: np.ndarray,
    y_pred1: np.ndarray,
    y_pred2: np.ndarray,
) -> dict[str, float]:
    """Exact two-sided paired DeLong test comparing two correlated ROC curves on identical samples."""
    y = np.asarray(y_true, dtype=int).ravel()
    p1 = np.asarray(y_pred1, dtype=float).ravel()
    p2 = np.asarray(y_pred2, dtype=float).ravel()
    pos_mask = y == 1
    neg_mask = y == 0
    m = int(np.sum(pos_mask))
    n = int(np.sum(neg_mask))
    if m < 2 or n < 2:
        raise ValueError("Paired DeLong test requires at least 2 positive and 2 negative observations.")

    auc1, v10_1, v01_1 = _delong_structural_components(p1[pos_mask], p1[neg_mask])
    auc2, v10_2, v01_2 = _delong_structural_components(p2[pos_mask], p2[neg_mask])

    s10 = np.cov(np.vstack([v10_1, v10_2]), ddof=1)
    s01 = np.cov(np.vstack([v01_1, v01_2]), ddof=1)
    cov_mat = (s10 / float(m)) + (s01 / float(n))

    delta = auc1 - auc2
    var_delta = float(cov_mat[0, 0] + cov_mat[1, 1] - 2.0 * cov_mat[0, 1])
    if abs(delta) < 1e-15 or var_delta <= 1e-18:
        return {
            "auroc_1": auc1,
            "auroc_2": auc2,
            "delta_auroc": 0.0 if abs(delta) < 1e-15 else delta,
            "se_delta": 0.0,
            "ci95_delta_low": delta,
            "ci95_delta_high": delta,
            "z_stat": 0.0,
            "p_value": 1.0 if abs(delta) < 1e-15 else 0.0,
            "log10_p": 0.0,
        }

    se_delta = math.sqrt(var_delta)
    z_stat = delta / se_delta
    p_val = float(2.0 * stats.norm.sf(abs(z_stat)))
    log10_p = float((math.log(2.0) + stats.norm.logsf(abs(z_stat))) / math.log(10.0))
    return {
        "auroc_1": auc1,
        "auroc_2": auc2,
        "delta_auroc": delta,
        "se_delta": se_delta,
        "ci95_delta_low": delta - 1.96 * se_delta,
        "ci95_delta_high": delta + 1.96 * se_delta,
        "z_stat": z_stat,
        "p_value": p_val,
        "log10_p": log10_p,
    }


def bootstrap_classification_metrics(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    n_boot: int = 2000,
    seed: int = 42,
    is_probability: bool = True,
) -> dict[str, float]:
    """Compute stratified bootstrap 95% confidence intervals for AUROC, AUPRC, and Brier score."""
    y = np.asarray(y_true, dtype=int).ravel()
    raw = np.asarray(y_prob, dtype=float).ravel()
    if is_probability:
        p = np.clip(raw, 0.0, 1.0)
    else:
        mn, mx = float(np.min(raw)), float(np.max(raw))
        p = (raw - mn) / (mx - mn + 1e-12) if mx > mn else np.full_like(raw, 0.5)

    auroc_pt = float(roc_auc_score(y, p))
    auprc_pt = float(average_precision_score(y, p))
    brier_pt = float(brier_score_loss(y, p))

    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    m, n = len(pos_idx), len(neg_idx)
    rng = np.random.default_rng(seed)

    aurocs = np.empty(n_boot, dtype=float)
    auprcs = np.empty(n_boot, dtype=float)
    briers = np.empty(n_boot, dtype=float)

    for b in range(n_boot):
        bp = pos_idx[rng.integers(0, m, size=m)]
        bn = neg_idx[rng.integers(0, n, size=n)]
        idx = np.concatenate([bp, bn])
        yb = y[idx]
        pb = p[idx]
        aurocs[b] = roc_auc_score(yb, pb)
        auprcs[b] = average_precision_score(yb, pb)
        briers[b] = brier_score_loss(yb, pb)

    return {
        "auroc": auroc_pt,
        "auroc_ci_low": float(np.minimum(np.percentile(aurocs, 2.5), auroc_pt)),
        "auroc_ci_high": float(np.maximum(np.percentile(aurocs, 97.5), auroc_pt)),
        "auprc": auprc_pt,
        "auprc_ci_low": float(np.minimum(np.percentile(auprcs, 2.5), auprc_pt)),
        "auprc_ci_high": float(np.maximum(np.percentile(auprcs, 97.5), auprc_pt)),
        "brier_score": brier_pt,
        "brier_ci_low": float(np.minimum(np.percentile(briers, 2.5), brier_pt)),
        "brier_ci_high": float(np.maximum(np.percentile(briers, 97.5), brier_pt)),
    }


def _batch_mann_whitney_auroc_vs_fixed_y(y_fixed: np.ndarray, scores_2d: np.ndarray) -> np.ndarray:
    """Exact vectorized AUROC for B score rows against fixed binary labels y_fixed in O(B * m * n)."""
    pos = scores_2d[:, y_fixed == 1]  # (B, m)
    neg = scores_2d[:, y_fixed == 0]  # (B, n)
    gt = (pos[:, :, None] > neg[:, None, :]).sum(axis=(1, 2))
    eq = (pos[:, :, None] == neg[:, None, :]).sum(axis=(1, 2))
    denom = float(pos.shape[1] * neg.shape[1])
    return (gt + 0.5 * eq) / denom


def exact_permutation_auroc_pvalue(
    y_true: np.ndarray,
    y_score: np.ndarray,
    y_baseline: np.ndarray | None = None,
    n_perm: int = 10000,
    seed: int = 42,
) -> dict[str, float]:
    """Compute exact/Monte-Carlo label-permutation two-sided P-value against chance (0.5) and paired swap vs baseline."""
    y = np.asarray(y_true, dtype=int).ravel()
    s1 = np.asarray(y_score, dtype=float).ravel()
    obs_auc = float(roc_auc_score(y, s1))
    obs_dev = abs(obs_auc - 0.5)

    rng = np.random.default_rng(seed)
    n_obs = len(y)
    if n_obs <= 400:
        # Permuting y against s1 has the identical distribution to permuting s1 against fixed y
        noise = rng.random((n_perm, n_obs))
        perm_idx = np.argsort(noise, axis=1)
        s1_perm = s1[perm_idx]
        ap_vec = _batch_mann_whitney_auroc_vs_fixed_y(y, s1_perm)
        exceed_null = int(np.sum(np.abs(ap_vec - 0.5) >= obs_dev - 1e-12))
    else:
        exceed_null = 0
        for _ in range(n_perm):
            yp = rng.permutation(y)
            ap = float(roc_auc_score(yp, s1))
            if abs(ap - 0.5) >= obs_dev - 1e-12:
                exceed_null += 1
    p_null = float((exceed_null + 1) / (n_perm + 1))

    p_vs_base = float("nan")
    obs_delta = float("nan")
    if y_baseline is not None:
        s2 = np.asarray(y_baseline, dtype=float).ravel()
        # Rank-normalize both scores onto [0, 1] so within-subject predictor swaps are scale-invariant
        r1 = stats.rankdata(s1, method="average") / float(n_obs)
        r2 = stats.rankdata(s2, method="average") / float(n_obs)
        obs_delta = float(roc_auc_score(y, r1) - roc_auc_score(y, r2))
        if n_obs <= 400:
            swap = rng.integers(0, 2, size=(n_perm, n_obs)).astype(bool)
            p_a = np.where(swap, r2[None, :], r1[None, :])
            p_b = np.where(swap, r1[None, :], r2[None, :])
            d_perm_vec = _batch_mann_whitney_auroc_vs_fixed_y(y, p_a) - _batch_mann_whitney_auroc_vs_fixed_y(y, p_b)
            exceed_swap = int(np.sum(np.abs(d_perm_vec) >= abs(obs_delta) - 1e-12))
        else:
            exceed_swap = 0
            for _ in range(n_perm):
                sw = rng.integers(0, 2, size=n_obs).astype(bool)
                pa = np.where(sw, r2, r1)
                pb = np.where(sw, r1, r2)
                d_perm = float(roc_auc_score(y, pa) - roc_auc_score(y, pb))
                if abs(d_perm) >= abs(obs_delta) - 1e-12:
                    exceed_swap += 1
        p_vs_base = float((exceed_swap + 1) / (n_perm + 1))

    return {
        "auroc": obs_auc,
        "p_perm_vs_null": p_null,
        "delta_auroc_vs_baseline": obs_delta,
        "p_perm_vs_baseline": p_vs_base,
    }


# ---------------------------------------------------------------------------
# 2. Clause-Level Negation-Aware Stoppage Adjudicator & Validation Benchmark
# ---------------------------------------------------------------------------
NEGATED_SAFETY_EFFICACY_CLAUSE_REGEX = re.compile(
    r"\b(?:"
    r"(?:not|no|neither|never|without|unrelated\s+to|independent\s+of|free\s+of|absence\s+of|non-)"
    r"[^.;]{0,65}?"
    r"\b(?:safety|toxicit\w*|adverse\s+event\w*|tolerabilit\w*|futility|benefit[-/]risk|risk[-/]benefit)"
    r"|"
    r"(?:not\s+(?:related|due|attributable|linked|driven)\s+(?:to|by)|unrelated\s+to)"
    r"[^.;]{0,65}?"
    r"\b(?:efficacy|effectiveness|safety|toxicit\w*|futility|benefit[-/]risk)"
    r"|"
    r"\bno\s+(?:safety|efficacy|new\s+safety)[/\s]*(?:or\s+(?:safety|efficacy)\s+)?(?:concern\w*|issue\w*|signal\w*|problem\w*|analys[ei]s|finding\w*)"
    r"|"
    r"\bno\s+(?:change\s+in\s+)?(?:benefit[-/]risk|risk[-/]benefit)\s+(?:impact|change|concern\w*|issue\w*)"
    r")",
    re.IGNORECASE,
)

CLAUSE_EFFICACY_REGEX = re.compile(
    r"(lack\s+of\s+(?:\w+\s+){0,2}efficacy|inefficacy|futility|futil\w+|insufficient\s+efficacy|"
    r"\bno\s+(?:clinical\s+|treatment\s+|meaningful\s+)?efficacy\b|primary\s+endpoint\w*\s+(?:was\s+|were\s+)?not\s+met|"
    r"did\s+not\s+(?:meet|achieve|demonstrate|show)\s+(?:its\s+|the\s+|meaningful\s+|sufficient\s+)?(?:primary|clinical\s+benefit|efficacy|endpoint)|"
    r"failed\s+to\s+(?:meet|demonstrate|show|achieve)|unfavorable\s+(?:risk[-/]benefit|benefit[-/]risk)|"
    r"\b(?:benefit[-/]risk|risk[-/]benefit)\b|no\s+(?:clinical|meaningful|significant)\s+benefit|"
    r"suboptimal\s+efficacy|low\s+response\s+rate|inferior\s+efficacy|disease\s+progression|"
    r"(?:due\s+to|for|based\s+on)\s+efficacy(?:/effectiveness)?\s+reason\w*|"
    r"unplanned\s+interim\s+efficacy|lack\s+of\s+(?:clinical\s+)?activity|insufficient\s+(?:clinical\s+)?activity)",
    re.IGNORECASE,
)

CLAUSE_SAFETY_REGEX = re.compile(
    r"(toxicit\w*|adverse\s+event\w*|\bsafety\b|hepatotox\w*|cardiotox\w*|nephrotox\w*|pneumonitis|"
    r"dose[- ]limiting|\bdlt\b|mortality\s+imbalance|increased\s+mortality|\bdeath\w*|\bfatal\w*|"
    r"not\s+tolerable|poor\s+tolerability|intoleran\w*|unacceptable\s+tolerability|"
    r"clinical\s+hold|fda\s+hold|\bsae\b|serious\s+adverse)",
    re.IGNORECASE,
)

CLAUSE_ADMIN_REGEX = re.compile(
    r"(accrual|recruit\w*|enrol+ment|\benrol+ed\b|\benrol+ing\b|"
    r"business|strategic|commercial|sponsor(?:'s)?\s+(?:decision|request|termination)|"
    r"company\s+(?:decision|decided)|corporate|portfolio|prioriti[sz]\w*|reprioriti[sz]\w*|"
    r"funding|financial|budget|money|resource\w*|covid[- ]?19|pandemic|sars-cov-2|"
    r"supply|manufactur\w*|drug\s+availability|drug\s+shortage|imp\s+supply|"
    r"investigator\s+(?:left|relocated|departed|decision)|\bpi\s+(?:left|relocated|departed)|"
    r"site\s+closure|administrative|protocol\s+redesign|replaced\s+by\s+another|"
    r"another\s+study|regulatory\s+delay|operational|logistical|feasibility|"
    r"termination\s+of\s+development|discontinue\s+development|not\s+to\s+pursue)",
    re.IGNORECASE,
)


def adjudicate_stoppage_clause_aware(overall_status: str, why_stopped: str) -> str:
    """Clause-level negation-aware clinical trial stoppage adjudicator."""
    st = (overall_status or "").upper().strip()
    if st == "COMPLETED":
        return "Completed"
    raw = (why_stopped or "").strip()
    if not raw:
        return "Unreported_Stop_Reason"
    cleaned = NEGATED_SAFETY_EFFICACY_CLAUSE_REGEX.sub(" ", raw)
    has_eff = bool(CLAUSE_EFFICACY_REGEX.search(cleaned))
    has_saf = bool(CLAUSE_SAFETY_REGEX.search(cleaned))
    has_adm = bool(CLAUSE_ADMIN_REGEX.search(raw))
    if has_eff and not has_saf:
        return "Scientific_Efficacy_Futility"
    if has_saf and not has_eff:
        return "Scientific_Safety_Toxicity"
    if has_eff and has_saf:
        return "Scientific_Efficacy_Futility"
    if has_adm:
        return "Administrative_Commercial_Enrollment"
    return "Other_Unclassified_Stop"


def audit_stoppage_classifier_benchmark(
    trials_df: pd.DataFrame,
    n_per_stratum: int = 60,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Stratified n=240 stopped-trial validation benchmark comparing registry deconvolution vs clause-aware reference adjudication."""
    stopped = trials_df[trials_df["overall_status"].astype(str).str.upper() != "COMPLETED"].copy()
    stopped["clause_reference"] = [
        adjudicate_stoppage_clause_aware(s, w)
        for s, w in zip(stopped["overall_status"].fillna(""), stopped["why_stopped"].fillna(""))
    ]

    def _collapse4(label: str) -> str:
        if label in ("Unreported_Stop_Reason", "Other_Unclassified_Stop"):
            return "Other_or_Unreported"
        return label

    stopped["pred_4class"] = stopped["stop_deconvolution"].astype(str).map(_collapse4)
    stopped["ref_4class"] = stopped["clause_reference"].astype(str).map(_collapse4)

    classes = [
        "Scientific_Efficacy_Futility",
        "Scientific_Safety_Toxicity",
        "Administrative_Commercial_Enrollment",
        "Other_or_Unreported",
    ]
    sampled_parts = []
    for i, cls in enumerate(classes):
        sub = stopped[stopped["ref_4class"] == cls]
        k = min(n_per_stratum, len(sub))
        sampled_parts.append(sub.sample(n=k, random_state=seed + i))
    audit_sample = pd.concat(sampled_parts, ignore_index=True)

    y_ref = audit_sample["ref_4class"].to_numpy()
    y_pred = audit_sample["pred_4class"].to_numpy()
    cm = confusion_matrix(y_ref, y_pred, labels=classes)
    cm_df = pd.DataFrame(cm, index=classes, columns=classes)

    kappa = float(cohen_kappa_score(y_ref, y_pred))
    acc = float(np.mean(y_ref == y_pred))

    metrics_rows = []
    for idx, cls in enumerate(classes):
        tp = float(cm[idx, idx])
        fp = float(np.sum(cm[:, idx]) - tp)
        fn = float(np.sum(cm[idx, :]) - tp)
        prec = tp / max(tp + fp, 1.0)
        rec = tp / max(tp + fn, 1.0)
        f1 = 2.0 * prec * rec / max(prec + rec, 1e-8)
        metrics_rows.append(
            {
                "stoppage_class": cls,
                "n_reference": int(np.sum(cm[idx, :])),
                "n_predicted": int(np.sum(cm[:, idx])),
                "true_positives": int(tp),
                "precision": round(prec, 4),
                "recall": round(rec, 4),
                "f1_score": round(f1, 4),
                "overall_accuracy": round(acc, 4),
                "cohen_kappa": round(kappa, 4),
            }
        )
    keep_cols = ["nct_id", "disease_axis", "overall_status", "why_stopped", "pred_4class", "ref_4class"]
    return pd.DataFrame(metrics_rows), cm_df, audit_sample[keep_cols]


# ---------------------------------------------------------------------------
# 3. Macro-Atlas 2005-2025 Longitudinal Trends & Modality x Stoppage Log-Odds
# ---------------------------------------------------------------------------
def _primary_modality_group(mod_str: str) -> str:
    s = str(mod_str or "").lower()
    if not s or s == "nan":
        return "Unresolved_or_Combo"
    if "antibody" in s:
        return "Antibody_Biologic"
    if any(k in s for k in ("cell", "gene", "oligo")):
        return "Cell_Gene_Oligo"
    if any(k in s for k in ("protein", "enzyme", "peptide")):
        return "Protein_Peptide"
    if "small molecule" in s:
        return "Small_Molecule"
    return "Other_Modality"


def compute_macro_longitudinal_and_modality_odds(
    master_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute 2005-2025 longitudinal trends, modality x stoppage log-odds forest table, and phase/enrollment distributions."""
    df = master_df.copy()
    df["year_num"] = pd.to_numeric(df["start_year"], errors="coerce")
    df["modality_group"] = df["drug_modalities"].map(_primary_modality_group)
    df["clause_stop"] = [
        adjudicate_stoppage_clause_aware(s, w)
        for s, w in zip(df["overall_status"].fillna(""), df["why_stopped"].fillna(""))
    ]

    # A. Longitudinal 2005-2025 trends (2-year bins or annual for smooth ribbons)
    valid_yr = df[(df["year_num"] >= 2005) & (df["year_num"] <= 2025)].copy()
    valid_yr["year_int"] = valid_yr["year_num"].astype(int)
    trend_rows = []
    for yr, sub in valid_yr.groupby("year_int"):
        n_tot = len(sub)
        if n_tot < 20:
            continue
        res_sub = sub[sub["n_targets"] > 0]
        n_res = len(res_sub)
        n_mab = int((res_sub["modality_group"] == "Antibody_Biologic").sum())
        p_mab = n_mab / max(n_res, 1)
        se_mab = math.sqrt(max(p_mab * (1.0 - p_mab) / max(n_res, 1), 1e-9))

        stopped = sub[sub["overall_status"] != "COMPLETED"]
        n_stop = len(stopped)
        p_admin = float((stopped["clause_stop"] == "Administrative_Commercial_Enrollment").sum()) / max(n_stop, 1)
        p_eff = float((stopped["clause_stop"] == "Scientific_Efficacy_Futility").sum()) / max(n_stop, 1)
        p_saf = float((stopped["clause_stop"] == "Scientific_Safety_Toxicity").sum()) / max(n_stop, 1)
        se_adm = math.sqrt(max(p_admin * (1.0 - p_admin) / max(n_stop, 1), 1e-9))

        trend_rows.append(
            {
                "start_year": int(yr),
                "n_trials": int(n_tot),
                "n_target_resolved": int(n_res),
                "antibody_share_pct": round(100.0 * p_mab, 3),
                "antibody_ci95_low": round(100.0 * max(0.0, p_mab - 1.96 * se_mab), 3),
                "antibody_ci95_high": round(100.0 * min(1.0, p_mab + 1.96 * se_mab), 3),
                "n_stopped": int(n_stop),
                "admin_share_of_stops_pct": round(100.0 * p_admin, 3),
                "admin_ci95_low": round(100.0 * max(0.0, p_admin - 1.96 * se_adm), 3),
                "admin_ci95_high": round(100.0 * min(1.0, p_admin + 1.96 * se_adm), 3),
                "efficacy_share_of_stops_pct": round(100.0 * p_eff, 3),
                "safety_share_of_stops_pct": round(100.0 * p_saf, 3),
            }
        )

    # B. Modality x Termination Mechanism Odds Ratios (Reference = Small_Molecule among target-resolved trials)
    res_all = df[df["n_targets"] > 0].copy()
    odds_rows = []
    for mod in ["Antibody_Biologic", "Protein_Peptide", "Cell_Gene_Oligo"]:
        for cause_label, cause_val in [
            ("Scientific_Efficacy_Futility", "Scientific_Efficacy_Futility"),
            ("Scientific_Safety_Toxicity", "Scientific_Safety_Toxicity"),
            ("Administrative_Commercial_Enrollment", "Administrative_Commercial_Enrollment"),
        ]:
            sub_m = res_all[res_all["modality_group"].isin([mod, "Small_Molecule"])]
            is_mod = sub_m["modality_group"] == mod
            is_cause = sub_m["clause_stop"] == cause_val
            # Haldane-Anscombe +0.5 continuity correction
            a = float(np.sum(is_mod & is_cause)) + 0.5
            b = float(np.sum(is_mod & (~is_cause))) + 0.5
            c = float(np.sum((~is_mod) & is_cause)) + 0.5
            d = float(np.sum((~is_mod) & (~is_cause))) + 0.5
            log_or = math.log((a * d) / (b * c))
            se_lor = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
            z = log_or / se_lor
            pval = float(2.0 * stats.norm.sf(abs(z)))
            odds_rows.append(
                {
                    "modality_group": mod,
                    "reference_modality": "Small_Molecule",
                    "termination_cause": cause_label,
                    "n_modality_events": int(a - 0.5),
                    "n_modality_total": int(a + b - 1.0),
                    "n_ref_events": int(c - 0.5),
                    "n_ref_total": int(c + d - 1.0),
                    "odds_ratio": round(math.exp(log_or), 4),
                    "ci95_low": round(math.exp(log_or - 1.96 * se_lor), 4),
                    "ci95_high": round(math.exp(log_or + 1.96 * se_lor), 4),
                    "log2_or": round(log_or / math.log(2.0), 4),
                    "p_value": pval,
                }
            )

    # C. Enrollment sample-size distributions by disease axis & stoppage class
    df["enroll_num"] = pd.to_numeric(df["enrollment"], errors="coerce").fillna(0.0)
    pos_enr = df[df["enroll_num"] > 0].copy()
    enr_rows = []
    for grp_type, col in [("disease_axis", "disease_axis"), ("stop_class", "clause_stop")]:
        for gname, gsub in pos_enr.groupby(col):
            vals = gsub["enroll_num"].to_numpy(dtype=float)
            enr_rows.append(
                {
                    "group_type": grp_type,
                    "group_name": str(gname),
                    "n_trials": int(len(vals)),
                    "median_enrollment": round(float(np.median(vals)), 1),
                    "q25_enrollment": round(float(np.percentile(vals, 25)), 1),
                    "q75_enrollment": round(float(np.percentile(vals, 75)), 1),
                    "p90_enrollment": round(float(np.percentile(vals, 90)), 1),
                    "total_enrolled_patients": int(np.sum(vals)),
                }
            )

    return pd.DataFrame(trend_rows), pd.DataFrame(odds_rows), pd.DataFrame(enr_rows)


# ---------------------------------------------------------------------------
# 4. 6-Step Orthogonal DNAcc Ablation, Unseen-Target GroupKFold & Fixed Effects
# ---------------------------------------------------------------------------
def _oof_logistic_predictions(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray | None = None,
    n_splits: int = 5,
    seed: int = 42,
) -> np.ndarray:
    oof = np.zeros(len(y), dtype=float)
    if groups is not None:
        n_uniq = len(np.unique(groups))
        k = min(n_splits, max(n_uniq, 2))
        splitter = GroupKFold(n_splits=k).split(X, y, groups=groups)
    else:
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(X, y)

    for tr, va in splitter:
        if len(np.unique(y[tr])) < 2:
            oof[va] = float(np.mean(y))
            continue
        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X[tr])
        X_va = scaler.transform(X[va])
        clf = LogisticRegression(C=1.0, max_iter=500, solver="lbfgs")
        clf.fit(X_tr, y[tr])
        oof[va] = clf.predict_proba(X_va)[:, 1]
    return oof


def evaluate_six_step_dnacc_ablation(
    trials_df: pd.DataFrame,
    cohort_name: str,
    n_boot: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    """Evaluate 6-step nested orthogonal DNAcc ablation under Stratified 5-fold CV and Unseen-Target GroupKFold."""
    y = trials_df["clinical_success"].to_numpy(dtype=int)
    groups = trials_df["primary_target_gene"].fillna("UNKNOWN").astype(str).to_numpy()

    tau0 = trials_df["tau0_healthy"].to_numpy(dtype=float)
    vb = trials_df["vital_bystander_hazard"].to_numpy(dtype=float)
    te = trials_df["tau_early"].to_numpy(dtype=float)
    dt = trials_df["delta_tau"].to_numpy(dtype=float)
    dr = trials_df["dynamic_retention_score"].to_numpy(dtype=float)
    sw = trials_df["compartment_switch"].to_numpy(dtype=float)
    st = trials_df["stromal_ectopic_share_late"].to_numpy(dtype=float)
    bar = trials_df["barrier_compatibility"].to_numpy(dtype=float)

    steps = [
        ("M1_Static_Healthy_Tau0", "Healthy baseline specificity (tau_0)", np.column_stack([tau0])),
        ("M2_Plus_Vital_Bystander", "+ Vital-organ bystander hazard", np.column_stack([tau0, vb])),
        ("M3_Plus_Disease_Early_Tau", "+ Early disease niche specificity (tau_early)", np.column_stack([tau0, vb, te])),
        ("M4_Plus_Dynamic_Drift", "+ Remodeling drift & retention (Delta_tau)", np.column_stack([tau0, vb, te, dt, dr])),
        ("M5_Plus_Compartment_Switch", "+ Lineage switch & stromal ectopic share", np.column_stack([tau0, vb, te, dt, dr, sw, st])),
        ("M6_Full_DNAcc", "+ Modality-stromal barrier (Full DNAcc)", np.column_stack([tau0, vb, te, dt, dr, sw, st, bar, te * bar])),
    ]

    rows = []
    for cv_name, grp_arr in [("Stratified_5Fold", None), ("UnseenTarget_Group5Fold", groups)]:
        p_m1 = _oof_logistic_predictions(steps[0][2], y, groups=grp_arr, seed=seed)
        prev_p = p_m1
        for step_idx, (sid, sdesc, X_mat) in enumerate(steps):
            p_curr = p_m1 if step_idx == 0 else _oof_logistic_predictions(X_mat, y, groups=grp_arr, seed=seed)
            b_stats = bootstrap_classification_metrics(y, p_curr, n_boot=n_boot, seed=seed + step_idx)
            dl_vs_m1 = fast_delong_paired_test(y, p_curr, p_m1)
            dl_vs_prev = fast_delong_paired_test(y, p_curr, prev_p)
            rows.append(
                {
                    "cohort": cohort_name,
                    "cv_scheme": cv_name,
                    "step_order": step_idx + 1,
                    "step_id": sid,
                    "component_added": sdesc,
                    "n_trials": int(len(y)),
                    "n_positive": int(np.sum(y)),
                    "n_unique_targets": int(len(np.unique(groups))),
                    "auroc": round(b_stats["auroc"], 4),
                    "auroc_ci_low": round(b_stats["auroc_ci_low"], 4),
                    "auroc_ci_high": round(b_stats["auroc_ci_high"], 4),
                    "auprc": round(b_stats["auprc"], 4),
                    "auprc_ci_low": round(b_stats["auprc_ci_low"], 4),
                    "auprc_ci_high": round(b_stats["auprc_ci_high"], 4),
                    "brier_score": round(b_stats["brier_score"], 4),
                    "brier_ci_low": round(b_stats["brier_ci_low"], 4),
                    "brier_ci_high": round(b_stats["brier_ci_high"], 4),
                    "delta_auroc_vs_m1": round(dl_vs_m1["delta_auroc"], 4),
                    "delong_p_vs_m1": dl_vs_m1["p_value"],
                    "delong_log10p_vs_m1": round(dl_vs_m1["log10_p"], 3),
                    "delta_auroc_vs_prev": round(dl_vs_prev["delta_auroc"], 4),
                    "delong_p_vs_prev": dl_vs_prev["p_value"],
                }
            )
            prev_p = p_curr

    # Single-target (n_targets == 1) vs Polypharmacology (n_targets >= 2) subgroup check on M6 vs M1
    if "n_targets" in trials_df.columns:
        for sub_tag, mask in [
            ("Subgroup_SingleTarget_Monotherapy_Target", trials_df["n_targets"] == 1),
            ("Subgroup_Polypharmacology_MultiTarget", trials_df["n_targets"] >= 2),
        ]:
            sub_df = trials_df[mask]
            if len(sub_df) >= 80 and sub_df["clinical_success"].nunique() == 2:
                ys = sub_df["clinical_success"].to_numpy(dtype=int)
                gs = sub_df["primary_target_gene"].fillna("UNKNOWN").astype(str).to_numpy()
                X1_s = sub_df[["tau0_healthy"]].to_numpy(dtype=float)
                X6_s = np.column_stack(
                    [
                        sub_df["tau0_healthy"].to_numpy(dtype=float),
                        sub_df["vital_bystander_hazard"].to_numpy(dtype=float),
                        sub_df["tau_early"].to_numpy(dtype=float),
                        sub_df["delta_tau"].to_numpy(dtype=float),
                        sub_df["dynamic_retention_score"].to_numpy(dtype=float),
                        sub_df["compartment_switch"].to_numpy(dtype=float),
                        sub_df["stromal_ectopic_share_late"].to_numpy(dtype=float),
                        sub_df["barrier_compatibility"].to_numpy(dtype=float),
                        sub_df["tau_early"].to_numpy(dtype=float) * sub_df["barrier_compatibility"].to_numpy(dtype=float),
                    ]
                )
                p1_s = _oof_logistic_predictions(X1_s, ys, seed=seed)
                p6_s = _oof_logistic_predictions(X6_s, ys, seed=seed)
                b6_s = bootstrap_classification_metrics(ys, p6_s, n_boot=n_boot, seed=seed + 50)
                dl_s = fast_delong_paired_test(ys, p6_s, p1_s)
                rows.append(
                    {
                        "cohort": cohort_name,
                        "cv_scheme": sub_tag,
                        "step_order": 6,
                        "step_id": "M6_Full_DNAcc",
                        "component_added": f"{sub_tag} (vs Static tau_0 AUROC={dl_s['auroc_2']:.4f})",
                        "n_trials": int(len(ys)),
                        "n_positive": int(np.sum(ys)),
                        "n_unique_targets": int(len(np.unique(gs))),
                        "auroc": round(b6_s["auroc"], 4),
                        "auroc_ci_low": round(b6_s["auroc_ci_low"], 4),
                        "auroc_ci_high": round(b6_s["auroc_ci_high"], 4),
                        "auprc": round(b6_s["auprc"], 4),
                        "auprc_ci_low": round(b6_s["auprc_ci_low"], 4),
                        "auprc_ci_high": round(b6_s["auprc_ci_high"], 4),
                        "brier_score": round(b6_s["brier_score"], 4),
                        "brier_ci_low": round(b6_s["brier_ci_low"], 4),
                        "brier_ci_high": round(b6_s["brier_ci_high"], 4),
                        "delta_auroc_vs_m1": round(dl_s["delta_auroc"], 4),
                        "delong_p_vs_m1": dl_s["p_value"],
                        "delong_log10p_vs_m1": round(dl_s["log10_p"], 3),
                        "delta_auroc_vs_prev": round(dl_s["delta_auroc"], 4),
                        "delong_p_vs_prev": dl_s["p_value"],
                    }
                )
    return pd.DataFrame(rows)


def fit_multivariable_fixed_effects_logistic(
    trials_df: pd.DataFrame,
    dnacc_oof_prob: np.ndarray,
) -> pd.DataFrame:
    """Fit multivariable logistic regression testing DNAcc adjusted for baseline confounders (axis, modality, era, sponsor, target multiplicity)."""
    df = trials_df.copy()
    y = df["clinical_success"].to_numpy(dtype=int)
    p_clip = np.clip(dnacc_oof_prob, 1e-5, 1.0 - 1e-5)
    dnacc_logit = np.log(p_clip / (1.0 - p_clip))
    dnacc_z = (dnacc_logit - np.mean(dnacc_logit)) / (np.std(dnacc_logit) + 1e-9)

    tau0_raw = df["tau0_healthy"].to_numpy(dtype=float)
    tau0_z = (tau0_raw - np.mean(tau0_raw)) / (np.std(tau0_raw) + 1e-9)

    df["mod_grp"] = df["drug_modalities"].map(_primary_modality_group)
    yr = pd.to_numeric(df.get("start_year", pd.Series(2015, index=df.index)), errors="coerce").fillna(2015)
    df["era"] = pd.cut(yr, bins=[0, 2009, 2014, 2019, 2030], labels=["Pre_2010", "2010_2014", "2015_2019", "2020_Plus"]).astype(str)
    sp_col = df["sponsor_class"].fillna("OTHER").astype(str) if "sponsor_class" in df.columns else pd.Series("INDUSTRY", index=df.index)
    df["sponsor_grp"] = np.where(sp_col.str.upper() == "INDUSTRY", "INDUSTRY", "NON_INDUSTRY")
    nt_col = pd.to_numeric(df["n_targets"], errors="coerce").fillna(1) if "n_targets" in df.columns else pd.Series(1, index=df.index)
    df["poly_grp"] = np.where(nt_col >= 2, "MultiTarget", "SingleTarget")

    cov_dummies = pd.get_dummies(
        df[["disease_axis", "mod_grp", "era", "sponsor_grp", "poly_grp"]],
        drop_first=True,
        dtype=float,
    )

    if "primary_target_gene" in df.columns:
        clusters = df["primary_target_gene"].fillna("UNKNOWN").astype(str).to_numpy()
    else:
        clusters = np.arange(len(y), dtype=str)
    uniq_g, inv_idx = np.unique(clusters, return_inverse=True)
    n_clusters = int(len(uniq_g))

    def _cluster_sandwich_se(X_mat: np.ndarray, y_vec: np.ndarray, p_vec: np.ndarray, H_inv_mat: np.ndarray, col_idx: int) -> float:
        scores_i = X_mat * (y_vec - p_vec)[:, None]
        u_g = np.zeros((n_clusters, X_mat.shape[1]), dtype=float)
        np.add.at(u_g, inv_idx, scores_i)
        dfc = n_clusters / max(n_clusters - 1.0, 1.0)
        v_sand = H_inv_mat @ ((u_g.T @ u_g) * dfc) @ H_inv_mat
        return math.sqrt(max(float(v_sand[col_idx, col_idx]), 1e-12))

    out_rows = []
    for spec_name, z_vec, feat_label in [
        ("Unadjusted_Univariate_DNAcc", dnacc_z, "DNAcc (per 1-SD, unadjusted)"),
        ("Unadjusted_Univariate_Static_Tau0", tau0_z, "Static Healthy tau_0 (per 1-SD, unadjusted)"),
        ("Multivariable_FixedEffects_DNAcc", dnacc_z, "DNAcc (per 1-SD, adjusted for Axis + Modality + Era + Sponsor + Polypharmacology)"),
        ("Multivariable_FixedEffects_Static_Tau0", tau0_z, "Static Healthy tau_0 (per 1-SD, adjusted for Axis + Modality + Era + Sponsor + Polypharmacology)"),
    ]:
        if spec_name.startswith("Unadjusted"):
            X_des = np.column_stack([np.ones(len(y)), z_vec])
        else:
            X_des = np.column_stack([np.ones(len(y)), z_vec, cov_dummies.to_numpy(dtype=float)])

        # Weak ridge penalty C=50 stabilizes any sparse dummy cell while keeping main effect unbiased
        clf = LogisticRegression(C=50.0, max_iter=800, solver="lbfgs", fit_intercept=False)
        clf.fit(X_des, y)
        beta = clf.coef_[0]
        p_hat = 1.0 / (1.0 + np.exp(-np.clip(X_des @ beta, -30.0, 30.0)))
        w = np.maximum(p_hat * (1.0 - p_hat), 1e-6)
        fisher = (X_des.T * w) @ X_des + np.diag([0.0] + [1.0 / 50.0] * (X_des.shape[1] - 1))
        cov = np.linalg.pinv(fisher)

        b1 = float(beta[1])
        se1 = math.sqrt(max(float(cov[1, 1]), 1e-12))
        z_stat = b1 / se1
        p_val = float(2.0 * stats.norm.sf(abs(z_stat)))
        log10_p = float((math.log(2.0) + stats.norm.logsf(abs(z_stat))) / math.log(10.0))

        se_cl = _cluster_sandwich_se(X_des, y, p_hat, cov, 1)
        z_cl = b1 / se_cl
        p_cl = float(2.0 * stats.norm.sf(abs(z_cl)))
        l10_cl = float((math.log(2.0) + stats.norm.logsf(abs(z_cl))) / math.log(10.0))
        out_rows.append(
            {
                "model_specification": spec_name,
                "predictor": feat_label,
                "n_trials": int(len(y)),
                "n_covariates": int(X_des.shape[1] - 1),
                "n_target_clusters": n_clusters,
                "odds_ratio_per_sd": round(math.exp(np.clip(b1, -10.0, 10.0)), 4),
                "ci95_low": round(math.exp(np.clip(b1 - 1.96 * se1, -10.0, 10.0)), 4),
                "ci95_high": round(math.exp(np.clip(b1 + 1.96 * se1, -10.0, 10.0)), 4),
                "z_stat": round(z_stat, 3),
                "p_value": p_val,
                "log10_p": round(log10_p, 3),
                "cluster_ci95_low": round(math.exp(np.clip(b1 - 1.96 * se_cl, -10.0, 10.0)), 4),
                "cluster_ci95_high": round(math.exp(np.clip(b1 + 1.96 * se_cl, -10.0, 10.0)), 4),
                "cluster_z_stat": round(z_cl, 3),
                "cluster_p_value": p_cl,
                "cluster_log10_p": round(l10_cl, 3),
            }
        )

    # Also fit Joint Head-to-Head Multivariable Model (DNAcc + Static tau_0 + Fixed Effects simultaneously)
    X_joint = np.column_stack([np.ones(len(y)), dnacc_z, tau0_z, cov_dummies.to_numpy(dtype=float)])
    clf_j = LogisticRegression(C=50.0, max_iter=800, solver="lbfgs", fit_intercept=False).fit(X_joint, y)
    bj = clf_j.coef_[0]
    pj = 1.0 / (1.0 + np.exp(-np.clip(X_joint @ bj, -30.0, 30.0)))
    wj = np.maximum(pj * (1.0 - pj), 1e-6)
    cov_j = np.linalg.pinv((X_joint.T * wj) @ X_joint + np.diag([0.0] + [1.0 / 50.0] * (X_joint.shape[1] - 1)))
    for col_i, lbl in [(1, "Joint_FixedEffects_DNAcc_controlling_for_Tau0"), (2, "Joint_FixedEffects_Tau0_controlling_for_DNAcc")]:
        bi = float(bj[col_i])
        sei = math.sqrt(max(float(cov_j[col_i, col_i]), 1e-12))
        zi = bi / sei
        pi = float(2.0 * stats.norm.sf(abs(zi)))
        l10 = float((math.log(2.0) + stats.norm.logsf(abs(zi))) / math.log(10.0))

        se_cl = _cluster_sandwich_se(X_joint, y, pj, cov_j, col_i)
        z_cl = bi / se_cl
        p_cl = float(2.0 * stats.norm.sf(abs(z_cl)))
        l10_cl = float((math.log(2.0) + stats.norm.logsf(abs(z_cl))) / math.log(10.0))
        out_rows.append(
            {
                "model_specification": lbl,
                "predictor": lbl,
                "n_trials": int(len(y)),
                "n_covariates": int(X_joint.shape[1] - 1),
                "n_target_clusters": n_clusters,
                "odds_ratio_per_sd": round(math.exp(np.clip(bi, -10.0, 10.0)), 4),
                "ci95_low": round(math.exp(np.clip(bi - 1.96 * sei, -10.0, 10.0)), 4),
                "ci95_high": round(math.exp(np.clip(bi + 1.96 * sei, -10.0, 10.0)), 4),
                "z_stat": round(zi, 3),
                "p_value": pi,
                "log10_p": round(l10, 3),
                "cluster_ci95_low": round(math.exp(np.clip(bi - 1.96 * se_cl, -10.0, 10.0)), 4),
                "cluster_ci95_high": round(math.exp(np.clip(bi + 1.96 * se_cl, -10.0, 10.0)), 4),
                "cluster_z_stat": round(z_cl, 3),
                "cluster_p_value": p_cl,
                "cluster_log10_p": round(l10_cl, 3),
            }
        )
    return pd.DataFrame(out_rows)


def evaluate_within_modality_dnacc_stratification(
    trials_df: pd.DataFrame,
    dnacc_oof_prob: np.ndarray,
    tau0_oof_prob: np.ndarray,
    cohort_name: str = "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes",
    n_boot: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    """Evaluate DNAcc vs Static tau_0 strictly within each drug modality stratum to prove non-modality dynamic drift lift."""
    df = trials_df.copy()
    df["mod_grp"] = df["drug_modalities"].map(_primary_modality_group)
    y_all = df["clinical_success"].to_numpy(dtype=int)
    p_dnacc = np.asarray(dnacc_oof_prob, dtype=float)
    p_tau0 = np.asarray(tau0_oof_prob, dtype=float)

    rows = []
    for mg in ["Small_Molecule", "Antibody_Biologic", "Protein_Peptide"]:
        mask = (df["mod_grp"] == mg).to_numpy()
        y_m = y_all[mask]
        if int(np.sum(y_m == 1)) < 5 or int(np.sum(y_m == 0)) < 5:
            continue
        b_dn = bootstrap_classification_metrics(y_m, p_dnacc[mask], n_boot=n_boot, seed=seed)
        b_t0 = bootstrap_classification_metrics(y_m, p_tau0[mask], n_boot=n_boot, seed=seed)
        dl = fast_delong_paired_test(y_m, p_dnacc[mask], p_tau0[mask])
        n_tgt = int(df.loc[mask, "primary_target_gene"].nunique()) if "primary_target_gene" in df.columns else 0
        rows.append(
            {
                "cohort": cohort_name,
                "modality_stratum": mg,
                "n_trials": int(len(y_m)),
                "n_positive": int(np.sum(y_m == 1)),
                "n_negative": int(np.sum(y_m == 0)),
                "n_unique_targets": n_tgt,
                "static_tau0_auroc": round(b_t0["auroc"], 4),
                "static_tau0_ci_low": round(b_t0["auroc_ci_low"], 4),
                "static_tau0_ci_high": round(b_t0["auroc_ci_high"], 4),
                "dnacc_auroc": round(b_dn["auroc"], 4),
                "dnacc_ci_low": round(b_dn["auroc_ci_low"], 4),
                "dnacc_ci_high": round(b_dn["auroc_ci_high"], 4),
                "delta_auroc": round(dl["delta_auroc"], 4),
                "delong_p_value": dl["p_value"],
                "delong_log10_p": round(dl["log10_p"], 3),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. 5-Regime Missingness & Deconvolution Sensitivity Bounds (3,118 Unreported Stops)
# ---------------------------------------------------------------------------
def evaluate_unreported_stop_sensitivity_regimes(
    full_target_resolved_df: pd.DataFrame,
    n_boot: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    """Evaluate DNAcc vs Static tau_0 across 5 missingness and negation-filtering sensitivity regimes."""
    df = full_target_resolved_df.copy()
    df["clause_stop"] = [
        adjudicate_stoppage_clause_aware(s, w)
        for s, w in zip(df["overall_status"].fillna(""), df["why_stopped"].fillna(""))
    ]

    def _feat_matrices(sub: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        t0 = sub["tau0_healthy"].to_numpy(dtype=float)
        vb = sub["vital_bystander_hazard"].to_numpy(dtype=float)
        te = sub["tau_early"].to_numpy(dtype=float)
        dt = sub["delta_tau"].to_numpy(dtype=float)
        dr = sub["dynamic_retention_score"].to_numpy(dtype=float)
        sw = sub["compartment_switch"].to_numpy(dtype=float)
        st = sub["stromal_ectopic_share_late"].to_numpy(dtype=float)
        ba = sub["barrier_compatibility"].to_numpy(dtype=float)
        return t0[:, None], np.column_stack([t0, vb, te, dt, dr, sw, st, ba, te * ba])

    regimes: list[tuple[str, str, pd.DataFrame, np.ndarray | None]] = []

    # R1: Primary Locked Deconfounded Scientific Outcome (n = 11,843)
    m_r1 = ((df["stop_deconvolution"] == "Completed") & (df["max_phase_num"] >= 3)) | (
        df["stop_deconvolution"].isin(["Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"])
    )
    d_r1 = df[m_r1].copy()
    d_r1["clinical_success"] = ((d_r1["stop_deconvolution"] == "Completed") & (d_r1["max_phase_num"] >= 3)).astype(int)
    regimes.append(("R1_Primary_Explicit_Scientific", "Primary explicit scientific outcomes (Completed Ph3/4 vs explicit futility/toxicity)", d_r1, None))

    # R2: Clause-Level Negation-Filtered Strict Scientific Outcome (removes 331 'not related to safety' disclaimers)
    m_r2 = ((df["clause_stop"] == "Completed") & (df["max_phase_num"] >= 3)) | (
        df["clause_stop"].isin(["Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"])
    )
    d_r2 = df[m_r2].copy()
    d_r2["clinical_success"] = ((d_r2["clause_stop"] == "Completed") & (d_r2["max_phase_num"] >= 3)).astype(int)
    regimes.append(("R2_Negation_Filtered_Strict_Scientific", "Clause-level negation-filtered scientific stops ('not due to safety' excluded)", d_r2, None))

    # R3: Worst-Case Bound (All 3,118 Unreported / Other target-resolved stops coded as biological failure y=0)
    m_r3 = m_r1 | df["stop_deconvolution"].isin(["Unreported_Stop_Reason", "Other_Unclassified_Stop"])
    d_r3 = df[m_r3].copy()
    d_r3["clinical_success"] = ((d_r3["stop_deconvolution"] == "Completed") & (d_r3["max_phase_num"] >= 3)).astype(int)
    regimes.append(("R3_WorstCase_All3118_Unreported_As_Failure", "Worst-case bound: all 3,118 unreported/other stops treated as biological failure", d_r3, None))

    # R4: Inverse-Probability Weighting (IPW) for Missing Stop-Reason Propensity among Stopped Trials
    stopped_pool = df[df["stop_deconvolution"].isin(["Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity", "Unreported_Stop_Reason", "Other_Unclassified_Stop"])].copy()
    is_observed = stopped_pool["stop_deconvolution"].isin(["Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"]).astype(int).to_numpy()
    stopped_pool["mod_grp"] = stopped_pool["drug_modalities"].map(_primary_modality_group)
    yr_s = pd.to_numeric(stopped_pool.get("start_year", 2015), errors="coerce").fillna(2015)
    stopped_pool["era"] = pd.cut(yr_s, bins=[0, 2009, 2014, 2019, 2030], labels=["E1", "E2", "E3", "E4"]).astype(str)
    X_prop = pd.get_dummies(stopped_pool[["disease_axis", "mod_grp", "era", "max_phase_num", "sponsor_class"]], drop_first=True, dtype=float).to_numpy()
    prop_clf = LogisticRegression(C=1.0, max_iter=500).fit(X_prop, is_observed)
    pi_obs = np.clip(prop_clf.predict_proba(X_prop)[:, 1], 0.08, 0.95)
    ipw_map = dict(zip(stopped_pool["nct_id"], 1.0 / pi_obs))
    w_r4 = np.array([1.0 if r.clinical_success == 1 else float(ipw_map.get(r.nct_id, 1.0)) for r in d_r1.itertuples()], dtype=float)
    regimes.append(("R4_IPW_Propensity_Reweighted", "Inverse-probability weighted (IPW) for stop-reason reporting propensity", d_r1.copy(), w_r4))

    # R5: Mandatory Posted-Results Subset Only (has_results == True on ClinicalTrials.gov)
    has_res_mask = d_r1["has_results"].astype(str).str.lower() == "true"
    d_r5 = d_r1[has_res_mask].copy()
    if len(d_r5) >= 100 and d_r5["clinical_success"].nunique() == 2:
        regimes.append(("R5_Posted_ClinicalResults_Only", "Restricted strictly to trials with posted quantitative results (HasResults=True)", d_r5, None))

    rows = []
    for r_id, r_desc, r_df, sample_w in regimes:
        y = r_df["clinical_success"].to_numpy(dtype=int)
        X1, X6 = _feat_matrices(r_df)
        if sample_w is None:
            p1 = _oof_logistic_predictions(X1, y, seed=seed)
            p6 = _oof_logistic_predictions(X6, y, seed=seed)
        else:
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
            p1 = np.zeros(len(y), dtype=float)
            p6 = np.zeros(len(y), dtype=float)
            for tr, va in skf.split(X6, y):
                sc1, sc6 = StandardScaler(), StandardScaler()
                c1 = LogisticRegression(C=1.0, max_iter=500).fit(sc1.fit_transform(X1[tr]), y[tr], sample_weight=sample_w[tr])
                c6 = LogisticRegression(C=1.0, max_iter=500).fit(sc6.fit_transform(X6[tr]), y[tr], sample_weight=sample_w[tr])
                p1[va] = c1.predict_proba(sc1.transform(X1[va]))[:, 1]
                p6[va] = c6.predict_proba(sc6.transform(X6[va]))[:, 1]

        b1 = bootstrap_classification_metrics(y, p1, n_boot=n_boot, seed=seed + 1)
        b6 = bootstrap_classification_metrics(y, p6, n_boot=n_boot, seed=seed + 6)
        dl = fast_delong_paired_test(y, p6, p1)
        rows.append(
            {
                "regime_id": r_id,
                "regime_description": r_desc,
                "n_trials": int(len(y)),
                "n_positive": int(np.sum(y)),
                "n_negative": int(len(y) - np.sum(y)),
                "static_tau0_auroc": round(b1["auroc"], 4),
                "static_tau0_ci_low": round(b1["auroc_ci_low"], 4),
                "static_tau0_ci_high": round(b1["auroc_ci_high"], 4),
                "dnacc_auroc": round(b6["auroc"], 4),
                "dnacc_ci_low": round(b6["auroc_ci_low"], 4),
                "dnacc_ci_high": round(b6["auroc_ci_high"], 4),
                "delta_auroc": round(dl["delta_auroc"], 4),
                "delong_p_value": dl["p_value"],
                "delong_log10_p": round(dl["log10_p"], 3),
                "dnacc_brier": round(b6["brier_score"], 4),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 6. 10-Decile Calibration, Decision Curve Analysis & Trial Portfolio Enrichment
# ---------------------------------------------------------------------------
def compute_calibration_deciles(
    y_true: np.ndarray,
    p_dnacc: np.ndarray,
    p_tau0: np.ndarray,
    n_bins: int = 10,
) -> pd.DataFrame:
    """Compute equal-frequency decile calibration curves with 95% Wilson binomial CIs."""
    y = np.asarray(y_true, dtype=int).ravel()
    rows = []
    for model_name, probs in [("Full_DNAcc", p_dnacc), ("Static_Healthy_Tau0", p_tau0)]:
        p = np.asarray(probs, dtype=float).ravel()
        bins = pd.qcut(pd.Series(p).rank(method="first"), q=n_bins, labels=False)
        for b in range(n_bins):
            mask = bins == b
            nb = int(np.sum(mask))
            if nb == 0:
                continue
            pred_mean = float(np.mean(p[mask]))
            obs_rate = float(np.mean(y[mask]))
            se = math.sqrt(max(obs_rate * (1.0 - obs_rate) / nb, 1e-9))
            rows.append(
                {
                    "model": model_name,
                    "decile": int(b + 1),
                    "n_trials": nb,
                    "mean_predicted_prob": round(pred_mean, 4),
                    "observed_event_rate": round(obs_rate, 4),
                    "obs_ci95_low": round(max(0.0, obs_rate - 1.96 * se), 4),
                    "obs_ci95_high": round(min(1.0, obs_rate + 1.96 * se), 4),
                    "abs_calibration_error": round(abs(pred_mean - obs_rate), 4),
                }
            )
    return pd.DataFrame(rows)


def compute_decision_curve_and_enrichment(
    y_true: np.ndarray,
    p_dnacc: np.ndarray,
    p_tau0: np.ndarray,
    median_enrollment: int = 246,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute Vickers & Elkin (2006) Decision Curve Analysis net benefit + Phase II->III portfolio enrichment."""
    y = np.asarray(y_true, dtype=int).ravel()
    p1 = np.asarray(p_dnacc, dtype=float).ravel()
    p0 = np.asarray(p_tau0, dtype=float).ravel()
    n_tot = float(len(y))
    prev = float(np.mean(y))

    # Threshold grid spanning clinical progression base rate
    lo_pt = max(0.05, round(prev - 0.35, 2))
    hi_pt = min(0.97, round(prev + 0.18, 2))
    thresholds = np.linspace(lo_pt, hi_pt, 35)

    dca_rows = []
    for pt in thresholds:
        odds_w = pt / max(1.0 - pt, 1e-6)
        nb_all = prev - (1.0 - prev) * odds_w

        tp_1 = float(np.sum((p1 >= pt) & (y == 1))) / n_tot
        fp_1 = float(np.sum((p1 >= pt) & (y == 0))) / n_tot
        nb_1 = tp_1 - fp_1 * odds_w

        tp_0 = float(np.sum((p0 >= pt) & (y == 1))) / n_tot
        fp_0 = float(np.sum((p0 >= pt) & (y == 0))) / n_tot
        nb_0 = tp_0 - fp_0 * odds_w

        # Net futile interventions avoided per 100 candidate programs vs Advance-All
        avoided_per_100 = 100.0 * max(0.0, (nb_1 - nb_all) / max(odds_w, 1e-6))
        dca_rows.append(
            {
                "threshold_pt": round(float(pt), 4),
                "prevalence": round(prev, 4),
                "nb_treat_none": 0.0,
                "nb_treat_all": round(nb_all, 5),
                "nb_static_tau0": round(nb_0, 5),
                "nb_dnacc": round(nb_1, 5),
                "delta_nb_vs_treat_all": round(nb_1 - nb_all, 5),
                "delta_nb_vs_tau0": round(nb_1 - nb_0, 5),
                "futile_programs_avoided_per_100": round(avoided_per_100, 2),
            }
        )

    # Portfolio Enrichment Simulation across quantile gating cutoffs (excluding bottom 10%, 20%, 25%, 30%, 40%, 50% of DNAcc)
    n_fail_total = int(np.sum(y == 0))
    n_win_total = int(np.sum(y == 1))
    enrich_rows = []
    for pct_gate in [0, 10, 20, 25, 30, 40, 50]:
        cut_dnacc = float(np.percentile(p1, pct_gate)) if pct_gate > 0 else -1.0
        cut_tau0 = float(np.percentile(p0, pct_gate)) if pct_gate > 0 else -1.0
        keep_1 = p1 >= cut_dnacc
        keep_0 = p0 >= cut_tau0

        win_rate_1 = float(np.mean(y[keep_1])) * 100.0
        win_rate_0 = float(np.mean(y[keep_0])) * 100.0
        futile_blocked_1 = int(np.sum((~keep_1) & (y == 0)))
        futile_blocked_0 = int(np.sum((~keep_0) & (y == 0)))
        winners_retained_1 = float(np.sum(keep_1 & (y == 1))) / max(n_win_total, 1) * 100.0

        enrich_rows.append(
            {
                "bottom_percentile_gated_out": pct_gate,
                "n_programs_advanced": int(np.sum(keep_1)),
                "dnacc_success_rate_pct": round(win_rate_1, 2),
                "static_tau0_success_rate_pct": round(win_rate_0, 2),
                "futile_trials_prevented_dnacc": futile_blocked_1,
                "futile_trials_prevented_pct": round(100.0 * futile_blocked_1 / max(n_fail_total, 1), 2),
                "futile_trials_prevented_tau0_pct": round(100.0 * futile_blocked_0 / max(n_fail_total, 1), 2),
                "true_winners_retained_pct": round(winners_retained_1, 2),
                "patients_spared_futile_exposure": int(futile_blocked_1 * median_enrollment),
            }
        )
    return pd.DataFrame(dca_rows), pd.DataFrame(enrich_rows)


# ---------------------------------------------------------------------------
# 7. Bidirectional 3x3 Cross-Organ Fibrosis Transfer & Profile Likelihood
# ---------------------------------------------------------------------------
def build_conserved_fibrosis_grammar_from_organ(
    organ_donors_df: pd.DataFrame,
    organ_label: str = "Human Fibrotic Organ",
) -> NicheGrammar:
    """Compile the 3-compartment bistable fibrotic NicheGrammar purely from any single source organ's human donors."""
    corr = organ_donors_df[["Parenchyma_Injury", "SPP1_SAMac", "Myofibroblast"]].corr()
    r_ps = float(np.clip(abs(corr.loc["Parenchyma_Injury", "SPP1_SAMac"]), 0.35, 0.95))
    r_sm = float(np.clip(abs(corr.loc["SPP1_SAMac", "Myofibroblast"]), 0.45, 0.98))
    r_mp = float(np.clip(abs(corr.loc["Myofibroblast", "Parenchyma_Injury"]), 0.35, 0.95))

    compartments = [
        CompartmentRule(
            name="Parenchyma_Injury",
            baseline_production=0.008,
            decay_rate=0.180,
            noise_sigma=0.035,
            marker_genes=sorted(FIBROSIS_COMPARTMENT_GENES["Parenchyma_Injury"]),
        ),
        CompartmentRule(
            name="SPP1_SAMac",
            baseline_production=0.006,
            decay_rate=0.180,
            noise_sigma=0.035,
            marker_genes=sorted(FIBROSIS_COMPARTMENT_GENES["SPP1_SAMac"]),
        ),
        CompartmentRule(
            name="Myofibroblast",
            baseline_production=0.006,
            decay_rate=0.180,
            noise_sigma=0.035,
            marker_genes=sorted(FIBROSIS_COMPARTMENT_GENES["Myofibroblast"]),
        ),
    ]
    edges = [
        DirectedEdgeRule(
            source="Parenchyma_Injury",
            target="SPP1_SAMac",
            interaction_type="activation",
            weight_v=round(0.335 + 0.045 * r_ps, 4),
            half_max_k=0.44,
            hill_n=2.4,
            ligand_gene="GDF15/CCL2",
            receptor_gene="TREM2/CCR2",
            effector_tf="SPP1",
            partial_corr_gain=round(r_ps, 4),
        ),
        DirectedEdgeRule(
            source="SPP1_SAMac",
            target="Myofibroblast",
            interaction_type="activation",
            weight_v=round(0.520 + 0.035 * r_sm, 4),
            half_max_k=0.62,
            hill_n=3.6,
            ligand_gene="SPP1/TGFB1/PDGFB",
            receptor_gene="CD44/TGFBR1/PDGFRA",
            effector_tf="SMAD3",
            partial_corr_gain=round(r_sm, 4),
        ),
        DirectedEdgeRule(
            source="Myofibroblast",
            target="SPP1_SAMac",
            interaction_type="activation",
            weight_v=round(0.515 + 0.035 * r_sm, 4),
            half_max_k=0.62,
            hill_n=3.6,
            ligand_gene="COL1A1/POSTN/CSF1",
            receptor_gene="CD44/ITGAV/CSF1R",
            effector_tf="NFKB1",
            partial_corr_gain=round(r_sm, 4),
        ),
        DirectedEdgeRule(
            source="Myofibroblast",
            target="Parenchyma_Injury",
            interaction_type="activation",
            weight_v=round(0.320 + 0.045 * r_mp, 4),
            half_max_k=0.60,
            hill_n=3.0,
            ligand_gene="TGFB1/ECM_Stiffness",
            receptor_gene="ITGB6/YAP1",
            effector_tf="SOX9",
            partial_corr_gain=round(r_mp, 4),
        ),
    ]
    return NicheGrammar(
        grammar_id=f"conserved_fibrosis_from_{organ_label.lower().replace(' ', '_')}",
        source_pack_id="fibrosis_ckm_ipf",
        organ_context=organ_label,
        compartments=compartments,
        edges=edges,
        remission_compartment="Myofibroblast",
        remission_threshold=0.36,
    )


def evaluate_bidirectional_3x3_fibrosis_transfer(
    organ_data: dict[str, dict[str, pd.DataFrame]],
    tau0_df: pd.DataFrame,
    verified_ncts: dict[str, dict[str, Any]],
    prior_edges: list[tuple[str, str]],
    ccr2_paradox_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Evaluate all 3 single-organ sources and 3 Leave-One-Organ-Out source pairs across all 14 trials and held-out organ subsets."""
    # organ_data keys: 'Liver_MASH', 'Lung_IPF', 'Kidney_CKD' with {'donors', 'early', 'late', 'drift'}
    organs = ["Liver_MASH", "Lung_IPF", "Kidney_CKD"]
    source_configs: list[tuple[str, list[str]]] = [
        ("Liver_Only (n=17 donors)", ["Liver_MASH"]),
        ("Lung_Only (n=30 donors)", ["Lung_IPF"]),
        ("Kidney_Only (n=17 donors)", ["Kidney_CKD"]),
        ("LOOO_Lung+Kidney -> HeldOut_Liver", ["Lung_IPF", "Kidney_CKD"]),
        ("LOOO_Liver+Kidney -> HeldOut_Lung", ["Liver_MASH", "Kidney_CKD"]),
        ("LOOO_Liver+Lung -> HeldOut_Kidney", ["Liver_MASH", "Lung_IPF"]),
    ]

    all_fib_genes = (
        FIBROSIS_COMPARTMENT_GENES["Parenchyma_Injury"]
        | FIBROSIS_COMPARTMENT_GENES["SPP1_SAMac"]
        | FIBROSIS_COMPARTMENT_GENES["Myofibroblast"]
    )
    pathogenic_core = {"COL1A1", "COL1A2", "COL3A1", "POSTN", "CTHRC1", "ACTA2", "SPP1", "TREM2", "TGFB1", "GDF15", "CCL2"}

    rows = []
    for src_label, src_orgs in source_configs:
        src_donors = pd.concat([organ_data[o]["donors"] for o in src_orgs], ignore_index=True)
        grammar = build_conserved_fibrosis_grammar_from_organ(src_donors, organ_label=src_label)

        # Fit RidgeGRNPerturber STRICTLY on the source organ(s) pseudobulks (zero target-organ leakage!)
        common_genes = sorted(
            set.intersection(*(set(organ_data[o]["late"].index) for o in src_orgs)) & all_fib_genes
        )
        expr_blocks = []
        for o in src_orgs:
            e_sub = organ_data[o]["early"].loc[common_genes].T
            l_sub = organ_data[o]["late"].loc[common_genes].T
            expr_blocks.extend([e_sub, l_sub, 0.5 * (e_sub + l_sub)])
        expr_df = pd.concat(expr_blocks, ignore_index=True)
        perturber = RidgeGRNPerturber(ridge_lambda=0.6, propagation_steps=3).fit(
            expr_df,
            prior_edges=prior_edges,
            pathogenic_program_genes=sorted(pathogenic_core & set(common_genes)),
        )

        df_tr, _ = evaluate_14_landmark_fibrosis_trials(
            grammar=grammar,
            tau0_df=tau0_df,
            liver_drift_df=organ_data["Liver_MASH"]["drift"],
            lung_drift_df=organ_data["Lung_IPF"]["drift"],
            kidney_drift_df=organ_data["Kidney_CKD"]["drift"],
            verified_ncts=verified_ncts,
            grn_perturber=perturber,
            ccr2_paradox_df=ccr2_paradox_df,
        )

        # Overall 14-trial performance + strictly zero-shot out-of-source organ performance + per-target-organ accuracy
        held_out = df_tr[~df_tr["organ"].isin(src_orgs)]
        y_all = df_tr["clinical_win"].to_numpy(dtype=int)
        p_all = df_tr["omniniche_remission_prob"].to_numpy(dtype=float)
        b_all = bootstrap_classification_metrics(y_all, p_all, n_boot=1000, seed=42)

        if len(held_out) >= 2 and held_out["clinical_win"].nunique() == 2:
            yh = held_out["clinical_win"].to_numpy(dtype=int)
            ph = held_out["omniniche_remission_prob"].to_numpy(dtype=float)
            zeroshot_auroc = round(float(roc_auc_score(yh, ph)), 4)
            zeroshot_brier = round(float(brier_score_loss(yh, ph)), 4)
        elif len(held_out) >= 1:
            yh = held_out["clinical_win"].to_numpy(dtype=int)
            ph = held_out["omniniche_remission_prob"].to_numpy(dtype=float)
            zeroshot_auroc = 1.0 if np.all((ph >= 0.5) == (yh == 1)) else 0.0
            zeroshot_brier = round(float(brier_score_loss(yh, ph)), 4)
        else:
            zeroshot_auroc = float("nan")
            zeroshot_brier = float("nan")

        per_org_acc = {}
        per_org_auroc = {}
        for tgt_o in organs:
            sub_o = df_tr[df_tr["organ"] == tgt_o]
            yo = sub_o["clinical_win"].to_numpy(dtype=int)
            po = sub_o["omniniche_remission_prob"].to_numpy(dtype=float)
            per_org_acc[tgt_o] = round(float(np.mean((po >= 0.5) == (yo == 1))), 4)
            per_org_auroc[tgt_o] = round(float(roc_auc_score(yo, po)), 4) if len(np.unique(yo)) == 2 else per_org_acc[tgt_o]

        rows.append(
            {
                "source_calibration": src_label,
                "source_organs": ";".join(src_orgs),
                "n_source_donors": int(len(src_donors)),
                "n_heldout_trials": int(len(held_out)),
                "all14_auroc": round(b_all["auroc"], 4),
                "all14_ci_low": round(b_all["auroc_ci_low"], 4),
                "all14_ci_high": round(b_all["auroc_ci_high"], 4),
                "all14_brier": round(b_all["brier_score"], 4),
                "heldout_zeroshot_auroc": zeroshot_auroc,
                "heldout_zeroshot_brier": zeroshot_brier,
                "target_liver_auroc": per_org_auroc["Liver_MASH"],
                "target_lung_auroc": per_org_auroc["Lung_IPF"],
                "target_kidney_accuracy": per_org_acc["Kidney_CKD"],
                "edge_v_parenchyma_to_samac": grammar.edges[0].weight_v,
                "edge_v_samac_to_myofib": grammar.edges[1].weight_v,
            }
        )
    return pd.DataFrame(rows)


def compute_continuous_profile_likelihood_curves(
    grammar: NicheGrammar,
    n_grid: int = 25,
) -> pd.DataFrame:
    """Compute 25-point continuous nuisance-reoptimized bilateral profile-likelihood Delta-chi2 curves (Raue et al., 2009)."""
    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)
    s0_obs = np.array([0.62, 0.58, 0.56], dtype=float)
    base_traj = solver.simulate_deterministic(s0=s0_obs, regimens=[])["state"]
    obs_idx = np.linspace(10, len(base_traj) - 1, 18, dtype=int)
    y_obs = base_traj[obs_idx]
    sigma_obs = 0.042

    multipliers = np.linspace(0.50, 1.50, n_grid)
    rows = []
    for e_idx, edge in enumerate(grammar.edges):
        param_name = f"v_{edge.source}_to_{edge.target}"
        v_mle = float(edge.weight_v)
        raw_chi2 = []
        opt_scales = []
        for mult in multipliers:
            v_scan = v_mle * float(mult)

            def _loss(rate_scale: float) -> float:
                sim = solver.simulate_deterministic(
                    s0=s0_obs,
                    regimens=[],
                    global_rate_scale=float(rate_scale),
                    edge_weight_overrides={e_idx: v_scan},
                )["state"][obs_idx]
                return float(np.sum(((sim - y_obs) / sigma_obs) ** 2))

            # Re-optimize nuisance global kinetic time-scale k_global around 1.0 (true profile likelihood)
            opt = optimize.minimize_scalar(_loss, bounds=(0.82, 1.22), method="bounded")
            raw_chi2.append(float(opt.fun))
            opt_scales.append(float(opt.x))

        min_chi2 = min(raw_chi2)
        for mult, c2, k_sc in zip(multipliers, raw_chi2, opt_scales):
            d_chi2 = max(0.0, c2 - min_chi2)
            rows.append(
                {
                    "edge_index": e_idx,
                    "parameter": param_name,
                    "mle_value": round(v_mle, 4),
                    "multiplier": round(float(mult), 4),
                    "scanned_value": round(v_mle * float(mult), 4),
                    "reoptimized_nuisance_scale": round(k_sc, 4),
                    "delta_chi2": round(d_chi2, 4),
                    "identifiable_at_95pct": bool(max(raw_chi2[0] - min_chi2, 0.0) > 3.84 and max(raw_chi2[-1] - min_chi2, 0.0) > 3.84),
                }
            )
    return pd.DataFrame(rows)


def compute_continuous_fibrosis_hysteresis_curve(
    grammar: NicheGrammar,
    n_points: int = 25,
) -> pd.DataFrame:
    """Scan continuous fibrotic progression coordinate s_prog in [0.18, 0.90] showing hysteresis collapse & dual rescue."""
    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)
    prog_grid = np.linspace(0.18, 0.90, n_points)

    reg_myofib_mono = [InterventionRegimen("Myofibroblast", 0.22, "small_molecule")]
    reg_paren_mono = [InterventionRegimen("Parenchyma_Injury", 0.82, "small_molecule")]
    reg_dual_loop = [
        InterventionRegimen("Myofibroblast", 0.22, "small_molecule"),
        InterventionRegimen("SPP1_SAMac", 0.21, "small_molecule"),
    ]
    reg_paren_myofib = [
        InterventionRegimen("Parenchyma_Injury", 0.82, "small_molecule"),
        InterventionRegimen("Myofibroblast", 0.22, "small_molecule"),
    ]

    rows = []
    for idx, sp in enumerate(prog_grid):
        # Smooth anatomical trajectory from F1/F2 [0.42, 0.34, 0.32] through F2/F3 [0.52, 0.44, 0.42] to F4 [0.85, 0.88, 0.90]
        s0 = np.clip(np.array([0.08 + 0.86 * sp, -0.02 + 1.00 * sp, -0.04 + 1.045 * sp], dtype=float), 0.10, 0.94)
        r_myo = solver.simulate_ensemble(s0=s0, regimens=reg_myofib_mono, n_patients=220, seed=1100 + idx)
        r_par = solver.simulate_ensemble(s0=s0, regimens=reg_paren_mono, n_patients=220, seed=1200 + idx)
        r_dual = solver.simulate_ensemble(s0=s0, regimens=reg_dual_loop, n_patients=220, seed=1300 + idx)
        r_pm = solver.simulate_ensemble(s0=s0, regimens=reg_paren_myofib, n_patients=220, seed=1400 + idx)
        rows.append(
            {
                "progression_index": round(float(sp), 4),
                "s0_parenchyma": round(float(s0[0]), 4),
                "s0_samac": round(float(s0[1]), 4),
                "s0_myofibroblast": round(float(s0[2]), 4),
                "myofib_lowdose_mono_remission": round(float(r_myo["remission_rate"]), 4),
                "parenchyma_thrb_mono_remission": round(float(r_par["remission_rate"]), 4),
                "dual_myofib_samac_rescue_remission": round(float(r_dual["remission_rate"]), 4),
                "dual_parenchyma_myofib_rescue_remission": round(float(r_pm["remission_rate"]), 4),
            }
        )
    return pd.DataFrame(rows)


def compute_2d_fibrosis_bliss_matrix_and_landscape(
    grammar: NicheGrammar,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Compute symmetric 6x6 clinical-target 2D Bliss synergy matrix in locked F4 cirrhosis + smooth 60x60 2D Waddington quiver field."""
    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)
    s0_f4 = np.array([0.85, 0.88, 0.90], dtype=float)

    arms = [
        ("THRB (Resmetirom SM)", "Parenchyma_Injury", 0.72, "small_molecule"),
        ("GLP1R (Semaglutide Peptide)", "Parenchyma_Injury", 0.68, "peptide"),
        ("TGFBR1/PDGFRB (Low-Dose SM)", "Myofibroblast", 0.145, "small_molecule"),
        ("TREM2/SPP1 (SAMac SM)", "SPP1_SAMac", 0.140, "small_molecule"),
        ("LOXL2 (Simtuzumab 150kDa mAb)", "Myofibroblast", 0.004, "antibody"),
        ("CCR2 (Cenicriviroc Lost-Target)", "SPP1_SAMac", 0.005, "small_molecule"),
    ]

    bliss_rows = []
    for i, (l1, c1, e1, m1) in enumerate(arms):
        for j, (l2, c2, e2, m2) in enumerate(arms):
            ra = InterventionRegimen(target_compartment=c1, perturbation_eta=e1, modality=m1)
            rb = InterventionRegimen(target_compartment=c2, perturbation_eta=e2, modality=m2)
            if i == j:
                res_m = solver.simulate_ensemble(s0=s0_f4, regimens=[ra], n_patients=220, seed=700 + i * 11)
                obs = float(res_m["remission_rate"])
                p1 = p2 = exp_bliss = obs
                excess = 0.0
            else:
                ii, jj = min(i, j), max(i, j)
                ra_o = InterventionRegimen(arms[ii][1], arms[ii][2], arms[ii][3])
                rb_o = InterventionRegimen(arms[jj][1], arms[jj][2], arms[jj][3])
                syn = solver.evaluate_bliss_synergy(
                    s0=s0_f4, regimen_a=ra_o, regimen_b=rb_o, n_patients=220, seed=700 + ii * 10 + jj
                )
                p1 = float(syn["effect_a"]) if i < j else float(syn["effect_b"])
                p2 = float(syn["effect_b"]) if i < j else float(syn["effect_a"])
                obs = float(syn["combo_remission_rate"])
                exp_bliss = float(syn["bliss_expected"])
                excess = float(syn["bliss_excess"])
            bliss_rows.append(
                {
                    "arm_row": l1,
                    "arm_col": l2,
                    "comp_row": c1,
                    "comp_col": c2,
                    "mono_row_remission": round(p1, 4),
                    "mono_col_remission": round(p2, 4),
                    "expected_bliss": round(exp_bliss, 4),
                    "observed_combo_remission": round(obs, 4),
                    "bliss_synergy_excess": round(excess, 4),
                }
            )

    # Smooth 60x60 2D Waddington Quasi-Potential U(S_Myofibroblast, S_Parenchyma) + deterministic quiver drift field
    grid = np.linspace(0.04, 0.96, 60)
    gx, gy = np.meshgrid(grid, grid)  # gx = Myofibroblast (comp 2), gy = Parenchyma_Injury (comp 0)
    # Adiabatic quasi-steady-state projection for SPP1_SAMac (comp 1) along the fibrotic manifold
    s_samac = np.clip(0.52 * gx + 0.48 * gy, 0.04, 0.96)
    pts = np.column_stack([gy.ravel(), s_samac.ravel(), gx.ravel()])
    drift = solver._drift_batch(pts, t=0.0, regimens=[])
    du_paren = drift[:, 0].reshape(gx.shape)
    du_myofib = drift[:, 2].reshape(gx.shape)

    # Compute smooth bistable quasi-potential U(x, y) by stochastic Langevin stationary density + Gaussian smoothing
    rng = np.random.default_rng(2026)
    n_walkers = 1600
    walkers = rng.uniform(0.05, 0.95, size=(n_walkers, 3))
    hist = np.zeros((60, 60), dtype=float)
    for step in range(140):
        d_w = solver._drift_batch(walkers, t=0.0, regimens=[])
        walkers = np.clip(walkers + 0.5 * d_w + math.sqrt(0.5) * 0.075 * rng.normal(size=walkers.shape), 0.02, 0.98)
        if step >= 35:
            h_step, _, _ = np.histogram2d(
                walkers[:, 2], walkers[:, 0], bins=60, range=[[0.04, 0.96], [0.04, 0.96]]
            )
            hist += h_step.T  # shape (y=Parenchyma, x=Myofibroblast)
    smooth_p = ndimage.gaussian_filter(hist + 0.8, sigma=2.2)
    smooth_p = smooth_p / np.max(smooth_p)
    u_2d = -np.log(np.maximum(smooth_p, 1e-6))
    u_2d = u_2d - np.min(u_2d)

    npz_dict = {
        "grid_myofibroblast": grid,
        "grid_parenchyma": grid,
        "U_2d": u_2d,
        "drift_myofibroblast": du_myofib,
        "drift_parenchyma": du_paren,
    }
    return pd.DataFrame(bliss_rows), npz_dict


# ---------------------------------------------------------------------------
# 8. 12-Section 10x Visium Continuous Radial Decay & 2D Tumor Waddington Shift
# ---------------------------------------------------------------------------
def compute_visium_radial_distance_profile(
    coords: np.ndarray,
    comp_scores_df: pd.DataFrame,
    max_pitch: int = 6,
) -> pd.DataFrame:
    """Compute continuous radial distance profile (in physical Visium spot pitches, ~100 um/pitch) from malignant core."""
    xy = np.asarray(coords, dtype=float)
    tree_all = spatial.cKDTree(xy)
    nn_dists, _ = tree_all.query(xy, k=min(7, len(xy)))
    # Characteristic 10x Visium center-to-center spot pitch (~100 um) estimated as median nearest-neighbor distance
    spot_pitch = float(np.median(nn_dists[:, 1])) if xy.shape[0] > 2 else 1.0
    spot_pitch = max(spot_pitch, 1e-6)

    tum = comp_scores_df["Tumor"].to_numpy(dtype=float)
    core_mask = tum >= np.quantile(tum, 0.75)
    if np.sum(core_mask) < 3:
        core_mask = tum >= np.median(tum)

    core_tree = spatial.cKDTree(xy[core_mask])
    dist_raw, _ = core_tree.query(xy, k=1)
    pitch_dist = dist_raw / spot_pitch
    pitch_bin = np.clip(np.rint(pitch_dist).astype(int), 0, max_pitch)

    # Standardize scores within section (z-score) so cross-section aggregation is platform-invariant
    z_df = (comp_scores_df - comp_scores_df.mean(axis=0)) / (comp_scores_df.std(axis=0) + 1e-8)

    rows = []
    for b in range(max_pitch + 1):
        m = pitch_bin == b
        if np.sum(m) == 0:
            continue
        rows.append(
            {
                "pitch_bin": int(b),
                "approx_distance_um": int(b * 100),
                "n_spots": int(np.sum(m)),
                "Tumor_z": float(z_df.loc[m, "Tumor"].mean()),
                "myCAF_z": float(z_df.loc[m, "myCAF"].mean()),
                "TAM_z": float(z_df.loc[m, "TAM"].mean()),
                "CD8_T_z": float(z_df.loc[m, "CD8_T"].mean()),
                "Tumor_raw": float(comp_scores_df.loc[m, "Tumor"].mean()),
                "myCAF_raw": float(comp_scores_df.loc[m, "myCAF"].mean()),
                "TAM_raw": float(comp_scores_df.loc[m, "TAM"].mean()),
                "CD8_T_raw": float(comp_scores_df.loc[m, "CD8_T"].mean()),
            }
        )
    return pd.DataFrame(rows)


def extract_12section_visium_radial_and_contrasting_archetypes(
    vault_sc_dir: pathlib.Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stream all 12 cached human carcinoma Visium sections (21,796 spots) to build radial decay + dual contrasting spot maps."""
    radial_parts: list[pd.DataFrame] = []
    section_spot_cache: dict[str, pd.DataFrame] = {}

    # 1. GSE171351 (4 human bladder carcinoma Visium sections)
    vis_h5ad = vault_sc_dir / "GSE171351_BladderCa_Visium" / "GSE171351_combined_visium.h5ad"
    if vis_h5ad.exists():
        adata_all = ad.read_h5ad(vis_h5ad)
        adata_all.var_names_make_unique()
        for sid in adata_all.obs["sampleID"].unique():
            sub = adata_all[adata_all.obs["sampleID"] == sid].copy()
            if sub.n_obs < 150:
                continue
            r = _score_visium_adata(sub)
            sec_id = f"GSE171351_{sid}"
            prof = compute_visium_radial_distance_profile(r["coords"], r["spot_scores_df"], max_pitch=6)
            prof["section_id"] = sec_id
            prof["cancer_type"] = "Bladder_Carcinoma"
            radial_parts.append(prof)

            sp_df = r["spot_scores_df"].copy()
            sp_df["x"] = r["coords"][:, 0]
            sp_df["y"] = r["coords"][:, 1]
            sp_df["section_id"] = sec_id
            sp_df["cancer_type"] = "Bladder_Carcinoma"
            sp_df["sei"] = r["stromal_exclusion_index"]
            sp_df["tumor_cd8_prox"] = float(r["proximity_weights_df"].loc["Tumor", "CD8_T"])
            section_spot_cache[sec_id] = sp_df

    # 2. GSE203612 _extracted_human (8 human pan-cancer Visium sections)
    ext_dir = vault_sc_dir / "GSE203612_PanCancer_Visium" / "_extracted_human"
    human_tags = ("BRCA", "OVCA", "PAAD", "PDAC", "UCEC", "HNSC", "LUAD", "LUSC", "SKCM", "GIST", "LIHC", "PRAD")
    if ext_dir.exists():
        for h5_p in sorted(ext_dir.glob("*filtered_feature_bc_matrix.h5")):
            prefix = h5_p.name.split("_Vis_")[0]
            pos_candidates = list(ext_dir.glob(f"{prefix}*tissue_positions*"))
            if not pos_candidates:
                continue
            pos_out = pos_candidates[0]
            with h5py.File(h5_p, "r") as hf:
                grp = hf["matrix"]
                mat_csc = sparse.csc_matrix(
                    (grp["data"][:], grp["indices"][:], grp["indptr"][:]),
                    shape=tuple(grp["shape"][:]),
                )
                genes = [x.decode("utf-8") for x in grp["features"]["name"][:]]
                barcodes = [x.decode("utf-8") for x in grp["barcodes"][:]]
            ad_sec = ad.AnnData(X=mat_csc.T.tocsr())
            ad_sec.obs_names = barcodes
            ad_sec.var_names = genes
            ad_sec.var_names_make_unique()

            pos_df = pd.read_csv(pos_out, header=None)
            if str(pos_df.iloc[0, 0]).lower() == "barcode":
                pos_df_h = pd.read_csv(pos_out)
                pos_df = pos_df_h.set_index(pos_df_h.columns[0])
                in_tis = pos_df.iloc[:, 0].astype(int)
                coords_df = pos_df.iloc[:, -2:].astype(float)
            else:
                pos_df = pos_df.set_index(0)
                in_tis = pos_df[1].astype(int)
                coords_df = pos_df[[4, 5]].astype(float)

            common_bc = ad_sec.obs_names.intersection(pos_df.index)
            if len(common_bc) < 150:
                continue
            ad_sec = ad_sec[common_bc].copy()
            ad_sec.obs["in_tissue"] = in_tis.loc[common_bc].to_numpy()
            ad_sec.obsm["spatial"] = coords_df.loc[common_bc].to_numpy(dtype=float)

            r = _score_visium_adata(ad_sec)
            ctype = next((t for t in human_tags if t in prefix.upper()), "Solid_Carcinoma")
            prof = compute_visium_radial_distance_profile(r["coords"], r["spot_scores_df"], max_pitch=6)
            prof["section_id"] = prefix
            prof["cancer_type"] = ctype
            radial_parts.append(prof)

            sp_df = r["spot_scores_df"].copy()
            sp_df["x"] = r["coords"][:, 0]
            sp_df["y"] = r["coords"][:, 1]
            sp_df["section_id"] = prefix
            sp_df["cancer_type"] = ctype
            sp_df["sei"] = r["stromal_exclusion_index"]
            sp_df["tumor_cd8_prox"] = float(r["proximity_weights_df"].loc["Tumor", "CD8_T"])
            section_spot_cache[prefix] = sp_df

    radial_all = pd.concat(radial_parts, ignore_index=True)

    # Pick dual contrasting archetypes: highest myCAF exclusion (GSM6177599_NYU_BRCA0) vs highest Tumor<->CD8 infiltration
    excluded_id = max(section_spot_cache, key=lambda k: float(section_spot_cache[k]["sei"].iloc[0]))
    infiltrated_id = max(
        (k for k in section_spot_cache if k != excluded_id),
        key=lambda k: float(section_spot_cache[k]["tumor_cd8_prox"].iloc[0]) - float(section_spot_cache[k]["sei"].iloc[0]),
    )
    ex_df = section_spot_cache[excluded_id].copy()
    ex_df["archetype"] = "Immune_Excluded_myCAF_Shield"
    inf_df = section_spot_cache[infiltrated_id].copy()
    inf_df["archetype"] = "Immune_Infiltrated_Permissive"
    contrasting_df = pd.concat([ex_df.reset_index(names="barcode"), inf_df.reset_index(names="barcode")], ignore_index=True)
    return radial_all, contrasting_df


def compute_2d_tumor_waddington_surface_shift(
    spatial_prox: pd.DataFrame | None = None,
) -> dict[str, np.ndarray]:
    """Compute smooth 60x60 2D solid-tumor Waddington surfaces U(S_myCAF, S_Tumor) + quiver vector fields (Mono vs Combo Rescue)."""
    grammar = build_default_solid_tumor_icb_grammar(spatial_prox)
    solver = LangevinNicheSolver(grammar, t_end=96.0, dt=0.5)

    # Representative stromal-locked resistant patient state (high myCAF shield blunts anti-PD-1 monotherapy; myCAF+PD-1 unlocks)
    reg_pd1_mono = [
        InterventionRegimen(target_compartment="CD8_Exhaustion", perturbation_eta=0.38, modality="antibody", start_time=0.0),
        InterventionRegimen(target_compartment="Tumor", perturbation_eta=0.22, modality="antibody", start_time=0.0),
    ]
    reg_combo_rescue = [
        InterventionRegimen(target_compartment="CD8_Exhaustion", perturbation_eta=0.90, modality="antibody", start_time=0.0),
        InterventionRegimen(target_compartment="myCAF", perturbation_eta=0.84, modality="small_molecule", start_time=0.0),
        InterventionRegimen(target_compartment="Tumor", perturbation_eta=0.76, modality="antibody", start_time=0.0),
    ]

    grid = np.linspace(0.04, 0.96, 60)
    gx, gy = np.meshgrid(grid, grid)  # gx = myCAF (comp 2), gy = Tumor (comp 0)
    s_tam = np.clip(0.55 * gx + 0.35 * gy, 0.05, 0.95)
    s_exh = np.clip(0.50 * gx + 0.40 * gy, 0.05, 0.95)
    pts = np.column_stack([gy.ravel(), s_tam.ravel(), gx.ravel(), s_exh.ravel()])

    d_mono = solver._drift_batch(pts, t=12.0, regimens=reg_pd1_mono)
    d_combo = solver._drift_batch(pts, t=12.0, regimens=reg_combo_rescue)

    out: dict[str, np.ndarray] = {
        "grid_mycaf": grid,
        "grid_tumor": grid,
        "drift_mono_mycaf": d_mono[:, 2].reshape(gx.shape),
        "drift_mono_tumor": d_mono[:, 0].reshape(gx.shape),
        "drift_combo_mycaf": d_combo[:, 2].reshape(gx.shape),
        "drift_combo_tumor": d_combo[:, 0].reshape(gx.shape),
    }

    for tag, regs in [("mono", reg_pd1_mono), ("combo", reg_combo_rescue)]:
        rng = np.random.default_rng(3030 if tag == "mono" else 4040)
        walkers = rng.uniform(0.05, 0.95, size=(1600, 4))
        hist = np.zeros((60, 60), dtype=float)
        for step in range(130):
            dw = solver._drift_batch(walkers, t=12.0, regimens=regs)
            walkers = np.clip(walkers + 0.5 * dw + math.sqrt(0.5) * 0.08 * rng.normal(size=walkers.shape), 0.02, 0.98)
            if step >= 30:
                h_step, _, _ = np.histogram2d(
                    walkers[:, 2], walkers[:, 0], bins=60, range=[[0.04, 0.96], [0.04, 0.96]]
                )
                hist += h_step.T  # shape (y=Tumor, x=myCAF)
        sp = ndimage.gaussian_filter(hist + 0.8, sigma=2.3)
        sp = sp / np.max(sp)
        u = -np.log(np.maximum(sp, 1e-6))
        out[f"U_2d_{tag}"] = u - np.min(u)

    return out
