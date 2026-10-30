from __future__ import annotations

from typing import Dict, List
import numpy as np
import pandas as pd
from scipy import stats


class ChronoTauTracker:
    """Disease-agnostic longitudinal / K-stage single-cell specificity tracker.
    Computes Yanai cell-compartment specificity tau(g, t_k), dynamic specificity drift Delta_tau,
    Sarle's on/off bimodality coefficient BC, and donor-level pseudobulk differential expression.
    """

    @staticmethod
    def compute_yanai_tau(compartment_expression_matrix: pd.DataFrame) -> pd.Series:
        """Compute Yanai specificity index tau in [0, 1] for each gene column across N compartment rows."""
        mat = np.maximum(compartment_expression_matrix.to_numpy(dtype=float), 0.0)
        n_comp = mat.shape[0]
        if n_comp < 2:
            raise ValueError("Yanai tau requires at least 2 cell compartments.")

        max_vals = mat.max(axis=0)
        taus = np.zeros(mat.shape[1], dtype=float)
        nonzero = max_vals > 1e-12
        normed = np.zeros_like(mat)
        normed[:, nonzero] = mat[:, nonzero] / max_vals[nonzero]
        taus[nonzero] = np.sum(1.0 - normed[:, nonzero], axis=0) / float(n_comp - 1)
        return pd.Series(np.clip(taus, 0.0, 1.0), index=compartment_expression_matrix.columns)

    @classmethod
    def compute_dynamic_tau_trajectories(
        cls,
        stage_compartment_profiles: Dict[str, pd.DataFrame],
    ) -> pd.DataFrame:
        """Track tau(g, t_k) across ordered stages t_0, t_1, ..., t_{K-1} and quantify delocalization drift."""
        stage_names = list(stage_compartment_profiles.keys())
        if len(stage_names) < 2:
            raise ValueError("Dynamic tau trajectory tracking requires at least 2 ordered stages (t0, t1, ...).")

        tau_cols: Dict[str, pd.Series] = {}
        peak_comp_t_last: pd.Series | None = None
        for idx, st in enumerate(stage_names):
            prof = stage_compartment_profiles[st]
            tau_cols[f"tau_{st}"] = cls.compute_yanai_tau(prof)
            if idx == len(stage_names) - 1:
                peak_comp_t_last = prof.idxmax(axis=0)

        df = pd.DataFrame(tau_cols)
        t0_col = f"tau_{stage_names[0]}"
        disease_cols = [f"tau_{s}" for s in stage_names[1:]]

        df["tau_t0"] = df[t0_col]
        df["min_disease_tau"] = df[disease_cols].min(axis=1)
        df["mean_disease_tau"] = df[disease_cols].mean(axis=1)
        df["delta_tau_drift"] = df["min_disease_tau"] - df["tau_t0"]
        df["dominant_compartment"] = peak_comp_t_last
        return df

    @staticmethod
    def compute_sarle_bimodality(values: np.ndarray) -> float:
        """Compute Sarle's finite-sample bimodality coefficient BC (BC > 5/9 = 0.555 indicates bimodal/on-off)."""
        arr = np.asarray(values, dtype=float)
        arr = arr[np.isfinite(arr)]
        n = arr.size
        if n < 4 or np.std(arr) < 1e-12:
            return 0.0
        skew = float(stats.skew(arr, bias=False))
        excess_kurt = float(stats.kurtosis(arr, fisher=True, bias=False))
        denom = excess_kurt + (3.0 * (n - 1) ** 2) / float((n - 2) * (n - 3))
        if denom <= 1e-12:
            return 0.0
        return float(np.clip((skew**2 + 1.0) / denom, 0.0, 1.0))

    @staticmethod
    def compute_donor_pseudobulk_de(
        cell_df: pd.DataFrame,
        gene_cols: List[str],
        donor_col: str = "donor_id",
        group_col: str = "group",
        control_label: str = "control",
        disease_label: str = "disease",
    ) -> pd.DataFrame:
        """Aggregate cells to biological donor means before hypothesis testing.
        Hard-blocks cell-level pseudoreplication when <2 independent donors exist per group.
        """
        grouped = cell_df.groupby([donor_col, group_col], observed=True)[gene_cols].mean().reset_index()
        ctrl_donors = grouped[grouped[group_col] == control_label]
        dis_donors = grouped[grouped[group_col] == disease_label]

        if len(ctrl_donors) < 2 or len(dis_donors) < 2:
            raise ValueError(
                f"Anti-pseudoreplication gate triggered: found {len(ctrl_donors)} control donors and "
                f"{len(dis_donors)} disease donors. At least 2 independent biological donors per group are required."
            )

        records = []
        eps = 1e-3
        for g in gene_cols:
            c_vals = np.maximum(ctrl_donors[g].to_numpy(dtype=float), 0.0)
            d_vals = np.maximum(dis_donors[g].to_numpy(dtype=float), 0.0)
            c_mean = float(np.mean(c_vals))
            d_mean = float(np.mean(d_vals))
            log2fc = float(np.log2((d_mean + eps) / (c_mean + eps)))
            _, pval = stats.ttest_ind(d_vals, c_vals, equal_var=False)
            if not np.isfinite(pval):
                pval = 1.0
            records.append(
                {
                    "gene": g,
                    "control_donor_mean": c_mean,
                    "disease_donor_mean": d_mean,
                    "log2fc": log2fc,
                    "pvalue": float(pval),
                }
            )

        res = pd.DataFrame(records).set_index("gene")
        # Benjamini-Hochberg FDR across tested genes
        pvals = res["pvalue"].to_numpy(dtype=float)
        m = len(pvals)
        order = np.argsort(pvals)
        ranked_p = pvals[order]
        fdr_sorted = np.minimum.accumulate((ranked_p * m / np.arange(1, m + 1))[::-1])[::-1]
        fdr = np.empty_like(fdr_sorted)
        fdr[order] = np.clip(fdr_sorted, 0.0, 1.0)
        res["fdr"] = fdr
        return res
