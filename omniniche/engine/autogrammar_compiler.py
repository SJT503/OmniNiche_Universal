from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
import numpy as np
import pandas as pd
import yaml

from omniniche.core.schema import (
    CompartmentRule,
    DirectedEdgeRule,
    NicheGrammar,
)


class DeconfoundedAutoGrammarCompiler:
    """Data-driven multicellular grammar compiler (benchmarked against CauTrigger causal deconfounding).
    Extracts directed cross-compartment interactions using time-lagged partial correlation controlling for
    autoregressive self-history and systemic inflammation background Z(t), applies spatial neighborhood
    proximity weights, prunes spurious edges via forward stepwise BIC, and emits NicheGrammar YAML.
    """

    def __init__(self, min_partial_lag_gain: float = 0.20):
        self.min_partial_lag_gain = float(min_partial_lag_gain)

    @staticmethod
    def _partial_correlation(x: np.ndarray, y: np.ndarray, covariates: np.ndarray) -> float:
        """Correlation between residuals of x and y after linear projection out of covariates."""
        design = np.column_stack([np.ones(len(x)), covariates])
        beta_x, _, _, _ = np.linalg.lstsq(design, x, rcond=None)
        beta_y, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        res_x = x - design @ beta_x
        res_y = y - design @ beta_y
        denom = float(np.linalg.norm(res_x) * np.linalg.norm(res_y))
        if denom < 1e-12:
            return 0.0
        return float(np.clip(np.dot(res_x, res_y) / denom, -1.0, 1.0))

    @staticmethod
    def _compute_bic(y: np.ndarray, X_predictors: np.ndarray) -> float:
        m = len(y)
        design = np.column_stack([np.ones(m), X_predictors])
        beta, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        rss = max(float(np.sum((y - design @ beta) ** 2)), 1e-12)
        k = design.shape[1]
        return float(k * np.log(max(m, 2)) + m * np.log(rss / m))

    def compile_from_longitudinal_states(
        self,
        longitudinal_df: pd.DataFrame,
        compartment_cols: List[str],
        pack_id: str,
        organ_context: str,
        remission_compartment: str,
        donor_col: str = "donor_id",
        stage_col: str = "stage_idx",
        confounder_col: Optional[str] = None,
        spatial_enrichment: Optional[pd.DataFrame] = None,
        compartment_markers: Optional[Dict[str, List[str]]] = None,
    ) -> NicheGrammar:
        df = longitudinal_df.sort_values([donor_col, stage_col]).reset_index(drop=True)

        # Construct lagged transition dataset (t_{k-1} -> t_k within each donor, or stage-mean transitions)
        prev_rows = []
        curr_rows = []
        for _, gdf in df.groupby(donor_col, observed=True):
            if len(gdf) < 2:
                continue
            prev_rows.append(gdf.iloc[:-1])
            curr_rows.append(gdf.iloc[1:])

        if not prev_rows:
            raise ValueError("Longitudinal state compilation requires >=2 ordered stages per donor/trajectory.")

        prev_df = pd.concat(prev_rows, ignore_index=True)
        curr_df = pd.concat(curr_rows, ignore_index=True)

        compartments: List[CompartmentRule] = []
        for col in compartment_cols:
            vals = np.clip(df[col].to_numpy(dtype=float), 0.0, 1.0)
            base_prod = float(np.clip(np.percentile(vals, 10) * 0.25, 0.01, 0.15))
            sigma = float(np.clip(np.std(vals) * 0.18, 0.02, 0.12))
            compartments.append(
                CompartmentRule(
                    name=col,
                    baseline_production=base_prod,
                    decay_rate=0.35,
                    noise_sigma=sigma,
                    marker_genes=(compartment_markers or {}).get(col, []),
                )
            )

        selected_edges: List[DirectedEdgeRule] = []
        total_bic = 0.0

        for dst in compartment_cols:
            y_curr = curr_df[dst].to_numpy(dtype=float)
            y_prev = prev_df[dst].to_numpy(dtype=float)

            if confounder_col and confounder_col in curr_df.columns:
                z_curr = curr_df[confounder_col].to_numpy(dtype=float)
            else:
                other_cols = [c for c in compartment_cols if c != dst]
                z_curr = curr_df[other_cols].mean(axis=1).to_numpy(dtype=float)

            base_cov = np.column_stack([y_prev, z_curr])
            best_bic = self._compute_bic(y_curr, base_cov)

            # Score all candidate incoming causal triggers src -> dst
            candidates = []
            for src in compartment_cols:
                if src == dst:
                    continue
                x_prev = prev_df[src].to_numpy(dtype=float)
                pcorr = self._partial_correlation(x_prev, y_curr, base_cov)

                sp_w = 1.0
                if spatial_enrichment is not None and src in spatial_enrichment.index and dst in spatial_enrichment.columns:
                    omega = float(spatial_enrichment.loc[src, dst])
                    sp_w = float(2.0 / (1.0 + np.exp(-omega)))

                eff_gain = pcorr * sp_w
                if abs(eff_gain) >= self.min_partial_lag_gain and abs(pcorr) >= self.min_partial_lag_gain * 0.85:
                    candidates.append((abs(eff_gain), pcorr, sp_w, src, x_prev))

            # Forward stepwise BIC pruning
            candidates.sort(key=lambda item: item[0], reverse=True)
            active_cov = base_cov.copy()
            for _, pcorr, sp_w, src, x_prev in candidates:
                trial_cov = np.column_stack([active_cov, x_prev])
                trial_bic = self._compute_bic(y_curr, trial_cov)
                if trial_bic < best_bic - 1.5:
                    best_bic = trial_bic
                    active_cov = trial_cov
                    inter_type = "activation" if pcorr > 0 else "inhibition"
                    weight_v = float(np.clip(abs(pcorr) * sp_w * 0.95, 0.20, 2.50))
                    hill_n = 3.0 if abs(pcorr) > 0.45 else 2.0
                    selected_edges.append(
                        DirectedEdgeRule(
                            source=src,
                            target=dst,
                            interaction_type=inter_type,
                            weight_v=weight_v,
                            half_max_k=0.40,
                            hill_n=hill_n,
                            partial_corr_gain=round(float(pcorr), 4),
                            spatial_proximity_weight=round(float(sp_w), 4),
                        )
                    )

            total_bic += best_bic

        return NicheGrammar(
            grammar_id=f"{pack_id}_autocompiled_v2",
            source_pack_id=pack_id,
            organ_context=organ_context,
            compartments=compartments,
            edges=selected_edges,
            remission_compartment=remission_compartment,
            remission_threshold=0.30,
            bic_score=round(total_bic, 3),
        )

    @staticmethod
    def export_yaml(grammar: NicheGrammar, path: Path | str) -> None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        body = yaml.safe_dump(grammar.model_dump(), sort_keys=False, allow_unicode=True)
        p.write_text(body, encoding="utf-8")

    @staticmethod
    def load_yaml(path: Path | str) -> NicheGrammar:
        raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
        return NicheGrammar.model_validate(raw)
