from __future__ import annotations

from typing import List, Tuple
import numpy as np
import pandas as pd

from omniniche.core.schema import DiseasePack
from omniniche.disease_packs.loader import compute_barrier_access


class SixGateTournamentReviewer:
    """Independent 6-Gate Scientific Reviewer + Leave-One-Gate-Out Ablation +
    Bradley-Terry-Luce (BTL) / Elo Pairwise Tournament Ranker adapted to any DiseasePack.
    """

    GATE_COLUMNS: List[str] = [
        "G1_organ_purity",
        "G2_donor_pseudobulk_robustness",
        "G3_dynamic_tau_retention",
        "G4_barrier_modality_alignment",
        "G5_vital_organ_safety",
        "G6_cross_cohort_concordance",
    ]

    def __init__(
        self,
        pack: DiseasePack,
        min_disease_tau: float = 0.60,
        max_negative_drift: float = -0.20,
        min_abs_log2fc: float = 0.50,
        max_fdr: float = 0.10,
        min_barrier_access: float = 0.35,
        min_concordance: float = 0.50,
    ):
        self.pack = pack
        self.min_disease_tau = min_disease_tau
        self.max_negative_drift = max_negative_drift
        self.min_abs_log2fc = min_abs_log2fc
        self.max_fdr = max_fdr
        self.min_barrier_access = min_barrier_access
        self.min_concordance = min_concordance

    def _run_btl_elo_tournament(self, utility_scores: np.ndarray, n_sweeps: int = 40, k_factor: float = 24.0) -> np.ndarray:
        n = len(utility_scores)
        elos = np.full(n, 1500.0, dtype=float)
        if n <= 1:
            return elos

        for _ in range(n_sweeps):
            for i in range(n):
                for j in range(i + 1, n):
                    # True BTL win probability from standardized multi-metric utility difference
                    p_i_wins = 1.0 / (1.0 + np.exp(-3.0 * (utility_scores[i] - utility_scores[j])))
                    expected_i = 1.0 / (1.0 + 10.0 ** ((elos[j] - elos[i]) / 400.0))
                    delta = k_factor * (p_i_wins - expected_i)
                    elos[i] += delta
                    elos[j] -= delta
        return elos

    def audit_and_rank(self, candidates: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        df = candidates.copy()
        purity_blacklist = {
            g.upper()
            for g in (
                self.pack.organ_purity_spec.ambient_contamination_genes
                + self.pack.organ_purity_spec.dissociation_stress_genes
            )
        }
        toxicity_blacklist = {g.upper() for g in self.pack.organ_purity_spec.vital_organ_toxicity_genes}

        genes_upper = [str(g).upper() for g in df.index]

        df["G1_organ_purity"] = [g not in purity_blacklist for g in genes_upper]
        df["G2_donor_pseudobulk_robustness"] = (
            (df["log2fc"].abs() >= self.min_abs_log2fc) & (df["fdr"] <= self.max_fdr)
        )
        df["G3_dynamic_tau_retention"] = (
            (df["min_disease_tau"] >= self.min_disease_tau)
            & (df["delta_tau_drift"] >= self.max_negative_drift)
        )

        modalities = df.get("preferred_modality", pd.Series("small_molecule", index=df.index))
        stage_progs = df.get("target_stage_progress", pd.Series(0.5, index=df.index))
        barrier_scores = [
            compute_barrier_access(self.pack, str(m), float(sp))
            for m, sp in zip(modalities, stage_progs)
        ]
        df["barrier_access_score"] = barrier_scores
        df["G4_barrier_modality_alignment"] = df["barrier_access_score"] >= self.min_barrier_access
        df["G5_vital_organ_safety"] = [g not in toxicity_blacklist for g in genes_upper]

        concordance = df.get("external_cohort_concordance", pd.Series(0.75, index=df.index))
        df["G6_cross_cohort_concordance"] = concordance >= self.min_concordance

        gate_mat = df[self.GATE_COLUMNS].to_numpy(dtype=bool)
        df["gates_passed_count"] = gate_mat.sum(axis=1)
        df["all_six_gates_pass"] = gate_mat.all(axis=1)

        # Leave-one-gate-out ablation summary
        ablation_rows = []
        full_pass_n = int(df["all_six_gates_pass"].sum())
        for idx, gcol in enumerate(self.GATE_COLUMNS):
            other_cols = [c for j, c in enumerate(self.GATE_COLUMNS) if j != idx]
            pass_without_gate = int(df[other_cols].all(axis=1).sum())
            ablation_rows.append(
                {
                    "ablated_gate": gcol,
                    "survivors_with_all_6": full_pass_n,
                    "survivors_without_gate": pass_without_gate,
                    "uniquely_filtered_decoys": pass_without_gate - full_pass_n,
                }
            )
        ablation_df = pd.DataFrame(ablation_rows)

        # Composite continuous utility for BTL/Elo pairwise tournament
        bc_vals = df.get("bimodality_bc", pd.Series(0.5, index=df.index)).to_numpy(dtype=float)
        utility = (
            1.4 * df["mean_disease_tau"].to_numpy(dtype=float)
            + 0.8 * np.clip(df["delta_tau_drift"].to_numpy(dtype=float), -1.0, 0.5)
            + 0.6 * np.tanh(df["log2fc"].abs().to_numpy(dtype=float) / 2.0)
            + 0.7 * df["barrier_access_score"].to_numpy(dtype=float)
            + 0.5 * bc_vals
            + 0.8 * concordance.to_numpy(dtype=float)
            + 1.5 * (df["gates_passed_count"].to_numpy(dtype=float) / 6.0)
        )
        df["composite_utility"] = utility
        df["elo_score"] = self._run_btl_elo_tournament(utility)
        df["tier"] = np.where(
            df["all_six_gates_pass"],
            "TIER_1_APPROVED",
            np.where(df["gates_passed_count"] >= 5, "TIER_2_CONDITIONAL", "REJECTED_DECOY"),
        )
        df = df.sort_values(["all_six_gates_pass", "elo_score"], ascending=[False, False])
        return df, ablation_df
