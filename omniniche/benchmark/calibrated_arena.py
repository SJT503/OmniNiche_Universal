from __future__ import annotations

from typing import Dict
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from omniniche.core.schema import DiseasePack, NicheGrammar
from omniniche.engine.langevin_solver import InterventionRegimen, LangevinNicheSolver


class CalibratedBenchmarkArena:
    """Head-to-head 5-Baseline Calibrated Benchmark Arena (aligned with Nat. Biotechnol. 2026
    calibrated perturbation benchmarking standards) plus Cross-Organ Zero-Shot Grammar Transfer evaluator.
    """

    def evaluate_five_baselines(self, df: pd.DataFrame) -> pd.DataFrame:
        y_true = df["y_clinical_success"].to_numpy(dtype=int)
        prevalence = float(np.mean(y_true))

        # B1: Static healthy-atlas tau_0 baseline (ignores disease delocalization, barrier & lock-in window)
        p_b1 = np.clip(df["tau_t0"].to_numpy(dtype=float), 0.01, 0.99)

        # B2: Uninformative prevalence / linear mean baseline (null reference)
        p_b2 = np.full(len(df), prevalence, dtype=float)

        # B3: Intra-compartment single-cell GRN perturbation only (CellOracle-style eta without cross-cell lock-in or barrier)
        p_b3 = np.clip(df["grn_eta"].to_numpy(dtype=float), 0.01, 0.99)

        # B4: Static ligand-receptor network centrality only (LIANA+/CellChat-style static interaction score)
        p_b4 = np.clip(df["lr_centrality"].to_numpy(dtype=float), 0.01, 0.99)

        # OmniNiche v2.0 Full: dynamic disease tau * GRN eta * organ barrier access * chrono-window/unlock factor
        chrono_unlock = np.where(
            (df["dosed_before_t50"].to_numpy(dtype=int) == 1) | (df["is_unlock_combo"].to_numpy(dtype=int) == 1),
            1.0,
            0.18,
        )
        raw_omni = (
            df["min_disease_tau"].to_numpy(dtype=float)
            * df["grn_eta"].to_numpy(dtype=float)
            * df["barrier_access"].to_numpy(dtype=float)
            * chrono_unlock
        )
        # Calibrated logistic mapping centered around separation threshold 0.35
        p_omni = 1.0 / (1.0 + np.exp(-8.0 * (raw_omni - 0.35)))

        model_preds = [
            ("B1_Static_Healthy_Tau", p_b1),
            ("B2_Uninformative_Linear_Mean", p_b2),
            ("B3_SingleCell_GRN_Only", p_b3),
            ("B4_Static_LR_Centrality", p_b4),
            ("OmniNiche_v2_Full", p_omni),
        ]

        null_auroc = float(roc_auc_score(y_true, p_b2))
        rows = []
        for name, preds in model_preds:
            auroc = float(roc_auc_score(y_true, preds))
            auprc = float(average_precision_score(y_true, preds))
            brier = float(brier_score_loss(y_true, preds))
            rows.append(
                {
                    "model_name": name,
                    "auroc": round(auroc, 4),
                    "auprc": round(auprc, 4),
                    "brier_score": round(brier, 4),
                    "calibrated_delta_auroc_vs_null": round(auroc - null_auroc, 4),
                }
            )
        return pd.DataFrame(rows)

    @staticmethod
    def evaluate_cross_organ_fibrosis_transfer(
        transferred_ipf_grammar: NicheGrammar,
        fibrosis_pack: DiseasePack,
    ) -> Dict[str, object]:
        """Zero-shot test: apply a grammar compiled on Kidney AKI-to-CKD fibrosis directly to
        Lung IPF compartments and compare late dense-ECM mAb monotherapy vs dual oral small-molecule unlocking.
        """
        solver = LangevinNicheSolver(grammar=transferred_ipf_grammar, pack=fibrosis_pack, t_end=96.0, dt=0.5)
        comp_names = [c.name for c in transferred_ipf_grammar.compartments]
        s0 = np.zeros(len(comp_names), dtype=float)
        s0[0] = 0.82
        if len(s0) > 1:
            s0[1] = 0.08
        if len(s0) > 2:
            s0[2] = 0.04

        epithelium_comp = comp_names[0]
        myofib_comp = transferred_ipf_grammar.remission_compartment

        # Scenario 1: Late anti-ECM mAb monotherapy at t=40h (suffers both dense-ECM mAb exclusion & upstream feed)
        late_mab = [
            InterventionRegimen(
                target_compartment=myofib_comp,
                perturbation_eta=0.75,
                modality="mAb",
                start_time=40.0,
            )
        ]
        res_mab = solver.simulate_ensemble(s0=s0, regimens=late_mab, n_patients=160, seed=2026)

        # Scenario 2: Late dual oral small-molecule antifibrotic unlocking (epithelium + myofibroblast)
        dual_sm = [
            InterventionRegimen(
                target_compartment=epithelium_comp,
                perturbation_eta=0.88,
                modality="small_molecule",
                start_time=40.0,
            ),
            InterventionRegimen(
                target_compartment=myofib_comp,
                perturbation_eta=0.65,
                modality="small_molecule",
                start_time=40.0,
            ),
        ]
        res_dual_sm = solver.simulate_ensemble(s0=s0, regimens=dual_sm, n_patients=160, seed=2026)

        return {
            "transferred_organ": transferred_ipf_grammar.organ_context,
            "late_anti_ecm_mab_monotherapy_remission": float(res_mab["remission_rate"]),
            "dual_small_molecule_remission": float(res_dual_sm["remission_rate"]),
        }
