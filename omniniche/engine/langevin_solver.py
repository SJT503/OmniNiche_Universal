from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize_scalar

from omniniche.core.schema import DiseasePack, NicheGrammar
from omniniche.disease_packs.loader import compute_barrier_access


@dataclass
class InterventionRegimen:
    target_compartment: str
    perturbation_eta: float
    modality: str = "small_molecule"
    start_time: float = 0.0
    end_time: float = 1e9


class LangevinNicheSolver:
    """Dimensionless stochastic Langevin multicellular dynamical solver with Waddington
    quasi-potential landscape U(S), Kramers escape-barrier quantification, monotherapy window
    cliff t50, staged multi-drug Bliss synergy, LOCO-CV calibration, and Profile-Likelihood
    identifiability auditing.
    """

    def __init__(
        self,
        grammar: NicheGrammar,
        pack: Optional[DiseasePack] = None,
        t_end: float = 72.0,
        dt: float = 0.5,
    ):
        self.grammar = grammar
        self.pack = pack
        self.t_end = float(t_end)
        self.dt = float(dt)
        self.comp_names = [c.name for c in grammar.compartments]
        self.comp_idx = {name: i for i, name in enumerate(self.comp_names)}
        self.n_comp = len(self.comp_names)

        self.alpha = np.array([c.baseline_production for c in grammar.compartments], dtype=float)
        self.gamma = np.array([c.decay_rate for c in grammar.compartments], dtype=float)
        self.sigma = np.array([c.noise_sigma for c in grammar.compartments], dtype=float)

    def _active_remain_factors(self, t: float, regimens: List[InterventionRegimen]) -> np.ndarray:
        remain = np.ones(self.n_comp, dtype=float)
        if not regimens:
            return remain
        norm_prog = float(np.clip(t / max(self.t_end, 1e-6), 0.0, 1.0))
        for reg in regimens:
            if reg.start_time <= t <= reg.end_time and reg.target_compartment in self.comp_idx:
                idx = self.comp_idx[reg.target_compartment]
                b_acc = (
                    compute_barrier_access(self.pack, reg.modality, norm_prog)
                    if self.pack is not None
                    else 0.85
                )
                eff = float(np.clip(reg.perturbation_eta * b_acc, -0.95, 0.98))
                remain[idx] *= (1.0 - eff)
        return np.clip(remain, 0.01, 2.0)

    def _drift_batch(
        self,
        S: np.ndarray,
        t: float,
        regimens: List[InterventionRegimen],
        global_rate_scale: float = 1.0,
        edge_weight_overrides: Optional[Dict[int, float]] = None,
    ) -> np.ndarray:
        """Compute deterministic drift vector dS/dt for batch S of shape (M, N) in [0, 1]^N."""
        m_batch = S.shape[0]
        remain = self._active_remain_factors(t, regimens)
        # Drug inhibition directly scales active signaling output and incoming activation by remaining unblocked fraction
        S_emit = S * remain[None, :]

        act_sum = np.zeros((m_batch, self.n_comp), dtype=float)
        inh_prod = np.ones((m_batch, self.n_comp), dtype=float)

        for e_idx, edge in enumerate(self.grammar.edges):
            if edge.source not in self.comp_idx or edge.target not in self.comp_idx:
                continue
            i = self.comp_idx[edge.source]
            j = self.comp_idx[edge.target]
            v = edge_weight_overrides[e_idx] if (edge_weight_overrides and e_idx in edge_weight_overrides) else edge.weight_v
            k_pow = edge.half_max_k ** edge.hill_n
            s_pow = np.maximum(S_emit[:, i], 0.0) ** edge.hill_n

            if edge.interaction_type == "activation":
                act_sum[:, j] += v * (s_pow / (k_pow + s_pow + 1e-12))
            else:
                inh_prod[:, j] *= k_pow / (k_pow + s_pow + 1e-12)

        prod = self.alpha[None, :] * remain[None, :] + (1.0 - S) * act_sum * inh_prod * remain[None, :]
        decay = self.gamma[None, :] * S
        return global_rate_scale * (prod - decay)

    def simulate_deterministic(
        self,
        s0: np.ndarray,
        regimens: List[InterventionRegimen],
        global_rate_scale: float = 1.0,
        edge_weight_overrides: Optional[Dict[int, float]] = None,
    ) -> Dict[str, np.ndarray]:
        n_steps = int(np.round(self.t_end / self.dt)) + 1
        times = np.linspace(0.0, self.t_end, n_steps)
        states = np.zeros((n_steps, self.n_comp), dtype=float)
        curr = np.clip(np.asarray(s0, dtype=float), 0.0, 1.0)[None, :]
        states[0] = curr[0]

        for step in range(1, n_steps):
            t_prev = times[step - 1]
            # Midpoint RK2 deterministic integration
            k1 = self._drift_batch(curr, t_prev, regimens, global_rate_scale, edge_weight_overrides)
            mid = np.clip(curr + 0.5 * self.dt * k1, 0.0, 1.0)
            k2 = self._drift_batch(mid, t_prev + 0.5 * self.dt, regimens, global_rate_scale, edge_weight_overrides)
            curr = np.clip(curr + self.dt * k2, 0.0, 1.0)
            states[step] = curr[0]

        return {"time": times, "state": states}

    def simulate_ensemble(
        self,
        s0: np.ndarray,
        regimens: List[InterventionRegimen],
        n_patients: int = 200,
        seed: int = 42,
    ) -> Dict[str, object]:
        rng = np.random.default_rng(seed)
        n_steps = int(np.round(self.t_end / self.dt)) + 1
        times = np.linspace(0.0, self.t_end, n_steps)
        traj = np.zeros((n_patients, n_steps, self.n_comp), dtype=float)

        base_s0 = np.clip(np.asarray(s0, dtype=float), 0.0, 1.0)
        curr = np.clip(
            base_s0[None, :] + rng.normal(0.0, 0.025, size=(n_patients, self.n_comp)),
            0.0,
            1.0,
        )
        traj[:, 0, :] = curr
        sqrt_dt = np.sqrt(self.dt)

        for step in range(1, n_steps):
            t_prev = times[step - 1]
            drift = self._drift_batch(curr, t_prev, regimens)
            diffusion = self.sigma[None, :] * np.sqrt(np.maximum(curr * (1.0 - curr), 0.0))
            dW = rng.normal(0.0, sqrt_dt, size=(n_patients, self.n_comp))
            curr = np.clip(curr + self.dt * drift + diffusion * dW, 0.0, 1.0)
            traj[:, step, :] = curr

        rem_idx = self.comp_idx[self.grammar.remission_compartment]
        final_rem_scores = traj[:, -1, rem_idx]
        remission_rate = float(np.mean(final_rem_scores <= self.grammar.remission_threshold))

        return {
            "time": times,
            "trajectories": traj,
            "mean_trajectory": traj.mean(axis=0),
            "ci95_low": np.percentile(traj, 2.5, axis=0),
            "ci95_high": np.percentile(traj, 97.5, axis=0),
            "final_effector_scores": final_rem_scores,
            "remission_rate": remission_rate,
        }

    def compute_monotherapy_window_t50(
        self,
        s0: np.ndarray,
        target_compartment: str,
        perturbation_eta: float,
        modality: str = "small_molecule",
        onset_grid: Optional[np.ndarray] = None,
        n_patients: int = 150,
    ) -> Dict[str, object]:
        if onset_grid is None:
            onset_grid = np.linspace(0.0, self.t_end * 0.75, 13)
        rates = []
        for idx, t_on in enumerate(onset_grid):
            reg = [
                InterventionRegimen(
                    target_compartment=target_compartment,
                    perturbation_eta=perturbation_eta,
                    modality=modality,
                    start_time=float(t_on),
                )
            ]
            res = self.simulate_ensemble(s0=s0, regimens=reg, n_patients=n_patients, seed=100 + idx)
            rates.append(float(res["remission_rate"]))

        rates_arr = np.asarray(rates, dtype=float)
        half_target = 0.5 * (float(rates_arr[0]) + float(rates_arr[-1]))
        # Find crossing of half_target
        above = np.where(rates_arr >= half_target)[0]
        if above.size == 0:
            t50 = 0.0
        elif above[-1] == len(onset_grid) - 1:
            t50 = float(onset_grid[-1])
        else:
            k = int(above[-1])
            r0, r1 = rates_arr[k], rates_arr[k + 1]
            frac = (r0 - half_target) / max(abs(r0 - r1), 1e-9)
            t50 = float(onset_grid[k] + np.clip(frac, 0.0, 1.0) * (onset_grid[k + 1] - onset_grid[k]))

        return {
            "onset_times": onset_grid,
            "remission_rates": rates_arr,
            "t50_window_cliff": t50,
        }

    def evaluate_bliss_synergy(
        self,
        s0: np.ndarray,
        regimen_a: InterventionRegimen,
        regimen_b: InterventionRegimen,
        n_patients: int = 200,
        seed: int = 42,
    ) -> Dict[str, float]:
        res_a = self.simulate_ensemble(s0=s0, regimens=[regimen_a], n_patients=n_patients, seed=seed)
        res_b = self.simulate_ensemble(s0=s0, regimens=[regimen_b], n_patients=n_patients, seed=seed)
        res_ab = self.simulate_ensemble(s0=s0, regimens=[regimen_a, regimen_b], n_patients=n_patients, seed=seed)

        ea = float(res_a["remission_rate"])
        eb = float(res_b["remission_rate"])
        eab = float(res_ab["remission_rate"])
        bliss_expected = ea + eb - ea * eb
        hsa_expected = max(ea, eb)

        return {
            "effect_a": ea,
            "effect_b": eb,
            "combo_remission_rate": eab,
            "bliss_expected": bliss_expected,
            "bliss_excess": eab - bliss_expected,
            "hsa_excess": eab - hsa_expected,
        }

    def compute_quasi_potential_landscape(
        self,
        regimens: List[InterventionRegimen],
        comp_x: str,
        comp_y: str,
        n_grid: int = 30,
        n_trajectories: int = 200,
        seed: int = 7,
    ) -> Dict[str, object]:
        """Compute Waddington 2D quasi-potential landscape U(x, y) = -ln(P_ss(x, y) + eps)
        and quantify pathological attractor lock depth.
        """
        rng = np.random.default_rng(seed)
        n_steps = int(np.round(self.t_end / self.dt)) + 1
        times = np.linspace(0.0, self.t_end, n_steps)
        curr = rng.uniform(0.05, 0.92, size=(n_trajectories, self.n_comp))

        ix = self.comp_idx[comp_x]
        iy = self.comp_idx[comp_y]
        burn_in = n_steps // 2
        xs_samples = []
        ys_samples = []
        sqrt_dt = np.sqrt(self.dt)

        for step in range(1, n_steps):
            t_prev = times[step - 1]
            drift = self._drift_batch(curr, t_prev, regimens)
            diffusion = np.maximum(self.sigma[None, :], 0.03) * np.sqrt(np.maximum(curr * (1.0 - curr), 0.0))
            curr = np.clip(curr + self.dt * drift + diffusion * rng.normal(0.0, sqrt_dt, size=curr.shape), 0.0, 1.0)
            if step >= burn_in:
                xs_samples.append(curr[:, ix].copy())
                ys_samples.append(curr[:, iy].copy())

        x_all = np.concatenate(xs_samples)
        y_all = np.concatenate(ys_samples)

        edges = np.linspace(0.0, 1.0, n_grid + 1)
        hist2d, _, _ = np.histogram2d(x_all, y_all, bins=[edges, edges], density=True)
        smooth_p = gaussian_filter(hist2d, sigma=1.2)
        smooth_p = smooth_p / max(float(np.sum(smooth_p)), 1e-12)

        eps = 1e-5
        U_grid = -np.log(smooth_p + eps)
        U_grid = U_grid - float(np.min(U_grid))

        centers = 0.5 * (edges[:-1] + edges[1:])
        rem_th = self.grammar.remission_threshold
        # Marginal 1D quasi-potential along y (remission/effector axis)
        p_y = np.sum(smooth_p, axis=0)
        u_y = -np.log(p_y + eps)
        u_y = u_y - float(np.min(u_y))

        saddle_mask = (centers >= rem_th * 0.85) & (centers <= min(rem_th * 1.65, 0.60))
        path_mask = centers > min(rem_th * 1.65, 0.60)
        u_ridge = float(np.max(u_y[saddle_mask])) if np.any(saddle_mask) else float(np.median(u_y))
        u_path_min = float(np.min(u_y[path_mask])) if np.any(path_mask) else float(np.max(u_y))
        lock_depth = max(0.0, u_ridge - u_path_min)

        return {
            "grid_centers": centers,
            "U_grid": U_grid,
            "u_y_profile": u_y,
            "pathological_lock_depth": lock_depth,
        }

    def run_loco_cv(self, obs_df: pd.DataFrame, s0: np.ndarray) -> Dict[str, object]:
        """Leave-One-Cohort-Out Cross-Validation calibrating global kinetic rate scale on C-1 cohorts."""
        cohorts = sorted(obs_df["cohort_id"].unique())
        per_cohort_rmse: Dict[str, float] = {}

        for held_out in cohorts:
            train_df = obs_df[obs_df["cohort_id"] != held_out]
            test_df = obs_df[obs_df["cohort_id"] == held_out]

            def _train_loss(rate_scale: float) -> float:
                sim = self.simulate_deterministic(s0=s0, regimens=[], global_rate_scale=rate_scale)
                sq_errs = []
                for _, row in train_df.iterrows():
                    idx = int(np.argmin(np.abs(sim["time"] - float(row["time_hr"]))))
                    obs_vec = np.array([float(row[c]) for c in self.comp_names])
                    sq_errs.append(np.mean((sim["state"][idx] - obs_vec) ** 2))
                return float(np.mean(sq_errs))

            opt = minimize_scalar(_train_loss, bounds=(0.5, 1.8), method="bounded")
            best_scale = float(opt.x)

            sim_test = self.simulate_deterministic(s0=s0, regimens=[], global_rate_scale=best_scale)
            test_sq = []
            for _, row in test_df.iterrows():
                idx = int(np.argmin(np.abs(sim_test["time"] - float(row["time_hr"]))))
                obs_vec = np.array([float(row[c]) for c in self.comp_names])
                test_sq.append(np.mean((sim_test["state"][idx] - obs_vec) ** 2))
            per_cohort_rmse[str(held_out)] = float(np.sqrt(np.mean(test_sq)))

        return {
            "per_cohort_rmse": per_cohort_rmse,
            "mean_out_of_sample_rmse": float(np.mean(list(per_cohort_rmse.values()))),
        }

    def profile_likelihood_identifiability(
        self,
        obs_df: pd.DataFrame,
        s0: np.ndarray,
        obs_sigma: float = 0.04,
    ) -> pd.DataFrame:
        """Scan 1D Profile Likelihood chi^2 across each edge's weight_v to flag structurally
        or practically unidentifiable parameters (95% chi^2_1 threshold = 3.841).
        """
        def _chi2(overrides: Optional[Dict[int, float]] = None) -> float:
            sim = self.simulate_deterministic(s0=s0, regimens=[], edge_weight_overrides=overrides)
            total = 0.0
            for _, row in obs_df.iterrows():
                idx = int(np.argmin(np.abs(sim["time"] - float(row["time_hr"]))))
                obs_vec = np.array([float(row[c]) for c in self.comp_names])
                total += float(np.sum(((sim["state"][idx] - obs_vec) / obs_sigma) ** 2))
            return total

        base_chi2 = _chi2(None)
        multipliers_low = [0.05, 0.25, 0.55]
        multipliers_high = [1.50, 2.10, 3.00]
        chi2_crit = 3.841

        records = []
        for e_idx, edge in enumerate(self.grammar.edges):
            v0 = edge.weight_v
            low_deltas = [max(0.0, _chi2({e_idx: v0 * m}) - base_chi2) for m in multipliers_low]
            high_deltas = [max(0.0, _chi2({e_idx: v0 * m}) - base_chi2) for m in multipliers_high]
            max_low = float(max(low_deltas))
            max_high = float(max(high_deltas))
            identifiable = (max_low >= chi2_crit) and (max_high >= chi2_crit)
            records.append(
                {
                    "edge_index": e_idx,
                    "source": edge.source,
                    "target": edge.target,
                    "nominal_weight_v": v0,
                    "max_delta_chi2_lower": round(max_low, 3),
                    "max_delta_chi2_upper": round(max_high, 3),
                    "identifiable_95ci": identifiable,
                }
            )
        return pd.DataFrame(records)
