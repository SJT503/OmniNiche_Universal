from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd


class RidgeGRNPerturber:
    """Intra-compartment regularized Gene Regulatory Network (GRN) perturbation engine.
    Fits a TF/signaling-prior-masked L2 ridge operator W^(c) on single-cell expression data and
    propagates in silico gene knockouts or polypharmacological drug profiles to compute bounded
    compartment suppression efficacy eta in [-1, 1].
    """

    def __init__(self, ridge_lambda: float = 1.0, propagation_steps: int = 3):
        self.ridge_lambda = float(ridge_lambda)
        self.propagation_steps = int(propagation_steps)
        self.genes_: List[str] = []
        self.gene_to_idx_: Dict[str, int] = {}
        self.W_: Optional[np.ndarray] = None
        self.gene_means_: Optional[np.ndarray] = None
        self.gene_stds_: Optional[np.ndarray] = None
        self.pathogenic_unit_vec_: Optional[np.ndarray] = None
        self.program_scale_: float = 1.0

    def fit(
        self,
        expression_df: pd.DataFrame,
        prior_edges: Optional[Iterable[Tuple[str, str]]] = None,
        pathogenic_program_genes: Optional[List[str]] = None,
    ) -> "RidgeGRNPerturber":
        self.genes_ = list(expression_df.columns)
        self.gene_to_idx_ = {g: i for i, g in enumerate(self.genes_)}
        g_count = len(self.genes_)

        X_raw = expression_df.to_numpy(dtype=float)
        self.gene_means_ = X_raw.mean(axis=0)
        self.gene_stds_ = np.maximum(X_raw.std(axis=0), 1e-6)
        X = X_raw - self.gene_means_

        # Build allowed regulator->target mask M (M[i, j] = True means gene i can regulate gene j)
        if prior_edges is not None:
            mask = np.zeros((g_count, g_count), dtype=bool)
            for src, dst in prior_edges:
                if src in self.gene_to_idx_ and dst in self.gene_to_idx_:
                    i, j = self.gene_to_idx_[src], self.gene_to_idx_[dst]
                    if i != j:
                        mask[i, j] = True
        else:
            mask = ~np.eye(g_count, dtype=bool)

        W = np.zeros((g_count, g_count), dtype=float)
        for j in range(g_count):
            reg_idx = np.where(mask[:, j])[0]
            if reg_idx.size == 0:
                continue
            X_reg = X[:, reg_idx]
            y_target = X[:, j]
            gram = X_reg.T @ X_reg + self.ridge_lambda * np.eye(reg_idx.size)
            coef = np.linalg.solve(gram, X_reg.T @ y_target)
            W[reg_idx, j] = coef

        # Spectral radius stabilization so multi-step propagation stays well-conditioned
        eig_vals = np.linalg.eigvals(W)
        rho = float(np.max(np.abs(eig_vals))) if eig_vals.size else 0.0
        if rho > 0.95:
            W = W * (0.95 / rho)
        self.W_ = W

        prog_genes = pathogenic_program_genes or self.genes_
        p_vec = np.zeros(g_count, dtype=float)
        for pg in prog_genes:
            if pg in self.gene_to_idx_:
                p_vec[self.gene_to_idx_[pg]] = 1.0
        norm = float(np.linalg.norm(p_vec))
        if norm < 1e-12:
            raise ValueError("None of pathogenic_program_genes matched expression_df columns.")
        self.pathogenic_unit_vec_ = p_vec / norm

        proj_scores = X @ self.pathogenic_unit_vec_
        self.program_scale_ = max(float(np.std(proj_scores)) * 1.5, 1e-3)
        return self

    def _propagate_delta(self, delta0: np.ndarray) -> float:
        if self.W_ is None or self.pathogenic_unit_vec_ is None:
            raise RuntimeError("RidgeGRNPerturber must be fitted before predicting perturbation efficacy.")
        total_delta = delta0.copy()
        curr = delta0.copy()
        for _ in range(self.propagation_steps):
            curr = curr @ self.W_
            total_delta = total_delta + curr
        suppression_raw = -float(np.dot(total_delta, self.pathogenic_unit_vec_))
        return float(np.tanh(suppression_raw / self.program_scale_))

    def predict_knockout_efficacy(self, target_gene: str) -> float:
        """Predict compartment pathogenic suppression efficacy eta in [-1, 1] upon target_gene knockout."""
        if self.gene_means_ is None or self.gene_stds_ is None:
            raise RuntimeError("RidgeGRNPerturber not fitted.")
        if target_gene not in self.gene_to_idx_:
            return 0.0
        idx = self.gene_to_idx_[target_gene]
        delta0 = np.zeros(len(self.genes_), dtype=float)
        delta0[idx] = -max(abs(float(self.gene_means_[idx])), float(self.gene_stds_[idx]))
        return self._propagate_delta(delta0)

    def predict_drug_efficacy(self, target_pic50_map: Dict[str, float], ref_pconc: float = 6.0) -> float:
        """Predict polypharmacological drug suppression efficacy eta weighted by target pIC50 occupancy."""
        if self.gene_means_ is None or self.gene_stds_ is None:
            raise RuntimeError("RidgeGRNPerturber not fitted.")
        delta0 = np.zeros(len(self.genes_), dtype=float)
        for gene, pic50 in target_pic50_map.items():
            if gene in self.gene_to_idx_:
                idx = self.gene_to_idx_[gene]
                occupancy = 1.0 / (1.0 + 10.0 ** (ref_pconc - float(pic50)))
                delta0[idx] = -occupancy * max(abs(float(self.gene_means_[idx])), float(self.gene_stds_[idx]))
        return self._propagate_delta(delta0)
