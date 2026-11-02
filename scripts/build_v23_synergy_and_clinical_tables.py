"""Build v2.3 Continuous 2D Dose-Synergy Surface, Anti-PD-1 Monotherapy Survival (n=74) & Continuous ROC (n=107) Tables."""
from __future__ import annotations

import pathlib
import sys
import tracemalloc

import pandas as pd

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.stage2.fibrosis_mash_deconvolution import build_conserved_fibrosis_grammar_from_liver
from omniniche.stage2.v2_statistical_fortification import (
    compute_2d_continuous_dose_synergy_surface,
    compute_antipd1_monotherapy_survival_and_roc,
)
from scripts.fetch_v22_genetics_and_monotherapy_cohorts import download_monotherapy_icb_geo_cohorts


def main() -> None:
    tracemalloc.start()
    vault_sc = ROOT / "data_vault" / "raw_single_cell"
    out_tables = ROOT / "outputs" / "stage2_tables"

    print(">>> [1/3] Verifying & SHA-256 locking monotherapy ICB clinical survival cache...")
    download_monotherapy_icb_geo_cohorts(vault_sc)

    print(">>> [2/3] Computing 13x13 Continuous 2D Locked-F4 Dose-Response Bliss Synergy Surface...")
    s4a = pd.read_csv(out_tables / "Table_S4a_PanOrgan_Fibrosis_All_Human_Donors.csv")
    liver_donors = s4a[s4a["organ"].astype(str).str.startswith("Liver")].copy()
    grammar_liv = build_conserved_fibrosis_grammar_from_liver(liver_donors)
    dose_syn_df = compute_2d_continuous_dose_synergy_surface(grammar_liv, n_grid=13)
    dose_syn_df.to_csv(out_tables / "Table_S4h_Continuous_2D_Dose_Synergy_Surface.csv", index=False, encoding="utf-8")
    print(
        f"    Table_S4h exported ({len(dose_syn_df)} grid points; max Bliss excess = {dose_syn_df['bliss_synergy_excess'].max():.4f}, "
        f"nonzero synergy cells = {int((dose_syn_df['bliss_synergy_excess'] > 0.02).sum())}/{len(dose_syn_df)})"
    )

    print(">>> [3/3] Computing n=74 Pre-Treatment Anti-PD-1 Monotherapy Survival & n=107 Continuous ROC Curves...")
    surv_74_df, surv_sum_df, roc_107_df = compute_antipd1_monotherapy_survival_and_roc(vault_sc, out_tables)
    surv_74_df.to_csv(out_tables / "Table_S5l_AntiPD1_Monotherapy_Survival_74Patients.csv", index=False, encoding="utf-8")
    surv_sum_df.to_csv(out_tables / "Table_S5m_AntiPD1_Survival_Summary.csv", index=False, encoding="utf-8")
    roc_107_df.to_csv(out_tables / "Table_S5n_PreTreatment_107Patients_ROC_Curves.csv", index=False, encoding="utf-8")
    print(surv_sum_df.to_string(index=False))

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(f"\n=== V2.3 TABLES EXPORTED DETERMINISTICALLY (Peak RAM = {peak_bytes / 1e6:.1f} MB << 7,000 MB) ===")


if __name__ == "__main__":
    main()
