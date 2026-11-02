"""Generate all 5 Main Flagship Figures + 4 Extended Data Figures for OmniNiche Stage II v2.3 (183 mm CNS Print Standard)."""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.viz.stage2_flagship_figures import (
    build_extended_data_figures,
    build_fig1_vector_architecture,
    build_fig2_macro_clinical_atlas,
    build_fig3_pan_disease_drift_law,
    build_fig4_fibrosis_mash_deconvolution,
    build_fig5_solid_tumor_spatial_icb,
)


def main() -> None:
    tables_dir = ROOT / "outputs" / "stage2_tables"
    out_dir = ROOT / "outputs" / "stage2_figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Purge superseded figure files so outputs/stage2_figures holds strictly the clean v2.3 9-figure suite
    for old_f in list(out_dir.glob("*")):
        if old_f.is_file():
            old_f.unlink()

    print(">>> [1/6] Rendering Fig. 1 Pure-Vector Mathematical & Cohort-Scale Architecture (4 Panels a-d, 183 mm, 600 DPI PNG + PDF + Editable SVG)...")
    build_fig1_vector_architecture(tables_dir, out_dir, artifact_dir=None)

    print(">>> [2/6] Rendering Fig. 2 Macro-Atlas (41,396 Systemic & Parenchymal Organ Trials, 6 Panels a-f)...")
    build_fig2_macro_clinical_atlas(tables_dir, out_dir, artifact_dir=None)

    print(">>> [3/6] Rendering Fig. 3 Pan-Organ Dynamic Drift & DNAcc Law (8 Panels a-h)...")
    build_fig3_pan_disease_drift_law(tables_dir, out_dir, artifact_dir=None)

    print(">>> [4/6] Rendering Fig. 4 Pan-Organ Human Fibrosis Zero-Shot Benchmark, Hysteresis & 2D Dose-Synergy Surface (8 Panels a-h)...")
    build_fig4_fibrosis_mash_deconvolution(tables_dir, out_dir, artifact_dir=None)

    print(">>> [5/6] Rendering Fig. 5 Solid-Tumor 10x Visium Spatial Barrier, Pre-ICB Responder & Monotherapy OS (8 Panels a-h)...")
    build_fig5_solid_tumor_spatial_icb(tables_dir, out_dir, artifact_dir=None)

    print(">>> [6/7] Rendering Extended Data Figs. 1-4 (Native 183 mm Print Width)...")
    build_extended_data_figures(tables_dir, out_dir, artifact_dir=None)

    print(">>> [7/7] Rendering Extended Data Figs. 5-6 & Syncing Main/Supp Figure Directories...")
    from omniniche.viz.stage2_supp_figures import organize_and_render_r3_figure_ecosystem
    organize_and_render_r3_figure_ecosystem(ROOT)

    print("\n=== ALL 11 STAGE II v3.1 183-MM FLAGSHIP & EXTENDED DATA FIGURES GENERATED SUCCESSFULLY ===")


if __name__ == "__main__":
    main()
