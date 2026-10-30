"""Generate all 5 Main Flagship Figures + 3 Extended Data Figures for OmniNiche v2.0 Stage II."""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.viz.stage2_flagship_figures import (
    build_extended_data_figures,
    build_fig2_macro_clinical_atlas,
    build_fig3_pan_disease_drift_law,
    build_fig4_fibrosis_mash_deconvolution,
    build_fig5_solid_tumor_spatial_icb,
    export_fig1_illustration,
)


def main() -> None:
    tables_dir = ROOT / "outputs" / "stage2_tables"
    out_dir = ROOT / "outputs" / "stage2_figures"
    artifact_dir = pathlib.Path(r"C:\Users\sheng\.gemini\antigravity\brain\c14f9ca0-84c8-456f-a314-0bd06b6f8bd7")

    fig1_src = artifact_dir / "fig1_omniniche_architecture_and_framework_v2_1793349563115.jpg"
    print(">>> [1/6] Exporting Fig. 1 BioRender Architecture Illustration (600 DPI PNG + PDF + SVG)...")
    export_fig1_illustration(fig1_src, out_dir, artifact_dir)

    print(">>> [2/6] Rendering Fig. 2 Macro-Atlas (41,396 Non-CNS Phase II-IV Clinical Trials)...")
    build_fig2_macro_clinical_atlas(tables_dir, out_dir, artifact_dir)

    print(">>> [3/6] Rendering Fig. 3 Pan-Organ Dynamic Drift x Biophysical Barrier Law...")
    build_fig3_pan_disease_drift_law(tables_dir, out_dir, artifact_dir)

    print(">>> [4/6] Rendering Fig. 4 Pure Human MASH & Cross-Organ Zero-Shot Fibrosis Benchmark...")
    build_fig4_fibrosis_mash_deconvolution(tables_dir, out_dir, artifact_dir)

    print(">>> [5/6] Rendering Fig. 5 Solid Tumor 10x Visium Spatial Barrier & Blind ICB Responder Arena...")
    build_fig5_solid_tumor_spatial_icb(tables_dir, out_dir, artifact_dir)

    print(">>> [6/6] Rendering Extended Data Figs. 1-3...")
    build_extended_data_figures(tables_dir, out_dir, artifact_dir)

    print("\n=== ALL 8 STAGE II FLAGSHIP & EXTENDED DATA FIGURES GENERATED SUCCESSFULLY ===")


if __name__ == "__main__":
    main()
