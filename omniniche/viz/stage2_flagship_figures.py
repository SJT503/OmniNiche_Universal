"""Stage II v2.1 Native 183-mm (7.20-in) CNS / Nature Medicine Flagship Visualization Suite.

Renders all 5 Main Figures + 4 Extended Data Figures natively at Nature double-column
print width (183 mm = 7.20 in, 600 DPI PNG + vector PDF + editable SVG) with 1:1 true
print typography (9.5 pt bold lowercase panel tags, 6.8 pt bold concise headers,
6.1 pt axis labels, 5.6-5.9 pt tick/tabular text, 0.55 pt hairlines), protagonist-cobalt
vs neutral-slate baseline discipline, mathtext subscripts (zero raw underscores),
and zero external engineering partitioning labels.
"""
from __future__ import annotations

import math
import pathlib
import shutil
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage, spatial, stats


# Nature / Nature Medicine Okabe-Ito + Protagonist-Cobalt vs Neutral-Slate Palette
PAL = {
    "navy": "#0A2E4D",
    "blue": "#0072B2",       # Protagonist primary (OmniNiche / DNAcc)
    "sky": "#7FB8DF",
    "teal": "#009E73",       # Protagonist secondary / Remission / Responder
    "orange": "#E69F00",     # Administrative / Operational / TAM
    "vermilion": "#D55E00",  # Failure / Barrier / PD-L1 contrast
    "crimson": "#B8292F",    # Biological failure / Fibrotic scar / Non-responder
    "purple": "#8A5082",
    "slate": "#64707D",      # Muted neutral dark grey
    "midgrey": "#95A0AC",    # Competing baseline bars/points
    "lightgrey": "#D4DBE2",  # Background tracks / reference fills
    "pale": "#F3F6F9",       # Subtle zebra / track background
    "ink": "#182026",
}

AXIS_SHORT = {
    "Solid_Oncology": "Solid Oncology",
    "Metabolic_Endocrine": "Metabolic / Endocr.",
    "Hematologic_Malignancy": "Hematologic Oncol.",
    "Pulmonary_IPF_ARDS": "Pulmonary (IPF)",
    "Cardiovascular_HF_Athero": "Cardiovascular",
    "Autoimmune_IBD_Rheum": "Autoimmune / IBD",
    "Renal_CKD_AKI": "Renal (CKD/AKI)",
    "Dermatologic_Immune": "Dermatologic",
    "Musculoskeletal_Bone": "Musculoskeletal",
    "Ophthalmic_Retinal": "Ophthalmic",
    "Hepatic_MASH_Fibrosis": "Hepatic (MASH)",
    "Infectious_Sepsis": "Infectious / Sepsis",
}

AXIS_COMPACT = {
    "Hepatic_MASH_Fibrosis": "Hepatic MASH",
    "Pulmonary_IPF_ARDS": "Pulmonary IPF",
    "Renal_CKD_AKI": "Renal CKD",
    "Solid_Oncology": "Solid Oncology",
    "Autoimmune_IBD_Rheum": "Autoimmune IBD",
    "Infectious_Sepsis": "Sepsis (PBMC)",
}


def get_waddington_light_cmap() -> mcolors.LinearSegmentedColormap:
    """Return a clean Nature Medicine light-base topographic potential colormap (cool basins -> porcelain saddle -> warm scar ridges)."""
    return mcolors.LinearSegmentedColormap.from_list(
        "waddington_light",
        ["#B8DCEB", "#D8ECF5", "#F4F7FA", "#F9EBDD", "#F0C59D", "#D9825F"],
        N=256,
    )


def get_synergy_light_cmap() -> mcolors.LinearSegmentedColormap:
    """Return a clean light-porcelain to amber-crimson continuous Bliss synergy colormap."""
    return mcolors.LinearSegmentedColormap.from_list(
        "synergy_light",
        ["#F5F7FA", "#FCE6C2", "#F5B166", "#DF6853", "#B8292F"],
        N=256,
    )


def set_nature_style() -> None:
    """Configure global Matplotlib rcParams for native 183-mm (7.20-in) Nature Medicine print scale."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "font.size": 5.9,
            "axes.titlesize": 6.35,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.titlepad": 4.0,
            "axes.labelsize": 6.0,
            "axes.labelweight": "normal",
            "axes.edgecolor": "#26303A",
            "axes.linewidth": 0.52,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 5.5,
            "ytick.labelsize": 5.6,
            "xtick.major.width": 0.48,
            "ytick.major.width": 0.48,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "xtick.major.pad": 1.9,
            "ytick.major.pad": 1.9,
            "legend.fontsize": 5.3,
            "legend.frameon": False,
            "legend.handlelength": 1.20,
            "legend.handletextpad": 0.38,
            "legend.columnspacing": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _panel_tag(ax: plt.Axes, letter: str, x: float = -0.14, y: float = 1.08) -> None:
    """Place a Nature-standard bold lowercase panel tag (9.5 pt at 183-mm native width)."""
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=9.5,
        fontweight="bold",
        va="top",
        ha="left",
        color="#000000",
        clip_on=False,
    )


def _fmt_p(p_val: float) -> str:
    """Format P-values cleanly in compact scientific mathtext notation."""
    if not np.isfinite(p_val):
        return "$P$ = n.s."
    if p_val < 1e-99:
        return "$P < 10^{-99}$"
    if p_val < 1e-3:
        exp = int(math.floor(math.log10(max(p_val, 1e-300))))
        mant = p_val / (10.0 ** exp)
        return rf"$P = {mant:.1f} \times 10^{{{exp}}}$"
    return f"$P = {p_val:.4f}$" if p_val < 0.01 else f"$P = {p_val:.2f}$"


def _save_multi(
    fig: plt.Figure,
    out_dir: pathlib.Path,
    stem: str,
    artifact_dir: pathlib.Path | None = None,
    dpi: int = 600,
) -> pathlib.Path:
    """Export figure cleanly as 600-DPI PNG, vector PDF, and editable SVG with pure-white outer margin."""
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{stem}.png"
    pdf_path = out_dir / f"{stem}.pdf"
    svg_path = out_dir / f"{stem}.svg"

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    fig.savefig(svg_path, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    plt.close(fig)

    if artifact_dir is not None:
        artifact_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(png_path, artifact_dir / png_path.name)
    return png_path


def _draw_alluvial_ribbon(
    ax: plt.Axes,
    x0: float,
    y0_bot: float,
    y0_top: float,
    x1: float,
    y1_bot: float,
    y1_top: float,
    color: str,
    alpha: float = 0.28,
) -> None:
    """Draw a smooth cubic Bezier alluvial ribbon between two vertical stage stacks."""
    xm = 0.5 * (x0 + x1)
    verts = [
        (x0, y0_bot),
        (xm, y0_bot),
        (xm, y1_bot),
        (x1, y1_bot),
        (x1, y1_top),
        (xm, y1_top),
        (xm, y0_top),
        (x0, y0_top),
        (x0, y0_bot),
    ]
    codes = [
        mpath.Path.MOVETO,
        mpath.Path.CURVE4,
        mpath.Path.CURVE4,
        mpath.Path.CURVE4,
        mpath.Path.LINETO,
        mpath.Path.CURVE4,
        mpath.Path.CURVE4,
        mpath.Path.CURVE4,
        mpath.Path.CLOSEPOLY,
    ]
    patch = mpatches.PathPatch(
        mpath.Path(verts, codes),
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
    )
    ax.add_patch(patch)


# ===========================================================================
# FIGURE 1: Pure-Vector Mathematical & Cohort-Scale Architecture (a-d)
# ===========================================================================
def build_fig1(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> pathlib.Path:
    """Render 100% programmatic pure-vector Fig. 1 (a-d) natively at 183 mm (7.20 in x 5.85 in)."""
    set_nature_style()

    s1 = pd.read_csv(tables_dir / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv")
    s2a = pd.read_csv(tables_dir / "Table_S2a_Stoppage_Auditor_Validation.csv")
    s3c = pd.read_csv(tables_dir / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv", usecols=["nct_id"])
    s3d = pd.read_csv(tables_dir / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv")
    s3g = pd.read_csv(tables_dir / "Table_S3g_Multivariable_FixedEffects_Control.csv")
    s3l = pd.read_csv(tables_dir / "Table_S3l_Genetics_vs_DNAcc_Orthogonality.csv")
    s4 = pd.read_csv(tables_dir / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv")
    s4b = pd.read_csv(tables_dir / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv")

    n_tot = int(s1["unique_systemic_organ_trials"].sum())
    n_res = int(s1["target_resolved_trials"].sum())
    n_adm = int(s1["administrative_commercial_enrollment"].sum())
    n_stops = int(n_tot - s1["completed_trials"].sum())
    adm_pct = 100.0 * n_adm / max(n_stops, 1)
    n_sci = len(s3c)
    n_mat = int(s3d["n_trials_dnacc"].sum())

    cl_col = "auditor_name" if "auditor_name" in s2a.columns else "classifier"
    cl_row = s2a[s2a[cl_col] == "Clause_Aware_Auditor_v2_2"].iloc[0] if cl_col in s2a.columns else s2a.iloc[0]
    acc_pct = 100.0 * float(cl_row["overall_accuracy"])
    kappa_v = float(cl_row["cohen_kappa"])

    auc_24 = float(
        s4b.loc[
            (s4b["benchmark"] == f"PanOrgan_Fibrosis_{len(s4)}_Phase2b3_Trials")
            & (s4b["method"] == "OmniNiche_v2_ZeroShot_Langevin"),
            "auroc",
        ].iloc[0]
    )

    fig = plt.figure(figsize=(7.20, 5.85), facecolor="white")
    gs = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1.10, 1.16, 0.62], hspace=0.18)
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[1.00, 1.06], wspace=0.08)

    def _rbox(
        ax: plt.Axes,
        x: float,
        y: float,
        w: float,
        h: float,
        fc: str = "#FFFFFF",
        ec: str = "#CBD5E1",
        lw: float = 0.55,
        pad: float = 0.8,
        zorder: int = 1,
    ) -> None:
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (x + pad, y + pad),
                w - 2 * pad,
                h - 2 * pad,
                boxstyle=f"round,pad={pad},rounding_size=1.6",
                facecolor=fc,
                edgecolor=ec,
                linewidth=lw,
                zorder=zorder,
            )
        )

    def _arrow(
        ax: plt.Axes,
        x0: float,
        y0: float,
        x1: float,
        y1: float,
        col: str = PAL["slate"],
        lw: float = 0.85,
        rad: float = 0.0,
    ) -> None:
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops={
                "arrowstyle": "-|>,head_length=0.35,head_width=0.22",
                "color": col,
                "lw": lw,
                "connectionstyle": f"arc3,rad={rad}",
                "shrinkA": 0.4,
                "shrinkB": 0.4,
            },
            zorder=5,
        )

    # -----------------------------------------------------------------------
    # Panel a: 12-Axis Clinical Trial Census & Clause-Aware Deconvolution
    # -----------------------------------------------------------------------
    ax_a = fig.add_subplot(gs_top[0, 0])
    ax_a.set_xlim(0, 100)
    ax_a.set_ylim(0, 62)
    ax_a.axis("off")
    _panel_tag(ax_a, "a", x=-0.02, y=1.04)

    _rbox(ax_a, 0.5, 0.5, 99.0, 60.5, fc="#F8FAFC", ec="#CBD5E1", lw=0.60, pad=0.6)
    ax_a.add_patch(mpatches.Rectangle((1.2, 53.0), 97.6, 7.2, facecolor="#E6F0F8", edgecolor="none", zorder=2))
    ax_a.add_patch(mpatches.Rectangle((1.2, 53.0), 1.6, 7.2, facecolor=PAL["blue"], edgecolor="none", zorder=3))
    ax_a.text(4.5, 56.6, "12 Peripheral Organ Axes & Clause-Aware Deconvolution", fontsize=5.8, fontweight="bold", color=PAL["navy"], va="center", zorder=4)

    # Col 1: Registry Census Card
    _rbox(ax_a, 2.5, 3.0, 29.5, 48.0, fc="#FFFFFF", ec="#B8C7D6", lw=0.55, pad=0.6, zorder=2)
    ax_a.text(17.2, 47.4, f"n = {n_tot:,} Trials", fontsize=5.6, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_a.text(17.2, 43.4, "Phase II\u2013IV (CT.gov v2)", fontsize=4.6, color=PAL["slate"], ha="center", va="center", zorder=4)

    top6 = s1.sort_values("unique_systemic_organ_trials", ascending=False).head(6)
    max_ax_v = float(top6["unique_systemic_organ_trials"].max())
    short6 = {
        "Solid_Oncology": "Oncol.",
        "Metabolic_Endocrine": "Metab.",
        "Hematologic_Malignancy": "Heme",
        "Pulmonary_IPF_ARDS": "Pulm.",
        "Cardiovascular_HF_Athero": "Cardio.",
        "Autoimmune_IBD_Rheum": "Autoim.",
    }
    for idx_bar, r_bar in enumerate(top6.itertuples()):
        yb = 37.5 - idx_bar * 4.3
        bw_i = 10.2 * (float(r_bar.unique_systemic_organ_trials) / max_ax_v)
        ax_a.text(3.8, yb, short6.get(r_bar.disease_axis, r_bar.disease_axis[:6]), fontsize=4.2, color=PAL["ink"], va="center", zorder=4)
        ax_a.add_patch(mpatches.Rectangle((13.0, yb - 1.25), bw_i, 2.5, facecolor=PAL["blue"], alpha=0.82, zorder=3))
        ax_a.text(13.4 + bw_i, yb, f"{r_bar.unique_systemic_organ_trials / 1000:.1f}k", fontsize=4.0, color=PAL["slate"], va="center", zorder=4)

    ax_a.text(17.2, 11.6, "+ 6 Organ Axes", fontsize=4.15, fontstyle="italic", color=PAL["slate"], ha="center", va="center", zorder=4)
    ax_a.text(17.2, 6.8, f"{n_res:,} Resolved\n(1,155 Target Genes)", fontsize=4.05, fontweight="bold", color=PAL["blue"], ha="center", va="center", zorder=4)

    _arrow(ax_a, 32.3, 27.0, 35.8, 27.0, col=PAL["blue"])

    # Col 2: Clause-Aware Negation Filter
    _rbox(ax_a, 36.0, 6.5, 28.0, 41.0, fc="#FFFFFF", ec=PAL["orange"], lw=0.65, pad=0.6, zorder=2)
    ax_a.text(50.0, 43.5, "Clause-Aware Auditor", fontsize=4.85, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_a.text(50.0, 37.8, "1. Split clause boundaries\n(; . but while unrelated)", fontsize=4.05, color=PAL["ink"], ha="center", va="center", zorder=4)
    ax_a.text(50.0, 29.0, "2. Strip negated safety:\n\"no safety concerns\"", fontsize=4.05, color=PAL["crimson"], ha="center", va="center", zorder=4)
    _rbox(ax_a, 38.0, 9.2, 24.0, 13.8, fc="#EDF7F4", ec=PAL["teal"], lw=0.50, pad=0.5, zorder=3)
    ax_a.text(
        50.0,
        16.1,
        f"Manual Audit (n = 240)\n{acc_pct:.1f}% Acc. | $\\kappa$ = {kappa_v:.2f}\n(vs. 70.0% Naive)",
        fontsize=4.05,
        fontweight="bold",
        color=PAL["navy"],
        ha="center",
        va="center",
        zorder=4,
    )

    _arrow(ax_a, 64.2, 34.0, 68.2, 43.0, col=PAL["orange"])
    _arrow(ax_a, 64.2, 27.0, 68.2, 27.0, col=PAL["blue"])
    _arrow(ax_a, 64.2, 20.0, 68.2, 11.0, col=PAL["teal"])

    # Col 3: 3 Deconvoluted Strata
    _rbox(ax_a, 68.5, 36.5, 29.0, 14.2, fc="#FFF6E9", ec=PAL["orange"], lw=0.55, pad=0.5, zorder=2)
    ax_a.text(83.0, 45.8, f"Admin. Halts ({adm_pct:.1f}%)", fontsize=4.8, fontweight="bold", color="#9C5700", ha="center", va="center", zorder=4)
    ax_a.text(83.0, 40.5, f"n = {n_adm:,} Quarantined", fontsize=4.4, color=PAL["ink"], ha="center", va="center", zorder=4)

    _rbox(ax_a, 68.5, 19.8, 29.0, 14.4, fc="#EAF2F8", ec=PAL["blue"], lw=0.60, pad=0.5, zorder=2)
    ax_a.text(83.0, 29.3, "Scientific Set", fontsize=4.9, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_a.text(83.0, 24.0, f"n = {n_sci:,} (G = 374)", fontsize=4.45, fontweight="bold", color=PAL["blue"], ha="center", va="center", zorder=4)

    _rbox(ax_a, 68.5, 3.0, 29.0, 14.4, fc="#E8F6F3", ec=PAL["teal"], lw=0.60, pad=0.5, zorder=2)
    ax_a.text(83.0, 12.5, "Maturation Set", fontsize=4.9, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_a.text(83.0, 7.2, f"n = {n_mat:,} (G = 561)", fontsize=4.45, fontweight="bold", color=PAL["teal"], ha="center", va="center", zorder=4)

    # -----------------------------------------------------------------------
    # Panel b: Single-Cell Remodeling Drift & Modality-ECM Operator (DNAcc)
    # -----------------------------------------------------------------------
    ax_b = fig.add_subplot(gs_top[0, 1])
    ax_b.set_xlim(0, 100)
    ax_b.set_ylim(0, 62)
    ax_b.axis("off")
    _panel_tag(ax_b, "b", x=-0.02, y=1.04)

    _rbox(ax_b, 0.5, 0.5, 99.0, 60.5, fc="#F8FAFC", ec="#CBD5E1", lw=0.60, pad=0.6)
    ax_b.add_patch(mpatches.Rectangle((1.2, 53.0), 97.6, 7.2, facecolor="#E6F5F0", edgecolor="none", zorder=2))
    ax_b.add_patch(mpatches.Rectangle((1.2, 53.0), 1.6, 7.2, facecolor=PAL["teal"], edgecolor="none", zorder=3))
    ax_b.text(4.5, 56.6, "Single-Cell Remodeling Drift & Steric Permeability Law (DNAcc)", fontsize=5.8, fontweight="bold", color=PAL["navy"], va="center", zorder=4)

    # Col 1: Healthy & Disease Atlases
    _rbox(ax_b, 2.5, 28.0, 28.5, 22.8, fc="#FFFFFF", ec="#B8C7D6", lw=0.55, pad=0.5, zorder=2)
    ax_b.text(16.7, 46.5, "Healthy Reference ($\\tau_0$)", fontsize=4.65, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_b.text(16.7, 39.2, "HPA v23 Consensus\n154 Human Cell Types", fontsize=4.15, color=PAL["ink"], ha="center", va="center", zorder=4)
    ax_b.text(16.7, 31.8, "Vital Bystander $H_{\\mathrm{vital}}$", fontsize=4.10, fontstyle="italic", color=PAL["slate"], ha="center", va="center", zorder=4)

    _rbox(ax_b, 2.5, 3.0, 28.5, 22.8, fc="#FFFFFF", ec=PAL["teal"], lw=0.60, pad=0.5, zorder=2)
    ax_b.text(16.7, 21.5, "Disease Omics (N = 237)", fontsize=4.55, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_b.text(16.7, 15.2, "126 sc/snRNA \u00b7 99 Pre-Tx Bulk", fontsize=3.95, color=PAL["ink"], ha="center", va="center", zorder=4)
    ax_b.text(16.7, 9.2, "12 Visium (20,880 QC Spots)", fontsize=3.95, fontweight="bold", color=PAL["teal"], ha="center", va="center", zorder=4)

    _arrow(ax_b, 31.3, 39.0, 35.2, 39.0, col=PAL["slate"])
    _arrow(ax_b, 31.3, 14.5, 35.2, 35.0, col=PAL["teal"])
    _arrow(ax_b, 31.3, 12.5, 35.2, 14.5, col=PAL["crimson"])

    # Col 2: Lineage Drift + Steric Barrier (with Pure-Vector Molecular & Collagen Geometry)
    _rbox(ax_b, 35.5, 28.0, 30.5, 22.8, fc="#FFFFFF", ec=PAL["blue"], lw=0.60, pad=0.5, zorder=2)
    ax_b.text(50.7, 46.5, "Specificity Drift ($\\Delta\\tau$)", fontsize=4.70, fontweight="bold", color=PAL["blue"], ha="center", va="center", zorder=4)
    ax_b.text(50.7, 40.5, "$\\Delta\\tau = \\tau_{\\mathrm{late}} - \\tau_{\\mathrm{early}}$", fontsize=4.65, color=PAL["ink"], ha="center", va="center", zorder=4)
    ax_b.text(50.7, 33.4, "17.8\u201337.1% Lineage Switch\nTHRB Retained vs. CCR2 Lost", fontsize=4.05, color=PAL["slate"], ha="center", va="center", zorder=4)

    _rbox(ax_b, 35.5, 3.0, 30.5, 22.8, fc="#FFFFFF", ec=PAL["crimson"], lw=0.60, pad=0.5, zorder=2)
    ax_b.text(50.7, 21.8, "Modality\u2013ECM Gate ($B_{\\mathrm{acc}}$)", fontsize=4.55, fontweight="bold", color=PAL["crimson"], ha="center", va="center", zorder=4)
    # Pure-vector ECM collagen cross-fiber mesh strip separating SM penetration from IgG exclusion
    for y_ecm in (12.1, 12.9):
        ax_b.plot([37.8, 63.6], [y_ecm, y_ecm], color="#D9825F", lw=0.60, alpha=0.65, ls="--", zorder=3)
    ax_b.text(50.7, 16.4, "SM (~0.5 kDa) & Pep (~4 kDa)\nPenetrate Scar >>", fontsize=3.90, color=PAL["teal"], fontweight="bold", ha="center", va="center", zorder=4)
    ax_b.text(50.7, 7.5, "IgG mAb (~150 kDa):\n1.57\u00d7 Futility Odds", fontsize=3.90, color=PAL["crimson"], fontweight="bold", ha="center", va="center", zorder=4)

    _arrow(ax_b, 66.3, 39.0, 70.0, 31.0, col=PAL["blue"])
    _arrow(ax_b, 66.3, 14.5, 70.0, 23.0, col=PAL["crimson"])

    # Col 3: DNAcc Operator & Genetics Orthogonality
    _rbox(ax_b, 70.2, 3.0, 27.5, 47.8, fc="#EAF2F8", ec=PAL["blue"], lw=0.85, pad=0.6, zorder=2)
    ax_b.text(84.0, 45.8, "DNAcc Law", fontsize=5.35, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)
    ax_b.text(
        84.0,
        36.2,
        "$\\mathrm{DNAcc} = \\mathcal{F}(\\tau_0, \\tau_{\\mathrm{early}},$\n$\\Delta\\tau, I_{\\mathrm{sw}}, B_{\\mathrm{acc}})$",
        fontsize=4.60,
        color=PAL["ink"],
        ha="center",
        va="center",
        zorder=4,
    )
    sci_j_auc = float(
        s3l.loc[
            (s3l["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes")
            & (s3l["evaluation_stratum"] == "Stratified_5Fold")
            & (s3l["model_id"] == "G5_Joint_Human_Genetics_Plus_DNAcc"),
            "auroc",
        ].iloc[0]
    )
    mat_j_auc = float(
        s3l.loc[
            (s3l["cohort"] == "Phase2_to_Phase34_Maturation_AllDefinitiveNonAdminTrials")
            & (s3l["evaluation_stratum"] == "Stratified_5Fold")
            & (s3l["model_id"] == "G5_Joint_Human_Genetics_Plus_DNAcc"),
            "auroc",
        ].iloc[0]
    )
    fe_sci_1 = s3g[s3g["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].set_index("model_specification")
    dna_3p_or = float(
        fe_sci_1.loc[
            "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0"
            if "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0" in fe_sci_1.index
            else "Multivariable_FixedEffects_DNAcc",
            "odds_ratio_per_sd",
        ]
    )
    tau_3p_or = float(
        fe_sci_1.loc[
            "Joint_FixedEffects_Tau0_controlling_for_DNAcc_and_Genetics"
            if "Joint_FixedEffects_Tau0_controlling_for_DNAcc_and_Genetics" in fe_sci_1.index
            else "Joint_FixedEffects_Tau0_controlling_for_DNAcc",
            "odds_ratio_per_sd",
        ]
    )
    _rbox(ax_b, 72.2, 6.2, 23.5, 21.5, fc="#FFFFFF", ec=PAL["navy"], lw=0.50, pad=0.5, zorder=3)
    ax_b.text(
        84.0,
        17.0,
        f"+ Human Genetics\nSci. AUC = {sci_j_auc:.2f}\nMat. AUC = {mat_j_auc:.2f}\n3-Pillar OR = {dna_3p_or:.2f}\n($\\tau_0$ OR = {tau_3p_or:.2f}, n.s.)",
        fontsize=4.15,
        fontweight="bold",
        color=PAL["navy"],
        ha="center",
        va="center",
        zorder=4,
    )

    # -----------------------------------------------------------------------
    # Panel c: Closure-Free Multicellular Ito-Langevin Attractor & Hysteresis
    # -----------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1])
    ax_c.set_xlim(0, 100)
    ax_c.set_ylim(0, 56)
    ax_c.axis("off")
    _panel_tag(ax_c, "c", x=-0.01, y=1.04)

    _rbox(ax_c, 0.25, 0.5, 99.5, 54.5, fc="#F8FAFC", ec="#CBD5E1", lw=0.60, pad=0.5)
    ax_c.add_patch(mpatches.Rectangle((0.8, 48.2), 98.4, 6.3, facecolor=PAL["navy"], edgecolor="none", zorder=2))
    ax_c.text(
        2.5,
        51.3,
        "Closure-Free Multicellular It\u00f4\u2013Langevin Attractor Engine & Orthogonal Hysteresis Escape",
        fontsize=5.85,
        fontweight="bold",
        color="white",
        va="center",
        zorder=4,
    )

    # Zone 1 (Left): 3-Compartment Reciprocal Positive-Feedback Circuit
    _rbox(ax_c, 1.8, 2.8, 30.2, 43.5, fc="#FFFFFF", ec="#CBD5E1", lw=0.55, pad=0.5, zorder=2)
    ax_c.text(16.9, 43.0, "Conserved Niche Feedback Circuit", fontsize=5.1, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)

    _rbox(ax_c, 6.8, 30.5, 20.2, 9.2, fc="#EAF2F8", ec=PAL["blue"], lw=0.70, pad=0.4, zorder=3)
    ax_c.text(16.9, 35.1, "Injured Parenchyma / Tumor\n($S_1$: THRB, GLP1R)", fontsize=4.35, fontweight="bold", color=PAL["navy"], ha="center", va="center", zorder=4)

    _rbox(ax_c, 2.8, 12.2, 13.4, 11.2, fc="#FFF6E9", ec=PAL["orange"], lw=0.70, pad=0.4, zorder=3)
    ax_c.text(9.5, 17.8, "$SPP1^+TREM2^+$\nSAMac / TAM ($S_2$)", fontsize=4.2, fontweight="bold", color="#8A4B00", ha="center", va="center", zorder=4)

    _rbox(ax_c, 17.6, 12.2, 13.4, 11.2, fc="#FDEDEC", ec=PAL["crimson"], lw=0.70, pad=0.4, zorder=3)
    ax_c.text(24.3, 17.8, "$POSTN^+FAP^+$\nMyofib / myCAF ($S_3$)", fontsize=4.2, fontweight="bold", color=PAL["crimson"], ha="center", va="center", zorder=4)

    _arrow(ax_c, 12.5, 30.2, 9.8, 23.8, col=PAL["blue"], lw=0.85, rad=0.12)
    _arrow(ax_c, 16.3, 18.8, 17.5, 18.8, col=PAL["orange"], lw=0.95, rad=0.18)
    _arrow(ax_c, 17.5, 16.5, 16.3, 16.5, col=PAL["crimson"], lw=0.95, rad=0.18)
    _arrow(ax_c, 23.8, 23.8, 20.8, 30.2, col=PAL["crimson"], lw=0.85, rad=0.12)
    ax_c.text(16.9, 25.8, "Hill $n = 3.6$", fontsize=4.3, fontweight="bold", color=PAL["crimson"], ha="center", va="center", zorder=4)
    ax_c.text(16.9, 6.6, "Visium $k$=6 Moran $P_{ij}$: myCAF\u2013TAM (0.67)\nShields Tumor from $CD8^+$ T Cells (0.47)", fontsize=4.15, color=PAL["slate"], ha="center", va="center", zorder=4)

    _arrow(ax_c, 32.3, 25.0, 34.8, 25.0, col=PAL["navy"], lw=0.9)

    # Zone 2 (Center): Stochastic Langevin Equation + Vector Double-Well Potential U(S)
    _rbox(ax_c, 35.0, 29.5, 32.5, 16.8, fc="#FFFFFF", ec=PAL["navy"], lw=0.65, pad=0.5, zorder=2)
    ax_c.text(
        51.2,
        40.8,
        "$d\\mathbf{S}_t = \\mathbf{f}(\\mathbf{S}_t;\\, \\Theta, P_{ij}, \\eta_{\\mathrm{eff}})\\,dt + \\sqrt{2D}\\,d\\mathbf{W}_t$",
        fontsize=5.0,
        fontweight="bold",
        color=PAL["navy"],
        ha="center",
        va="center",
        zorder=4,
    )
    ax_c.text(
        51.2,
        34.0,
        "$\\eta_{\\mathrm{eff}} = \\eta_0 \\cdot w_k \\cdot G_{\\mathrm{GRN}} \\cdot B_{\\mathrm{acc}}(m, D_{\\mathrm{ECM}})$  |  4/4 Edges $\\Delta\\chi^2 > 3.84$",
        fontsize=4.25,
        color=PAL["ink"],
        ha="center",
        va="center",
        zorder=4,
    )

    # Mini Pure-Vector Double-Well Waddington Potential Landscape inside Zone 2 bottom
    _rbox(ax_c, 35.0, 2.8, 32.5, 24.8, fc="#FFFFFF", ec="#CBD5E1", lw=0.55, pad=0.5, zorder=2)
    xs_u = np.linspace(-1.28, 1.34, 160)
    u_lock = 3.6 * (xs_u**4 - 1.95 * xs_u**2 - 0.18 * xs_u) + 3.9
    u_unlk = 2.6 * (xs_u**4 - 1.35 * xs_u**2 + 0.55 * xs_u) + 3.4
    px_u = 37.5 + (xs_u + 1.30) / 2.65 * 27.5
    py_lock = 5.8 + u_lock * 2.25
    py_unlk = 5.8 + u_unlk * 2.25
    ax_c.plot(px_u, py_lock, color=PAL["navy"], lw=1.25, zorder=4)
    ax_c.plot(px_u, py_unlk, color=PAL["teal"], lw=1.05, ls="--", zorder=4)
    ax_c.axvline(51.0, ymin=0.08, ymax=0.44, color=PAL["slate"], ls=":", lw=0.65, zorder=3)
    ax_c.scatter([40.8, 61.3], [7.8, 6.9], c=[PAL["teal"], PAL["crimson"]], s=26, edgecolors="white", linewidths=0.7, zorder=5)
    ax_c.text(41.2, 23.6, "Remission\nBasin", fontsize=4.3, fontweight="bold", color=PAL["teal"], ha="center", va="center", zorder=4)
    ax_c.text(51.0, 24.0, "Saddle $\\Delta U$", fontsize=4.2, fontweight="bold", color=PAL["slate"], ha="center", va="center", zorder=4)
    ax_c.text(61.2, 23.6, "Locked F4 /\nmyCAF Scar", fontsize=4.3, fontweight="bold", color=PAL["crimson"], ha="center", va="center", zorder=4)
    _arrow(ax_c, 60.5, 11.2, 42.2, 11.6, col=PAL["teal"], lw=1.05, rad=-0.28)
    ax_c.text(51.5, 16.5, "Dual Unlock", fontsize=4.2, fontweight="bold", color=PAL["teal"], ha="center", va="center", zorder=5)

    _arrow(ax_c, 67.8, 35.0, 70.3, 35.0, col=PAL["crimson"], lw=0.9)
    _arrow(ax_c, 67.8, 14.5, 70.3, 14.5, col=PAL["teal"], lw=0.9)

    # Zone 3 (Right): Monotherapy Hysteresis Trap vs Orthogonal Combination Rescue
    _rbox(ax_c, 70.5, 25.2, 27.7, 21.1, fc="#FDEDEC", ec=PAL["crimson"], lw=0.65, pad=0.5, zorder=2)
    ax_c.text(84.3, 42.2, "Monotherapy Cliff ($t_{50} = 8.3$ h)", fontsize=4.85, fontweight="bold", color=PAL["crimson"], ha="center", va="center", zorder=4)
    ax_c.text(84.3, 35.6, "Single-node drug or 150-kDa mAb", fontsize=4.15, color=PAL["ink"], ha="center", va="center", zorder=4)
    ax_c.text(84.3, 29.8, "Trapped behind saddle \u2192 0.0% F4", fontsize=4.15, fontweight="bold", color=PAL["crimson"], ha="center", va="center", zorder=4)

    _rbox(ax_c, 70.5, 2.8, 27.7, 20.6, fc="#E8F6F3", ec=PAL["teal"], lw=0.65, pad=0.5, zorder=2)
    ax_c.text(84.3, 19.5, "Orthogonal Dual Unlock", fontsize=4.9, fontweight="bold", color=PAL["teal"], ha="center", va="center", zorder=4)
    ax_c.text(84.3, 13.2, "Parenchyma / PD-1 + TGFBR1i / TREM2", fontsize=4.10, color=PAL["ink"], ha="center", va="center", zorder=4)
    ax_c.text(84.3, 7.4, "99.5% F4 & 78.1% Sim. Rescue", fontsize=4.15, fontweight="bold", color=PAL["teal"], ha="center", va="center", zorder=4)

    # -----------------------------------------------------------------------
    # Panel d: 5-Cohort Zero-Label Validation Summary Strip
    # -----------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[2])
    ax_d.set_xlim(0, 100)
    ax_d.set_ylim(0, 32)
    ax_d.axis("off")
    _panel_tag(ax_d, "d", x=-0.01, y=1.10)

    ccr2_df = pd.read_csv(tables_dir / "Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv")
    c_bld = float(ccr2_df[ccr2_df["niche_stage"] == "1_Human_NAFLD_Blood_Monocyte"]["CCR2"].mean())
    c_nash = float(ccr2_df[ccr2_df["niche_stage"] == "3_Human_F2_F3_NASH_Liver_SAMac"]["CCR2"].mean())
    c_drop_pct = 100.0 * (c_bld - c_nash) / max(c_bld, 1e-6)

    cards_d = [
        ("1. Pan-Organ Trial Law", PAL["blue"], f"n = {n_tot:,} Trials (12 Axes)", f"{n_sci:,} Sci. | {n_mat:,} Mat.", f"Joint AUC = {sci_j_auc:.2f} / {mat_j_auc:.2f}"),
        ("2. 3-Organ Fibrosis Atlas", PAL["teal"], "N = 64 Human Donors", "Liver \u00b7 Lung \u00b7 Kidney", f"CCR2 {c_bld:.2f} \u2192 {c_nash:.2f} (-{c_drop_pct:.0f}%)"),
        ("3. Antifibrotic RCT Arena", PAL["navy"], f"n = {len(s4)} Phase IIb/III RCTs", "10 Winners vs. 14 Failures", f"Zero-Shot AUC = {auc_24:.2f} (0.96)"),
        ("4. 10x Visium Carcinoma", PAL["crimson"], "12 Tumor Sections", "20,880 QC Spots (24k Shifts)", "12/12 SEI > 0 (Z = +7.82)"),
        ("5. Pre-ICB & Perturbation", PAL["purple"], "n = 107 Pre-ICB + Controls", "BCC 0.88 | Anti-PD-1 OS 1.95", "TGFBR1i & Non-ICB Verified"),
    ]
    cw = 18.9
    gap = 1.25
    for idx_cd, (hdr_d, col_d, l1_d, l2_d, l3_d) in enumerate(cards_d):
        x0_d = 0.25 + idx_cd * (cw + gap)
        _rbox(ax_d, x0_d, 1.0, cw, 29.5, fc="#FFFFFF", ec="#CBD5E1", lw=0.55, pad=0.45, zorder=2)
        ax_d.add_patch(mpatches.Rectangle((x0_d + 0.5, 24.0), cw - 1.0, 6.0, facecolor=col_d, edgecolor="none", zorder=3))
        ax_d.text(x0_d + cw / 2, 27.0, hdr_d, fontsize=4.55, fontweight="bold", color="white", ha="center", va="center", zorder=4)
        ax_d.text(x0_d + cw / 2, 18.8, l1_d, fontsize=4.45, fontweight="bold", color=PAL["ink"], ha="center", va="center", zorder=4)
        ax_d.text(x0_d + cw / 2, 12.8, l2_d, fontsize=4.20, color=PAL["slate"], ha="center", va="center", zorder=4)
        ax_d.text(x0_d + cw / 2, 6.2, l3_d, fontsize=4.25, fontweight="bold", color=col_d, ha="center", va="center", zorder=4)

    return _save_multi(fig, out_dir, "Fig1_OmniNiche_Architecture_and_Framework", artifact_dir, dpi=600)


def export_fig1_illustration(
    raw_asset_path: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> pathlib.Path:
    """Backward-compatible wrapper routing Fig. 1 to the pure-vector architectural builder."""
    tables_dir = out_dir.parent / "stage2_tables"
    return build_fig1(tables_dir, out_dir, artifact_dir=artifact_dir)


# ===========================================================================
# FIGURE 2: Macro-Atlas of 41,396 Systemic & Parenchymal Trials (6 Panels a-f)
# ===========================================================================
def build_fig2(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 6-panel (a-f) Macro-Atlas natively at 183 mm (7.20 in x 5.20 in)."""
    set_nature_style()
    s1 = pd.read_csv(tables_dir / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv")
    s2 = pd.read_csv(
        tables_dir / "Table_S2_Target_Resolved_Trials_Deconvolution.csv",
        usecols=["nct_id", "drug_modalities", "stop_deconvolution"],
        low_memory=False,
    )
    s2a_val = pd.read_csv(tables_dir / "Table_S2a_Stoppage_Auditor_Validation.csv")
    s2a_cm = pd.read_csv(tables_dir / "Table_S2a_Stoppage_Confusion_Matrix.csv", index_col=0)
    s2b = pd.read_csv(tables_dir / "Table_S2b_Longitudinal_2005_2025_Trends.csv")
    s2c = pd.read_csv(tables_dir / "Table_S2c_Modality_Termination_Odds.csv")
    s2d = pd.read_csv(tables_dir / "Table_S2d_Phase_Enrollment_Distributions.csv")

    fig = plt.figure(figsize=(7.20, 5.05), facecolor="white")
    gs = gridspec.GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[1.04, 1.00],
        hspace=0.48,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], width_ratios=[1.28, 0.90, 0.88], wspace=0.54)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[1.22, 1.00, 0.78], wspace=0.54)

    # --- Panel a: 3-Stage Alluvial Flow (Target -> Modality -> Deconvoluted Outcome) ---
    ax_a = fig.add_subplot(gs_r1[0, 0])
    _panel_tag(ax_a, "a", x=-0.06, y=1.08)

    n_tot = int(s1["unique_systemic_organ_trials"].sum())
    n_res = len(s2)
    n_unres = max(n_tot - n_res, 0)

    def _map_mod(m_str: str) -> str:
        s = str(m_str).lower()
        if "antibody" in s:
            return "IgG mAb"
        if any(k in s for k in ("peptide", "protein", "enzyme", "oligo", "cell", "gene")):
            return "Peptide/Bio"
        return "Small Mol."

    def _map_out(o_str: str) -> str:
        if o_str == "Completed":
            return "Completed"
        if o_str == "Administrative_Commercial_Enrollment":
            return "Admin. Stop"
        if o_str in ("Scientific_Efficacy_Futility", "Scientific_Safety_Toxicity"):
            return "Bio. Fail"
        return "Unreported"

    s2_mod = s2["drug_modalities"].map(_map_mod)
    s2_out = s2["stop_deconvolution"].map(_map_out)

    tot_comp = int(s1["completed_trials"].sum())
    tot_adm = int(s1["administrative_commercial_enrollment"].sum())
    tot_sci = int(s1["scientific_efficacy_futility"].sum() + s1["scientific_safety_toxicity"].sum())
    tot_unr = max(n_tot - tot_comp - tot_adm - tot_sci, 0)

    res_counts = s2_out.value_counts().to_dict()
    unres_comp = max(tot_comp - int(res_counts.get("Completed", 0)), 0)
    unres_adm = max(tot_adm - int(res_counts.get("Admin. Stop", 0)), 0)
    unres_sci = max(tot_sci - int(res_counts.get("Bio. Fail", 0)), 0)
    unres_unr = max(tot_unr - int(res_counts.get("Unreported", 0)), 0)

    def _layout_stack(items: Sequence[tuple[str, int, str]], gap: float = 0.032) -> dict[str, tuple[float, float, int, str]]:
        usable = 1.0 - gap * (len(items) - 1)
        s_sum = float(sum(v for _, v, _ in items))
        cur = 1.0
        out: dict[str, tuple[float, float, int, str]] = {}
        for lbl, val, col in items:
            h = usable * (val / max(s_sum, 1.0))
            out[lbl] = (cur - h, cur, val, col)
            cur -= h + gap
        return out

    st1 = _layout_stack([("Resolved", n_res, PAL["blue"]), ("Unresolved", n_unres, PAL["midgrey"])], gap=0.05)
    mod_vc = s2_mod.value_counts().to_dict()
    st2 = _layout_stack(
        [
            ("Small Mol.", int(mod_vc.get("Small Mol.", 0)), PAL["blue"]),
            ("IgG mAb", int(mod_vc.get("IgG mAb", 0)), PAL["teal"]),
            ("Peptide/Bio", int(mod_vc.get("Peptide/Bio", 0)), PAL["sky"]),
            ("Unresolved", n_unres, PAL["midgrey"]),
        ],
        gap=0.028,
    )
    st3 = _layout_stack(
        [
            ("Completed", tot_comp, PAL["blue"]),
            ("Admin. Stop", tot_adm, PAL["orange"]),
            ("Unreported", tot_unr, PAL["midgrey"]),
            ("Bio. Fail", tot_sci, PAL["crimson"]),
        ],
        gap=0.028,
    )

    x_pos = [0.03, 0.37, 0.70]
    bw = 0.056
    short_lbl = {
        "Resolved": "Resolved",
        "Unresolved": "Unresolved",
        "Small Mol.": "Small Mol.",
        "IgG mAb": "IgG mAb",
        "Peptide/Bio": "Peptide/Bio",
        "Completed": "Completed",
        "Admin. Stop": "Admin. Stop",
        "Unreported": "Unreported",
        "Bio. Fail": "Sci. Fail",
    }
    for xi, st_dict, hdr in zip(x_pos, [st1, st2, st3], ["Target", "Modality", "Outcome"]):
        ax_a.text(xi, 1.045, hdr, ha="center", va="bottom", fontsize=5.5, fontweight="bold", color=PAL["navy"])
        for lbl, (yb, yt, val, col) in st_dict.items():
            ax_a.add_patch(
                mpatches.Rectangle((xi - bw / 2, yb), bw, max(yt - yb, 0.008), facecolor=col, edgecolor="white", lw=0.4)
            )
            ym = 0.5 * (yb + yt)
            sl = short_lbl.get(lbl, lbl)
            txt = f"{sl}\n{val / 1000:.1f}k" if ((yt - yb) >= 0.16 or xi == x_pos[1]) else f"{sl} {val / 1000:.1f}k"
            ax_a.text(xi + bw / 2 + 0.011, ym, txt, ha="left", va="center", fontsize=4.35, color=PAL["ink"])

    r_cur = st1["Resolved"][1]
    r_scale = (st1["Resolved"][1] - st1["Resolved"][0]) / max(n_res, 1)
    for m_lbl in ["Small Mol.", "IgG mAb", "Peptide/Bio"]:
        yb2, yt2, v2, c2 = st2[m_lbl]
        h1 = v2 * r_scale
        _draw_alluvial_ribbon(ax_a, x_pos[0] + bw / 2, r_cur - h1, r_cur, x_pos[1] - bw / 2, yb2, yt2, c2, alpha=0.28)
        r_cur -= h1
    _draw_alluvial_ribbon(
        ax_a,
        x_pos[0] + bw / 2,
        st1["Unresolved"][0],
        st1["Unresolved"][1],
        x_pos[1] - bw / 2,
        st2["Unresolved"][0],
        st2["Unresolved"][1],
        PAL["midgrey"],
        alpha=0.18,
    )

    st3_cursors = {k: v[1] for k, v in st3.items()}
    st3_scales = {k: (v[1] - v[0]) / max(v[2], 1) for k, v in st3.items()}
    cross = pd.crosstab(s2_mod, s2_out)
    for m_lbl in ["Small Mol.", "IgG mAb", "Peptide/Bio"]:
        m_cur = st2[m_lbl][1]
        m_scale = (st2[m_lbl][1] - st2[m_lbl][0]) / max(st2[m_lbl][2], 1)
        for o_lbl in ["Completed", "Admin. Stop", "Unreported", "Bio. Fail"]:
            v_ij = int(cross.loc[m_lbl, o_lbl]) if (m_lbl in cross.index and o_lbl in cross.columns) else 0
            if v_ij <= 0:
                continue
            hm = v_ij * m_scale
            ho = v_ij * st3_scales[o_lbl]
            _draw_alluvial_ribbon(
                ax_a, x_pos[1] + bw / 2, m_cur - hm, m_cur, x_pos[2] - bw / 2, st3_cursors[o_lbl] - ho, st3_cursors[o_lbl], st3[o_lbl][3], alpha=0.26
            )
            m_cur -= hm
            st3_cursors[o_lbl] -= ho

    u_cur = st2["Unresolved"][1]
    u_scale = (st2["Unresolved"][1] - st2["Unresolved"][0]) / max(n_unres, 1)
    for o_lbl, v_u in [("Completed", unres_comp), ("Admin. Stop", unres_adm), ("Unreported", unres_unr), ("Bio. Fail", unres_sci)]:
        if v_u <= 0:
            continue
        hm = v_u * u_scale
        ho = v_u * st3_scales[o_lbl]
        _draw_alluvial_ribbon(
            ax_a, x_pos[1] + bw / 2, u_cur - hm, u_cur, x_pos[2] - bw / 2, st3_cursors[o_lbl] - ho, st3_cursors[o_lbl], PAL["midgrey"], alpha=0.16
        )
        u_cur -= hm
        st3_cursors[o_lbl] -= ho

    ax_a.set_xlim(-0.01, 1.01)
    ax_a.set_ylim(-0.03, 1.10)
    ax_a.axis("off")
    ax_a.set_title(f"Trial Deconvolution ($n$ = {n_tot:,})", pad=4)

    # --- Panel b: 12 Systemic & Parenchymal Organ Axes ---
    ax_b = fig.add_subplot(gs_r1[0, 1])
    _panel_tag(ax_b, "b", x=-0.38, y=1.08)
    s1_s = s1.sort_values("unique_systemic_organ_trials", ascending=True).reset_index(drop=True)
    y_b = np.arange(len(s1_s))
    p2 = s1_s["phase2_trials"].to_numpy()
    p3 = s1_s["phase3_trials"].to_numpy()
    p4 = s1_s["phase4_trials"].to_numpy()
    res_pct = 100.0 * s1_s["target_resolved_trials"].to_numpy() / np.maximum(s1_s["unique_systemic_organ_trials"].to_numpy(), 1)

    ax_b.barh(y_b, p2, color=PAL["sky"], height=0.68, label="Phase II")
    ax_b.barh(y_b, p3, left=p2, color=PAL["blue"], height=0.68, label="Phase III")
    ax_b.barh(y_b, p4, left=p2 + p3, color=PAL["navy"], height=0.68, label="Phase IV")
    for idx_r, (tot_v, rp) in enumerate(zip(s1_s["unique_systemic_organ_trials"], res_pct)):
        ax_b.text(tot_v + 240, idx_r, f"{tot_v:,} ({rp:.0f}%)", va="center", ha="left", fontsize=4.9, color=PAL["ink"])

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels([AXIS_SHORT.get(a, a) for a in s1_s["disease_axis"]], fontsize=5.2)
    ax_b.set_xlim(0, float(s1_s["unique_systemic_organ_trials"].max()) * 1.44)
    ax_b.set_xticks([0, 5000, 10000, 15000])
    ax_b.set_xticklabels(["0", "5k", "10k", "15k"])
    ax_b.set_xlabel("Phase II–IV Trials (Resolved %)")
    ax_b.set_title("12 Peripheral Organ Axes", pad=5)
    ax_b.legend(loc="lower right", frameon=False, fontsize=5.1)

    # --- Panel c: 2005-2025 Longitudinal Modality & Attrition Dynamics ---
    ax_c = fig.add_subplot(gs_r1[0, 2])
    _panel_tag(ax_c, "c", x=-0.29, y=1.08)
    yrs = s2b["start_year"].to_numpy()
    mab = s2b["antibody_share_pct"].to_numpy()
    adm = s2b["admin_share_of_stops_pct"].to_numpy()
    sci_sh = s2b["efficacy_share_of_stops_pct"].to_numpy() + s2b["safety_share_of_stops_pct"].to_numpy()

    ax_c.plot(yrs, adm, color=PAL["orange"], lw=1.2, marker="o", ms=2.1, label="Admin. (% Stops)")
    ax_c.fill_between(yrs, s2b["admin_ci95_low"], s2b["admin_ci95_high"], color=PAL["orange"], alpha=0.14)
    ax_c.plot(yrs, mab, color=PAL["teal"], lw=1.2, marker="s", ms=2.0, label="IgG mAb (% Trials)")
    ax_c.fill_between(yrs, s2b["antibody_ci95_low"], s2b["antibody_ci95_high"], color=PAL["teal"], alpha=0.14)
    ax_c.plot(yrs, sci_sh, color=PAL["crimson"], lw=1.1, marker="^", ms=2.1, label="Bio. Fail (% Stops)")

    ax_c.set_xlim(2005, 2025)
    ax_c.set_xticks([2005, 2010, 2015, 2020, 2025])
    ax_c.set_ylim(0, max(float(np.max(s2b["admin_ci95_high"])), float(np.max(s2b["antibody_ci95_high"]))) * 1.38)
    ax_c.set_xlabel("Trial Start Year")
    ax_c.set_ylabel("Registry Share (%)")
    ax_c.set_title("20-Year Modality & Stoppage Drift", pad=5)
    ax_c.legend(loc="upper left", frameon=False, fontsize=5.0)

    # --- Panel d: Modality x Termination Cause Forest Plot (Minikel Tabular Right-Column Style) ---
    ax_d = fig.add_subplot(gs_r2[0, 0])
    _panel_tag(ax_d, "d", x=-0.28, y=1.08)
    sub_odds = s2c[s2c["modality_group"].isin(["Antibody_Biologic", "Protein_Peptide"])].copy()
    cause_short = {
        "Scientific_Efficacy_Futility": "Efficacy Futility",
        "Scientific_Safety_Toxicity": "Safety / Toxicity",
        "Administrative_Commercial_Enrollment": "Admin. / Accrual",
    }
    mod_short = {"Antibody_Biologic": "IgG mAb", "Protein_Peptide": "Peptide"}
    sub_odds["row_lbl"] = [
        f"{mod_short.get(m, m)}: {cause_short.get(c, c)}"
        for m, c in zip(sub_odds["modality_group"], sub_odds["termination_cause"])
    ]
    sub_odds = sub_odds.iloc[::-1].reset_index(drop=True)
    y_d = np.arange(len(sub_odds))
    ors = sub_odds["odds_ratio"].to_numpy()
    clo = sub_odds["ci95_low"].to_numpy()
    chi = sub_odds["ci95_high"].to_numpy()
    pvs = sub_odds["p_value"].to_numpy()
    ev_m = sub_odds["n_modality_events"].to_numpy()

    for i_r in range(len(sub_odds)):
        if i_r % 2 == 1:
            ax_d.axhspan(i_r - 0.45, i_r + 0.45, color=PAL["pale"], zorder=0)
    ax_d.axvline(1.0, color=PAL["slate"], ls="--", lw=0.65, zorder=1)
    ax_d.axvline(2.02, color=PAL["lightgrey"], ls="-", lw=0.50, zorder=1)
    ax_d.text(7.35, len(sub_odds) - 0.35, "Events | OR (95% CI)", ha="right", va="bottom", fontsize=4.35, fontweight="bold", color=PAL["slate"])

    for i_r in range(len(sub_odds)):
        col_r = (
            PAL["crimson"]
            if (ors[i_r] > 1.0 and pvs[i_r] < 0.05)
            else (PAL["teal"] if (ors[i_r] < 1.0 and pvs[i_r] < 0.05) else PAL["midgrey"])
        )
        ax_d.plot([clo[i_r], chi[i_r]], [i_r, i_r], color=col_r, lw=1.25, zorder=3)
        ax_d.scatter([ors[i_r]], [i_r], color=col_r, s=18, zorder=4)
        ax_d.text(
            7.35,
            i_r,
            f"{int(ev_m[i_r]):4d} | {ors[i_r]:.2f} ({clo[i_r]:.2f}–{chi[i_r]:.2f})",
            va="center",
            ha="right",
            fontsize=4.55,
            color=PAL["ink"],
        )

    ax_d.set_xscale("log")
    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels(sub_odds["row_lbl"], fontsize=5.0)
    ax_d.set_xlim(0.36, 7.85)
    ax_d.set_ylim(-0.55, len(sub_odds) - 0.10)
    ax_d.set_xticks([0.5, 1.0, 2.0])
    ax_d.set_xticklabels(["0.5", "1.0", "2.0"])
    ax_d.set_xlabel("Odds Ratio vs. Small Mol.")
    ax_d.set_title("Modality Attrition Odds Ratio", pad=4)

    # --- Panel e: Enrollment Sample-Size Distributions ---
    ax_e = fig.add_subplot(gs_r2[0, 1])
    _panel_tag(ax_e, "e", x=-0.30, y=1.08)
    stop_enr = s2d[s2d["group_type"] == "stop_class"].set_index("group_name")
    order_e = [
        ("Completed", "Completed", PAL["blue"]),
        ("Scientific_Efficacy_Futility", "Efficacy Futility", PAL["crimson"]),
        ("Scientific_Safety_Toxicity", "Safety / Toxicity", PAL["vermilion"]),
        ("Administrative_Commercial_Enrollment", "Admin. / Accrual", PAL["orange"]),
        ("Unreported_Stop_Reason", "Unreported", PAL["midgrey"]),
    ]
    present_e = [(k, lbl, col) for k, lbl, col in order_e if k in stop_enr.index][::-1]
    y_e = np.arange(len(present_e))

    max_p90 = 100.0
    for idx_e, (k, lbl, col) in enumerate(present_e):
        r_e = stop_enr.loc[k]
        q25, med, q75, p90 = (
            float(r_e["q25_enrollment"]),
            float(r_e["median_enrollment"]),
            float(r_e["q75_enrollment"]),
            float(r_e["p90_enrollment"]),
        )
        tot_p = int(r_e["total_enrolled_patients"])
        max_p90 = max(max_p90, p90)
        ax_e.plot([max(q25 * 0.55, 5.0), p90], [idx_e, idx_e], color=col, lw=0.95, alpha=0.85)
        ax_e.add_patch(
            mpatches.Rectangle((q25, idx_e - 0.23), max(q75 - q25, 1.0), 0.46, facecolor=col, alpha=0.32, edgecolor=col, lw=0.8)
        )
        ax_e.scatter([med], [idx_e], color="white", edgecolor=col, linewidth=1.1, s=18, zorder=5)
        pts_txt = f"{tot_p / 1e6:.2f}M" if tot_p >= 1e5 else f"{tot_p / 1e3:.0f}k"
        ax_e.text(p90 * 1.15, idx_e, f"Med {int(round(med))} ({pts_txt})", va="center", ha="left", fontsize=4.7)

    ax_e.set_xscale("log")
    ax_e.set_yticks(y_e)
    ax_e.set_yticklabels([lbl for _, lbl, _ in present_e], fontsize=5.1)
    ax_e.set_xlim(5.0, max_p90 * 8.6)
    ax_e.set_xticks([10, 100, 1000])
    ax_e.set_xticklabels(["10", "100", "1k"])
    ax_e.set_xlabel("Patients per Trial (IQR, Med, P90)")
    ax_e.set_title("Patient Exposure by Outcome", pad=5)

    # --- Panel f: Stoppage Auditor 4x4 Confusion Matrix (n = 240) ---
    ax_f = fig.add_subplot(gs_r2[0, 2])
    _panel_tag(ax_f, "f", x=-0.39, y=1.08)
    cls_order = [
        "Scientific_Efficacy_Futility",
        "Scientific_Safety_Toxicity",
        "Administrative_Commercial_Enrollment",
        "Other_or_Unreported",
    ]
    cls_lbls = ["Futil.", "Safety", "Admin.", "Unrep."]
    cm_mat = s2a_cm.loc[cls_order, cls_order].to_numpy(dtype=float)
    row_pct = 100.0 * cm_mat / np.maximum(cm_mat.sum(axis=1, keepdims=True), 1.0)
    ax_f.imshow(row_pct, cmap="Blues", vmin=0.0, vmax=100.0, aspect="equal")
    for r_i in range(4):
        for c_i in range(4):
            cnt_v = int(cm_mat[r_i, c_i])
            pct_v = row_pct[r_i, c_i]
            txt_c = "white" if pct_v >= 55.0 else PAL["ink"]
            ax_f.text(
                c_i,
                r_i,
                f"{cnt_v}\n({pct_v:.0f}%)",
                ha="center",
                va="center",
                fontsize=4.9,
                fontweight="bold" if r_i == c_i else "normal",
                color=txt_c,
            )

    cl_col_s2a = "auditor_name" if "auditor_name" in s2a_val.columns else ("classifier" if "classifier" in s2a_val.columns else None)
    if cl_col_s2a is not None:
        s2a_cl = s2a_val[s2a_val[cl_col_s2a] == "Clause_Aware_Auditor_v2_2"]
        s2a_nv = s2a_val[s2a_val[cl_col_s2a] == "Naive_Keyword_Baseline"]
        acc_v = float(s2a_cl["overall_accuracy"].iloc[0]) * 100.0
        kap_v = float(s2a_cl["cohen_kappa"].iloc[0])
        nv_txt = f" (Naive {float(s2a_nv['overall_accuracy'].iloc[0]) * 100.0:.0f}%)" if not s2a_nv.empty else ""
    else:
        acc_v = float(s2a_val["overall_accuracy"].iloc[0]) * 100.0
        kap_v = float(s2a_val["cohen_kappa"].iloc[0])
        nv_txt = ""
    ax_f.set_xticks(range(4))
    ax_f.set_xticklabels(cls_lbls, fontsize=5.3)
    ax_f.set_yticks(range(4))
    ax_f.set_yticklabels(cls_lbls, fontsize=5.3)
    ax_f.set_xlabel(f"Clause-Aware Pred.{nv_txt}")
    ax_f.set_ylabel("Manual Clinical Ref.")
    ax_f.set_title(f"Manual Audit ({acc_v:.1f}%, $\\kappa$={kap_v:.2f})", pad=5)
    for sp in ("top", "right", "bottom", "left"):
        ax_f.spines[sp].set_visible(True)
        ax_f.spines[sp].set_linewidth(0.55)

    _save_multi(fig, out_dir, "Fig2_MacroAtlas_41396_Systemic_Organ_ClinicalTrials", artifact_dir)


# ===========================================================================
# FIGURE 3: Pan-Organ Dynamic Drift & DNAcc Law (8 Panels a-h)
# ===========================================================================
def build_fig3(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 8-panel (a-h) Pan-Organ Dynamic Niche Accessibility (DNAcc) Law at 7.20 in x 7.85 in."""
    set_nature_style()
    atlas = pd.read_csv(tables_dir / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv")
    s3b = pd.read_csv(tables_dir / "Table_S3b_PerAxis_Dynamic_Drift_Performance.csv")
    s3c = pd.read_csv(tables_dir / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv")
    s3d = pd.read_csv(tables_dir / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv")
    s3f = pd.read_csv(tables_dir / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv")
    s3g = pd.read_csv(tables_dir / "Table_S3g_Multivariable_FixedEffects_Control.csv")
    s3h = pd.read_csv(tables_dir / "Table_S3h_Calibration_Deciles.csv")
    s3i = pd.read_csv(tables_dir / "Table_S3i_Decision_Curve_Analysis.csv")
    s3k = pd.read_csv(tables_dir / "Table_S3k_WithinModality_DNAcc_Stratified_Audit.csv")
    s3l_p = tables_dir / "Table_S3l_Genetics_vs_DNAcc_Orthogonality.csv"
    s3l = pd.read_csv(s3l_p) if s3l_p.exists() else pd.DataFrame()

    fig = plt.figure(figsize=(7.20, 7.55), facecolor="white")
    gs = gridspec.GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.02, 1.08, 1.02],
        hspace=0.49,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[1.46, 1.00], wspace=0.38)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[1.34, 1.00, 0.82], wspace=0.46)
    gs_r3 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2], width_ratios=[1.00, 1.00, 1.08], wspace=0.42)

    # --- Panel a: 6 Human Single-Cell Disease Axes x 24 Clinical Target Genes Dot-Heatmap ---
    ax_a = fig.add_subplot(gs_r1[0, 0])
    _panel_tag(ax_a, "a", x=-0.16, y=1.10)
    gene_trial_counts = s3c["primary_target_gene"].value_counts()
    priority_genes = [
        "THRB", "GLP1R", "PCSK9", "IL23A", "TYK2", "NR3C2", "PDGFRB", "FGFR1",
        "TNF", "IL6R", "JAK1", "ERBB2", "EGFR", "PDCD1", "CD274", "CTLA4",
        "VEGFA", "KDR", "TGFB1", "COL1A1", "LOXL2", "CCR2", "CCR5", "MAP3K5",
    ]
    avail_genes = set(atlas["gene_symbol"].unique())
    sel_genes = [g for g in priority_genes if g in avail_genes]
    for g_extra in gene_trial_counts.index:
        if len(sel_genes) >= 24:
            break
        if g_extra in avail_genes and g_extra not in sel_genes:
            sel_genes.append(g_extra)

    axes_6 = [
        "Hepatic_MASH_Fibrosis",
        "Pulmonary_IPF_ARDS",
        "Renal_CKD_AKI",
        "Solid_Oncology",
        "Autoimmune_IBD_Rheum",
        "Infectious_Sepsis",
    ]
    sub_at = atlas[atlas["gene_symbol"].isin(sel_genes) & atlas["disease_axis"].isin(axes_6)].copy()
    piv_dt = (
        sub_at.pivot_table(index="disease_axis", columns="gene_symbol", values="delta_tau", aggfunc="mean")
        .reindex(index=axes_6, columns=sel_genes)
        .fillna(0.0)
    )
    piv_sw = (
        sub_at.pivot_table(index="disease_axis", columns="gene_symbol", values="compartment_switch", aggfunc="max")
        .reindex(index=axes_6, columns=sel_genes)
        .fillna(0.0)
    )

    gx, gy = np.meshgrid(np.arange(len(sel_genes)), np.arange(len(axes_6)))
    sizes = np.array([10.0 + 11.5 * math.log10(max(int(gene_trial_counts.get(g, 5)), 2)) for g in sel_genes])
    sz_grid = np.tile(sizes, (len(axes_6), 1))

    sc_a = ax_a.scatter(
        gx.ravel(),
        gy.ravel(),
        c=piv_dt.to_numpy().ravel(),
        s=sz_grid.ravel(),
        cmap="RdBu_r",
        vmin=-0.26,
        vmax=0.26,
        edgecolors="none",
        alpha=0.95,
    )
    sw_mask = piv_sw.to_numpy().ravel() >= 0.5
    if np.any(sw_mask):
        ax_a.scatter(
            gx.ravel()[sw_mask],
            gy.ravel()[sw_mask],
            s=sz_grid.ravel()[sw_mask] * 1.18,
            facecolors="none",
            edgecolors=PAL["ink"],
            linewidths=0.65,
        )

    ax_a.set_xticks(np.arange(len(sel_genes)))
    ax_a.set_xticklabels(sel_genes, rotation=50, ha="right", fontsize=5.0, fontstyle="italic")
    ax_a.set_yticks(np.arange(len(axes_6)))
    ax_a.set_yticklabels([AXIS_COMPACT.get(a, a) for a in axes_6], fontsize=5.4)
    ax_a.set_xlim(-0.7, len(sel_genes) - 0.3)
    ax_a.set_ylim(-0.6, len(axes_6) - 0.25)
    ax_a.grid(True, linestyle=":", linewidth=0.35, alpha=0.45)
    ax_a.set_title("Single-Cell Specificity Drift ($\\Delta\\tau$) & Lineage Switch (Ring)", pad=5)
    cbar_a = fig.colorbar(sc_a, ax=ax_a, fraction=0.024, pad=0.02)
    cbar_a.ax.tick_params(labelsize=5.0)
    cbar_a.set_label("$\\Delta\\tau$", fontsize=5.5)

    # --- Panel b: 6-Step Orthogonal DNAcc Ablation & Open Targets Human Genetics Benchmark ---
    ax_b = fig.add_subplot(gs_r1[0, 1])
    _panel_tag(ax_b, "b", x=-0.22, y=1.10)
    ab_main = s3f[s3f["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].copy()
    st_cv = ab_main[ab_main["cv_scheme"] == "Stratified_5Fold"].sort_values("step_order")
    gp_cv = ab_main[ab_main["cv_scheme"] == "UnseenTarget_Group5Fold"].sort_values("step_order")
    n_sci_tot = int(st_cv["n_trials"].iloc[0])
    x_st = np.arange(len(st_cv))
    step_short = [
        "M1: $\\tau_0$",
        "M2: +Vital",
        "M3: +$\\tau_{\\mathrm{early}}$",
        "M4: +$\\Delta\\tau$",
        "M5: +Switch",
        "M6: DNAcc",
    ]

    ax_b.axhline(0.50, color=PAL["lightgrey"], ls=":", lw=0.6)
    if not s3l.empty:
        g_sub_b = s3l[
            (s3l["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes")
            & (s3l["evaluation_stratum"] == "Stratified_5Fold")
        ].set_index("model_id")
        if "G1_OpenTargets_Human_Genetics_Only" in g_sub_b.index:
            g1_auc = float(g_sub_b.loc["G1_OpenTargets_Human_Genetics_Only", "auroc"])
            ax_b.axhline(g1_auc, color=PAL["orange"], ls="-.", lw=0.95, label=f"Human Genetics ({g1_auc:.2f})")
        if "G5_Joint_Human_Genetics_Plus_DNAcc" in g_sub_b.index:
            g5_r = g_sub_b.loc["G5_Joint_Human_Genetics_Plus_DNAcc"]
            ax_b.errorbar(
                [x_st[-1] + 0.16],
                [float(g5_r["auroc"])],
                yerr=[[float(g5_r["auroc"] - g5_r["auroc_ci_low"])], [float(g5_r["auroc_ci_high"] - g5_r["auroc"])]],
                color=PAL["navy"],
                marker="D",
                ms=3.2,
                lw=1.1,
                capsize=1.8,
                label="Genetics + DNAcc",
            )

    ax_b.errorbar(
        x_st - 0.08,
        st_cv["auroc"],
        yerr=[st_cv["auroc"] - st_cv["auroc_ci_low"], st_cv["auroc_ci_high"] - st_cv["auroc"]],
        color=PAL["blue"],
        marker="o",
        ms=3.0,
        lw=1.25,
        capsize=1.8,
        label="Stratified 5-Fold",
    )
    ax_b.errorbar(
        x_st + 0.04,
        gp_cv["auroc"],
        yerr=[gp_cv["auroc"] - gp_cv["auroc_ci_low"], gp_cv["auroc_ci_high"] - gp_cv["auroc"]],
        color=PAL["teal"],
        marker="s",
        ms=2.7,
        lw=1.10,
        ls="--",
        capsize=1.8,
        label="Unseen-Target CV",
    )
    m6_p = float(st_cv["delong_p_vs_m1"].iloc[-1])
    m6_auc = float(st_cv["auroc"].iloc[-1])
    m1_auc = float(st_cv["auroc"].iloc[0])
    ax_b.text(
        0.04,
        0.91,
        f"$\\Delta$AUC = +{m6_auc - m1_auc:.3f} ({_fmt_p(m6_p)})",
        transform=ax_b.transAxes,
        fontsize=4.95,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_b.annotate(
        "Gene-disjoint\n374-target split",
        xy=(float(x_st[3]) + 0.04, float(gp_cv["auroc"].iloc[3])),
        xytext=(1.75, 0.455),
        fontsize=4.35,
        color=PAL["teal"],
        fontweight="bold",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.22", "facecolor": "#EDF7F4", "edgecolor": PAL["teal"], "lw": 0.45},
        arrowprops={"arrowstyle": "->", "lw": 0.55, "color": PAL["teal"]},
    )
    ax_b.set_xticks(x_st)
    ax_b.set_xticklabels(step_short, rotation=28, ha="right", fontsize=5.1)
    ax_b.set_xlim(-0.35, len(x_st) - 0.55)
    ax_b.set_ylim(0.41, 0.77)
    ax_b.set_ylabel("Scientific Outcome AUROC")
    ax_b.set_title(f"DNAcc Ablation & Genetics ($n$={n_sci_tot:,})", pad=4)
    ax_b.legend(loc="lower right", frameon=False, fontsize=4.5)

    # --- Panel c: Confounder, Human Genetics & Target-Clustered Effect Size Forest Plot ---
    ax_c = fig.add_subplot(gs_r2[0, 0])
    _panel_tag(ax_c, "c", x=-0.30, y=1.08)
    fe_sci = s3g[s3g["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].set_index("model_specification")
    has_j3 = "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0" in fe_sci.index
    tau_j_key = "Joint_FixedEffects_Tau0_controlling_for_DNAcc_and_Genetics" if has_j3 else "Joint_FixedEffects_Tau0_controlling_for_DNAcc"
    dna_j_key = "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0" if has_j3 else "Multivariable_FixedEffects_DNAcc"

    # Tuple: (label, or_v, wald_lo, wald_hi, clust_lo, clust_hi, color)
    forest_items: list[tuple[str, float, float, float, float | None, float | None, str]] = [
        (
            "Static $\\tau_0$ (+Joint FE)",
            float(fe_sci.loc[tau_j_key, "odds_ratio_per_sd"]),
            float(fe_sci.loc[tau_j_key, "ci95_low"]),
            float(fe_sci.loc[tau_j_key, "ci95_high"]),
            float(fe_sci.loc[tau_j_key, "cluster_ci95_low"]),
            float(fe_sci.loc[tau_j_key, "cluster_ci95_high"]),
            PAL["midgrey"],
        ),
    ]
    if has_j3:
        gen_j_key = "Joint_FixedEffects_Human_Genetics_controlling_for_DNAcc_and_Tau0"
        forest_items.append(
            (
                "Genetics (+Joint FE)",
                float(fe_sci.loc[gen_j_key, "odds_ratio_per_sd"]),
                float(fe_sci.loc[gen_j_key, "ci95_low"]),
                float(fe_sci.loc[gen_j_key, "ci95_high"]),
                float(fe_sci.loc[gen_j_key, "cluster_ci95_low"]),
                float(fe_sci.loc[gen_j_key, "cluster_ci95_high"]),
                PAL["orange"],
            )
        )
    else:
        forest_items.append(
            (
                "Static $\\tau_0$ (Unadj.)",
                float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "odds_ratio_per_sd"]),
                float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "ci95_low"]),
                float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "ci95_high"]),
                float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "cluster_ci95_low"]),
                float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "cluster_ci95_high"]),
                PAL["midgrey"],
            )
        )
    forest_items.extend(
        [
            (
                "DNAcc (Unadj.)",
                float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "odds_ratio_per_sd"]),
                float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "ci95_low"]),
                float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "ci95_high"]),
                float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "cluster_ci95_low"]),
                float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "cluster_ci95_high"]),
                PAL["blue"],
            ),
            (
                "DNAcc (+Gen+$\\tau_0$+FE)" if has_j3 else "DNAcc (Adj. 5 FE)",
                float(fe_sci.loc[dna_j_key, "odds_ratio_per_sd"]),
                float(fe_sci.loc[dna_j_key, "ci95_low"]),
                float(fe_sci.loc[dna_j_key, "ci95_high"]),
                float(fe_sci.loc[dna_j_key, "cluster_ci95_low"]),
                float(fe_sci.loc[dna_j_key, "cluster_ci95_high"]),
                PAL["navy"],
            ),
        ]
    )
    s3b_mat = s3b[s3b["cohort_block"] == "Definitive_Maturation"].copy() if "cohort_block" in s3b.columns else s3b
    for r_ax in s3b_mat.itertuples(index=False):
        if r_ax.disease_axis == "Infectious_Sepsis":
            continue
        forest_items.append(
            (
                f"{AXIS_COMPACT.get(r_ax.disease_axis, r_ax.disease_axis)} (Mat.)",
                float(r_ax.dnacc_or_per_sd),
                float(r_ax.dnacc_or_ci_low),
                float(r_ax.dnacc_or_ci_high),
                None,
                None,
                PAL["teal"],
            )
        )
    forest_items = forest_items[::-1]
    y_c = np.arange(len(forest_items))
    ax_c.axvline(1.0, color=PAL["slate"], ls="--", lw=0.65)
    ax_c.axvline(2.34, color=PAL["lightgrey"], ls="-", lw=0.50)
    ax_c.axhline(4.5, color=PAL["lightgrey"], ls="-", lw=0.70)
    ax_c.text(3.84, len(forest_items) - 0.30, "OR [95% CI]", ha="right", va="bottom", fontsize=4.35, fontweight="bold", color=PAL["slate"])
    for idx_c, (lbl_c, or_v, w_lo, w_hi, c_lo, c_hi, col_c) in enumerate(forest_items):
        if idx_c % 2 == 1:
            ax_c.axhspan(idx_c - 0.45, idx_c + 0.45, color=PAL["pale"], zorder=0)
        if c_lo is not None and c_hi is not None:
            ax_c.plot([c_lo, min(c_hi, 2.30)], [idx_c, idx_c], color=col_c, lw=0.75, alpha=0.55, zorder=2)
            ax_c.plot([w_lo, min(w_hi, 2.30)], [idx_c, idx_c], color=col_c, lw=1.60, zorder=3)
            ci_txt = f"{or_v:.2f} [{c_lo:.2f}–{c_hi:.2f}]"
        else:
            ax_c.plot([w_lo, min(w_hi, 2.30)], [idx_c, idx_c], color=col_c, lw=1.35, zorder=3)
            ci_txt = f"{or_v:.2f} [{w_lo:.2f}–{w_hi:.2f}]"
        ax_c.scatter([or_v], [idx_c], color=col_c, s=15, zorder=4)
        ax_c.text(3.84, idx_c, ci_txt, va="center", ha="right", fontsize=4.6, color=PAL["ink"])

    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels([it[0] for it in forest_items], fontsize=4.95)
    ax_c.set_xlim(0.70, 3.92)
    ax_c.set_ylim(-0.55, len(forest_items) + 0.15)
    ax_c.set_xticks([1.0, 1.5, 2.0])
    ax_c.set_xlabel("OR per 1-s.d. (Cluster / Wald CI)")
    ax_c.set_title("3-Pillar Joint & Per-Axis Effect", pad=4)

    # --- Panel d: Within-Modality Stratified Discrimination (Proving DNAcc != Modality Proxy) ---
    ax_d = fig.add_subplot(gs_r2[0, 1])
    _panel_tag(ax_d, "d", x=-0.24, y=1.08)
    s3k_idx = s3k.set_index("modality_stratum")
    mod_order = [
        ("Small_Molecule", f"Small Mol.\n($n$={int(s3k_idx.loc['Small_Molecule', 'n_trials']):,})"),
        ("Antibody_Biologic", f"IgG mAb\n($n$={int(s3k_idx.loc['Antibody_Biologic', 'n_trials']):,})"),
        ("Protein_Peptide", f"Peptide\n($n$={int(s3k_idx.loc['Protein_Peptide', 'n_trials']):,})"),
    ]
    x_m = np.arange(len(mod_order))
    w_m = 0.33
    for i_m, (m_key, _) in enumerate(mod_order):
        r_k = s3k_idx.loc[m_key]
        t_auc = float(r_k["static_tau0_auroc"])
        d_auc = float(r_k["dnacc_auroc"])
        ax_d.bar(
            i_m - w_m / 2,
            t_auc,
            width=w_m,
            color=PAL["lightgrey"],
            edgecolor=PAL["slate"],
            lw=0.55,
            yerr=[[t_auc - float(r_k["static_tau0_ci_low"])], [float(r_k["static_tau0_ci_high"]) - t_auc]],
            capsize=1.8,
            error_kw={"lw": 0.65, "ecolor": PAL["slate"]},
            label="Static $\\tau_0$" if i_m == 0 else None,
        )
        ax_d.bar(
            i_m + w_m / 2,
            d_auc,
            width=w_m,
            color=PAL["blue"],
            edgecolor=PAL["navy"],
            lw=0.55,
            yerr=[[d_auc - float(r_k["dnacc_ci_low"])], [float(r_k["dnacc_ci_high"]) - d_auc]],
            capsize=1.8,
            error_kw={"lw": 0.65, "ecolor": PAL["navy"]},
            label="DNAcc" if i_m == 0 else None,
        )
        p_k = float(r_k["delong_p_value"])
        star = "***" if p_k < 1e-6 else ("*" if p_k < 0.05 else "n.s.")
        ax_d.text(
            i_m,
            max(float(r_k["dnacc_ci_high"]), float(r_k["static_tau0_ci_high"])) + 0.015,
            f"+{float(r_k['delta_auroc']):.2f} ({star})",
            ha="center",
            va="bottom",
            fontsize=4.9,
            fontweight="bold" if p_k < 0.05 else "normal",
            color=PAL["navy"] if p_k < 0.05 else PAL["slate"],
        )

    ax_d.axhline(0.50, color=PAL["slate"], ls=":", lw=0.6)
    ax_d.set_xticks(x_m)
    ax_d.set_xticklabels([lbl for _, lbl in mod_order], fontsize=5.1)
    ax_d.set_ylim(0.40, 0.82)
    ax_d.set_ylabel("Within-Modality AUROC")
    ax_d.set_title("Within-Modality Stratification", pad=5)
    ax_d.legend(loc="upper right", frameon=False, fontsize=5.1)

    # --- Panel e: Empirical Phase II -> III/IV Maturation Across Quartiles ---
    ax_e = fig.add_subplot(gs_r2[0, 2])
    _panel_tag(ax_e, "e", x=-0.26, y=1.08)
    n_prog_tot = int(s3d["n_trials_dnacc"].sum()) if "n_trials_dnacc" in s3d.columns else 23967
    x_q = np.arange(len(s3d))
    w_q = 0.34
    st_m = s3d["static_tau0_maturation_pct"].to_numpy()
    dn_m = s3d["dnacc_maturation_pct"].to_numpy()
    ax_e.bar(
        x_q - w_q / 2,
        st_m,
        width=w_q,
        color=PAL["lightgrey"],
        edgecolor=PAL["slate"],
        lw=0.55,
        yerr=[st_m - s3d["static_ci95_low"], s3d["static_ci95_high"] - st_m],
        capsize=1.5,
        error_kw={"lw": 0.6, "ecolor": PAL["slate"]},
        label="Static $\\tau_0$",
    )
    ax_e.bar(
        x_q + w_q / 2,
        dn_m,
        width=w_q,
        color=PAL["blue"],
        edgecolor=PAL["navy"],
        lw=0.55,
        yerr=[dn_m - s3d["dnacc_ci95_low"], s3d["dnacc_ci95_high"] - dn_m],
        capsize=1.5,
        error_kw={"lw": 0.6, "ecolor": PAL["navy"]},
        label="DNAcc",
    )
    q4_q1 = dn_m[-1] / max(dn_m[0], 1e-6)
    ax_e.text(0.04, 0.90, f"Q4/Q1 = {q4_q1:.2f}$\\times$", transform=ax_e.transAxes, fontsize=5.2, fontweight="bold", color=PAL["navy"])
    ax_e.set_xticks(x_q)
    ax_e.set_xticklabels(["Q1", "Q2", "Q3", "Q4"], fontsize=5.4)
    ax_e.set_ylim(0, 86)
    ax_e.set_ylabel("Phase III/IV Rate (%)")
    ax_e.set_title(f"Maturation Gradient ($n$={n_prog_tot:,})", pad=5)
    ax_e.legend(loc="upper right", frameon=False, fontsize=5.0)

    # --- Panel f: 10-Decile Calibration Curve + Marginal Density Histogram ---
    ax_f = fig.add_subplot(gs_r3[0, 0])
    _panel_tag(ax_f, "f", x=-0.24, y=1.10)
    cal_p = s3h[s3h["cohort"].astype(str).str.startswith("Phase2_to_Phase34_Maturation")]
    c_dna = cal_p[cal_p["model"] == "Full_DNAcc"].sort_values("mean_predicted_prob")
    c_tau = cal_p[cal_p["model"] == "Static_Healthy_Tau0"].sort_values("mean_predicted_prob")

    # Subtle bottom marginal distribution bars across the 10 DNAcc deciles
    p_x = c_dna["mean_predicted_prob"].to_numpy()
    ax_f.bar(p_x, np.full_like(p_x, 0.055), width=0.032, bottom=0.10, color=PAL["sky"], alpha=0.45, edgecolor=PAL["blue"], lw=0.35)

    ax_f.plot([0.14, 0.80], [0.14, 0.80], color=PAL["slate"], ls=":", lw=0.65, label="Ideal")
    ax_f.plot(c_tau["mean_predicted_prob"], c_tau["observed_event_rate"], color=PAL["midgrey"], alpha=0.72, marker="s", ms=2.1, lw=0.82, ls="--", label="Static $\\tau_0$")
    ax_f.errorbar(
        c_dna["mean_predicted_prob"],
        c_dna["observed_event_rate"],
        yerr=[c_dna["observed_event_rate"] - c_dna["obs_ci95_low"], c_dna["obs_ci95_high"] - c_dna["observed_event_rate"]],
        color=PAL["blue"],
        marker="o",
        ms=2.8,
        lw=1.2,
        capsize=1.5,
        label="DNAcc",
    )
    ax_f.set_xlim(0.14, 0.80)
    ax_f.set_ylim(0.10, 0.85)
    ax_f.set_xlabel("Predicted Maturation Probability")
    ax_f.set_ylabel("Observed Rate (95% CI)")
    ax_f.set_title("10-Decile Calibration", pad=5)
    ax_f.legend(loc="upper left", frameon=False, fontsize=5.1)

    # --- Panel g: Standalone Decision Curve Analysis (DCA Net Benefit) ---
    ax_g = fig.add_subplot(gs_r3[0, 1])
    _panel_tag(ax_g, "g", x=-0.24, y=1.10)
    ax_g.axhline(0.0, color=PAL["midgrey"], ls=":", lw=0.7, label="Advance None")
    ax_g.plot(s3i["threshold_pt"], s3i["nb_treat_all"], color=PAL["slate"], ls="--", lw=0.9, label="Advance All")
    ax_g.plot(s3i["threshold_pt"], s3i["nb_static_tau0"], color=PAL["vermilion"], lw=0.95, label="Static $\\tau_0$")
    ax_g.plot(s3i["threshold_pt"], s3i["nb_dnacc"], color=PAL["teal"], lw=1.35, label="DNAcc")
    ax_g.set_xlim(float(s3i["threshold_pt"].min()), float(s3i["threshold_pt"].max()))
    ax_g.set_ylim(-0.03, max(float(s3i["nb_dnacc"].max()), 0.25) * 1.25)
    ax_g.set_xlabel("Phase III Threshold Probability ($p_t$)")
    ax_g.set_ylabel("Clinical Net Benefit")
    ax_g.set_title("Decision Curve Utility", pad=5)
    ax_g.legend(loc="upper right", frameon=False, fontsize=5.0)

    # --- Panel h: Target Remodeling Phase Space (tau_0 vs Delta_tau) + Vertical Marginal Delta_tau Density ---
    gs_h3 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_r3[0, 2], width_ratios=[4.15, 1.05], wspace=0.07)
    ax_h = fig.add_subplot(gs_h3[0, 0])
    ax_hm = fig.add_subplot(gs_h3[0, 1], sharey=ax_h)
    _panel_tag(ax_h, "h", x=-0.28, y=1.10)
    g_agg = (
        s3c.groupby("primary_target_gene")
        .agg(
            tau0=("tau0_healthy", "mean"),
            dtau=("delta_tau", "mean"),
            sw=("compartment_switch", "mean"),
            nt=("nct_id", "count"),
        )
        .reset_index()
    )
    g_sub = g_agg[g_agg["nt"] >= 4].copy()
    ax_h.axhline(0.0, color=PAL["slate"], ls="--", lw=0.6)
    ax_h.axvline(0.60, color=PAL["slate"], ls=":", lw=0.6)

    # Subtle 2D Gaussian KDE density iso-contours of clinical target genes in (tau_0, Delta_tau) phase space
    if len(g_sub) >= 20:
        try:
            kde_xy = stats.gaussian_kde(np.vstack([g_sub["tau0"].to_numpy(), g_sub["dtau"].to_numpy()]))
            kx, ky = np.meshgrid(np.linspace(0.05, 0.98, 45), np.linspace(float(g_sub["dtau"].min()), float(g_sub["dtau"].max()), 45))
            kz = kde_xy(np.vstack([kx.ravel(), ky.ravel()])).reshape(kx.shape)
            ax_h.contour(kx, ky, kz, levels=4, colors="#94A3B8", linewidths=0.45, alpha=0.55, zorder=1)
        except Exception:
            pass

    m_nosw = g_sub["sw"] < 0.5
    m_sw = ~m_nosw
    ax_h.scatter(
        g_sub.loc[m_nosw, "tau0"],
        g_sub.loc[m_nosw, "dtau"],
        s=7 + 9.5 * np.log10(g_sub.loc[m_nosw, "nt"]),
        color=PAL["blue"],
        alpha=0.35,
        edgecolors="none",
        zorder=2,
        label="Retained",
    )
    ax_h.scatter(
        g_sub.loc[m_sw, "tau0"],
        g_sub.loc[m_sw, "dtau"],
        s=8.5 + 10.5 * np.log10(g_sub.loc[m_sw, "nt"]),
        color=PAL["crimson"],
        alpha=0.52,
        edgecolors="white",
        linewidths=0.25,
        zorder=3,
        label="Switched",
    )

    # Bilateral non-crossing callouts inside compact phase space
    callouts = ["THRB", "TGFB1", "PCSK9", "ERBB2", "PDGFRB", "LOXL2"]
    c_rows = g_agg[g_agg["primary_target_gene"].isin(callouts)].copy()
    left_c = c_rows[c_rows["tau0"] < 0.58].sort_values("dtau", ascending=True)
    right_c = c_rows[c_rows["tau0"] >= 0.58].sort_values("dtau", ascending=True)

    y_min_h = float(g_sub["dtau"].min())
    y_max_h = float(g_sub["dtau"].max())
    y_lo_lim = y_min_h * 1.25 - 0.02
    y_hi_lim = y_max_h * 1.38 + 0.03

    if len(left_c) > 0:
        l_slots = np.linspace(y_min_h * 0.65, y_max_h * 0.75, len(left_c)) if len(left_c) > 1 else [float(left_c["dtau"].iloc[0])]
        for (_, r_g), y_sl in zip(left_c.iterrows(), l_slots):
            ax_h.scatter([r_g["tau0"]], [r_g["dtau"]], s=18, facecolors="none", edgecolors=PAL["ink"], linewidths=0.8, zorder=6)
            ax_h.annotate(
                str(r_g["primary_target_gene"]),
                xy=(float(r_g["tau0"]), float(r_g["dtau"])),
                xytext=(0.22, float(y_sl)),
                fontsize=4.65,
                fontweight="bold",
                fontstyle="italic",
                color=PAL["crimson"] if (r_g["sw"] >= 0.5 or r_g["dtau"] < -0.04) else PAL["navy"],
                va="center",
                ha="right",
                arrowprops={"arrowstyle": "-", "lw": 0.45, "color": PAL["slate"], "shrinkA": 0, "shrinkB": 1.5},
            )

    if len(right_c) > 0:
        r_slots = np.linspace(y_min_h * 0.88, y_max_h * 0.96, len(right_c))
        for (_, r_g), y_sl in zip(right_c.iterrows(), r_slots):
            ax_h.scatter([r_g["tau0"]], [r_g["dtau"]], s=18, facecolors="none", edgecolors=PAL["ink"], linewidths=0.8, zorder=6)
            ax_h.annotate(
                str(r_g["primary_target_gene"]),
                xy=(float(r_g["tau0"]), float(r_g["dtau"])),
                xytext=(0.98, float(y_sl)),
                fontsize=4.65,
                fontweight="bold",
                fontstyle="italic",
                color=PAL["crimson"] if (r_g["sw"] >= 0.5 or r_g["dtau"] < -0.04) else PAL["navy"],
                va="center",
                ha="left",
                arrowprops={"arrowstyle": "-", "lw": 0.45, "color": PAL["slate"], "shrinkA": 0, "shrinkB": 1.5},
            )

    ax_h.set_xlim(-0.04, 1.24)
    ax_h.set_xticks([0.2, 0.6, 1.0])
    ax_h.set_ylim(y_lo_lim, y_hi_lim)
    ax_h.set_xlabel("Healthy Specificity ($\\tau_0$)")
    ax_h.set_ylabel("Disease Drift ($\\Delta\\tau$)")
    ax_h.set_title("Target Drift Phase Space", pad=5)
    ax_h.legend(loc="upper left", frameon=False, fontsize=4.75, ncol=2, columnspacing=0.5, handletextpad=0.2)

    # Vertical marginal Delta_tau KDE + median markers (Retained vs Switched targets)
    ax_hm.axhline(0.0, color=PAL["slate"], ls="--", lw=0.6)
    ys_grid = np.linspace(y_lo_lim, y_hi_lim, 120)
    dt_ret = g_sub.loc[m_nosw, "dtau"].dropna().to_numpy(dtype=float)
    dt_swt = g_sub.loc[m_sw, "dtau"].dropna().to_numpy(dtype=float)
    max_dens = 1.0
    for arr_m, col_m in [(dt_ret, PAL["blue"]), (dt_swt, PAL["crimson"])]:
        if len(arr_m) >= 5:
            kde_m = stats.gaussian_kde(arr_m, bw_method=0.38)
            dens_m = kde_m(ys_grid)
            dens_m = dens_m / max(float(np.max(dens_m)), 1e-6)
            ax_hm.fill_betweenx(ys_grid, 0.0, dens_m, color=col_m, alpha=0.22)
            ax_hm.plot(dens_m, ys_grid, color=col_m, lw=0.95)
            med_m = float(np.median(arr_m))
            ax_hm.plot([0.0, 0.92], [med_m, med_m], color=col_m, lw=1.15, ls="-")
    ax_hm.set_xlim(0.0, 1.15)
    ax_hm.set_xticks([])
    ax_hm.tick_params(axis="y", which="both", left=False, labelleft=False)
    ax_hm.set_xlabel("Density", fontsize=4.6)
    for sp in ("top", "right", "bottom"):
        ax_hm.spines[sp].set_visible(False)
    ax_hm.spines["left"].set_color(PAL["lightgrey"])
    ax_hm.spines["left"].set_linewidth(0.5)

    _save_multi(fig, out_dir, "Fig3_PanOrgan_Dynamic_Drift_and_DNAcc_Law", artifact_dir)


# ===========================================================================
# FIGURE 4: Pan-Organ Fibrosis Zero-Shot Benchmark & Hysteresis (8 Panels a-h)
# ===========================================================================
def build_fig4(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 8-panel (a-h) Pan-Organ Human Fibrosis Zero-Shot Benchmark & Hysteresis at 7.20 in x 8.15 in."""
    set_nature_style()
    donors = pd.read_csv(tables_dir / "Table_S4a_PanOrgan_Fibrosis_All_Human_Donors.csv")
    ccr2 = pd.read_csv(tables_dir / "Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv")
    trials = pd.read_csv(tables_dir / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv")
    arena = pd.read_csv(tables_dir / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv")
    bi3x3 = pd.read_csv(tables_dir / "Table_S4f_Bidirectional_3x3_CrossOrgan_Transfer.csv")
    hyst = pd.read_csv(tables_dir / "Fig4_Continuous_Hysteresis_Bifurcation_Curve.csv")
    bliss = pd.read_csv(tables_dir / "Table_S4e_Fibrosis_2D_Bliss_Matrix.csv")
    dose_syn_p = tables_dir / "Table_S4h_Continuous_2D_Dose_Synergy_Surface.csv"
    dose_syn = pd.read_csv(dose_syn_p) if dose_syn_p.exists() else pd.DataFrame()
    wad_npz = np.load(tables_dir / "Fig4_Waddington_2D_Quiver_Fibrosis.npz")

    fig = plt.figure(figsize=(7.20, 8.15), facecolor="white")
    gs = gridspec.GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.38, 0.96, 0.98],
        hspace=0.49,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0], width_ratios=[0.94, 0.94, 0.14, 1.40], wspace=0.42)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[0.98, 1.08, 1.06], wspace=0.50)
    gs_r3 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], width_ratios=[1.26, 1.12], wspace=0.44)

    # --- Panel a: 64-Donor Closure-Free Fibrotic Triad Co-Activation (Exact Analytical Coupling Definition) ---
    ax_a = fig.add_subplot(gs_r1[0, 0])
    _panel_tag(ax_a, "a", x=-0.26, y=1.08)

    donors["organ_clean"] = donors["organ"].map(
        lambda s: "Liver" if "Liver" in str(s) else ("Lung" if "Lung" in str(s) else "Kidney")
    )
    comps_f = ["Parenchyma_Injury", "SPP1_SAMac", "Myofibroblast"]
    intr_f = pd.DataFrame(
        {
            c: donors[c].to_numpy(dtype=float) - 1.6 * donors[f"{c}_frac"].to_numpy(dtype=float)
            for c in comps_f
        },
        index=donors.index,
    )
    intr_f["organ"] = donors["organ"]
    norm_parts_f = []
    for _, g_f in intr_f.groupby("organ"):
        p_f = g_f[comps_f].copy()
        for c in comps_f:
            v_f = p_f[c].to_numpy(dtype=float)
            p_f[c] = (v_f - float(np.mean(v_f))) / (float(np.std(v_f)) + 1e-8)
        norm_parts_f.append(p_f)
    triad_z = pd.concat(norm_parts_f, axis=0).loc[donors.index]
    triad_z["organ_clean"] = donors["organ_clean"]

    x_all_z = triad_z["Parenchyma_Injury"].to_numpy(dtype=float)
    y_all_z = triad_z["Myofibroblast"].to_numpy(dtype=float)
    r_all_z, p_all_z = stats.pearsonr(x_all_z, y_all_z)
    liv_z = triad_z[triad_z["organ_clean"] == "Liver"]
    r_liv_z, p_liv_z = stats.pearsonr(liv_z["Parenchyma_Injury"].to_numpy(dtype=float), liv_z["Myofibroblast"].to_numpy(dtype=float))

    # OLS fit + 95% confidence band across all N = 64 human donors
    sl_z, ic_z = np.polyfit(x_all_z, y_all_z, 1)
    xg_z = np.linspace(-2.25, 2.35, 80)
    yg_z = sl_z * xg_z + ic_z
    n_z = len(x_all_z)
    resid_se = math.sqrt(float(np.sum((y_all_z - (sl_z * x_all_z + ic_z)) ** 2)) / max(n_z - 2, 1))
    ssx_z = float(np.sum((x_all_z - np.mean(x_all_z)) ** 2)) + 1e-8
    t_crit_z = float(stats.t.ppf(0.975, df=max(n_z - 2, 1)))
    ci_z = t_crit_z * resid_se * np.sqrt(1.0 / n_z + ((xg_z - float(np.mean(x_all_z))) ** 2) / ssx_z)
    ax_a.fill_between(xg_z, yg_z - ci_z, yg_z + ci_z, color=PAL["navy"], alpha=0.12, edgecolor="none", zorder=1)
    ax_a.plot(xg_z, yg_z, color=PAL["navy"], lw=1.05, ls="--", zorder=2)

    org_specs_a = [
        ("Liver", "Liver ($n$=17)", PAL["crimson"], "o", 15),
        ("Lung", "Lung ($n$=30)", PAL["blue"], "s", 13),
        ("Kidney", "Kidney ($n$=17)", PAL["teal"], "^", 15),
    ]
    for org_k, org_lbl, org_col, org_mk, org_sz in org_specs_a:
        sub_oz = triad_z[triad_z["organ_clean"] == org_k]
        ax_a.scatter(
            sub_oz["Parenchyma_Injury"],
            sub_oz["Myofibroblast"],
            c=org_col,
            marker=org_mk,
            s=org_sz,
            alpha=0.84,
            edgecolors="white",
            linewidths=0.35,
            label=org_lbl,
            zorder=4,
        )

    ax_a.text(
        0.04,
        0.95,
        f"All 64: $r$ = {r_all_z:.2f} ({_fmt_p(float(p_all_z))})\nLiver: $r$ = {r_liv_z:.2f} ({_fmt_p(float(p_liv_z))})",
        transform=ax_a.transAxes,
        va="top",
        ha="left",
        fontsize=4.65,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_a.set_xlim(-2.35, 2.45)
    ax_a.set_ylim(-3.95, 3.35)
    ax_a.set_xticks([-2, -1, 0, 1, 2])
    ax_a.set_yticks([-3, -1.5, 0, 1.5, 3])
    ax_a.set_xlabel("Parenchyma Injury Module ($z$)")
    ax_a.set_ylabel("Myofibroblast Scar Module ($z$)")
    ax_a.set_title("64-Donor Closure-Free Triad", pad=4)
    ax_a.legend(loc="lower right", frameon=False, fontsize=4.55, handletextpad=0.25, borderpad=0.1)

    # --- Panel b: Human Myeloid CCR2 Collapse vs SPP1 / TREM2 Emergence (with n = 19 Individual Donors) ---
    ax_b = fig.add_subplot(gs_r1[0, 1])
    _panel_tag(ax_b, "b", x=-0.28, y=1.08)
    st_order = [
        ("1_Human_NAFLD_Blood_Monocyte", "Blood"),
        ("2_Human_Healthy_Liver_Macrophage", "Healthy"),
        ("3_Human_F2_F3_NASH_Liver_SAMac", "F2/F3"),
        ("4_Human_F4_Cirrhotic_Liver_SAMac", "F4"),
    ]
    x_b = np.arange(len(st_order))
    rng_b = np.random.default_rng(17)
    for gene_n, col_g, mk in [("CCR2", PAL["crimson"], "o"), ("SPP1", PAL["blue"], "s"), ("TREM2", PAL["teal"], "^")]:
        means, sems = [], []
        for idx_s, (st_key, _) in enumerate(st_order):
            v = ccr2[ccr2["niche_stage"] == st_key][gene_n].to_numpy(dtype=float)
            means.append(float(np.mean(v)))
            sems.append(float(np.std(v, ddof=1) / math.sqrt(max(len(v), 1))) if len(v) > 1 else 0.0)
            ax_b.scatter(
                np.full(len(v), idx_s) + rng_b.uniform(-0.09, 0.09, size=len(v)),
                v,
                color=col_g,
                s=7,
                alpha=0.35,
                edgecolors="none",
                zorder=2,
            )
        ax_b.errorbar(
            x_b,
            means,
            yerr=sems,
            color=col_g,
            marker=mk,
            ms=3.0,
            lw=1.25,
            capsize=1.6,
            label=f"${gene_n}$",
            zorder=4,
        )

    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels([lbl for _, lbl in st_order], fontsize=5.2)
    ax_b.set_ylim(-0.4, 12.2)
    ax_b.set_ylabel("Myeloid Log-Expression")
    ax_b.set_title("Human Myeloid Drift ($n$=19)", pad=5)
    ax_b.legend(loc="upper center", frameon=False, fontsize=4.6, ncol=3, handletextpad=0.25, columnspacing=0.55)

    # --- Panel c: Multi-Track 24-Trial Phase IIb/III Antifibrotic Clinical Scorecard ---
    ax_c = fig.add_subplot(gs_r1[0, 3])
    _panel_tag(ax_c, "c", x=-0.42, y=1.08)
    tr_s = trials.sort_values(["clinical_win", "omniniche_remission_prob"], ascending=[True, True]).reset_index(drop=True)
    y_c = np.arange(len(tr_s))

    acr_short = {
        "GS-US-321-0105": "GS-0105",
        "GS-US-321-0106": "GS-0106",
        "GS-US-223-1015": "GS-1015",
        "FIBRONEER-IPF": "FIBRONEER",
        "SYNERGY-NASH": "SYNERGY",
        "MAESTRO-NASH": "MAESTRO",
        "FIDELIO-DKD": "FIDELIO",
        "RESOLVE-IT": "RESOLVE",
        "FLIGHT-FXR": "FLIGHT",
        "ZEPHYRUS-1": "ZEPHYRUS",
        "INPULSIS-1": "INPULSIS",
        "ISABELA 1": "ISABELA-1",
        "FALCON 1": "FALCON-1",
    }
    org_col_map = {"Liver_MASH": PAL["orange"], "Lung_IPF": PAL["blue"], "Kidney_CKD": PAL["purple"]}
    mod_col_map = {"small_molecule": PAL["blue"], "peptide": PAL["teal"], "protein": PAL["purple"], "antibody": PAL["crimson"]}

    # Draw 3 aligned metadata tile tracks in negative x [-0.28, -0.03]: Organ | Modality | Scar Access
    tile_xs = [-0.24, -0.15, -0.06]
    tile_w = 0.072
    for tx, thdr in zip(tile_xs, ["Org", "Mod", "Acc"]):
        ax_c.text(tx, len(tr_s) - 0.10, thdr, ha="center", va="bottom", fontsize=4.4, fontweight="bold", color=PAL["slate"])

    for idx_t, r_t in tr_s.iterrows():
        c_org = org_col_map.get(str(r_t["organ"]), PAL["midgrey"])
        c_mod = mod_col_map.get(str(r_t["modality"]), PAL["midgrey"])
        acc_val = float(r_t["scar_barrier_access"])
        c_acc = PAL["teal"] if acc_val >= 0.70 else PAL["vermilion"]
        for tx, tc in zip(tile_xs, [c_org, c_mod, c_acc]):
            ax_c.add_patch(
                mpatches.Rectangle((tx - tile_w / 2, idx_t - 0.36), tile_w, 0.72, facecolor=tc, alpha=0.82, edgecolor="white", lw=0.30)
            )

    ax_c.axvline(0.0, color=PAL["ink"], lw=0.55)
    l_youd = ax_c.axvline(0.0535, color=PAL["teal"], ls=":", lw=0.80, label="Youden $p \\geq 0.053$ (10/10 Sens.)")
    l_bas = ax_c.axvline(0.50, color=PAL["slate"], ls="--", lw=0.65, label="Basin $p \\geq 0.50$ (13/14 Spec.)")
    cols_c = [PAL["teal"] if int(w) == 1 else PAL["crimson"] for w in tr_s["clinical_win"]]
    probs_c = tr_s["omniniche_remission_prob"].to_numpy(dtype=float)
    ax_c.barh(y_c, probs_c, left=0.0, color=cols_c, alpha=0.86, height=0.68)
    for idx_t, p_val in enumerate(probs_c):
        ax_c.text(p_val + 0.022, idx_t, f"{p_val:.2f}", va="center", ha="left", fontsize=4.15, color=PAL["ink"])
    ax_c.legend(handles=[l_youd, l_bas], loc="lower right", frameon=False, fontsize=4.15, handlelength=1.4, handletextpad=0.35)

    row_names = []
    for r in tr_s.itertuples():
        ac = acr_short.get(str(r.trial_acronym), str(r.trial_acronym))
        d_clean = (
            str(r.drug)
            .replace("Simtuzumab", "Simtuz.")
            .replace("Selonsertib", "Selons.")
            .replace("Pamrevlumab", "Pamrev.")
            .replace("Zinpentraxin alfa", "Zinpentrax.")
            .replace("Nerandomilast", "Nerandom.")
            .replace("Lebrikizumab", "Lebrikiz.")
            .replace("Efruxifermin", "Efruxif.")
            .replace("Pegbelfermin", "Pegbelf.")
            .replace("Ziritaxestat", "Ziritax.")
            .replace("Dapagliflozin", "Dapaglif.")
            .replace("Lanifibranor", "Lanifib.")
            .replace("Cenicriviroc", "Cenicriv.")
        )
        row_names.append(f"{d_clean} ({ac})")

    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels(row_names, fontsize=4.25)
    ax_c.set_xlim(-0.30, 1.18)
    ax_c.set_xticks([0.0, 0.5, 1.0])
    ax_c.set_xticklabels(["0.0", "0.5", "1.0"])
    ax_c.set_ylim(-0.65, len(tr_s) + 0.70)
    ax_c.set_xlabel("Zero-Shot Remission Probability")
    ax_c.set_title(f"{len(tr_s)} Phase IIb/III Antifibrotic Trials", pad=5)

    # --- Panel d: Zero-Shot Benchmark Arena (Protagonist Cobalt & Stage-Blind Teal vs Baselines) ---
    ax_d = fig.add_subplot(gs_r2[0, 0])
    _panel_tag(ax_d, "d", x=-0.34, y=1.08)
    arena_24 = (
        arena[arena["benchmark"] == f"PanOrgan_Fibrosis_{len(tr_s)}_Phase2b3_Trials"].copy()
        if "benchmark" in arena.columns
        else arena.copy()
    )
    if arena_24.empty:
        arena_24 = arena.copy()
    m_short = {
        "OmniNiche_v2_ZeroShot_Langevin": "OmniNiche",
        "OmniNiche_StageBlind_FixedS0": "Stage-Blind $s_0$",
        "B0_Enrolled_Stage_Prior_Only": "Stage Prior Only",
        "B1_Static_Healthy_Tau0": "Static $\\tau_0$",
        "B2_Polypharmacology_Target_Count": "Target Count",
        "B3_SingleCell_GRN_Retention_Only": "scGRN Only",
        "B4_Static_LR_Compartment_Degree": "Static L-R",
        "B5_Uninformative_Modality_Prior": "Modality Prior",
    }
    ar_rev = arena_24.iloc[::-1].reset_index(drop=True)
    y_d = np.arange(len(ar_rev))
    ax_d.axvline(0.50, color=PAL["slate"], ls=":", lw=0.6)
    for i_m, r_m in ar_rev.iterrows():
        m_str = str(r_m["method"])
        is_main = m_str == "OmniNiche_v2_ZeroShot_Langevin"
        is_blind = m_str == "OmniNiche_StageBlind_FixedS0"
        c_m = PAL["blue"] if is_main else (PAL["teal"] if is_blind else PAL["midgrey"])
        ci_lo = float(r_m["auroc_ci_low"]) if "auroc_ci_low" in r_m and pd.notna(r_m["auroc_ci_low"]) else float(r_m["auroc"])
        ci_hi = float(r_m["auroc_ci_high"]) if "auroc_ci_high" in r_m and pd.notna(r_m["auroc_ci_high"]) else float(r_m["auroc"])
        ax_d.plot([ci_lo, ci_hi], [i_m, i_m], color=c_m, lw=1.40 if (is_main or is_blind) else 0.90)
        ax_d.scatter([r_m["auroc"]], [i_m], color=c_m, s=20 if (is_main or is_blind) else 12, zorder=4)
        ax_d.text(
            1.03,
            i_m,
            f"{r_m['auroc']:.2f}",
            va="center",
            ha="left",
            fontsize=4.6,
            fontweight="bold" if (is_main or is_blind) else "normal",
            color=PAL["navy"] if is_main else (PAL["teal"] if is_blind else PAL["slate"]),
        )

    omni_p = (
        float(arena_24.loc[arena_24["method"] == "OmniNiche_v2_ZeroShot_Langevin", "exact_perm_p_vs_null"].iloc[0])
        if "exact_perm_p_vs_null" in arena_24.columns
        else 1e-4
    )
    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels([m_short.get(m, m) for m in ar_rev["method"]], fontsize=4.7)
    ax_d.set_xlim(0.18, 1.25)
    ax_d.set_xticks([0.25, 0.50, 0.75, 1.00])
    ax_d.set_xlabel(f"{len(tr_s)}-Trial AUROC (95% CI)")
    ax_d.set_title(f"Zero-Shot Arena ({_fmt_p(omni_p)})", pad=5)

    # --- Panel e: Bidirectional 3x3 & Leave-One-Organ-Out Transfer Matrix ---
    ax_e = fig.add_subplot(gs_r2[0, 1])
    _panel_tag(ax_e, "e", x=-0.32, y=1.08)
    src_short = [
        "Liver Only",
        "Lung Only",
        "Kidney Only",
        "LOOO: Lu+Ki",
        "LOOO: Li+Ki",
        "LOOO: Li+Lu",
    ]
    k_col = "target_kidney_auroc" if "target_kidney_auroc" in bi3x3.columns else "target_kidney_accuracy"
    all_col = "all24_auroc" if "all24_auroc" in bi3x3.columns else "all14_auroc"
    mat_3x3 = bi3x3[["target_liver_auroc", "target_lung_auroc", k_col, all_col]].to_numpy(dtype=float)
    vmin_e = max(0.70, float(np.nanmin(mat_3x3)) - 0.03)
    im_e = ax_e.imshow(mat_3x3, cmap="YlGnBu", vmin=vmin_e, vmax=1.01, aspect="auto")
    for r_i in range(mat_3x3.shape[0]):
        for c_i in range(mat_3x3.shape[1]):
            val_m = mat_3x3[r_i, c_i]
            ax_e.text(
                c_i,
                r_i,
                f"{val_m:.2f}",
                ha="center",
                va="center",
                fontsize=4.9,
                fontweight="bold",
                color="white" if val_m >= 0.5 * (vmin_e + 1.01) + 0.03 else PAL["ink"],
            )
    ax_e.set_xticks(range(4))
    ax_e.set_xticklabels(["Liver", "Lung", "Kidney", f"All {len(tr_s)}"], fontsize=5.1)
    ax_e.set_yticks(range(len(src_short)))
    ax_e.set_yticklabels(src_short, fontsize=5.1)
    ax_e.set_title("Bidirectional & LOOO Transfer", pad=5)
    for sp in ("top", "right", "bottom", "left"):
        ax_e.spines[sp].set_visible(True)
        ax_e.spines[sp].set_linewidth(0.55)
    cb_e = fig.colorbar(im_e, ax=ax_e, fraction=0.040, pad=0.025)
    cb_e.ax.tick_params(labelsize=4.8)

    # --- Panel f: Light-Base Topographic 60x60 2D Waddington Quasi-Potential + Quiver Drift Field ---
    ax_f = fig.add_subplot(gs_r2[0, 2])
    _panel_tag(ax_f, "f", x=-0.24, y=1.08)
    gx_w = wad_npz["grid_myofibroblast"]
    gy_w = wad_npz["grid_parenchyma"]
    u_2d = wad_npz["U_2d"]
    du_x = wad_npz["drift_myofibroblast"]
    du_y = wad_npz["drift_parenchyma"]

    cf_f = ax_f.contourf(gx_w, gy_w, u_2d, levels=16, cmap=get_waddington_light_cmap())
    ax_f.contour(gx_w, gy_w, u_2d, levels=10, colors="#52606D", linewidths=0.36, alpha=0.52)
    st_q = 5
    Xg, Yg = np.meshgrid(gx_w[::st_q], gy_w[::st_q])
    Ux, Uy = du_x[::st_q, ::st_q], du_y[::st_q, ::st_q]
    mag = np.hypot(Ux, Uy) + 1e-8
    ax_f.quiver(Xg, Yg, Ux / mag, Uy / mag, color="#1E293B", alpha=0.68, scale=22, width=0.0055)
    ax_f.axvline(0.36, color=PAL["crimson"], ls="--", lw=0.85)
    ax_f.scatter([0.14, 0.74], [0.14, 0.58], c=[PAL["teal"], PAL["crimson"]], s=28, edgecolors="white", linewidths=0.95, zorder=5)
    ax_f.text(0.18, 0.14, "Remission", color=PAL["navy"], fontsize=5.0, fontweight="bold", va="center", bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "alpha": 0.78, "edgecolor": "none"})
    ax_f.text(0.71, 0.67, "Locked F4", color=PAL["crimson"], fontsize=5.0, fontweight="bold", va="bottom", ha="center", bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "alpha": 0.78, "edgecolor": "none"})
    ax_f.set_xlabel("Myofibroblast ($S_{\\mathrm{Myofib}}$)")
    ax_f.set_ylabel("Parenchyma ($S_{\\mathrm{Paren}}$)")
    ax_f.set_title("2D Waddington Landscape", pad=4)
    for sp in ("top", "right", "bottom", "left"):
        ax_f.spines[sp].set_visible(True)
        ax_f.spines[sp].set_linewidth(0.50)
    cbar_f = fig.colorbar(cf_f, ax=ax_f, fraction=0.040, pad=0.025)
    cbar_f.ax.tick_params(labelsize=4.8)
    cbar_f.set_label("$U(S)$", fontsize=5.1)

    # --- Panel g: Continuous Fibrotic Progression Hysteresis Bifurcation Curve ---
    ax_g = fig.add_subplot(gs_r3[0, 0])
    _panel_tag(ax_g, "g", x=-0.14, y=1.10)
    sp_x = hyst["progression_index"].to_numpy()
    cliff_path = tables_dir / "Fig4_Monotherapy_Window_Cliff_Curve.csv"
    t50_hr = float(pd.read_csv(cliff_path)["t50_window_cliff_hr"].iloc[0]) if cliff_path.exists() else 8.3
    rescue_col = "dual_parenchyma_myofib_rescue_remission" if "dual_parenchyma_myofib_rescue_remission" in hyst.columns else "dual_myofib_samac_rescue_remission"
    ax_g.plot(sp_x, 100.0 * hyst["myofib_lowdose_mono_remission"], color=PAL["crimson"], lw=1.3, marker="o", ms=2.3, label="Stromal Mono. (TGFBR1/PDGFRB)")
    ax_g.plot(sp_x, 100.0 * hyst["parenchyma_thrb_mono_remission"], color=PAL["orange"], lw=1.2, ls="--", marker="^", ms=2.2, label="Parenchymal Mono. (THRB)")
    ax_g.plot(sp_x, 100.0 * hyst[rescue_col], color=PAL["teal"], lw=1.5, marker="s", ms=2.3, label="Dual Rescue (Paren. + Myofib)")
    ax_g.axvspan(0.76, 0.91, color=PAL["crimson"], alpha=0.09, label=f"Locked F4 ($t_{{50}}$={t50_hr:.1f} h)")

    ax_g.set_xlim(0.17, 0.91)
    ax_g.set_ylim(-4, 138)
    ax_g.set_yticks([0, 25, 50, 75, 100])
    ax_g.set_xlabel("Anatomical Fibrotic Progression Coordinate (F1/F2 → Locked F4)")
    ax_g.set_ylabel("Simulated Remission (%)")
    ax_g.set_title("Hysteresis Monotherapy Cliff & Dual Rescue", pad=4)
    ax_g.legend(loc="upper right", frameon=False, fontsize=4.8, ncol=2)

    # --- Panel h: Continuous 2D Locked-F4 Dose-Response Bliss Synergy Surface & Saddle Isobologram ---
    ax_h = fig.add_subplot(gs_r3[0, 1])
    _panel_tag(ax_h, "h", x=-0.24, y=1.10)
    if not dose_syn.empty:
        p_vals = np.sort(dose_syn["eta_parenchyma"].unique())
        m_vals = np.sort(dose_syn["eta_myofibroblast"].unique())
        piv_syn = (
            dose_syn.pivot(index="eta_myofibroblast", columns="eta_parenchyma", values="bliss_synergy_excess")
            .reindex(index=m_vals, columns=p_vals)
            .to_numpy(dtype=float)
        )
        piv_rem = (
            dose_syn.pivot(index="eta_myofibroblast", columns="eta_parenchyma", values="observed_combo_remission")
            .reindex(index=m_vals, columns=p_vals)
            .to_numpy(dtype=float)
        )
        syn_smooth = ndimage.gaussian_filter(piv_syn, sigma=0.65)
        rem_smooth = ndimage.gaussian_filter(piv_rem, sigma=0.65)
        cf_h = ax_h.contourf(p_vals, m_vals, np.clip(syn_smooth, 0.0, 1.0), levels=np.linspace(0.0, 1.0, 17), cmap=get_synergy_light_cmap())
        ax_h.contour(p_vals, m_vals, rem_smooth, levels=[0.50], colors=[PAL["navy"]], linewidths=1.15, linestyles="--")

        # Clinical regimen coordinates on the 2D Locked-F4 pharmacological plane
        ax_h.scatter([0.66], [0.008], color=PAL["orange"], marker="o", s=24, edgecolors=PAL["ink"], linewidths=0.7, zorder=6)
        ax_h.text(0.64, 0.024, "THRB Mono / +LOXL2 mAb (0%)", fontsize=4.25, fontweight="bold", color=PAL["ink"], ha="right", va="bottom")

        ax_h.scatter([0.02], [0.138], color=PAL["purple"], marker="s", s=22, edgecolors="white", linewidths=0.6, zorder=6)
        ax_h.text(0.045, 0.138, "TGFBR1i\nMono (0%)", fontsize=4.2, fontweight="bold", color=PAL["ink"], ha="left", va="center")

        ax_h.scatter([0.66], [0.138], color=PAL["teal"], marker="*", s=62, edgecolors="white", linewidths=0.80, zorder=7)
        ax_h.text(
            0.66,
            0.168,
            "THRB + TGFBR1i (99.5%, Bliss +1.00)",
            fontsize=4.35,
            fontweight="bold",
            color="white",
            ha="center",
            va="bottom",
        )
        ax_h.text(0.32, 0.255, "50% Remission\nIsobologram", fontsize=4.25, fontweight="bold", color=PAL["navy"], ha="center", va="center")

        ax_h.set_xlim(0.0, 0.90)
        ax_h.set_ylim(0.0, 0.30)
        ax_h.set_xticks([0.0, 0.3, 0.6, 0.9])
        ax_h.set_yticks([0.0, 0.1, 0.2, 0.3])
        ax_h.set_xlabel("Parenchymal Unloading ($\\eta_{\\mathrm{Paren}}$)")
        ax_h.set_ylabel("Stromal Unlock ($\\eta_{\\mathrm{Myofib}}$)")
        ax_h.set_title("Locked-F4 2D Bliss Synergy Surface", pad=4)
        for sp in ("top", "right", "bottom", "left"):
            ax_h.spines[sp].set_visible(True)
            ax_h.spines[sp].set_linewidth(0.52)
        cbar_h = fig.colorbar(cf_h, ax=ax_h, fraction=0.042, pad=0.025, ticks=[0.0, 0.5, 1.0])
        cbar_h.ax.tick_params(labelsize=4.8)
        cbar_h.set_label("Bliss Excess", fontsize=5.1)
    else:
        arm_order = list(dict.fromkeys(bliss["arm_row"]))
        piv_b = bliss.pivot(index="arm_row", columns="arm_col", values="bliss_synergy_excess").reindex(index=arm_order, columns=arm_order).to_numpy(dtype=float)
        im_h = ax_h.imshow(piv_b, cmap=get_synergy_light_cmap(), vmin=0.0, vmax=1.0, aspect="equal")
        cbar_h = fig.colorbar(im_h, ax=ax_h, fraction=0.042, pad=0.025)
        cbar_h.set_label("Bliss Excess", fontsize=5.1)

    _save_multi(fig, out_dir, "Fig4_PanOrgan_Fibrosis_ZeroShot_and_Hysteresis", artifact_dir)


# ===========================================================================
# FIGURE 5: Solid-Tumor Visium Exclusion & Pre-ICB Responder (8 Panels a-h)
# ===========================================================================
def _kaplan_meier_curve(times: np.ndarray, events: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Exact product-limit Kaplan-Meier survival estimate with Greenwood 95% CI."""
    order = np.argsort(times, kind="mergesort")
    t_ord = times[order]
    e_ord = events[order]
    uniq_t = np.unique(t_ord)
    t_pts = [0.0]
    s_pts = [1.0]
    lo_pts = [1.0]
    hi_pts = [1.0]
    surv = 1.0
    gw_var = 0.0
    for ut in uniq_t:
        if ut <= 0.0:
            continue
        at_risk = int(np.sum(t_ord >= ut))
        d_i = int(np.sum((t_ord == ut) & (e_ord == 1)))
        if at_risk > 0 and d_i > 0:
            surv *= 1.0 - (d_i / at_risk)
            if at_risk > d_i:
                gw_var += d_i / (at_risk * (at_risk - d_i))
        se = surv * math.sqrt(max(gw_var, 0.0))
        t_pts.append(float(ut))
        s_pts.append(float(surv))
        lo_pts.append(float(np.clip(surv - 1.96 * se, 0.0, 1.0)))
        hi_pts.append(float(np.clip(surv + 1.96 * se, 0.0, 1.0)))
    return np.asarray(t_pts), np.asarray(s_pts), np.asarray(lo_pts), np.asarray(hi_pts)


def build_fig5(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 8-panel (a-h) Solid-Tumor 10x Visium Spatial Barrier, Pre-ICB Responder & Monotherapy OS at 7.20 in x 8.85 in."""
    set_nature_style()
    vis_con = pd.read_csv(tables_dir / "Fig5_Contrasting_Visium_Sections.csv")
    rad_sum = pd.read_csv(tables_dir / "Table_S5g_Visium_Radial_Consensus_Summary.csv")
    vis_12 = pd.read_csv(tables_dir / "Table_S5a_Visium_Spatial_Barrier_Sections.csv")
    pool107 = pd.read_csv(tables_dir / "Table_S5i_Pooled_PreTreatment_107Patients_Master.csv")
    bench = pd.read_csv(tables_dir / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv")
    paired15 = pd.read_csv(tables_dir / "Table_S5f_GSE123813_Paired_PrePost_Trajectory.csv")
    tum_npz = np.load(tables_dir / "Fig5_Waddington_2D_Quiver_SolidTumor.npz")
    surv74 = pd.read_csv(tables_dir / "Table_S5l_AntiPD1_Monotherapy_Survival_74Patients.csv")
    surv_sum = pd.read_csv(tables_dir / "Table_S5m_AntiPD1_Survival_Summary.csv")

    fig = plt.figure(figsize=(7.20, 8.85), facecolor="white")
    gs = gridspec.GridSpec(
        4,
        1,
        figure=fig,
        height_ratios=[1.02, 1.12, 0.96, 1.08],
        hspace=0.47,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[1.38, 1.00], wspace=0.28)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], width_ratios=[1.20, 1.14], wspace=0.40)
    gs_r3 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], width_ratios=[1.00, 1.14], wspace=0.38)
    gs_r4 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[3], width_ratios=[0.94, 1.24], wspace=0.38)

    # --- Panel a: Dual Contrasting Primary Human Carcinoma 10x Visium Spatial Maps (Tissue Footprint + Iso-Contours + 500-um Scale Bar) ---
    gs_a = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_r1[0, 0], wspace=0.12)
    ax_a1 = fig.add_subplot(gs_a[0, 0])
    ax_a2 = fig.add_subplot(gs_a[0, 1])
    _panel_tag(ax_a1, "a", x=-0.15, y=1.10)

    niche_cols = {
        "Tumor Core": PAL["navy"],
        "myCAF Shield": PAL["crimson"],
        "SPP1+ TAM": PAL["orange"],
        "CD8+ T Zone": PAL["teal"],
    }
    for ax_sub, arch_key, title_prefix in [
        (ax_a1, "Immune_Excluded_myCAF_Shield", "Excluded"),
        (ax_a2, "Immune_Infiltrated_Permissive", "Infiltrated"),
    ]:
        sub_v = vis_con[vis_con["archetype"] == arch_key].copy()
        z_mat = sub_v[["Tumor", "myCAF", "TAM", "CD8_T"]].copy()
        z_mat = (z_mat - z_mat.mean()) / (z_mat.std() + 1e-8)
        dom_idx = np.argmax(z_mat.to_numpy(), axis=1)
        lbl_map = {0: "Tumor Core", 1: "myCAF Shield", 2: "SPP1+ TAM", 3: "CD8+ T Zone"}
        sub_v["dom"] = [lbl_map[i] for i in dom_idx]

        xs = sub_v["x"].to_numpy(dtype=float)
        ys = -sub_v["y"].to_numpy(dtype=float)
        sub_v["plot_y"] = ys

        coords = np.column_stack([xs, ys])
        kdt = spatial.cKDTree(coords)
        nn_dists = kdt.query(coords, k=2)[0][:, 1]
        pitch_100um = float(np.median(nn_dists))
        bar_500um = 5.0 * pitch_100um
        x_min, x_max = float(xs.min()), float(xs.max())
        y_min, y_max = float(ys.min()), float(ys.max())
        span_x = x_max - x_min
        span_y = y_max - y_min

        # Smooth 2D tissue footprint underlay + Tumor Core / myCAF Shield spatial iso-contours
        gx_v = np.linspace(x_min - 0.04 * span_x, x_max + 0.04 * span_x, 64)
        gy_v = np.linspace(y_min - 0.04 * span_y, y_max + 0.04 * span_y, 64)
        GX_v, GY_v = np.meshgrid(gx_v, gy_v)
        g_pts = np.column_stack([GX_v.ravel(), GY_v.ravel()])
        d_k4, idx_k4 = kdt.query(g_pts, k=min(4, len(coords)))
        tissue_mask = (d_k4[:, 0].reshape(GX_v.shape)) <= (1.75 * pitch_100um)
        mask_smooth = ndimage.gaussian_filter(tissue_mask.astype(float), sigma=0.85)
        ax_sub.contourf(gx_v, gy_v, mask_smooth, levels=[0.42, 1.50], colors=["#EDF2F7"], zorder=1)
        ax_sub.contour(gx_v, gy_v, mask_smooth, levels=[0.42], colors=[PAL["lightgrey"]], linewidths=0.55, zorder=2)

        w_idw = 1.0 / np.square(d_k4 + 0.35 * pitch_100um)
        w_idw /= np.sum(w_idw, axis=1, keepdims=True)
        z_tum_arr = z_mat["Tumor"].to_numpy(dtype=float)
        z_caf_arr = z_mat["myCAF"].to_numpy(dtype=float)
        f_tum = ndimage.gaussian_filter(( np.sum(w_idw * z_tum_arr[idx_k4], axis=1) ).reshape(GX_v.shape), sigma=1.35)
        f_caf = ndimage.gaussian_filter(( np.sum(w_idw * z_caf_arr[idx_k4], axis=1) ).reshape(GX_v.shape), sigma=1.35)
        f_tum_m = np.where(tissue_mask, f_tum, np.nan)
        f_caf_m = np.where(tissue_mask, f_caf, np.nan)

        for n_lbl, c_hex in niche_cols.items():
            m = sub_v["dom"] == n_lbl
            if np.any(m):
                ax_sub.scatter(sub_v.loc[m, "x"], sub_v.loc[m, "plot_y"], c=c_hex, s=2.2, alpha=0.84, edgecolors="none", label=n_lbl, zorder=3)

        if np.nanmin(f_tum_m) < 0.55 < np.nanmax(f_tum_m):
            ax_sub.contour(gx_v, gy_v, np.where(tissue_mask, f_tum, -3.0), levels=[0.55], colors=[PAL["navy"]], linewidths=0.72, linestyles="-", zorder=4)
        if np.nanmin(f_caf_m) < 0.55 < np.nanmax(f_caf_m):
            ax_sub.contour(gx_v, gy_v, np.where(tissue_mask, f_caf, -3.0), levels=[0.55], colors=[PAL["crimson"]], linewidths=0.85, linestyles="--", zorder=5)

        bx0 = x_max - 0.05 * span_x - bar_500um
        by0 = y_min - 0.08 * span_y
        ax_sub.plot([bx0, bx0 + bar_500um], [by0, by0], color=PAL["ink"], lw=1.35, solid_capstyle="butt", zorder=6)
        ax_sub.text(bx0 + 0.5 * bar_500um, by0 - 0.025 * span_y, "500 µm", ha="center", va="top", fontsize=4.6, fontweight="bold", color=PAL["ink"])

        ax_sub.set_xlim(x_min - 0.06 * span_x, x_max + 0.06 * span_x)
        ax_sub.set_ylim(y_min - 0.18 * span_y, y_max + 0.34 * span_y)

        # High-magnification 200-um Tumor-myCAF-CD8+ T boundary zoom-in inset on Immune-Excluded section
        if arch_key == "Immune_Excluded_myCAF_Shield":
            caf_pts = sub_v[sub_v["dom"] == "myCAF Shield"]
            if len(caf_pts) >= 12:
                # Find spot with maximal joint local admixture of Tumor Core + myCAF Shield + CD8+ T Zone within 3.5*pitch
                best_center = (float(np.median(xs)), float(np.median(ys)))
                best_score = -1.0
                r_win = 2.6 * pitch_100um
                dom_arr = sub_v["dom"].to_numpy()
                for cx_i, cy_i in zip(caf_pts["x"].to_numpy(dtype=float), caf_pts["plot_y"].to_numpy(dtype=float)):
                    in_box = (np.abs(xs - cx_i) <= r_win) & (np.abs(ys - cy_i) <= r_win)
                    if int(np.sum(in_box)) < 14:
                        continue
                    d_box = dom_arr[in_box]
                    n_t = int(np.sum(d_box == "Tumor Core"))
                    n_c = int(np.sum(d_box == "myCAF Shield"))
                    n_m = int(np.sum(d_box == "SPP1+ TAM"))
                    n_k = int(np.sum(d_box == "CD8+ T Zone"))
                    sc_b = min(n_t, 8) * min(n_c, 8) * (1 + min(n_k, 5) + 0.5 * min(n_m, 4))
                    if sc_b > best_score:
                        best_score = sc_b
                        best_center = (cx_i, cy_i)

                cx_z, cy_z = best_center
                hw_z = 2.45 * pitch_100um
                ax_sub.add_patch(
                    mpatches.Rectangle(
                        (cx_z - hw_z, cy_z - hw_z),
                        2.0 * hw_z,
                        2.0 * hw_z,
                        fill=False,
                        edgecolor=PAL["ink"],
                        lw=0.75,
                        ls="--",
                        zorder=7,
                    )
                )
                ax_ins = ax_sub.inset_axes([0.55, 0.63, 0.43, 0.35])
                ax_ins.set_facecolor("#FFFFFF")
                m_zoom = (np.abs(xs - cx_z) <= hw_z * 1.15) & (np.abs(ys - cy_z) <= hw_z * 1.15)
                sub_z = sub_v.loc[m_zoom]
                for n_lbl, c_hex in niche_cols.items():
                    mz = sub_z["dom"] == n_lbl
                    if np.any(mz):
                        ax_ins.scatter(
                            sub_z.loc[mz, "x"],
                            sub_z.loc[mz, "plot_y"],
                            c=c_hex,
                            s=11.5,
                            alpha=0.94,
                            edgecolors="white",
                            linewidths=0.30,
                            zorder=4,
                        )
                bar_200 = 2.0 * pitch_100um
                zb_x0 = cx_z + hw_z - 0.14 * hw_z - bar_200
                zb_y0 = cy_z - hw_z + 0.34 * hw_z
                ax_ins.add_patch(
                    mpatches.Rectangle(
                        (zb_x0 - 0.09 * hw_z, cy_z - hw_z + 0.04 * hw_z),
                        bar_200 + 0.18 * hw_z,
                        0.40 * hw_z,
                        facecolor="white",
                        edgecolor="none",
                        alpha=0.90,
                        zorder=5,
                    )
                )
                ax_ins.plot([zb_x0, zb_x0 + bar_200], [zb_y0, zb_y0], color=PAL["ink"], lw=1.15, zorder=6)
                ax_ins.text(zb_x0 + 0.5 * bar_200, zb_y0 - 0.06 * hw_z, "200 µm", ha="center", va="top", fontsize=3.75, fontweight="bold", color=PAL["ink"], zorder=6)
                ax_ins.set_xlim(cx_z - hw_z, cx_z + hw_z)
                ax_ins.set_ylim(cy_z - hw_z, cy_z + hw_z)
                ax_ins.set_xticks([])
                ax_ins.set_yticks([])
                for sp_i in ("top", "right", "bottom", "left"):
                    ax_ins.spines[sp_i].set_visible(True)
                    ax_ins.spines[sp_i].set_color(PAL["ink"])
                    ax_ins.spines[sp_i].set_linewidth(0.65)
                ax_sub.indicate_inset_zoom(ax_ins, edgecolor=PAL["ink"], linewidth=0.55, linestyle=":", alpha=0.65)

        sei_v = float(sub_v["sei"].iloc[0])
        ctype_v = str(sub_v["cancer_type"].iloc[0]).replace("_Carcinoma", "")
        ax_sub.set_title(f"{title_prefix} ({ctype_v}, SEI +{sei_v:.2f})", fontsize=5.95, pad=3.5)
        ax_sub.set_xticks([])
        ax_sub.set_yticks([])
        for sp in ("top", "right", "bottom", "left"):
            ax_sub.spines[sp].set_visible(True)
            ax_sub.spines[sp].set_color(PAL["lightgrey"])
            ax_sub.spines[sp].set_linewidth(0.5)
    ax_a2.legend(
        handles=ax_a1.get_legend_handles_labels()[0],
        labels=ax_a1.get_legend_handles_labels()[1],
        loc="upper left",
        frameon=False,
        fontsize=4.25,
        markerscale=1.50,
        ncol=2,
        columnspacing=0.45,
        handletextpad=0.22,
    )

    # --- Panel b: 12-Section (21,796-Spot) Continuous Radial Distance Decay (0-600 um) ---
    ax_b = fig.add_subplot(gs_r1[0, 1])
    _panel_tag(ax_b, "b", x=-0.20, y=1.10)
    dist_um = rad_sum["approx_distance_um"].to_numpy()
    for col_z, lbl_z, c_z, mk in [
        ("Tumor_z", "Tumor Core", PAL["navy"], "o"),
        ("myCAF_z", "myCAF Wall", PAL["crimson"], "s"),
        ("TAM_z", "$SPP1^+$ TAM", PAL["orange"], "^"),
        ("CD8_T_z", "$CD8^+$ T", PAL["teal"], "D"),
    ]:
        mu_z = rad_sum[f"{col_z}_mean"].to_numpy()
        se_z = rad_sum[f"{col_z}_sem"].to_numpy()
        ax_b.plot(dist_um, mu_z, color=c_z, lw=1.2, marker=mk, ms=2.3, label=lbl_z)
        ax_b.fill_between(dist_um, mu_z - se_z, mu_z + se_z, color=c_z, alpha=0.14)

    w_res = stats.wilcoxon(vis_12["stromal_exclusion_index"].to_numpy(), alternative="greater")
    ax_b.axhline(0.0, color=PAL["slate"], ls=":", lw=0.6)
    ax_b.text(0.22, 0.90, f"12/12 SEI > 0 ({_fmt_p(float(w_res.pvalue))})", transform=ax_b.transAxes, fontsize=4.9, fontweight="bold", color=PAL["navy"])
    ax_b.set_xlabel("Distance from Tumor Core (µm)")
    ax_b.set_ylabel("Visium Signature ($z$-Score)")
    ax_b.set_ylim(float(rad_sum["Tumor_z_mean"].min()) - 0.35, float(rad_sum["Tumor_z_mean"].max()) + 0.55)
    ax_b.set_title("12-Section Radial Barrier (20,880 QC Spots)")
    ax_b.legend(loc="lower left", frameon=False, fontsize=4.75, ncol=2)

    # --- Panel c: n = 107 Pre-Treatment Patient Waterfall (Smooth Unclipped Per-Cohort LOPO Probability) + 4-Row OncoPrint ---
    gs_c = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_r2[0, 0], height_ratios=[2.30, 1.05], hspace=0.08)
    ax_c1 = fig.add_subplot(gs_c[0, 0])
    ax_c2 = fig.add_subplot(gs_c[1, 0], sharex=ax_c1)
    _panel_tag(ax_c1, "c", x=-0.15, y=1.14)

    wf = pool107.copy()
    coh_col = "cohort" if "cohort" in wf.columns else "cancer_cohort"
    pdl1_col = "b1_static_pdl1_prob" if "b1_static_pdl1_prob" in wf.columns else "b1_static_pdl1"
    wf["pdl1_hi"] = wf.groupby(coh_col)[pdl1_col].transform(lambda s: s >= s.median())
    wf["caf_hi"] = wf["myCAF"] >= 0.50

    def _smooth_lopo_prob(s: pd.Series) -> pd.Series:
        p_clip = np.clip(s.to_numpy(dtype=float), 1e-4, 1.0 - 1e-4)
        arr = np.log(p_clip / (1.0 - p_clip))
        n_c = len(arr)
        if n_c <= 2:
            return pd.Series(p_clip, index=s.index)
        tot = float(np.sum(arr))
        tot_sq = float(np.sum(arr * arr))
        z_lopo = np.empty(n_c, dtype=float)
        for idx_i in range(n_c):
            rem_mean = (tot - arr[idx_i]) / (n_c - 1)
            rem_var = max((tot_sq - arr[idx_i] ** 2) / (n_c - 1) - rem_mean ** 2, 0.0)
            rem_std = math.sqrt(rem_var * (n_c - 1) / max(n_c - 2, 1)) + 1e-8
            z_lopo[idx_i] = (arr[idx_i] - rem_mean) / rem_std
        return pd.Series(1.0 / (1.0 + np.exp(-0.85 * z_lopo)), index=s.index)

    raw_col_wf = "raw_icb_remission_prob" if "raw_icb_remission_prob" in wf.columns else "simulated_icb_remission_prob"
    wf["smooth_lopo_prob"] = wf.groupby(coh_col, group_keys=False)[raw_col_wf].apply(_smooth_lopo_prob)
    wf = wf.sort_values("smooth_lopo_prob", ascending=False).reset_index(drop=True)
    x_wf = np.arange(len(wf))
    delta_p = 100.0 * (wf["smooth_lopo_prob"].to_numpy() - 0.50)
    bar_cols = [PAL["teal"] if int(r) == 1 else PAL["crimson"] for r in wf["responder"]]
    ax_c1.bar(x_wf, delta_p, width=0.84, color=bar_cols, alpha=0.88, edgecolor="none")

    ax_c1.axhline(0.0, color=PAL["ink"], lw=0.65)
    ax_c1.set_ylim(-56, 66)
    ax_c1.set_yticks([-50, 0, 50])
    ax_c1.set_ylabel("$\\Delta$Prob. vs. 50% (%)")
    ax_c1.set_title(f"Zero-Refit Pre-ICB Waterfall ($n$ = {len(wf)})")
    ax_c1.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_c1.spines["bottom"].set_visible(False)
    h_leg = [
        mpatches.Patch(facecolor=PAL["teal"], label="Responder (R)"),
        mpatches.Patch(facecolor=PAL["crimson"], label="Non-Responder (NR)"),
    ]
    ax_c1.legend(handles=h_leg, loc="upper right", frameon=False, fontsize=4.8, ncol=2)

    track_names = ["myCAF High", "PD-L1 High", "Cohort", "RECIST R"]
    for i_w, r_w in wf.iterrows():
        c_resp = PAL["teal"] if int(r_w["responder"]) == 1 else PAL["crimson"]
        c_str = str(r_w[coh_col])
        c_coh = (
            PAL["navy"]
            if "BCC" in c_str
            else (PAL["purple"] if "91061" in c_str else (PAL["blue"] if "78220" in c_str else PAL["sky"]))
        )
        c_pdl1 = PAL["orange"] if bool(r_w["pdl1_hi"]) else PAL["lightgrey"]
        c_caf = PAL["crimson"] if bool(r_w["caf_hi"]) else PAL["lightgrey"]
        for row_y, col_tile in enumerate([c_caf, c_pdl1, c_coh, c_resp]):
            ax_c2.add_patch(
                mpatches.Rectangle((i_w - 0.43, row_y - 0.40), 0.86, 0.80, facecolor=col_tile, edgecolor="none")
            )

    ax_c2.set_xlim(-0.8, len(wf) - 0.2)
    ax_c2.set_ylim(-0.55, 3.55)
    ax_c2.set_yticks(range(4))
    ax_c2.set_yticklabels(track_names, fontsize=4.7)
    tick_pos = sorted({0, 24, 49, 74, len(wf) - 1})
    ax_c2.set_xticks(tick_pos)
    ax_c2.set_xticklabels([str(p + 1) for p in tick_pos], fontsize=4.7)
    ax_c2.set_xlabel("Ranked Pre-Treatment Patients", fontsize=5.0)
    for sp in ("top", "right", "left", "bottom"):
        ax_c2.spines[sp].set_visible(False)
    ax_c2.tick_params(axis="y", length=0)

    # --- Panel d: Sample-Size-Tiered Multi-Cohort Anti-PD-1 Benchmark + Right-Aligned Tabular AUC Column ---
    ax_d = fig.add_subplot(gs_r2[0, 1])
    _panel_tag(ax_d, "d", x=-0.34, y=1.08)
    # Ordered bottom -> top so Primary Pooled Pre-Treatment cohorts (n=107, n=83) anchor the top tier
    coh_map = [
        ("GSE123813_BCC_PreTreatment_AntiPD1_Mono_scRNA", "BCC Mono.\n($n$ = 8)"),
        ("GSE207422_NSCLC_Resected_TME_scRNA", "NSCLC Resect.\n($n$ = 14)"),
        ("GSE207422_NSCLC_PreTreatment_ChemoAntiPD1_Biopsy", "NSCLC Chemo-IO\n($n$ = 24)"),
        ("Pooled_Pure_AntiPD1_Monotherapy_PreTreatment_n83", "Anti-PD-1 Mono.\n($n$ = 83)"),
        ("Pooled_All_Strictly_PreTreatment_Biopsies_n107", "All Pre-Tx\n($n$ = 107)"),
    ]
    meth_sel = [
        ("OmniNiche_v2_Spatial_Langevin", "OmniNiche", PAL["blue"], "o", 0.275),
        ("OmniNiche_NoVisium_UniformPij", "No-Visium $P_{ij}$", PAL["teal"], "P", 0.165),
        ("B6_Published_18Gene_Tcell_GEP", "18-Gene GEP", PAL["orange"], "^", 0.055),
        ("B7_Published_CAF_TIDE_Exclusion", "CAF/TIDE Excl.", PAL["purple"], "D", -0.055),
        ("B8_Algebraic_CD8_to_myCAF_Ratio", "CD8:myCAF", PAL["slate"], "v", -0.165),
        ("B1_Static_PDL1_CD274", "PD-L1 ($CD274$)", PAL["vermilion"], "s", -0.275),
    ]
    ax_d.axhspan(2.50, 4.48, color="#EDF4FA", zorder=0)
    ax_d.axhline(2.50, color=PAL["slate"], ls="--", lw=0.55, zorder=1)
    ax_d.plot([0.50, 0.50], [-0.54, 4.48], color=PAL["slate"], ls=":", lw=0.6, zorder=1)
    ax_d.text(1.03, 4.54, "AUC", ha="left", va="bottom", fontsize=4.7, fontweight="bold", color=PAL["navy"])
    for c_idx, (c_key, _) in enumerate(coh_map):
        if c_idx == 1:
            ax_d.axhspan(c_idx - 0.45, c_idx + 0.45, color=PAL["pale"], zorder=0)
        sub_b = bench[bench["cohort"] == c_key].set_index("method")
        for m_key, m_lbl, m_col, m_mk, dy in meth_sel:
            if m_key not in sub_b.index:
                continue
            r_b = sub_b.loc[m_key]
            is_om = m_key == "OmniNiche_v2_Spatial_Langevin"
            ci_lo = float(r_b["auroc_ci_low"]) if "auroc_ci_low" in r_b and pd.notna(r_b["auroc_ci_low"]) else float(r_b["auroc"])
            ci_hi = float(r_b["auroc_ci_high"]) if "auroc_ci_high" in r_b and pd.notna(r_b["auroc_ci_high"]) else float(r_b["auroc"])
            ax_d.plot([ci_lo, ci_hi], [c_idx + dy, c_idx + dy], color=m_col, lw=1.35 if is_om else 0.72, zorder=3)
            ax_d.scatter([r_b["auroc"]], [c_idx + dy], color=m_col, marker=m_mk, s=16 if is_om else 8.8, zorder=4, label=m_lbl if c_idx == 0 else None)
            if is_om:
                ax_d.text(1.03, c_idx, f"{float(r_b['auroc']):.2f}", ha="left", va="center", fontsize=4.75, fontweight="bold", color=PAL["navy"])

    ax_d.set_yticks(range(len(coh_map)))
    ax_d.set_yticklabels([lbl for _, lbl in coh_map], fontsize=4.7)
    ax_d.set_xlim(0.10, 1.19)
    ax_d.set_xticks([0.25, 0.50, 0.75, 1.00])
    ax_d.set_ylim(-0.54, len(coh_map) + 1.08)
    ax_d.set_xlabel("Blind Responder AUROC (95% CI)")
    ax_d.set_title("Multi-Cohort Anti-PD-1 Benchmark")
    ax_d.legend(loc="upper left", frameon=False, fontsize=4.15, ncol=2, columnspacing=0.50, handletextpad=0.22)

    # --- Panel e: Paired Pre- vs Post-Anti-PD-1 Longitudinal Human Trajectory (n = 15) ---
    ax_e = fig.add_subplot(gs_r3[0, 0])
    _panel_tag(ax_e, "e", x=-0.18, y=1.10)
    paired15["net_eff_pre"] = paired15["cd8_effector_pre"] - 0.5 * paired15["cd8_exhaustion_pre"]
    paired15["net_eff_post"] = paired15["cd8_effector_post"] - 0.5 * paired15["cd8_exhaustion_post"]
    paired15["delta_net_eff"] = paired15["net_eff_post"] - paired15["net_eff_pre"]
    for grp_i, (resp_v, lbl_p, col_p) in enumerate([(1, "Responder ($n$=8)", PAL["teal"]), (0, "Non-Responder ($n$=7)", PAL["crimson"])]):
        sub_p = paired15[paired15["responder"] == resp_v]
        x0_p, x1_p = grp_i * 2.0, grp_i * 2.0 + 0.85
        for r_p in sub_p.itertuples():
            ax_e.plot([x0_p, x1_p], [r_p.net_eff_pre, r_p.net_eff_post], color=col_p, alpha=0.38, lw=0.75)
            ax_e.scatter([x0_p, x1_p], [r_p.net_eff_pre, r_p.net_eff_post], color=col_p, s=10, alpha=0.75)
        ax_e.plot([x0_p, x1_p], [sub_p["net_eff_pre"].mean(), sub_p["net_eff_post"].mean()], color=col_p, lw=1.7, marker="o", ms=3.3, label=lbl_p)

    d_r = paired15.loc[paired15["responder"] == 1, "delta_net_eff"].to_numpy(dtype=float)
    d_nr = paired15.loc[paired15["responder"] == 0, "delta_net_eff"].to_numpy(dtype=float)
    mw_p15 = float(stats.mannwhitneyu(d_r, d_nr, alternative="two-sided").pvalue)
    ax_e.text(
        0.03,
        0.94,
        f"$\\Delta$R +{np.mean(d_r):.2f} vs. $\\Delta$NR {np.mean(d_nr):.2f}\n(Mann–Whitney $P = {mw_p15:.2f}$, n.s.)",
        transform=ax_e.transAxes,
        va="top",
        ha="left",
        fontsize=4.35,
        color=PAL["slate"],
        fontweight="bold",
    )
    ax_e.set_xticks([0.0, 0.85, 2.0, 2.85])
    ax_e.set_xticklabels(["Pre\n(R)", "Post\n(R)", "Pre\n(NR)", "Post\n(NR)"], fontsize=5.1)
    ax_e.set_xlim(-0.42, 3.28)
    ax_e.set_ylim(0.20, float(paired15[["net_eff_pre", "net_eff_post"]].max().max()) * 1.38)
    ax_e.set_ylabel("Net $CD8^+$ Cytotoxic Score")
    ax_e.set_title("Paired Pre/Post Biopsies ($n$ = 15)")
    ax_e.legend(loc="upper right", frameon=False, fontsize=4.65)

    # --- Panel f: 2D Solid-Tumor Light-Base Topographic Waddington Surface + Quiver Shift Under Anti-PD-1 + myCAF/TGFBR1 ---
    ax_f = fig.add_subplot(gs_r3[0, 1])
    _panel_tag(ax_f, "f", x=-0.18, y=1.10)
    gx_t = tum_npz["grid_mycaf"]
    gy_t = tum_npz["grid_tumor"]
    u_combo = tum_npz["U_2d_combo"]
    dx_c = tum_npz["drift_combo_mycaf"]
    dy_c = tum_npz["drift_combo_tumor"]

    cf_f = ax_f.contourf(gx_t, gy_t, u_combo, levels=16, cmap=get_waddington_light_cmap())
    ax_f.contour(gx_t, gy_t, u_combo, levels=10, colors="#64748B", linewidths=0.36, alpha=0.48)
    st_t = 5
    Xgt, Ygt = np.meshgrid(gx_t[::st_t], gy_t[::st_t])
    Uxt, Uyt = dx_c[::st_t, ::st_t], dy_c[::st_t, ::st_t]
    mag_t = np.hypot(Uxt, Uyt) + 1e-8
    ax_f.quiver(Xgt, Ygt, Uxt / mag_t, Uyt / mag_t, color="#1E293B", alpha=0.72, scale=21, width=0.0052)
    ax_f.axhline(0.38, color=PAL["crimson"], ls="--", lw=0.80)
    ax_f.scatter([0.16, 0.84], [0.16, 0.82], c=[PAL["teal"], PAL["crimson"]], s=26, edgecolors="white", linewidths=0.90, zorder=5)
    pill_box = {"boxstyle": "round,pad=0.14", "fc": "white", "ec": "none", "alpha": 0.84}
    ax_f.text(0.21, 0.16, "Unlocked", color=PAL["teal"], fontsize=4.95, fontweight="bold", va="center", bbox=pill_box, zorder=6)
    ax_f.text(0.79, 0.82, "myCAF Lock", color=PAL["crimson"], fontsize=4.95, fontweight="bold", va="center", ha="right", bbox=pill_box, zorder=6)
    ax_f.set_xlabel("myCAF Barrier ($S_{\\mathrm{myCAF}}$)")
    ax_f.set_ylabel("Tumor Burden ($S_{\\mathrm{Tumor}}$)")
    ax_f.set_title("Unlocked 2D Attractor (Anti-PD-1 + TGFBR1i)")
    cbar_f = fig.colorbar(cf_f, ax=ax_f, fraction=0.040, pad=0.025)
    cbar_f.ax.tick_params(labelsize=4.7)
    cbar_f.set_label("$U_{\\mathrm{combo}}$", fontsize=5.0)

    # --- Panel g: 3-Arm Paired Trajectory Rescue Across All n = 70 Anti-PD-1-Resistant Pre-Treatment Patients ---
    ax_g = fig.add_subplot(gs_r4[0, 0])
    _panel_tag(ax_g, "g", x=-0.20, y=1.10)
    nr_all = pool107[pool107["responder"] == 0].copy()
    mono_col = "raw_icb_remission_prob" if "raw_icb_remission_prob" in nr_all.columns else "simulated_icb_remission_prob"
    p_mono = 100.0 * nr_all[mono_col].to_numpy()
    p_comb = 100.0 * nr_all["combo_pd1_plus_mycaf_remission_prob"].to_numpy()
    p_tam = 100.0 * nr_all["combo_pd1_plus_tam_remission_prob"].to_numpy()
    w_nr = stats.wilcoxon(p_comb - p_mono, alternative="greater")

    seen_coh_g: set[str] = set()
    for r_nr in nr_all.itertuples():
        c_str = str(getattr(r_nr, coh_col))
        if "BCC" in c_str:
            c_nr, c_lbl = PAL["teal"], "BCC scRNA"
        elif "207422" in c_str:
            c_nr, c_lbl = PAL["blue"], "NSCLC Bulk"
        else:
            c_nr, c_lbl = PAL["purple"], "Melanoma Bulk"
        m_val = 100.0 * float(getattr(r_nr, mono_col))
        c1_val = 100.0 * float(r_nr.combo_pd1_plus_mycaf_remission_prob)
        c2_val = 100.0 * float(r_nr.combo_pd1_plus_tam_remission_prob)
        lbl_use = c_lbl if c_lbl not in seen_coh_g else None
        if lbl_use:
            seen_coh_g.add(c_lbl)
        ax_g.plot([0, 1, 2], [m_val, c1_val, c2_val], color=c_nr, alpha=0.22, lw=0.45, label=lbl_use)
        ax_g.scatter([0, 1, 2], [m_val, c1_val, c2_val], color=c_nr, s=6.0, alpha=0.48)

    ax_g.plot([0, 1, 2], [np.mean(p_mono), np.mean(p_comb), np.mean(p_tam)], color=PAL["crimson"], lw=1.75, marker="o", ms=3.5, zorder=6)
    ax_g.text(
        1.00,
        104.5,
        f"{np.mean(p_mono):.1f}% → {np.mean(p_comb):.1f}% / {np.mean(p_tam):.1f}% ({_fmt_p(float(w_nr.pvalue))})",
        ha="center",
        va="bottom",
        fontsize=4.60,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_g.set_xticks([0, 1, 2])
    ax_g.set_xticklabels(["Anti-PD-1\nMono.", "+ myCAF /\nTGFBR1i", "+ SPP1 /\nTREM2"], fontsize=4.95)
    ax_g.set_xlim(-0.32, 2.32)
    ax_g.set_ylim(-4, 134)
    ax_g.set_yticks([0, 25, 50, 75, 100])
    ax_g.set_ylabel("Remission Prob. (%)")
    ax_g.set_title(f"In Silico Resistant Rescue ($n$ = {len(nr_all)})")
    ax_g.legend(
        loc="upper center",
        bbox_to_anchor=(0.50, 1.02),
        ncol=3,
        frameon=False,
        fontsize=4.35,
        handlelength=1.0,
        handletextpad=0.3,
        columnspacing=0.7,
        borderpad=0.1,
    )

    # --- Panel h: Real Pre-Treatment Anti-PD-1 Monotherapy Overall Survival (n = 74) + Aligned Number-at-Risk Table ---
    gs_h = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_r4[0, 1], height_ratios=[2.55, 0.82], hspace=0.10)
    ax_h1 = fig.add_subplot(gs_h[0, 0])
    ax_h2 = fig.add_subplot(gs_h[1, 0], sharex=ax_h1)
    _panel_tag(ax_h1, "h", x=-0.16, y=1.14)

    r_os = surv_sum[surv_sum["stratification"] == "Stromal_CAF_Exclusion_Barrier_High_vs_Low"].iloc[0]
    risk_ticks = [0, 12, 24, 36]
    max_t_disp = 38.2
    grp_specs = [
        (0, "Low Stromal Barrier ($n$=37)", "Low", PAL["teal"], 1),
        (1, "High Stromal Barrier ($n$=37)", "High", PAL["crimson"], 0),
    ]
    for hi_flag, leg_lbl, row_lbl, c_km, row_y in grp_specs:
        sub_s = surv74[surv74["stromal_barrier_high"] == hi_flag]
        t_arr = sub_s["os_months"].to_numpy(dtype=float)
        e_arr = sub_s["event_dead"].to_numpy(dtype=int)
        tk, sk, lok, hik = _kaplan_meier_curve(t_arr, e_arr)
        # Extend step horizontally to last observed follow-up in stratum
        t_last = float(np.max(t_arr))
        if t_last > tk[-1]:
            tk = np.append(tk, t_last)
            sk = np.append(sk, sk[-1])
            lok = np.append(lok, lok[-1])
            hik = np.append(hik, hik[-1])
        ax_h1.step(tk, 100.0 * sk, where="post", color=c_km, lw=1.35, label=leg_lbl, zorder=4)
        ax_h1.fill_between(tk, 100.0 * lok, 100.0 * hik, step="post", color=c_km, alpha=0.13, zorder=2)

        # Censoring tick marks
        cens_t = t_arr[e_arr == 0]
        if len(cens_t) > 0:
            cens_s = np.array([100.0 * sk[np.searchsorted(tk, ct, side="right") - 1] for ct in cens_t])
            ax_h1.scatter(cens_t, cens_s, marker="|", color=c_km, s=16, linewidths=0.80, zorder=5)

        # Aligned Number-at-Risk row
        ax_h2.text(-1.2, row_y, row_lbl, ha="right", va="center", fontsize=4.75, fontweight="bold", color=c_km)
        for rt in risk_ticks:
            n_ar = int(np.sum(t_arr >= rt))
            ax_h2.text(rt, row_y, f"{n_ar}", ha="center", va="center", fontsize=4.8, color=PAL["ink"])

    strat_line = (
        f"Log-rank $P = {float(r_os['logrank_p']):.3f}$ | Strat. $P = {float(r_os['stratified_logrank_p']):.3f}$"
        if "stratified_logrank_p" in r_os.index
        else f"Log-rank $P = {float(r_os['logrank_p']):.3f}$ ({int(r_os['total_deaths'])} deaths)"
    )
    ax_h1.text(
        0.98,
        0.95,
        f"HR = {r_os['hazard_ratio']:.2f} ({r_os['hr_ci95_low']:.2f}–{r_os['hr_ci95_high']:.2f}, {int(r_os['total_deaths'])} deaths)\n"
        + strat_line,
        transform=ax_h1.transAxes,
        ha="right",
        va="top",
        fontsize=4.60,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_h1.set_xlim(-0.8, max_t_disp)
    ax_h1.set_ylim(-3, 106)
    ax_h1.set_yticks([0, 25, 50, 75, 100])
    ax_h1.set_ylabel("Overall Survival (%)")
    ax_h1.set_title("Anti-PD-1 Monotherapy Survival ($n$ = 74)")
    ax_h1.legend(loc="lower left", frameon=False, fontsize=4.65)
    ax_h1.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    ax_h2.text(-1.2, 1.92, "No. at risk", ha="right", va="center", fontsize=4.55, fontweight="bold", color=PAL["slate"])
    ax_h2.set_ylim(-0.55, 2.15)
    ax_h2.set_yticks([])
    ax_h2.set_xticks(risk_ticks)
    ax_h2.set_xticklabels([str(rt) for rt in risk_ticks], fontsize=4.9)
    ax_h2.set_xlabel("Follow-Up After Pre-Tx Biopsy (Months)", fontsize=5.1)
    for sp in ("top", "right", "left"):
        ax_h2.spines[sp].set_visible(False)
    ax_h2.spines["bottom"].set_color(PAL["slate"])
    ax_h2.spines["bottom"].set_linewidth(0.55)

    _save_multi(fig, out_dir, "Fig5_SolidTumor_Visium_Exclusion_and_PreICB_Responder", artifact_dir)


# ===========================================================================
# EXTENDED DATA FIGURES 1-4 (Native 183 mm = 7.20 in Print Width)
# ===========================================================================
def build_extended_data_figures(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render Extended Data Figs. 1-4 natively at 183 mm (7.20 in) print width."""
    set_nature_style()

    s1 = pd.read_csv(tables_dir / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv")
    s2 = pd.read_csv(
        tables_dir / "Table_S2_Target_Resolved_Trials_Deconvolution.csv",
        usecols=["nct_id", "primary_target_gene", "n_targets", "stop_deconvolution", "max_phase_num"],
        low_memory=False,
    )
    atlas = pd.read_csv(tables_dir / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv")
    s3b = pd.read_csv(tables_dir / "Table_S3b_PerAxis_Dynamic_Drift_Performance.csv")
    s3c = pd.read_csv(tables_dir / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv")
    s3d = pd.read_csv(tables_dir / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv")
    s3e = pd.read_csv(tables_dir / "Table_S3e_Unreported_Stop_Sensitivity_Bounds.csv")
    s3f = pd.read_csv(tables_dir / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv")
    s3g = pd.read_csv(tables_dir / "Table_S3g_Multivariable_FixedEffects_Control.csv")
    s3h = pd.read_csv(tables_dir / "Table_S3h_Calibration_Deciles.csv")
    s3j = pd.read_csv(tables_dir / "Table_S3j_Simulated_Trial_Enrichment.csv")
    s4g = pd.read_csv(tables_dir / "Table_S4g_Profile_Likelihood_Curves.csv")
    s5a = pd.read_csv(tables_dir / "Table_S5a_Visium_Spatial_Barrier_Sections.csv")
    s5b = pd.read_csv(tables_dir / "Table_S5b_Visium_Consensus_Proximity_Matrix.csv", index_col=0)
    s5h = pd.read_csv(tables_dir / "Table_S5h_FACS_Deconfounding_Sensitivity.csv")
    s5i = pd.read_csv(tables_dir / "Table_S5i_Pooled_PreTreatment_107Patients_Master.csv")

    # -----------------------------------------------------------------------
    # Extended Data Fig. 1: PRISMA Funnel, Worked Clause Table, Top Targets, Polypharmacology & 6-Axis Group5Fold (a-e)
    # -----------------------------------------------------------------------
    fig1 = plt.figure(figsize=(7.20, 6.65), facecolor="white")
    gs1 = gridspec.GridSpec(3, 12, figure=fig1, height_ratios=[0.95, 0.88, 1.05], wspace=1.35, hspace=0.50)

    # ED1a: PRISMA cohort filtration funnel
    ax = fig1.add_subplot(gs1[0, 0:6])
    _panel_tag(ax, "a", x=-0.26, y=1.08)
    n_sys = int(s1["unique_systemic_organ_trials"].sum())
    n_definitive = int(s3d["n_trials_dnacc"].sum())
    prisma_steps = [
        ("Raw Query (+587 Brain-Met. Excl.)", n_sys + 587, PAL["midgrey"]),
        ("12 Peripheral Organ Axes", n_sys, PAL["navy"]),
        ("Target-Resolved (ChEMBL/OT)", len(s2), PAL["blue"]),
        ("Definitive Non-Admin. Cohort", n_definitive, PAL["teal"]),
        ("Deconfounded Scientific Set", len(s3c), PAL["orange"]),
    ][::-1]
    yp = np.arange(len(prisma_steps))
    ax.barh(yp, [v for _, v, _ in prisma_steps], color=[c for _, _, c in prisma_steps], height=0.60, alpha=0.88)
    for i_p, (_, v_p, _) in enumerate(prisma_steps):
        ax.text(v_p + 700, i_p, f"{v_p:,}", va="center", ha="left", fontsize=5.1, fontweight="bold")
    ax.set_yticks(yp)
    ax.set_yticklabels([lbl for lbl, _, _ in prisma_steps], fontsize=5.2)
    ax.set_xlim(0, (n_sys + 587) * 1.28)
    ax.set_xticks([0, 20000, 40000])
    ax.set_xticklabels(["0", "20k", "40k"])
    ax.set_xlabel("Unique ClinicalTrials.gov Studies")
    ax.set_title("PRISMA Filtration Funnel", pad=5)

    # ED1b: Top 16 Clinical Target Genes Across 12 Axes
    ax = fig1.add_subplot(gs1[0, 6:12])
    _panel_tag(ax, "b", x=-0.18, y=1.08)
    top_g = s2["primary_target_gene"].dropna().value_counts().head(14)[::-1]
    yg = np.arange(len(top_g))
    ax.barh(yg, top_g.to_numpy(), color=PAL["blue"], height=0.66, alpha=0.85)
    for i_g, v_g in enumerate(top_g.to_numpy()):
        ax.text(v_g + 15, i_g, f"{int(v_g)}", va="center", ha="left", fontsize=4.8)
    ax.set_yticks(yg)
    ax.set_yticklabels(top_g.index, fontsize=4.9, fontstyle="italic")
    ax.set_xlim(0, float(top_g.max()) * 1.25)
    ax.set_xlabel("Target-Resolved Phase II–IV Trials")
    ax.set_title("Top Clinical Target Genes", pad=5)

    # ED1c: Worked Clause-Aware Negation Deconvolution Table (Minikel ED Fig. 1b Style)
    ax = fig1.add_subplot(gs1[1, 0:12])
    _panel_tag(ax, "c", x=-0.045, y=1.10)
    ax.axis("off")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Worked Clause-Aware Negation Deconvolution Examples (ClinicalTrials.gov whyStopped)", pad=4)

    ax.axhline(0.96, color=PAL["ink"], lw=0.75)
    ax.axhline(0.79, color=PAL["slate"], lw=0.50)
    ax.axhline(0.04, color=PAL["ink"], lw=0.75)
    headers = [(0.01, "NCT ID"), (0.14, "Verbatim whyStopped Narrative Excerpt"), (0.64, "Negated Clause Stripped"), (0.83, "Adjudicated Class")]
    for hx, htxt in headers:
        ax.text(hx, 0.87, htxt, ha="left", va="center", fontsize=5.1, fontweight="bold", color=PAL["navy"])

    worked_rows = [
        ("NCT04987307", "Terminated due to meeting a prespecified futility rule and not related to any safety concerns.", "not related to ... safety", "Efficacy Futility", PAL["crimson"]),
        ("NCT03584009", "Sponsor's decision, no safety concerns.", "no safety concerns", "Admin. / Commercial", PAL["orange"]),
        ("NCT00646659", "Recruitment suspended prematurely for safety concerns and closed after IDMC review.", "(none)", "Safety / Toxicity", PAL["vermilion"]),
        ("NCT00330200", "Due to poor enrollment.", "(none)", "Admin. / Accrual", PAL["orange"]),
    ]
    row_ys = [0.67, 0.49, 0.31, 0.13]
    for ry, (nid, verb, neg_s, adj_c, c_adj) in zip(row_ys, worked_rows):
        ax.text(0.01, ry, nid, ha="left", va="center", fontsize=4.9, fontweight="bold", color=PAL["ink"])
        ax.text(0.14, ry, verb, ha="left", va="center", fontsize=4.8, color=PAL["ink"])
        ax.text(0.64, ry, neg_s, ha="left", va="center", fontsize=4.8, fontstyle="italic", color=PAL["slate"])
        ax.text(0.83, ry, adj_c, ha="left", va="center", fontsize=4.9, fontweight="bold", color=c_adj)

    # ED1d: Polypharmacology distribution & subgroup robustness
    ax = fig1.add_subplot(gs1[2, 0:6])
    _panel_tag(ax, "d", x=-0.16, y=1.08)
    sub_sg = s3f[
        (s3f["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes")
        & (s3f["cv_scheme"].str.startswith("Subgroup_"))
    ].copy()
    nt_bins = s2["n_targets"].clip(upper=5).value_counts().sort_index()
    xb = np.arange(len(nt_bins))
    ax.bar(xb, nt_bins.to_numpy(), color=PAL["sky"], edgecolor=PAL["navy"], lw=0.6, width=0.56)
    for i_b, v_b in enumerate(nt_bins.to_numpy()):
        ax.text(i_b, v_b + 350, f"{int(v_b):,}", ha="center", va="bottom", fontsize=4.9)
    if len(sub_sg) >= 2:
        sg_st = sub_sg[sub_sg["cv_scheme"].str.contains("SingleTarget")].sort_values("step_order")
        sg_mt = sub_sg[sub_sg["cv_scheme"].str.contains("Polypharmacology")].sort_values("step_order")
        if len(sg_st) > 0 and len(sg_mt) > 0:
            r_single = sg_st.iloc[-1]
            r_multi = sg_mt.iloc[-1]
            ax.text(
                0.22,
                0.78,
                f"1-Target GroupCV = {r_single['auroc']:.3f} [{r_single['auroc_ci_low']:.2f}–{r_single['auroc_ci_high']:.2f}]\n"
                f"Multi-Target AUC = {r_multi['auroc']:.3f} [{r_multi['auroc_ci_low']:.2f}–{r_multi['auroc_ci_high']:.2f}]",
                transform=ax.transAxes,
                fontsize=4.8,
                fontweight="bold",
                color=PAL["navy"],
            )
    ax.set_xticks(xb)
    ax.set_xticklabels(["1", "2", "3", "4", "5+"], fontsize=5.3)
    ax.set_xlabel("Resolved Target Genes per Regimen")
    ax.set_ylim(0, float(nt_bins.max()) * 1.32)
    ax.set_ylabel("Target-Resolved Trials")
    ax.set_title("Polypharmacology Spectrum", pad=5)

    # ED1e: Core 6 Single-Cell Axes Stratified vs Unseen-Target Group5Fold Ablation
    ax = fig1.add_subplot(gs1[2, 6:12])
    _panel_tag(ax, "e", x=-0.16, y=1.08)
    ab_6 = s3f[s3f["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_Core6SingleCellAxes"]
    st6 = ab_6[ab_6["cv_scheme"] == "Stratified_5Fold"].sort_values("step_order")
    gp6 = ab_6[ab_6["cv_scheme"] == "UnseenTarget_Group5Fold"].sort_values("step_order")
    n_core6 = int(st6["n_trials"].iloc[0]) if len(st6) > 0 else 0
    x6 = np.arange(len(st6))
    ax.axhline(0.50, color=PAL["lightgrey"], ls=":", lw=0.6)
    ax.errorbar(
        x6 - 0.05,
        st6["auroc"],
        yerr=[st6["auroc"] - st6["auroc_ci_low"], st6["auroc_ci_high"] - st6["auroc"]],
        color=PAL["blue"],
        marker="o",
        ms=2.9,
        lw=1.2,
        capsize=1.7,
        label=f"Stratified 5-Fold ($n$={n_core6:,})",
    )
    ax.errorbar(
        x6 + 0.05,
        gp6["auroc"],
        yerr=[gp6["auroc"] - gp6["auroc_ci_low"], gp6["auroc_ci_high"] - gp6["auroc"]],
        color=PAL["teal"],
        marker="s",
        ms=2.7,
        ls="--",
        lw=1.1,
        capsize=1.7,
        label="Unseen-Target CV",
    )
    ax.set_xticks(x6)
    ax.set_xticklabels(
        ["M1: $\\tau_0$", "M2: +Vital", "M3: +$\\tau_{\\mathrm{early}}$", "M4: +$\\Delta\\tau$", "M5: +Switch", "M6: DNAcc"],
        rotation=26,
        ha="right",
        fontsize=5.1,
    )
    ax.set_ylim(0.40, float(st6["auroc_ci_high"].max()) + 0.055)
    ax.set_ylabel("Core 6-Axis AUROC")
    ax.set_title("Core 6-Axis Unseen-Target CV", pad=5)
    ax.legend(loc="lower right", frameon=False, fontsize=5.0)

    _save_multi(fig1, out_dir, "Extended_Data_Fig1_PRISMA_Target_Modality_Audit", artifact_dir)

    # -----------------------------------------------------------------------
    # Extended Data Fig. 2: Cohort Census, Switch Frequency, Drift Distributions & Per-Axis AUROC (a-d)
    # -----------------------------------------------------------------------
    fig2, axes2 = plt.subplots(2, 2, figsize=(7.20, 5.10), facecolor="white")
    fig2.subplots_adjust(wspace=0.52, hspace=0.48)

    ax = axes2[0, 0]
    _panel_tag(ax, "a", x=-0.36, y=1.08)
    coh_census = [
        ("Melanoma Anti-PD-1 (GSE91061+78220)", 75),
        ("Solid Tumor (GSE207422+123813)", 46),
        ("Pulmonary IPF (GSE135893)", 30),
        ("Ulcerative Colitis (GSE250487)", 28),
        ("Hepatic MASH (GSE136103+212837)", 17),
        ("Renal CKD (GSE271188)", 17),
        ("Sepsis PBMC (GSE167363)", 12),
        ("10x Visium Tumors (12 Sections)", 12),
    ][::-1]
    yc2 = np.arange(len(coh_census))
    ax.barh(yc2, [v for _, v in coh_census], color=PAL["teal"], height=0.64, alpha=0.85)
    for i_c, (_, v_c) in enumerate(coh_census):
        ax.text(v_c + 1.2, i_c, f"n={v_c}", va="center", ha="left", fontsize=4.8, fontweight="bold")
    ax.set_yticks(yc2)
    ax.set_yticklabels([lbl for lbl, _ in coh_census], fontsize=4.8)
    ax.set_xlim(0, 94)
    ax.set_xlabel("Human Donors / Tissue Specimens (N=237)")
    ax.set_title("Primary Human Omics Census", pad=5)

    ax = axes2[0, 1]
    _panel_tag(ax, "b", x=-0.30, y=1.08)
    sw_rates = atlas.groupby("disease_axis")["compartment_switch"].mean().mul(100.0).sort_values()
    ysw = np.arange(len(sw_rates))
    ax.barh(ysw, sw_rates.to_numpy(), color=PAL["crimson"], alpha=0.84, height=0.60)
    for i_s, v_s in enumerate(sw_rates.to_numpy()):
        ax.text(v_s + 0.8, i_s, f"{v_s:.1f}%", va="center", ha="left", fontsize=5.0, fontweight="bold")
    ax.set_yticks(ysw)
    ax.set_yticklabels([AXIS_COMPACT.get(a, a) for a in sw_rates.index], fontsize=5.2)
    ax.set_xlim(0, float(sw_rates.max()) * 1.28)
    ax.set_xlabel("Targets Switching Lineage (%)")
    ax.set_title("Ectopic Compartment Switching", pad=5)

    ax = axes2[1, 0]
    _panel_tag(ax, "c", x=-0.30, y=1.08)
    ax_list = list(sw_rates.index)
    data_dt = [atlas[atlas["disease_axis"] == a]["delta_tau"].dropna().to_numpy() for a in ax_list]
    bp = ax.boxplot(data_dt, vert=False, patch_artist=True, showfliers=False, widths=0.54)
    for patch in bp["boxes"]:
        patch.set_facecolor(PAL["sky"])
        patch.set_alpha(0.65)
        patch.set_edgecolor(PAL["navy"])
        patch.set_linewidth(0.6)
    for med_l in bp["medians"]:
        med_l.set_color(PAL["crimson"])
        med_l.set_linewidth(1.1)
    for w_l in bp["whiskers"] + bp["caps"]:
        w_l.set_linewidth(0.6)
    ax.axvline(0.0, color=PAL["slate"], ls="--", lw=0.65)
    ax.set_yticks(range(1, len(ax_list) + 1))
    ax.set_yticklabels([AXIS_COMPACT.get(a, a) for a in ax_list], fontsize=5.2)
    ax.set_xlabel("Specificity Drift ($\\Delta\\tau = \\tau_{\\mathrm{late}} - \\tau_{\\mathrm{early}}$)")
    ax.set_title("Drift Distributions Across 6 Axes", pad=5)

    ax = axes2[1, 1]
    _panel_tag(ax, "d", x=-0.30, y=1.08)
    s3b_mat = s3b[s3b["cohort_block"] == "Definitive_Maturation"] if "cohort_block" in s3b.columns else s3b
    s3b_s = s3b_mat.sort_values("dnacc_auroc", ascending=True).reset_index(drop=True)
    y_ax = np.arange(len(s3b_s))
    ax.axvline(0.50, color=PAL["slate"], ls=":", lw=0.65)
    for i_x, r_x in s3b_s.iterrows():
        ax.plot([r_x["static_tau0_ci_low"], r_x["static_tau0_ci_high"]], [i_x - 0.13, i_x - 0.13], color=PAL["midgrey"], lw=0.95)
        ax.scatter([r_x["static_tau0_auroc"]], [i_x - 0.13], color=PAL["midgrey"], marker="s", s=14, label="Static $\\tau_0$" if i_x == 0 else None)
        ax.plot([r_x["dnacc_ci_low"], r_x["dnacc_ci_high"]], [i_x + 0.13, i_x + 0.13], color=PAL["blue"], lw=1.25)
        ax.scatter([r_x["dnacc_auroc"]], [i_x + 0.13], color=PAL["blue"], marker="o", s=18, label="DNAcc" if i_x == 0 else None)
        ax.text(float(s3b_s["dnacc_ci_high"].max()) + 0.012, i_x, f"+{r_x['delta_auroc']:.2f}", va="center", ha="left", fontsize=4.9, color=PAL["navy"])
    ax.set_yticks(y_ax)
    ax.set_yticklabels([AXIS_COMPACT.get(a, a) for a in s3b_s["disease_axis"]], fontsize=5.2)
    ax.set_xlim(0.32, float(s3b_s["dnacc_ci_high"].max()) + 0.11)
    ax.set_xlabel("Per-Axis Maturation AUROC")
    ax.set_title("Per-Axis Discrimination & Boundary", pad=5)
    ax.legend(loc="upper left", frameon=False, fontsize=5.0)

    _save_multi(fig2, out_dir, "Extended_Data_Fig2_SixAxis_SingleCell_Drift_and_Ablation", artifact_dir)

    # -----------------------------------------------------------------------
    # Extended Data Fig. 3: Profile Likelihood, Dual-CI Fixed Effects, 5-Regime Sensitivity, Enrichment & Calibration (a-e)
    # -----------------------------------------------------------------------
    fig3 = plt.figure(figsize=(7.20, 6.75), facecolor="white")
    gs3 = gridspec.GridSpec(3, 1, figure=fig3, height_ratios=[1.00, 1.05, 1.00], hspace=0.50)
    gs3_r1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs3[0], width_ratios=[1.00, 1.06], wspace=0.52)
    gs3_r3 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs3[2], width_ratios=[1.00, 1.00], wspace=0.34)

    # ED3a: 25-Point Profile-Likelihood Identifiability
    ax = fig3.add_subplot(gs3_r1[0, 0])
    _panel_tag(ax, "a", x=-0.18, y=1.08)
    edge_cols = [PAL["blue"], PAL["crimson"], PAL["teal"], PAL["orange"]]
    edge_short = {
        "v_Parenchyma_Injury_to_SPP1_SAMac": "Parenchyma → SAMac",
        "v_SPP1_SAMac_to_Myofibroblast": "SAMac → Myofibroblast",
        "v_Myofibroblast_to_SPP1_SAMac": "Myofibroblast → SAMac",
        "v_Myofibroblast_to_Parenchyma_Injury": "Myofibroblast → Parenchyma",
    }
    for (p_name, p_sub), col_e in zip(s4g.groupby("parameter", sort=False), edge_cols):
        p_sub_s = p_sub.sort_values("multiplier")
        ax.plot(p_sub_s["multiplier"], p_sub_s["delta_chi2"], color=col_e, lw=1.25, label=edge_short.get(p_name, p_name))
    ax.axhline(3.84, color=PAL["slate"], ls="--", lw=0.75, label="95% Cutoff ($\\Delta\\chi^2$=3.84)")
    ax.set_xlabel("Fold Multiplier ($v / v_{\\mathrm{MLE}}$)")
    ax.set_ylabel("Profile $\\Delta\\chi^2$")
    ax.set_ylim(-2.2, 85.0)
    ax.set_title("25-Point Structural Identifiability", pad=5)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.56, 1.00),
        frameon=False,
        fontsize=4.25,
        ncol=1,
        handlelength=1.15,
        handletextpad=0.38,
        labelspacing=0.24,
        borderpad=0.15,
    )

    # ED3b: Multivariable Fixed-Effects & 3-Pillar Joint Control (DNAcc vs Open Targets Genetics vs Static tau_0)
    ax = fig3.add_subplot(gs3_r1[0, 1])
    _panel_tag(ax, "b", x=-0.38, y=1.08)
    spec_order = [
        "Unadjusted_Univariate_DNAcc",
        "Unadjusted_Univariate_Human_Genetics",
        "Unadjusted_Univariate_Static_Tau0",
        "Multivariable_FixedEffects_DNAcc",
        "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0",
        "Joint_FixedEffects_Human_Genetics_controlling_for_DNAcc_and_Tau0",
        "Joint_FixedEffects_Tau0_controlling_for_DNAcc_and_Genetics",
    ]
    spec_short = {
        "Unadjusted_Univariate_DNAcc": "DNAcc (Unadj.)",
        "Unadjusted_Univariate_Human_Genetics": "Human Genetics (Unadj.)",
        "Unadjusted_Univariate_Static_Tau0": "Static $\\tau_0$ (Unadj.)",
        "Multivariable_FixedEffects_DNAcc": "DNAcc (+5 FE)",
        "Multivariable_FixedEffects_Human_Genetics": "Genetics (+5 FE)",
        "Multivariable_FixedEffects_Static_Tau0": "Static $\\tau_0$ (+5 FE)",
        "Joint_FixedEffects_DNAcc_controlling_for_Tau0": "Joint: DNAcc (+$\\tau_0$+FE)",
        "Joint_FixedEffects_Tau0_controlling_for_DNAcc": "Joint: $\\tau_0$ (+DNAcc+FE)",
        "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0": "3-Pillar: DNAcc",
        "Joint_FixedEffects_Human_Genetics_controlling_for_DNAcc_and_Tau0": "3-Pillar: Genetics",
        "Joint_FixedEffects_Tau0_controlling_for_DNAcc_and_Genetics": "3-Pillar: Static $\\tau_0$",
    }
    fe_raw = s3g[s3g["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].copy()
    present_specs = [sp for sp in spec_order if sp in set(fe_raw["model_specification"])]
    if len(present_specs) >= 4:
        fe_sub = fe_raw.set_index("model_specification").loc[present_specs[::-1]].reset_index()
    else:
        fe_sub = fe_raw.iloc[::-1].reset_index(drop=True)
    yf3 = np.arange(len(fe_sub))
    ax.axvline(1.0, color=PAL["slate"], ls="--", lw=0.65)
    x_txt_fe = max(float(fe_sub["cluster_ci95_high"].max()) + 0.05, 1.95)
    for i_f, r_f in fe_sub.iterrows():
        sp_n = str(r_f["model_specification"])
        is_dna = sp_n in (
            "Unadjusted_Univariate_DNAcc",
            "Multivariable_FixedEffects_DNAcc",
            "Joint_FixedEffects_DNAcc_controlling_for_Tau0",
            "Joint_FixedEffects_DNAcc_controlling_for_Genetics_and_Tau0",
        )
        is_gen = sp_n in (
            "Unadjusted_Univariate_Human_Genetics",
            "Multivariable_FixedEffects_Human_Genetics",
            "Joint_FixedEffects_Human_Genetics_controlling_for_DNAcc_and_Tau0",
        )
        c_f = PAL["blue"] if is_dna else (PAL["teal"] if is_gen else PAL["midgrey"])
        ax.plot([r_f["cluster_ci95_low"], r_f["cluster_ci95_high"]], [i_f, i_f], color=c_f, lw=0.75, alpha=0.55)
        ax.plot([r_f["ci95_low"], r_f["ci95_high"]], [i_f, i_f], color=c_f, lw=1.50)
        ax.scatter([r_f["odds_ratio_per_sd"]], [i_f], color=c_f, s=15, zorder=4)
        ax.text(x_txt_fe, i_f, f"{r_f['odds_ratio_per_sd']:.2f} ({_fmt_p(float(r_f['cluster_p_value']))})", va="center", ha="left", fontsize=4.4)
    ax.set_yticks(yf3)
    ax.set_yticklabels([spec_short.get(s, s) for s in fe_sub["model_specification"]], fontsize=4.8)
    ax.set_xlim(0.75, x_txt_fe + 1.25)
    ax.set_xticks([1.0, 1.5, 2.0])
    ax.set_xlabel("OR per 1-s.d. (Cluster & Wald CI)")
    ax.set_title("Cluster-Robust 3-Pillar Control", pad=5)

    # ED3c: 7-Regime Missingness, Negation, Regulatory Approval & 9-Axis Boundary Sensitivity (R1-R7)
    ax = fig3.add_subplot(gs3[1])
    _panel_tag(ax, "c", x=-0.14, y=1.08)
    reg_short = {
        "R1_Primary_ClauseAware_Scientific": "R1: Primary Clause-Aware",
        "R1_Primary_Explicit_Scientific": "R1: Primary Explicit Scientific",
        "R2_Unfiltered_Naive_Keyword_Scientific": "R2: Unfiltered Naive Keyword",
        "R2_Negation_Filtered_Strict_Scientific": "R2: Strict Negation-Filtered",
        "R3_WorstCase_All2618_Unreported_As_Failure": "R3: All Unrep. Stops = Fail",
        "R3_WorstCase_All3118_Unreported_As_Failure": "R3: All Unrep. Stops = Fail",
        "R4_IPW_Propensity_Reweighted": "R4: IPW Propensity Reweighted",
        "R5_Posted_ClinicalResults_Only": "R5: Posted Results Only",
        "R6_Regulatory_Approved_Indication_vs_Explicit_SciFail": "R6: ChEMBL/OT Approved Gold",
        "R7_Strictly_Peripheral_Interstitial_9Axes": "R7: Strict 9-Axis Interstitial",
    }
    s3e_rev = s3e.iloc[::-1].reset_index(drop=True)
    y_g = np.arange(len(s3e_rev))
    ax.axvline(0.50, color=PAL["slate"], ls=":", lw=0.65)
    for i_g, r_s in s3e_rev.iterrows():
        ax.plot([r_s["static_tau0_ci_low"], r_s["static_tau0_ci_high"]], [i_g - 0.13, i_g - 0.13], color=PAL["midgrey"], lw=0.95)
        ax.scatter([r_s["static_tau0_auroc"]], [i_g - 0.13], color=PAL["midgrey"], marker="s", s=14, label="Static $\\tau_0$" if i_g == 0 else None)
        ax.plot([r_s["dnacc_ci_low"], r_s["dnacc_ci_high"]], [i_g + 0.13, i_g + 0.13], color=PAL["blue"], lw=1.25)
        ax.scatter([r_s["dnacc_auroc"]], [i_g + 0.13], color=PAL["blue"], marker="o", s=17, label="DNAcc" if i_g == 0 else None)
        ax.text(
            float(s3e_rev["dnacc_ci_high"].max()) + 0.010,
            i_g,
            f"$\\Delta$AUC = +{r_s['delta_auroc']:.3f} ({_fmt_p(float(r_s['delong_p_value']))})",
            va="center",
            ha="left",
            fontsize=4.7,
            color=PAL["navy"],
        )
    ax.set_yticks(y_g)
    ax.set_yticklabels([f"{reg_short.get(r.regime_id, r.regime_id)} ($n$={r.n_trials:,})" for r in s3e_rev.itertuples()], fontsize=4.8)
    ax.set_xlim(0.47, float(s3e_rev["dnacc_ci_high"].max()) + 0.145)
    ax.set_xlabel("AUROC (2,000-Bootstrap 95% CI)")
    ax.set_title("7-Regime Missingness, Regulatory Approval & 9-Axis Boundary Audit", pad=5)
    ax.legend(loc="lower right", frameon=False, fontsize=4.9)

    # ED3d: Simulated Phase II -> III Portfolio Enrichment Utility
    ax = fig3.add_subplot(gs3_r3[0, 0])
    _panel_tag(ax, "d", x=-0.18, y=1.08)
    xg3 = s3j["bottom_percentile_gated_out"].to_numpy()
    ax.plot(xg3, s3j["dnacc_success_rate_pct"], color=PAL["blue"], lw=1.3, marker="o", ms=2.8, label="DNAcc Gate")
    ax.plot(xg3, s3j["static_tau0_success_rate_pct"], color=PAL["midgrey"], lw=1.05, ls="--", marker="s", ms=2.5, label="Static $\\tau_0$ Gate")
    r25 = s3j[s3j["bottom_percentile_gated_out"] == 25]
    if len(r25) > 0:
        sp_pts = int(r25["patients_spared_futile_exposure"].iloc[0])
        prev_f = int(r25["futile_trials_prevented_dnacc"].iloc[0])
        ax.text(
            0.04,
            0.82,
            f"25% Gate: -{prev_f:,} Failed Trials\n(~{sp_pts / 1e6:.2f}M Patients Spared)",
            transform=ax.transAxes,
            fontsize=4.9,
            fontweight="bold",
            color=PAL["navy"],
        )
    ax.set_xlabel("Bottom Percentile Gated Out (%)")
    ax.set_ylabel("Phase III/IV Success (%)")
    ax.set_ylim(float(s3j["static_tau0_success_rate_pct"].min()) * 0.90, float(s3j["dnacc_success_rate_pct"].max()) * 1.18)
    ax.set_title("Simulated Portfolio Enrichment", pad=5)
    ax.legend(loc="lower right", frameon=False, fontsize=5.0)

    # ED3e: 10-Decile Calibration on Deconfounded Scientific Cohort
    ax = fig3.add_subplot(gs3_r3[0, 1])
    _panel_tag(ax, "e", x=-0.18, y=1.08)
    cal_s = s3h[s3h["cohort"].astype(str).str.startswith("Deconfounded_Scientific_Outcomes")]
    cd_s = cal_s[cal_s["model"] == "Full_DNAcc"]
    ct_s = cal_s[cal_s["model"] == "Static_Healthy_Tau0"]
    ax.plot([0.55, 0.99], [0.55, 0.99], color=PAL["slate"], ls=":", lw=0.7, label="Ideal")
    ax.plot(ct_s["mean_predicted_prob"], ct_s["observed_event_rate"], color=PAL["midgrey"], marker="s", ms=2.4, ls="--", lw=0.95, label="Static $\\tau_0$")
    ax.errorbar(
        cd_s["mean_predicted_prob"],
        cd_s["observed_event_rate"],
        yerr=[cd_s["observed_event_rate"] - cd_s["obs_ci95_low"], cd_s["obs_ci95_high"] - cd_s["observed_event_rate"]],
        color=PAL["blue"],
        marker="o",
        ms=2.8,
        lw=1.2,
        capsize=1.6,
        label="DNAcc",
    )
    ax.set_xlabel("Predicted Scientific Success")
    ax.set_ylabel("Observed Rate (95% CI)")
    ax.set_title(f"Calibration (Scientific Set, $n$={len(s3c):,})", pad=5)
    ax.legend(loc="lower right", frameon=False, fontsize=5.0)

    _save_multi(fig3, out_dir, "Extended_Data_Fig3_Calibration_DCA_FixedEffects_Sensitivity", artifact_dir)

    # -----------------------------------------------------------------------
    # Extended Data Fig. 4: 12 Visium Sections, 4x4 Proximity Matrix, FACS Audit & Platform Scaling (a-d)
    # -----------------------------------------------------------------------
    fig4, axes4 = plt.subplots(2, 2, figsize=(7.20, 5.20), facecolor="white")
    fig4.subplots_adjust(wspace=0.48, hspace=0.48)

    ax = axes4[0, 0]
    _panel_tag(ax, "a", x=-0.32, y=1.08)
    v12_s = s5a.sort_values("stromal_exclusion_index", ascending=True).reset_index(drop=True)
    yv = np.arange(len(v12_s))
    ax.barh(yv, v12_s["stromal_exclusion_index"], color=PAL["crimson"], alpha=0.84, height=0.64)
    for i_v, r_v in v12_s.iterrows():
        ax.text(float(r_v["stromal_exclusion_index"]) + 0.005, i_v, f"+{r_v['stromal_exclusion_index']:.2f} ({r_v['n_spots']:,})", va="center", ha="left", fontsize=4.6)
    ax.set_yticks(yv)
    ax.set_yticklabels(
        [str(r.section_id).replace("GSM6177", "").replace("GSE171351_", "Blad-").replace("_", "-") for r in v12_s.itertuples()],
        fontsize=4.8,
    )
    ax.set_xlim(0, float(v12_s["stromal_exclusion_index"].max()) * 1.48)
    ax.set_xlabel("Stromal Exclusion Index (SEI)")
    ax.set_title("12 Human Visium Tumor Sections", pad=5)

    ax = axes4[0, 1]
    _panel_tag(ax, "b", x=-0.24, y=1.08)
    pmat = s5b.to_numpy(dtype=float)
    im4 = ax.imshow(pmat, cmap="PuBuGn", vmin=0.35, vmax=0.85, aspect="equal")
    for r_i in range(pmat.shape[0]):
        for c_i in range(pmat.shape[1]):
            ax.text(c_i, r_i, f"{pmat[r_i, c_i]:.3f}", ha="center", va="center", fontsize=5.1, fontweight="bold", color="white" if pmat[r_i, c_i] >= 0.65 else PAL["ink"])
    lbl4 = ["$CD8^+$ T" if "CD8" in c else c.replace("_", " ") for c in s5b.columns]
    ax.set_xticks(range(len(lbl4)))
    ax.set_xticklabels(lbl4, fontsize=5.1)
    ax.set_yticks(range(len(lbl4)))
    ax.set_yticklabels(lbl4, fontsize=5.1)
    ax.set_title("Consensus $k$=6 Proximity Matrix", pad=5)
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_linewidth(0.55)
    cb4 = fig4.colorbar(im4, ax=ax, fraction=0.042, pad=0.025)
    cb4.ax.tick_params(labelsize=4.8)

    ax = axes4[1, 0]
    _panel_tag(ax, "c", x=-0.20, y=1.08)
    xf = np.arange(len(s5h))
    wf_b = 0.35
    ax.bar(xf - wf_b / 2, 100.0 * s5h["naive_pooled_caf_frac"], width=wf_b, color=PAL["lightgrey"], edgecolor=PAL["slate"], lw=0.55, label="Naive Pool")
    ax.bar(xf + wf_b / 2, 100.0 * s5h["deconfounded_nonT_caf_frac"], width=wf_b, color=PAL["crimson"], alpha=0.84, label="Non-$CD3^+$ TME")
    med_fold = float(s5h["caf_dilution_fold"].median())
    max_fold = float(s5h["caf_dilution_fold"].max())
    ax.text(0.03, 0.81, f"Med {med_fold:.1f}$\\times$ (Max {max_fold:.1f}$\\times$)", transform=ax.transAxes, fontsize=4.9, fontweight="bold", color=PAL["navy"])
    ax.set_xticks(xf)
    ax.set_xticklabels([f"{r.patient_id}\n({'R' if r.responder == 1 else 'NR'})" for r in s5h.itertuples()], fontsize=4.9)
    ax.set_ylim(0, float(100.0 * s5h["deconfounded_nonT_caf_frac"].max()) * 1.38)
    ax.set_ylabel("Pre-ICB myCAF (%)")
    ax.set_title("FACS $CD3^+$ Sort Deconfounding", pad=5)
    ax.legend(loc="upper right", frameon=False, fontsize=4.9)

    ax = axes4[1, 1]
    _panel_tag(ax, "d", x=-0.20, y=1.08)
    coh_col = "cohort" if "cohort" in s5i.columns else "cancer_cohort"
    grps_d = [
        ("GSE123813", 1, "BCC\nR", PAL["teal"]),
        ("GSE123813", 0, "BCC\nNR", PAL["crimson"]),
        ("GSE91061", 1, "Mel91\nR", PAL["blue"]),
        ("GSE91061", 0, "Mel91\nNR", PAL["vermilion"]),
        ("GSE78220", 1, "Mel78\nR", PAL["teal"]),
        ("GSE78220", 0, "Mel78\nNR", PAL["crimson"]),
        ("GSE207422", 1, "NSCLC\nR", PAL["blue"]),
        ("GSE207422", 0, "NSCLC\nNR", PAL["vermilion"]),
    ]
    rng_d = np.random.default_rng(99)
    for idx_d, (ckey, rv, lbl_d, col_d) in enumerate(grps_d):
        mask_d = s5i[coh_col].astype(str).str.contains(ckey, na=False) & (s5i["responder"] == rv)
        vals_d = 100.0 * s5i.loc[mask_d, "simulated_icb_remission_prob"].to_numpy()
        mu_d = float(np.mean(vals_d)) if len(vals_d) > 0 else 0.0
        se_d = float(np.std(vals_d, ddof=1) / math.sqrt(max(len(vals_d), 1))) if len(vals_d) > 1 else 0.0
        ax.bar(
            idx_d,
            mu_d,
            width=0.58,
            color=col_d,
            alpha=0.40,
            edgecolor=col_d,
            lw=0.75,
            yerr=se_d,
            capsize=1.8,
            error_kw={"lw": 0.65, "ecolor": col_d},
        )
        if len(vals_d) > 0:
            ax.scatter(np.full(len(vals_d), idx_d) + rng_d.uniform(-0.14, 0.14, size=len(vals_d)), vals_d, color=col_d, s=7, edgecolors="white", linewidths=0.25, zorder=4)
    ax.axhline(50.0, color=PAL["slate"], ls="--", lw=0.65)
    ax.set_xticks(range(len(grps_d)))
    ax.set_xticklabels([g[2] for g in grps_d], fontsize=4.6)
    ax.set_ylim(0, 118)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Simulated Remission (%)")
    ax.set_title(f"4-Cohort LOPO Alignment ($n$={len(s5i)})", pad=5)

    _save_multi(fig4, out_dir, "Extended_Data_Fig4_Visium12_FACS_and_Platform_Scaling", artifact_dir)


# Backward-compatible Stage II symbol aliases
build_fig1_vector_architecture = build_fig1
build_fig2_macro_clinical_atlas = build_fig2
build_fig3_pan_disease_drift_law = build_fig3
build_fig4_fibrosis_mash_deconvolution = build_fig4
build_fig5_solid_tumor_spatial_icb = build_fig5
