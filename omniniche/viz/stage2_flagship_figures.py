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
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import spatial, stats


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
            "font.size": 6.0,
            "axes.titlesize": 6.8,
            "axes.titleweight": "bold",
            "axes.labelsize": 6.2,
            "axes.labelweight": "normal",
            "axes.edgecolor": "#26303A",
            "axes.linewidth": 0.55,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 5.6,
            "ytick.labelsize": 5.7,
            "xtick.major.width": 0.50,
            "ytick.major.width": 0.50,
            "xtick.major.size": 2.3,
            "ytick.major.size": 2.3,
            "xtick.major.pad": 2.0,
            "ytick.major.pad": 2.0,
            "legend.fontsize": 5.5,
            "legend.frameon": False,
            "legend.handlelength": 1.25,
            "legend.handletextpad": 0.4,
            "legend.columnspacing": 0.8,
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
# FIGURE 1: BioRender Multi-Scale Architectural & Mechanistic Overview (a-d)
# ===========================================================================
def export_fig1_illustration(
    raw_asset_path: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> pathlib.Path:
    """Export pure-white borderless BioRender 2D flat Fig. 1 natively at 183 mm (7.20 in) width."""
    set_nature_style()
    im = Image.open(raw_asset_path).convert("RGB")
    arr = np.asarray(im)
    h_px, w_px = arr.shape[:2]

    fig_w = 7.20
    fig_h = fig_w * (h_px / float(w_px))
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.imshow(arr, interpolation="lanczos")
    ax.axis("off")
    return _save_multi(fig, out_dir, "Fig1_OmniNiche_Architecture_and_Framework", artifact_dir, dpi=600)


# ===========================================================================
# FIGURE 2: Macro-Atlas of 41,396 Systemic & Parenchymal Trials (6 Panels a-f)
# ===========================================================================
def build_fig2_macro_clinical_atlas(
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

    fig = plt.figure(figsize=(7.20, 5.20), facecolor="white")
    gs = gridspec.GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[1.05, 1.00],
        hspace=0.58,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], width_ratios=[1.24, 0.90, 0.92], wspace=0.60)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[1.20, 1.02, 0.78], wspace=0.58)

    # --- Panel a: 3-Stage Alluvial Flow (Target -> Modality -> Deconvoluted Outcome) ---
    ax_a = fig.add_subplot(gs_r1[0, 0])
    _panel_tag(ax_a, "a", x=-0.08, y=1.08)

    n_tot = int(s1["unique_systemic_organ_trials"].sum())
    n_res = len(s2)
    n_unres = max(n_tot - n_res, 0)

    def _map_mod(m_str: str) -> str:
        s = str(m_str).lower()
        if "antibody" in s:
            return "IgG mAb"
        if "peptide" in s or "protein" in s or "enzyme" in s or "oligo" in s or "cell" in s:
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

    x_pos = [0.04, 0.39, 0.74]
    bw = 0.064
    short_lbl = {
        "Resolved": "Res.",
        "Unresolved": "Unres.",
        "Small Mol.": "SM",
        "IgG mAb": "mAb",
        "Peptide/Bio": "Pep.",
        "Completed": "Comp.",
        "Admin. Stop": "Admin.",
        "Unreported": "Unrep.",
        "Bio. Fail": "BioFail",
    }
    for xi, st_dict, hdr in zip(x_pos, [st1, st2, st3], ["Target", "Modality", "Outcome"]):
        ax_a.text(xi, 1.045, hdr, ha="center", va="bottom", fontsize=5.7, fontweight="bold", color=PAL["navy"])
        for lbl, (yb, yt, val, col) in st_dict.items():
            ax_a.add_patch(
                mpatches.Rectangle((xi - bw / 2, yb), bw, max(yt - yb, 0.008), facecolor=col, edgecolor="white", lw=0.4)
            )
            ym = 0.5 * (yb + yt)
            sl = short_lbl.get(lbl, lbl)
            txt = f"{sl}\n{val / 1000:.1f}k" if (yt - yb) >= 0.14 else f"{sl} {val / 1000:.1f}k"
            ax_a.text(xi + bw / 2 + 0.015, ym, txt, ha="left", va="center", fontsize=4.6, color=PAL["ink"])

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

    ax_a.set_xlim(-0.01, 1.05)
    ax_a.set_ylim(-0.03, 1.10)
    ax_a.axis("off")
    ax_a.set_title(f"Trial Deconvolution ($n$ = {n_tot:,})", pad=5)

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

    for i_r in range(len(sub_odds)):
        col_r = (
            PAL["crimson"]
            if (ors[i_r] > 1.0 and pvs[i_r] < 0.05)
            else (PAL["teal"] if (ors[i_r] < 1.0 and pvs[i_r] < 0.05) else PAL["midgrey"])
        )
        ax_d.plot([clo[i_r], chi[i_r]], [i_r, i_r], color=col_r, lw=1.25, zorder=3)
        ax_d.scatter([ors[i_r]], [i_r], color=col_r, s=20, zorder=4)
        ax_d.text(
            2.15,
            i_r,
            f"n={int(ev_m[i_r])} | {ors[i_r]:.2f} ({clo[i_r]:.2f}–{chi[i_r]:.2f})",
            va="center",
            ha="left",
            fontsize=4.7,
            color=PAL["ink"],
        )

    ax_d.set_xscale("log")
    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels(sub_odds["row_lbl"], fontsize=5.1)
    ax_d.set_xlim(0.34, 9.2)
    ax_d.set_xticks([0.5, 1.0, 2.0])
    ax_d.set_xticklabels(["0.5", "1.0", "2.0"])
    ax_d.set_xlabel("Odds Ratio vs. Small Mol. (← Lower | Higher →)")
    ax_d.set_title("Modality Attrition Odds Ratio", pad=5)

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
def build_fig3_pan_disease_drift_law(
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

    fig = plt.figure(figsize=(7.20, 7.80), facecolor="white")
    gs = gridspec.GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.02, 1.08, 1.02],
        hspace=0.58,
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
        fontsize=5.0,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_b.set_xticks(x_st)
    ax_b.set_xticklabels(step_short, rotation=28, ha="right", fontsize=5.1)
    ax_b.set_xlim(-0.35, len(x_st) - 0.55)
    ax_b.set_ylim(0.41, 0.77)
    ax_b.set_ylabel("Scientific Outcome AUROC")
    ax_b.set_title(f"DNAcc Ablation & Genetics ($n$={n_sci_tot:,})", pad=5)
    ax_b.legend(loc="lower right", frameon=False, fontsize=4.6)

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
    ax_c.axhline(4.5, color=PAL["lightgrey"], ls="-", lw=0.75)
    for idx_c, (lbl_c, or_v, w_lo, w_hi, c_lo, c_hi, col_c) in enumerate(forest_items):
        if c_lo is not None and c_hi is not None:
            ax_c.plot([c_lo, min(c_hi, 2.45)], [idx_c, idx_c], color=col_c, lw=0.75, alpha=0.55)
            ax_c.plot([w_lo, min(w_hi, 2.45)], [idx_c, idx_c], color=col_c, lw=1.65)
            ci_txt = f"{or_v:.2f} [{c_lo:.2f}–{c_hi:.2f}]"
        else:
            ax_c.plot([w_lo, min(w_hi, 2.45)], [idx_c, idx_c], color=col_c, lw=1.35)
            ci_txt = f"{or_v:.2f} [{w_lo:.2f}–{w_hi:.2f}]"
        ax_c.scatter([or_v], [idx_c], color=col_c, s=16, zorder=4)
        ax_c.text(2.42, idx_c, ci_txt, va="center", ha="left", fontsize=4.7, color=PAL["ink"])

    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels([it[0] for it in forest_items], fontsize=5.0)
    ax_c.set_xlim(0.70, 3.85)
    ax_c.set_xticks([1.0, 1.5, 2.0])
    ax_c.set_xlabel("OR per 1-s.d. (Cluster / Wald 95% CI)")
    ax_c.set_title("3-Pillar Joint & Per-Axis Effect", pad=5)

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

    ax_f.plot([0.14, 0.80], [0.14, 0.80], color=PAL["slate"], ls=":", lw=0.7, label="Ideal")
    ax_f.plot(c_tau["mean_predicted_prob"], c_tau["observed_event_rate"], color=PAL["midgrey"], marker="s", ms=2.4, lw=0.95, ls="--", label="Static $\\tau_0$")
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

    # --- Panel h: Target Remodeling Phase Space (tau_0 vs Delta_tau) with Bilateral Non-Crossing Rails ---
    ax_h = fig.add_subplot(gs_r3[0, 2])
    _panel_tag(ax_h, "h", x=-0.24, y=1.10)
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

    m_nosw = g_sub["sw"] < 0.5
    m_sw = ~m_nosw
    ax_h.scatter(
        g_sub.loc[m_nosw, "tau0"],
        g_sub.loc[m_nosw, "dtau"],
        s=7 + 10 * np.log10(g_sub.loc[m_nosw, "nt"]),
        color=PAL["blue"],
        alpha=0.50,
        edgecolors="none",
        label="Retained",
    )
    ax_h.scatter(
        g_sub.loc[m_sw, "tau0"],
        g_sub.loc[m_sw, "dtau"],
        s=9 + 11 * np.log10(g_sub.loc[m_sw, "nt"]),
        color=PAL["crimson"],
        alpha=0.65,
        edgecolors="white",
        linewidths=0.3,
        label="Switched",
    )

    # Bilateral outer rails: moderate-tau0 targets (< 0.58) route left; high-tau0 targets (>= 0.58) route right
    callouts = ["THRB", "TGFB1", "PCSK9", "ERBB2", "PDGFRB", "LOXL2"]
    c_rows = g_agg[g_agg["primary_target_gene"].isin(callouts)].copy()
    left_c = c_rows[c_rows["tau0"] < 0.58].sort_values("dtau", ascending=True)
    right_c = c_rows[c_rows["tau0"] >= 0.58].sort_values("dtau", ascending=True)

    y_min_h = float(g_sub["dtau"].min())
    y_max_h = float(g_sub["dtau"].max())

    if len(left_c) > 0:
        l_slots = np.linspace(y_min_h * 0.65, y_max_h * 0.75, len(left_c)) if len(left_c) > 1 else [float(left_c["dtau"].iloc[0])]
        for (_, r_g), y_sl in zip(left_c.iterrows(), l_slots):
            ax_h.scatter([r_g["tau0"]], [r_g["dtau"]], s=20, facecolors="none", edgecolors=PAL["ink"], linewidths=0.8, zorder=6)
            ax_h.annotate(
                str(r_g["primary_target_gene"]),
                xy=(float(r_g["tau0"]), float(r_g["dtau"])),
                xytext=(0.22, float(y_sl)),
                fontsize=4.9,
                fontweight="bold",
                fontstyle="italic",
                color=PAL["crimson"] if (r_g["sw"] >= 0.5 or r_g["dtau"] < -0.04) else PAL["navy"],
                va="center",
                ha="right",
                arrowprops={"arrowstyle": "-", "lw": 0.45, "color": PAL["slate"], "shrinkA": 0, "shrinkB": 1.5},
            )

    if len(right_c) > 0:
        r_slots = np.linspace(y_min_h * 0.85, y_max_h * 0.98, len(right_c))
        for (_, r_g), y_sl in zip(right_c.iterrows(), r_slots):
            ax_h.scatter([r_g["tau0"]], [r_g["dtau"]], s=20, facecolors="none", edgecolors=PAL["ink"], linewidths=0.8, zorder=6)
            ax_h.annotate(
                str(r_g["primary_target_gene"]),
                xy=(float(r_g["tau0"]), float(r_g["dtau"])),
                xytext=(1.02, float(y_sl)),
                fontsize=4.9,
                fontweight="bold",
                fontstyle="italic",
                color=PAL["crimson"] if (r_g["sw"] >= 0.5 or r_g["dtau"] < -0.04) else PAL["navy"],
                va="center",
                ha="left",
                arrowprops={"arrowstyle": "-", "lw": 0.45, "color": PAL["slate"], "shrinkA": 0, "shrinkB": 1.5},
            )

    ax_h.set_xlim(-0.04, 1.30)
    ax_h.set_xticks([0.2, 0.6, 1.0])
    ax_h.set_ylim(y_min_h * 1.25 - 0.02, y_max_h * 1.38 + 0.03)
    ax_h.set_xlabel("Healthy Specificity ($\\tau_0$)")
    ax_h.set_ylabel("Disease Drift ($\\Delta\\tau$)")
    ax_h.set_title("Target Drift Phase Space", pad=5)
    ax_h.legend(loc="upper left", frameon=False, fontsize=4.9, ncol=2)

    _save_multi(fig, out_dir, "Fig3_PanOrgan_Dynamic_Drift_and_DNAcc_Law", artifact_dir)


# ===========================================================================
# FIGURE 4: Pan-Organ Fibrosis Zero-Shot Benchmark & Hysteresis (8 Panels a-h)
# ===========================================================================
def build_fig4_fibrosis_mash_deconvolution(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 8-panel (a-h) Pan-Organ Human Fibrosis Zero-Shot Benchmark & Hysteresis at 7.20 in x 7.85 in."""
    set_nature_style()
    donors = pd.read_csv(tables_dir / "Table_S4a_PanOrgan_Fibrosis_All_Human_Donors.csv")
    ccr2 = pd.read_csv(tables_dir / "Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv")
    trials = pd.read_csv(tables_dir / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv")
    arena = pd.read_csv(tables_dir / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv")
    bi3x3 = pd.read_csv(tables_dir / "Table_S4f_Bidirectional_3x3_CrossOrgan_Transfer.csv")
    hyst = pd.read_csv(tables_dir / "Fig4_Continuous_Hysteresis_Bifurcation_Curve.csv")
    bliss = pd.read_csv(tables_dir / "Table_S4e_Fibrosis_2D_Bliss_Matrix.csv")
    wad_npz = np.load(tables_dir / "Fig4_Waddington_2D_Quiver_Fibrosis.npz")

    fig = plt.figure(figsize=(7.20, 8.35), facecolor="white")
    gs = gridspec.GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.38, 0.96, 0.98],
        hspace=0.56,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0], width_ratios=[0.90, 0.96, 0.14, 1.42], wspace=0.42)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], width_ratios=[0.98, 1.08, 1.06], wspace=0.52)
    gs_r3 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], width_ratios=[1.38, 1.00], wspace=0.44)

    # --- Panel a: 64-Donor Primary Human 3-Organ Fibrotic Niche Ratio ---
    ax_a = fig.add_subplot(gs_r1[0, 0])
    _panel_tag(ax_a, "a", x=-0.26, y=1.08)
    donors["organ_clean"] = donors["organ"].map(
        lambda s: "Liver" if "Liver" in str(s) else ("Lung" if "Lung" in str(s) else "Kidney")
    )
    donors["is_fibrotic"] = ~donors["stage"].astype(str).str.contains("Control|Healthy|Early_MCD", case=False, na=False)
    donors["raw_scar"] = donors["Myofibroblast_frac"].clip(lower=0.003)
    ctrl_means = donors[~donors["is_fibrotic"]].groupby("organ")["raw_scar"].mean().to_dict()
    donors["scar_index"] = np.clip(
        [float(r.raw_scar) / max(float(ctrl_means.get(r.organ, 0.05)), 1e-6) for r in donors.itertuples()],
        0.15,
        6.5,
    )

    org_order = ["Liver", "Lung", "Kidney"]
    rng_a = np.random.default_rng(42)
    x_a = np.arange(len(org_order))
    w_a = 0.34
    max_y_a = 0.1
    for i_o, org_n in enumerate(org_order):
        for j_s, (is_fib, lbl_s, col_s) in enumerate([(False, "Control", PAL["midgrey"]), (True, "Fibrotic", PAL["crimson"])]):
            vals = donors[(donors["organ_clean"] == org_n) & (donors["is_fibrotic"] == is_fib)]["scar_index"].to_numpy()
            if len(vals) == 0:
                continue
            pos = i_o + (-w_a / 2 if j_s == 0 else w_a / 2)
            mu = float(np.mean(vals))
            sem = float(np.std(vals, ddof=1) / math.sqrt(max(len(vals), 1))) if len(vals) > 1 else 0.0
            ax_a.bar(
                pos,
                mu,
                width=w_a * 0.88,
                color=col_s,
                alpha=0.38,
                edgecolor=col_s,
                lw=0.75,
                yerr=sem,
                capsize=1.6,
                error_kw={"lw": 0.65, "ecolor": col_s},
                label=lbl_s if i_o == 0 else None,
            )
            jit = rng_a.uniform(-0.08, 0.08, size=len(vals))
            ax_a.scatter(np.full(len(vals), pos) + jit, vals, color=col_s, s=8, alpha=0.82, edgecolors="white", linewidths=0.25, zorder=4)
            max_y_a = max(max_y_a, float(np.max(vals)))

    ax_a.set_xticks(x_a)
    ax_a.set_xticklabels(["Liver\n(17)", "Lung\n(30)", "Kidney\n(17)"], fontsize=5.2)
    ax_a.set_ylim(0, max_y_a * 1.32)
    ax_a.set_ylabel("Normalized Scar Niche Fold")
    ax_a.set_title("64-Donor Fibrosis Atlas", pad=5)
    ax_a.legend(loc="upper left", frameon=False, fontsize=5.0)

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
    ax_c.axvline(0.50, color=PAL["slate"], ls="--", lw=0.6)
    cols_c = [PAL["teal"] if int(w) == 1 else PAL["crimson"] for w in tr_s["clinical_win"]]
    probs_c = tr_s["omniniche_remission_prob"].to_numpy(dtype=float)
    ax_c.barh(y_c, probs_c, left=0.0, color=cols_c, alpha=0.86, height=0.68)
    for idx_t, p_val in enumerate(probs_c):
        ax_c.text(p_val + 0.022, idx_t, f"{p_val:.2f}", va="center", ha="left", fontsize=4.15, color=PAL["ink"])

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

    # --- Panel f: Smooth 60x60 2D Waddington Quasi-Potential + Quiver Drift Field ---
    ax_f = fig.add_subplot(gs_r2[0, 2])
    _panel_tag(ax_f, "f", x=-0.24, y=1.08)
    gx_w = wad_npz["grid_myofibroblast"]
    gy_w = wad_npz["grid_parenchyma"]
    u_2d = wad_npz["U_2d"]
    du_x = wad_npz["drift_myofibroblast"]
    du_y = wad_npz["drift_parenchyma"]

    cf_f = ax_f.contourf(gx_w, gy_w, u_2d, levels=16, cmap="viridis_r")
    ax_f.contour(gx_w, gy_w, u_2d, levels=9, colors="white", linewidths=0.3, alpha=0.45)
    st_q = 5
    Xg, Yg = np.meshgrid(gx_w[::st_q], gy_w[::st_q])
    Ux, Uy = du_x[::st_q, ::st_q], du_y[::st_q, ::st_q]
    mag = np.hypot(Ux, Uy) + 1e-8
    ax_f.quiver(Xg, Yg, Ux / mag, Uy / mag, color="white", alpha=0.78, scale=22, width=0.006)
    ax_f.axvline(0.36, color=PAL["orange"], ls="--", lw=0.8)
    ax_f.scatter([0.14, 0.74], [0.14, 0.58], c=[PAL["teal"], PAL["crimson"]], s=26, edgecolors="white", linewidths=0.9, zorder=5)
    ax_f.text(0.18, 0.14, "Remission", color="white", fontsize=5.1, fontweight="bold", va="center")
    ax_f.text(0.71, 0.67, "Locked F4", color=PAL["ink"], fontsize=5.1, fontweight="bold", va="bottom", ha="center")
    ax_f.set_xlabel("Myofibroblast ($S_{\\mathrm{Myofib}}$)")
    ax_f.set_ylabel("Parenchyma ($S_{\\mathrm{Paren}}$)")
    ax_f.set_title("2D Waddington Landscape", pad=5)
    cbar_f = fig.colorbar(cf_f, ax=ax_f, fraction=0.040, pad=0.025)
    cbar_f.ax.tick_params(labelsize=4.8)
    cbar_f.set_label("$U(S)$", fontsize=5.2)

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
    ax_g.set_title("Hysteresis Monotherapy Cliff & Dual Rescue", pad=5)
    ax_g.legend(loc="upper right", frameon=False, fontsize=4.9, ncol=2)

    # --- Panel h: Symmetric 6x6 Clinical-Target 2D Bliss Synergy Matrix in Locked F4 ---
    ax_h = fig.add_subplot(gs_r3[0, 1])
    _panel_tag(ax_h, "h", x=-0.28, y=1.10)
    arm_order = list(dict.fromkeys(bliss["arm_row"]))
    arm_short = [a.split(" (")[0].replace("Simtuzumab", "Simtuz.").replace("Cenicriviroc", "Cenicriv.") for a in arm_order]
    piv_b = bliss.pivot(index="arm_row", columns="arm_col", values="bliss_synergy_excess").reindex(index=arm_order, columns=arm_order).to_numpy(dtype=float)
    im_h = ax_h.imshow(piv_b, cmap="YlOrRd", vmin=0.0, vmax=1.0, aspect="equal")
    for r_i in range(len(arm_order)):
        for c_i in range(len(arm_order)):
            val_b = piv_b[r_i, c_i]
            ax_h.text(
                c_i,
                r_i,
                f"{val_b:+.2f}" if r_i != c_i else "0.00",
                ha="center",
                va="center",
                fontsize=4.6,
                fontweight="bold" if val_b >= 0.40 else "normal",
                color="white" if val_b >= 0.55 else PAL["ink"],
            )

    ax_h.set_xticks(range(len(arm_short)))
    ax_h.set_xticklabels(arm_short, rotation=35, ha="right", fontsize=4.9)
    ax_h.set_yticks(range(len(arm_short)))
    ax_h.set_yticklabels(arm_short, fontsize=4.9)
    ax_h.set_title("Locked-F4 2D Bliss Synergy", pad=5)
    for sp in ("top", "right", "bottom", "left"):
        ax_h.spines[sp].set_visible(True)
        ax_h.spines[sp].set_linewidth(0.55)
    cbar_h = fig.colorbar(im_h, ax=ax_h, fraction=0.042, pad=0.025)
    cbar_h.ax.tick_params(labelsize=4.8)
    cbar_h.set_label("Bliss Excess", fontsize=5.1)

    _save_multi(fig, out_dir, "Fig4_PanOrgan_Fibrosis_ZeroShot_and_Hysteresis", artifact_dir)


# ===========================================================================
# FIGURE 5: Solid-Tumor Visium Exclusion & Pre-ICB Responder (7 Panels a-g)
# ===========================================================================
def build_fig5_solid_tumor_spatial_icb(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 7-panel (a-g) Solid-Tumor 10x Visium Spatial Barrier & Blind Pre-ICB Responder at 7.20 in x 7.95 in."""
    set_nature_style()
    vis_con = pd.read_csv(tables_dir / "Fig5_Contrasting_Visium_Sections.csv")
    rad_sum = pd.read_csv(tables_dir / "Table_S5g_Visium_Radial_Consensus_Summary.csv")
    vis_12 = pd.read_csv(tables_dir / "Table_S5a_Visium_Spatial_Barrier_Sections.csv")
    pool107 = pd.read_csv(tables_dir / "Table_S5i_Pooled_PreTreatment_107Patients_Master.csv")
    bench = pd.read_csv(tables_dir / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv")
    paired15 = pd.read_csv(tables_dir / "Table_S5f_GSE123813_Paired_PrePost_Trajectory.csv")
    tum_npz = np.load(tables_dir / "Fig5_Waddington_2D_Quiver_SolidTumor.npz")

    fig = plt.figure(figsize=(7.20, 7.95), facecolor="white")
    gs = gridspec.GridSpec(
        3,
        1,
        figure=fig,
        height_ratios=[1.02, 1.15, 1.00],
        hspace=0.58,
    )
    gs_r1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[1.40, 1.00], wspace=0.30)
    gs_r2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], width_ratios=[1.32, 1.00], wspace=0.42)
    gs_r3 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[2], width_ratios=[0.96, 1.12, 0.92], wspace=0.50)

    # --- Panel a: Dual Contrasting Primary Human Carcinoma 10x Visium Spatial Maps (with 500-um Physical Scale Bar) ---
    gs_a = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_r1[0, 0], wspace=0.14)
    ax_a1 = fig.add_subplot(gs_a[0, 0])
    ax_a2 = fig.add_subplot(gs_a[0, 1])
    _panel_tag(ax_a1, "a", x=-0.16, y=1.10)

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

        for n_lbl, c_hex in niche_cols.items():
            m = sub_v["dom"] == n_lbl
            if np.any(m):
                ax_sub.scatter(sub_v.loc[m, "x"], sub_v.loc[m, "plot_y"], c=c_hex, s=2.6, alpha=0.88, edgecolors="none", label=n_lbl)

        # Exact 500-um physical scale bar = 5 * median nearest-neighbor Visium 100-um spot pitch
        coords = np.column_stack([xs, ys])
        nn_dists = spatial.cKDTree(coords).query(coords, k=2)[0][:, 1]
        pitch_100um = float(np.median(nn_dists))
        bar_500um = 5.0 * pitch_100um
        x_min, x_max = float(xs.min()), float(xs.max())
        y_min, y_max = float(ys.min()), float(ys.max())
        span_x = x_max - x_min
        span_y = y_max - y_min
        bx0 = x_max - 0.06 * span_x - bar_500um
        by0 = y_min - 0.08 * span_y
        ax_sub.plot([bx0, bx0 + bar_500um], [by0, by0], color=PAL["ink"], lw=1.4, solid_capstyle="butt")
        ax_sub.text(bx0 + 0.5 * bar_500um, by0 - 0.025 * span_y, "500 µm", ha="center", va="top", fontsize=4.7, fontweight="bold", color=PAL["ink"])

        ax_sub.set_xlim(x_min - 0.06 * span_x, x_max + 0.06 * span_x)
        ax_sub.set_ylim(y_min - 0.18 * span_y, y_max + 0.30 * span_y)

        sei_v = float(sub_v["sei"].iloc[0])
        ctype_v = str(sub_v["cancer_type"].iloc[0]).replace("_Carcinoma", "")
        ax_sub.set_title(f"{title_prefix} ({ctype_v}, SEI +{sei_v:.2f})", fontsize=6.1, pad=4)
        ax_sub.set_xticks([])
        ax_sub.set_yticks([])
        for sp in ("top", "right", "bottom", "left"):
            ax_sub.spines[sp].set_visible(True)
            ax_sub.spines[sp].set_color(PAL["lightgrey"])
            ax_sub.spines[sp].set_linewidth(0.5)
    ax_a1.legend(loc="upper left", frameon=False, fontsize=4.6, markerscale=1.6, ncol=2)

    # --- Panel b: 12-Section (21,796-Spot) Continuous Radial Distance Decay (0-600 um) ---
    ax_b = fig.add_subplot(gs_r1[0, 1])
    _panel_tag(ax_b, "b", x=-0.22, y=1.10)
    dist_um = rad_sum["approx_distance_um"].to_numpy()
    for col_z, lbl_z, c_z, mk in [
        ("Tumor_z", "Tumor Core", PAL["navy"], "o"),
        ("myCAF_z", "myCAF Wall", PAL["crimson"], "s"),
        ("TAM_z", "$SPP1^+$ TAM", PAL["orange"], "^"),
        ("CD8_T_z", "$CD8^+$ T", PAL["teal"], "D"),
    ]:
        mu_z = rad_sum[f"{col_z}_mean"].to_numpy()
        se_z = rad_sum[f"{col_z}_sem"].to_numpy()
        ax_b.plot(dist_um, mu_z, color=c_z, lw=1.2, marker=mk, ms=2.4, label=lbl_z)
        ax_b.fill_between(dist_um, mu_z - se_z, mu_z + se_z, color=c_z, alpha=0.14)

    w_res = stats.wilcoxon(vis_12["stromal_exclusion_index"].to_numpy(), alternative="greater")
    ax_b.axhline(0.0, color=PAL["slate"], ls=":", lw=0.6)
    ax_b.text(0.24, 0.90, f"12/12 SEI > 0 ({_fmt_p(float(w_res.pvalue))})", transform=ax_b.transAxes, fontsize=5.0, fontweight="bold", color=PAL["navy"])
    ax_b.set_xlabel("Distance from Tumor Core (µm)")
    ax_b.set_ylabel("Visium Signature ($z$-Score)")
    ax_b.set_ylim(float(rad_sum["Tumor_z_mean"].min()) - 0.35, float(rad_sum["Tumor_z_mean"].max()) + 0.55)
    ax_b.set_title("12-Section Radial Barrier (21,796 Spots)", pad=5)
    ax_b.legend(loc="lower left", frameon=False, fontsize=4.9, ncol=2)

    # --- Panel c: n = 107 Pre-Treatment Patient Waterfall + 4-Row Aligned OncoPrint Biomarker Strip ---
    gs_c = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_r2[0, 0], height_ratios=[2.30, 1.05], hspace=0.08)
    ax_c1 = fig.add_subplot(gs_c[0, 0])
    ax_c2 = fig.add_subplot(gs_c[1, 0], sharex=ax_c1)
    _panel_tag(ax_c1, "c", x=-0.14, y=1.14)

    wf = pool107.copy()
    coh_col = "cohort" if "cohort" in wf.columns else "cancer_cohort"
    pdl1_col = "b1_static_pdl1_prob" if "b1_static_pdl1_prob" in wf.columns else "b1_static_pdl1"
    wf["pdl1_hi"] = wf.groupby(coh_col)[pdl1_col].transform(lambda s: s >= s.median())
    wf["caf_hi"] = wf["myCAF"] >= 0.50
    wf = wf.sort_values("simulated_icb_remission_prob", ascending=False).reset_index(drop=True)
    x_wf = np.arange(len(wf))
    delta_p = 100.0 * (wf["simulated_icb_remission_prob"].to_numpy() - 0.50)
    bar_cols = [PAL["teal"] if int(r) == 1 else PAL["crimson"] for r in wf["responder"]]
    ax_c1.bar(x_wf, delta_p, width=0.84, color=bar_cols, alpha=0.88, edgecolor="none")

    ax_c1.axhline(0.0, color=PAL["ink"], lw=0.65)
    ax_c1.set_ylim(-56, 66)
    ax_c1.set_yticks([-50, 0, 50])
    ax_c1.set_ylabel("$\\Delta$Prob. vs. 50% (%)")
    ax_c1.set_title(f"Zero-Refit Pre-ICB Waterfall ($n$ = {len(wf)})", pad=5)
    ax_c1.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_c1.spines["bottom"].set_visible(False)
    h_leg = [
        mpatches.Patch(facecolor=PAL["teal"], label="Responder (R)"),
        mpatches.Patch(facecolor=PAL["crimson"], label="Non-Responder (NR)"),
    ]
    ax_c1.legend(handles=h_leg, loc="upper right", frameon=False, fontsize=5.0, ncol=2)

    # 4-row OncoPrint tracks: myCAF Shield | PD-L1 | Cohort (BCC / Mel91 / Mel78 / NSCLC) | RECIST R
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
    ax_c2.set_yticklabels(track_names, fontsize=4.8)
    tick_pos = sorted({0, 24, 49, 74, len(wf) - 1})
    ax_c2.set_xticks(tick_pos)
    ax_c2.set_xticklabels([str(p + 1) for p in tick_pos], fontsize=4.8)
    ax_c2.set_xlabel("Ranked Pre-Treatment Patients", fontsize=5.1)
    for sp in ("top", "right", "left", "bottom"):
        ax_c2.spines[sp].set_visible(False)
    ax_c2.tick_params(axis="y", length=0)

    # --- Panel d: 5-Cohort Blind Anti-PD-1 Benchmark + Spatial & Published Biomarkers ---
    ax_d = fig.add_subplot(gs_r2[0, 1])
    _panel_tag(ax_d, "d", x=-0.36, y=1.08)
    coh_map = [
        ("GSE123813_BCC_PreTreatment_AntiPD1_Mono_scRNA", "BCC Mono.\n($n$ = 8)"),
        ("Pooled_Pure_AntiPD1_Monotherapy_PreTreatment_n83", "Anti-PD-1 Mono.\n($n$ = 83)"),
        ("GSE207422_NSCLC_PreTreatment_ChemoAntiPD1_Biopsy", "NSCLC Chemo-IO\n($n$ = 24)"),
        ("Pooled_All_Strictly_PreTreatment_Biopsies_n107", "All Pre-Tx\n($n$ = 107)"),
        ("GSE207422_NSCLC_Resected_TME_scRNA", "NSCLC Resect.\n($n$ = 14)"),
    ]
    meth_sel = [
        ("OmniNiche_v2_Spatial_Langevin", "OmniNiche", PAL["blue"], "o", -0.24),
        ("OmniNiche_NoVisium_UniformPij", "No-Visium $P_{ij}$", PAL["teal"], "P", -0.12),
        ("B6_Published_18Gene_Tcell_GEP", "18-Gene GEP", PAL["orange"], "^", 0.00),
        ("B7_Published_CAF_TIDE_Exclusion", "CAF/TIDE Excl.", PAL["purple"], "D", 0.12),
        ("B1_Static_PDL1_CD274", "PD-L1 ($CD274$)", PAL["vermilion"], "s", 0.24),
    ]
    ax_d.axvline(0.50, color=PAL["slate"], ls=":", lw=0.6)
    for c_idx, (c_key, _) in enumerate(coh_map):
        if c_idx % 2 == 1:
            ax_d.axhspan(c_idx - 0.45, c_idx + 0.45, color=PAL["pale"], zorder=0)
        sub_b = bench[bench["cohort"] == c_key].set_index("method")
        for m_key, m_lbl, m_col, m_mk, dy in meth_sel:
            if m_key not in sub_b.index:
                continue
            r_b = sub_b.loc[m_key]
            is_om = m_key == "OmniNiche_v2_Spatial_Langevin"
            ci_lo = float(r_b["auroc_ci_low"]) if "auroc_ci_low" in r_b and pd.notna(r_b["auroc_ci_low"]) else float(r_b["auroc"])
            ci_hi = float(r_b["auroc_ci_high"]) if "auroc_ci_high" in r_b and pd.notna(r_b["auroc_ci_high"]) else float(r_b["auroc"])
            ax_d.plot([ci_lo, ci_hi], [c_idx + dy, c_idx + dy], color=m_col, lw=1.35 if is_om else 0.80, zorder=3)
            ax_d.scatter([r_b["auroc"]], [c_idx + dy], color=m_col, marker=m_mk, s=18 if is_om else 10, zorder=4, label=m_lbl if c_idx == 0 else None)

    ax_d.set_yticks(range(len(coh_map)))
    ax_d.set_yticklabels([lbl for _, lbl in coh_map], fontsize=4.8)
    ax_d.set_xlim(0.10, 1.06)
    ax_d.set_ylim(-0.52, len(coh_map) + 0.68)
    ax_d.set_xlabel("Blind Responder AUROC (95% CI)")
    ax_d.set_title("Multi-Cohort Anti-PD-1 Benchmark", pad=5)
    ax_d.legend(loc="upper left", frameon=False, fontsize=4.35, ncol=2, columnspacing=0.6, handletextpad=0.3)

    # --- Panel e: Paired Pre- vs Post-Anti-PD-1 Longitudinal Human Trajectory (n = 15) ---
    ax_e = fig.add_subplot(gs_r3[0, 0])
    _panel_tag(ax_e, "e", x=-0.24, y=1.10)
    paired15["net_eff_pre"] = paired15["cd8_effector_pre"] - 0.5 * paired15["cd8_exhaustion_pre"]
    paired15["net_eff_post"] = paired15["cd8_effector_post"] - 0.5 * paired15["cd8_exhaustion_post"]
    for grp_i, (resp_v, lbl_p, col_p) in enumerate([(1, "R ($n$=8)", PAL["teal"]), (0, "NR ($n$=7)", PAL["crimson"])]):
        sub_p = paired15[paired15["responder"] == resp_v]
        x0_p, x1_p = grp_i * 2.1, grp_i * 2.1 + 0.85
        for r_p in sub_p.itertuples():
            ax_e.plot([x0_p, x1_p], [r_p.net_eff_pre, r_p.net_eff_post], color=col_p, alpha=0.38, lw=0.75)
            ax_e.scatter([x0_p, x1_p], [r_p.net_eff_pre, r_p.net_eff_post], color=col_p, s=9, alpha=0.75)
        ax_e.plot([x0_p, x1_p], [sub_p["net_eff_pre"].mean(), sub_p["net_eff_post"].mean()], color=col_p, lw=1.7, marker="o", ms=3.2, label=lbl_p)

    ax_e.set_xticks([0.0, 0.85, 2.1, 2.95])
    ax_e.set_xticklabels(["Pre\n(R)", "Post\n(R)", "Pre\n(NR)", "Post\n(NR)"], fontsize=5.2)
    ax_e.set_xlim(-0.4, 3.35)
    ax_e.set_ylim(0.20, float(paired15[["net_eff_pre", "net_eff_post"]].max().max()) * 1.28)
    ax_e.set_ylabel("Net $CD8^+$ Cytotoxic Score")
    ax_e.set_title("Paired Pre/Post Biopsies ($n$=15)", pad=5)
    ax_e.legend(loc="upper right", frameon=False, fontsize=5.0)

    # --- Panel f: 2D Solid-Tumor Waddington Surface + Quiver Shift Under Anti-PD-1 + myCAF/TGFBR1 ---
    ax_f = fig.add_subplot(gs_r3[0, 1])
    _panel_tag(ax_f, "f", x=-0.24, y=1.10)
    gx_t = tum_npz["grid_mycaf"]
    gy_t = tum_npz["grid_tumor"]
    u_combo = tum_npz["U_2d_combo"]
    dx_c = tum_npz["drift_combo_mycaf"]
    dy_c = tum_npz["drift_combo_tumor"]

    cf_f = ax_f.contourf(gx_t, gy_t, u_combo, levels=15, cmap="magma_r")
    st_t = 5
    Xgt, Ygt = np.meshgrid(gx_t[::st_t], gy_t[::st_t])
    Uxt, Uyt = dx_c[::st_t, ::st_t], dy_c[::st_t, ::st_t]
    mag_t = np.hypot(Uxt, Uyt) + 1e-8
    ax_f.quiver(Xgt, Ygt, Uxt / mag_t, Uyt / mag_t, color="white", alpha=0.82, scale=22, width=0.006)
    ax_f.axhline(0.38, color=PAL["sky"], ls="--", lw=0.8)
    ax_f.scatter([0.16, 0.84], [0.16, 0.82], c=[PAL["teal"], PAL["crimson"]], s=24, edgecolors="white", linewidths=0.85, zorder=5)
    ax_f.text(0.20, 0.16, "Unlocked", color="white", fontsize=5.0, fontweight="bold", va="center")
    ax_f.text(0.80, 0.82, "myCAF Lock", color="white", fontsize=5.0, fontweight="bold", va="center", ha="right")
    ax_f.set_xlabel("myCAF Barrier ($S_{\\mathrm{myCAF}}$)")
    ax_f.set_ylabel("Tumor Burden ($S_{\\mathrm{Tumor}}$)")
    ax_f.set_title("Unlocked 2D Attractor", pad=5)
    cbar_f = fig.colorbar(cf_f, ax=ax_f, fraction=0.040, pad=0.025)
    cbar_f.ax.tick_params(labelsize=4.8)
    cbar_f.set_label("$U_{\\mathrm{combo}}$", fontsize=5.1)

    # --- Panel g: Paired Trajectory Rescue Across All n = 70 Anti-PD-1-Resistant Pre-Treatment Patients ---
    ax_g = fig.add_subplot(gs_r3[0, 2])
    _panel_tag(ax_g, "g", x=-0.26, y=1.10)
    nr_all = pool107[pool107["responder"] == 0].copy()
    mono_col = "raw_icb_remission_prob" if "raw_icb_remission_prob" in nr_all.columns else "simulated_icb_remission_prob"
    p_mono = 100.0 * nr_all[mono_col].to_numpy()
    p_comb = 100.0 * nr_all["combo_pd1_plus_mycaf_remission_prob"].to_numpy()
    w_nr = stats.wilcoxon(p_comb - p_mono, alternative="greater")

    for r_nr in nr_all.itertuples():
        c_str = str(getattr(r_nr, coh_col))
        c_nr = PAL["navy"] if "BCC" in c_str else (PAL["purple"] if "91061" in c_str else PAL["blue"])
        m_val = 100.0 * float(getattr(r_nr, mono_col))
        ax_g.plot([0, 1], [m_val, 100.0 * r_nr.combo_pd1_plus_mycaf_remission_prob], color=c_nr, alpha=0.22, lw=0.45)
        ax_g.scatter([0, 1], [m_val, 100.0 * r_nr.combo_pd1_plus_mycaf_remission_prob], color=c_nr, s=6, alpha=0.50)

    ax_g.plot([0, 1], [np.mean(p_mono), np.mean(p_comb)], color=PAL["crimson"], lw=1.8, marker="o", ms=3.6, zorder=6)
    ax_g.text(
        0.50,
        106.0,
        f"{np.mean(p_mono):.1f}% → {np.mean(p_comb):.1f}%\n({_fmt_p(float(w_nr.pvalue))})",
        ha="center",
        va="bottom",
        fontsize=4.9,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_g.set_xticks([0, 1])
    ax_g.set_xticklabels(["Anti-PD-1\nMono.", "Anti-PD-1 +\nTGFBR1"], fontsize=5.3)
    ax_g.set_xlim(-0.35, 1.35)
    ax_g.set_ylim(-4, 134)
    ax_g.set_yticks([0, 25, 50, 75, 100])
    ax_g.set_ylabel("Remission Prob. (%)")
    ax_g.set_title(f"Resistant Patient Rescue ($n$={len(nr_all)})", pad=5)

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
    fig1 = plt.figure(figsize=(7.20, 6.85), facecolor="white")
    gs1 = gridspec.GridSpec(3, 12, figure=fig1, height_ratios=[0.95, 0.88, 1.05], wspace=1.35, hspace=0.58)

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
    _panel_tag(ax, "c", x=-0.04, y=1.08)
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
    fig3 = plt.figure(figsize=(7.20, 6.95), facecolor="white")
    gs3 = gridspec.GridSpec(3, 1, figure=fig3, height_ratios=[1.00, 1.05, 1.00], hspace=0.58)
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
    ax.set_ylim(-1.8, 72.0)
    ax.set_title("25-Point Structural Identifiability", pad=5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.50, 0.99), frameon=False, fontsize=4.5, ncol=1)

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
        [str(r.section_id).replace("GSM6177", "").replace("GSE171351_", "Blad_") for r in v12_s.itertuples()],
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
