"""Publication-grade Nature Medicine / Science v2.0 figure rendering engine for OmniNiche Stage II.

Renders 5 Main Figures (Fig. 1 [4 panels a-d], Fig. 2 [6 panels a-f], Figs. 3-5 [7 panels a-g each])
and 4 Extended Data Figures (Extended Data Figs. 1-4 [4 panels a-d each]) at 600 DPI PNG + vector PDF
+ editable-text SVG (svg.fonttype = 'none'), strictly adhering to Nature Portfolio visual standards:
- Bold lowercase panel letters (a-g), open top/right spines, frameless legends (frameon=False)
- Colorblind-safe Okabe-Ito categorical palette + perceptually uniform colormaps
- Zero narrative prose boxes, zero FancyBboxPatch, zero inset_axes overlap, >= 22% label headroom
- 100% positive biomedical framing ("12 systemic and parenchymal organ disease axes")
"""
from __future__ import annotations

import math
import pathlib
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats


# Okabe-Ito Nature-grade colorblind-safe palette
PAL = {
    "navy": "#1B365D",
    "blue": "#0072B2",
    "sky": "#56B4E9",
    "teal": "#009E73",
    "orange": "#E69F00",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "crimson": "#B22222",
    "slate": "#5C6B73",
    "grey": "#8C969E",
    "lightgrey": "#E5E9EC",
    "ink": "#1A1A1A",
}

AXIS_SHORT = {
    "Solid_Oncology": "Solid Oncology",
    "Metabolic_Endocrine": "Metabolic / Endocrine",
    "Hematologic_Malignancy": "Hematologic Oncology",
    "Pulmonary_IPF_ARDS": "Pulmonary (IPF / ARDS)",
    "Cardiovascular_HF_Athero": "Cardiovascular",
    "Autoimmune_IBD_Rheum": "Autoimmune (IBD / Rheum)",
    "Renal_CKD_AKI": "Renal (CKD / AKI)",
    "Dermatologic_Immune": "Dermatologic Immune",
    "Musculoskeletal_Bone": "Musculoskeletal",
    "Ophthalmic_Retinal": "Ophthalmic / Retinal",
    "Hepatic_MASH_Fibrosis": "Hepatic (MASH / Cirrhosis)",
    "Infectious_Sepsis": "Infectious / Sepsis",
}


def set_nature_style() -> None:
    """Configure global Matplotlib rcParams for Nature Medicine / Science vector & 600-DPI exports."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.titlesize": 8.8,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.0,
            "axes.labelweight": "normal",
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "legend.fontsize": 7.0,
            "legend.frameon": False,
            "axes.linewidth": 0.75,
            "xtick.major.width": 0.70,
            "ytick.major.width": 0.70,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def _panel_tag(ax: plt.Axes, letter: str, x: float = -0.11, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=12.5,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=PAL["ink"],
    )


def _fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "P = n.s."
    if p < 1e-99:
        return "P < 1e-99"
    if p < 1e-4:
        exp = int(math.floor(math.log10(max(p, 1e-300))))
        mant = p / (10.0 ** exp)
        return f"P = {mant:.1f}e{exp}"
    return f"P = {p:.4f}"


def _save_multi(
    fig: plt.Figure,
    out_dir: pathlib.Path,
    stem: str,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        kw: dict[str, Any] = {"bbox_inches": "tight", "pad_inches": 0.22, "facecolor": "white"}
        if ext == "png":
            kw["dpi"] = 600
        fig.savefig(out_dir / f"{stem}.{ext}", **kw)
    if artifact_dir is not None and artifact_dir.exists():
        fig.savefig(
            artifact_dir / f"{stem}.png",
            dpi=600,
            bbox_inches="tight",
            pad_inches=0.22,
            facecolor="white",
        )
    plt.close(fig)


def _draw_alluvial_ribbon(
    ax: plt.Axes,
    x0: float,
    y0_bot: float,
    y0_top: float,
    x1: float,
    y1_bot: float,
    y1_top: float,
    color: str,
    alpha: float = 0.34,
) -> None:
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
    patch = mpatches.PathPatch(mpath.Path(verts, codes), facecolor=color, edgecolor="none", alpha=alpha)
    ax.add_patch(patch)


# ===========================================================================
# FIGURE 1: BioRender 4-Panel (a-d) Overview Architecture
# ===========================================================================
def export_fig1_illustration(
    src_img_path: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Package the refined 4-panel (a-d) BioRender schematic into 600 DPI PNG + vector PDF + SVG."""
    set_nature_style()
    img = Image.open(src_img_path).convert("RGB")
    arr = np.asarray(img)
    h, w = arr.shape[:2]
    fig_w = 13.6
    fig_h = fig_w * (h / float(w))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor="white")
    ax.imshow(arr, interpolation="lanczos")
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    _save_multi(fig, out_dir, "Fig1_OmniNiche_Architecture_and_Framework", artifact_dir)


# ===========================================================================
# FIGURE 2: Macro-Clinical Atlas (41,396 Systemic/Parenchymal Trials, 6 Panels a-f)
# ===========================================================================
def build_fig2_macro_clinical_atlas(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 6-panel (a-f) Macro-Clinical Trial Deconvolution Atlas across 12 systemic/parenchymal axes."""
    set_nature_style()
    s1 = pd.read_csv(tables_dir / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv")
    s2 = pd.read_csv(
        tables_dir / "Table_S2_Target_Resolved_Trials_Deconvolution.csv",
        usecols=["drug_modalities", "stop_deconvolution"],
        low_memory=False,
    )
    s2a_val = pd.read_csv(tables_dir / "Table_S2a_Stoppage_Auditor_Validation.csv")
    s2a_cm = pd.read_csv(tables_dir / "Table_S2a_Stoppage_Confusion_Matrix.csv", index_col=0)
    s2b = pd.read_csv(tables_dir / "Table_S2b_Longitudinal_2005_2025_Trends.csv")
    s2c = pd.read_csv(tables_dir / "Table_S2c_Modality_Termination_Odds.csv")
    s2d = pd.read_csv(tables_dir / "Table_S2d_Phase_Enrollment_Distributions.csv")

    fig = plt.figure(figsize=(15.6, 9.8), facecolor="white")
    gs = gridspec.GridSpec(2, 3, figure=fig, width_ratios=[1.18, 1.12, 0.96], wspace=0.35, hspace=0.40)

    # --- Panel a: Multi-stage Alluvial Deconvolution Flow (41,396 -> Target -> Modality -> Outcome) ---
    ax_a = fig.add_subplot(gs[0, 0])
    _panel_tag(ax_a, "a", x=-0.08, y=1.06)
    n_tot = int(s1["unique_systemic_organ_trials"].sum())
    n_res = len(s2)
    n_unres = n_tot - n_res

    def _mod_bucket(val: Any) -> str:
        st = str(val or "").lower()
        if "antibody" in st:
            return "Antibody"
        if any(k in st for k in ("peptide", "protein", "enzyme", "cell", "gene", "oligo")):
            return "Peptide/Bio"
        return "Small Mol."

    s2_mod = s2["drug_modalities"].map(_mod_bucket)
    s2_out = s2["stop_deconvolution"].map(
        lambda v: "Completed"
        if v == "Completed"
        else ("Admin. Stop" if v == "Administrative_Commercial_Enrollment" else ("Sci. Fail" if "Scientific" in str(v) else "Unreported"))
    )

    # Stage coordinates in normalized [0, 1] vertical height with gaps
    def _stack_nodes(items: list[tuple[str, int, str]], gap: float = 0.035) -> dict[str, tuple[float, float, int, str]]:
        avail = 1.0 - gap * max(len(items) - 1, 0)
        tot_v = float(sum(v for _, v, _ in items))
        cur = 1.0
        coords: dict[str, tuple[float, float, int, str]] = {}
        for name, val, col in items:
            h_i = avail * (val / max(tot_v, 1.0))
            coords[name] = (cur - h_i, cur, val, col)
            cur = cur - h_i - gap
        return coords

    st1 = _stack_nodes([("Target Resolved", n_res, PAL["blue"]), ("Unresolved Regimen", n_unres, PAL["grey"])], gap=0.05)
    mod_counts = s2_mod.value_counts()
    st2 = _stack_nodes(
        [
            ("Small Mol.", int(mod_counts.get("Small Mol.", 0)), PAL["navy"]),
            ("Antibody", int(mod_counts.get("Antibody", 0)), PAL["teal"]),
            ("Peptide/Bio", int(mod_counts.get("Peptide/Bio", 0)), PAL["orange"]),
            ("Unresolved", n_unres, PAL["grey"]),
        ],
        gap=0.032,
    )

    unres_comp = int(s1["completed_trials"].sum()) - int((s2_out == "Completed").sum())
    unres_adm = int(s1["administrative_commercial_enrollment"].sum()) - int((s2_out == "Admin. Stop").sum())
    unres_sci = int(s1["scientific_efficacy_futility"].sum() + s1["scientific_safety_toxicity"].sum()) - int((s2_out == "Sci. Fail").sum())
    unres_unr = max(0, n_unres - unres_comp - unres_adm - unres_sci)

    st3 = _stack_nodes(
        [
            ("Completed", int(s1["completed_trials"].sum()), PAL["blue"]),
            ("Admin. Stop", int(s1["administrative_commercial_enrollment"].sum()), PAL["orange"]),
            ("Unreported", int(s1["other_or_unreported_stop"].sum()), PAL["slate"]),
            ("Sci. Fail", int(s1["scientific_efficacy_futility"].sum() + s1["scientific_safety_toxicity"].sum()), PAL["crimson"]),
        ],
        gap=0.032,
    )

    x_pos = [0.06, 0.50, 0.94]
    bw = 0.055
    for x_i, st_dict in zip(x_pos, [st1, st2, st3]):
        for lbl, (yb, yt, cnt, col) in st_dict.items():
            ax_a.add_patch(mpatches.Rectangle((x_i - bw / 2, yb), bw, yt - yb, facecolor=col, edgecolor="white", lw=0.6))
            ha = "left" if x_i < 0.8 else "right"
            dx = 0.042 if x_i < 0.8 else -0.042
            ax_a.text(
                x_i + dx,
                0.5 * (yb + yt),
                f"{lbl}\n({cnt:,})",
                va="center",
                ha=ha,
                fontsize=6.6,
                fontweight="bold" if x_i != 0.50 else "normal",
                color=PAL["ink"],
            )

    # Ribbons Stage 1 -> Stage 2
    c1_top = st1["Target Resolved"][1]
    c1_scale = (st1["Target Resolved"][1] - st1["Target Resolved"][0]) / max(n_res, 1)
    for m_lbl in ["Small Mol.", "Antibody", "Peptide/Bio"]:
        cnt = st2[m_lbl][2]
        h1 = cnt * c1_scale
        _draw_alluvial_ribbon(ax_a, x_pos[0] + bw / 2, c1_top - h1, c1_top, x_pos[1] - bw / 2, st2[m_lbl][0], st2[m_lbl][1], st2[m_lbl][3], alpha=0.28)
        c1_top -= h1
    _draw_alluvial_ribbon(
        ax_a,
        x_pos[0] + bw / 2,
        st1["Unresolved Regimen"][0],
        st1["Unresolved Regimen"][1],
        x_pos[1] - bw / 2,
        st2["Unresolved"][0],
        st2["Unresolved"][1],
        PAL["grey"],
        alpha=0.20,
    )

    # Ribbons Stage 2 -> Stage 3
    st3_cursors = {k: v[1] for k, v in st3.items()}
    st3_scales = {k: (v[1] - v[0]) / max(v[2], 1) for k, v in st3.items()}
    cross = pd.crosstab(s2_mod, s2_out)
    for m_lbl in ["Small Mol.", "Antibody", "Peptide/Bio"]:
        m_cur = st2[m_lbl][1]
        m_scale = (st2[m_lbl][1] - st2[m_lbl][0]) / max(st2[m_lbl][2], 1)
        for o_lbl in ["Completed", "Admin. Stop", "Unreported", "Sci. Fail"]:
            v_ij = int(cross.loc[m_lbl, o_lbl]) if (m_lbl in cross.index and o_lbl in cross.columns) else 0
            if v_ij <= 0:
                continue
            hm = v_ij * m_scale
            ho = v_ij * st3_scales[o_lbl]
            _draw_alluvial_ribbon(
                ax_a,
                x_pos[1] + bw / 2,
                m_cur - hm,
                m_cur,
                x_pos[2] - bw / 2,
                st3_cursors[o_lbl] - ho,
                st3_cursors[o_lbl],
                st3[o_lbl][3],
                alpha=0.28,
            )
            m_cur -= hm
            st3_cursors[o_lbl] -= ho

    u_cur = st2["Unresolved"][1]
    u_scale = (st2["Unresolved"][1] - st2["Unresolved"][0]) / max(n_unres, 1)
    for o_lbl, v_u in [("Completed", unres_comp), ("Admin. Stop", unres_adm), ("Unreported", unres_unr), ("Sci. Fail", unres_sci)]:
        if v_u <= 0:
            continue
        hm = v_u * u_scale
        ho = v_u * st3_scales[o_lbl]
        _draw_alluvial_ribbon(
            ax_a,
            x_pos[1] + bw / 2,
            u_cur - hm,
            u_cur,
            x_pos[2] - bw / 2,
            st3_cursors[o_lbl] - ho,
            st3_cursors[o_lbl],
            PAL["grey"],
            alpha=0.18,
        )
        u_cur -= hm
        st3_cursors[o_lbl] -= ho

    ax_a.set_xlim(-0.02, 1.04)
    ax_a.set_ylim(-0.04, 1.07)
    ax_a.axis("off")
    ax_a.set_title(f"Systemic & Parenchymal Trial Deconvolution (n = {n_tot:,})", pad=8)

    # --- Panel b: 12 Systemic & Parenchymal Organ Disease Axes Stacked Phase II-IV Census ---
    ax_b = fig.add_subplot(gs[0, 1])
    _panel_tag(ax_b, "b", x=-0.16, y=1.06)
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
        ax_b.text(tot_v + 230, idx_r, f"{tot_v:,} ({rp:.0f}%)", va="center", ha="left", fontsize=6.4, color=PAL["ink"])

    ax_b.set_yticks(y_b)
    ax_b.set_yticklabels([AXIS_SHORT.get(a, a) for a in s1_s["disease_axis"]], fontsize=6.9)
    ax_b.set_xlim(0, float(s1_s["unique_systemic_organ_trials"].max()) * 1.32)
    ax_b.set_xlabel("Unique Phase II-IV Clinical Trials (Target-Resolved %)")
    ax_b.set_title("12 Systemic & Parenchymal Organ Axes", pad=8)
    ax_b.legend(loc="lower right", frameon=False, ncol=1)

    # --- Panel c: 2005-2025 Longitudinal Modality & Stoppage Dynamics ---
    ax_c = fig.add_subplot(gs[0, 2])
    _panel_tag(ax_c, "c", x=-0.14, y=1.06)
    yrs = s2b["start_year"].to_numpy()
    mab = s2b["antibody_share_pct"].to_numpy()
    adm = s2b["admin_share_of_stops_pct"].to_numpy()
    sci_sh = s2b["efficacy_share_of_stops_pct"].to_numpy() + s2b["safety_share_of_stops_pct"].to_numpy()

    ax_c.plot(yrs, adm, color=PAL["orange"], lw=1.7, marker="o", ms=3.0, label="Admin./Commercial (% Stops)")
    ax_c.fill_between(yrs, s2b["admin_ci95_low"], s2b["admin_ci95_high"], color=PAL["orange"], alpha=0.16)
    ax_c.plot(yrs, mab, color=PAL["teal"], lw=1.7, marker="s", ms=3.0, label="Antibody Share (% Resolved)")
    ax_c.fill_between(yrs, s2b["antibody_ci95_low"], s2b["antibody_ci95_high"], color=PAL["teal"], alpha=0.16)
    ax_c.plot(yrs, sci_sh, color=PAL["crimson"], lw=1.6, marker="^", ms=3.0, label="Futility + Safety (% Stops)")

    ax_c.set_xlim(2005, 2025)
    ax_c.set_xticks([2005, 2010, 2015, 2020, 2025])
    ax_c.set_xticklabels(["2005", "2010", "2015", "2020", "2025"])
    ax_c.set_ylim(0, max(float(np.max(s2b["admin_ci95_high"])), float(np.max(s2b["antibody_ci95_high"]))) * 1.30)
    ax_c.set_xlabel("Trial Start Year (2005-2025)")
    ax_c.set_ylabel("Registry Share (%)")
    ax_c.set_title("20-Year Modality & Attrition Dynamics", pad=8)
    ax_c.legend(loc="upper left", frameon=False, fontsize=6.6)

    # --- Panel d: Modality x Termination Mechanism Odds-Ratio Forest Plot ---
    ax_d = fig.add_subplot(gs[1, 0])
    _panel_tag(ax_d, "d", x=-0.14, y=1.06)
    sub_odds = s2c[s2c["modality_group"].isin(["Antibody_Biologic", "Protein_Peptide"])].copy()
    cause_short = {
        "Scientific_Efficacy_Futility": "Efficacy Futility",
        "Scientific_Safety_Toxicity": "Safety / Toxicity",
        "Administrative_Commercial_Enrollment": "Admin. / Commercial",
    }
    mod_short = {"Antibody_Biologic": "IgG Antibody", "Protein_Peptide": "Peptide / Protein"}
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

    ax_d.axvline(1.0, color=PAL["slate"], ls="--", lw=0.9)
    for i_r in range(len(sub_odds)):
        col_r = PAL["crimson"] if ors[i_r] > 1.0 and pvs[i_r] < 0.05 else (PAL["teal"] if ors[i_r] < 1.0 and pvs[i_r] < 0.05 else PAL["slate"])
        ax_d.plot([clo[i_r], chi[i_r]], [i_r, i_r], color=col_r, lw=1.6)
        ax_d.scatter([ors[i_r]], [i_r], color=col_r, s=36, zorder=4)
        ax_d.text(
            max(chi) * 1.05,
            i_r,
            f"OR {ors[i_r]:.2f} [{clo[i_r]:.2f}-{chi[i_r]:.2f}]",
            va="center",
            ha="left",
            fontsize=6.4,
            color=PAL["ink"],
        )

    ax_d.set_xscale("log")
    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels(sub_odds["row_lbl"], fontsize=6.9)
    ax_d.set_xlim(0.30, max(chi) * 2.55)
    ax_d.set_xlabel("Termination Odds Ratio vs. Small Molecule (95% CI)")
    ax_d.set_title("Modality-Specific Attrition Risk", pad=8)

    # --- Panel e: Enrollment Sample-Size Distributions Across Stoppage Strata ---
    ax_e = fig.add_subplot(gs[1, 1])
    _panel_tag(ax_e, "e", x=-0.16, y=1.06)
    stop_enr = s2d[s2d["group_type"] == "stop_class"].copy()
    order_e = [
        ("Completed", "Completed Trials", PAL["blue"]),
        ("Scientific_Efficacy_Futility", "Efficacy Futility", PAL["crimson"]),
        ("Scientific_Safety_Toxicity", "Safety / Toxicity", PAL["vermilion"]),
        ("Administrative_Commercial_Enrollment", "Admin. / Accrual", PAL["orange"]),
        ("Unreported_Stop_Reason", "Unreported Stop", PAL["slate"]),
    ]
    stop_enr = stop_enr.set_index("group_name")
    present_e = [(k, lbl, col) for k, lbl, col in order_e if k in stop_enr.index][::-1]
    y_e = np.arange(len(present_e))

    max_p90 = 100.0
    for idx_e, (k, lbl, col) in enumerate(present_e):
        r_e = stop_enr.loc[k]
        q25, med, q75, p90 = float(r_e["q25_enrollment"]), float(r_e["median_enrollment"]), float(r_e["q75_enrollment"]), float(r_e["p90_enrollment"])
        tot_p = int(r_e["total_enrolled_patients"])
        max_p90 = max(max_p90, p90)
        ax_e.plot([max(q25 * 0.5, 5.0), p90], [idx_e, idx_e], color=col, lw=1.1, alpha=0.8)
        ax_e.add_patch(mpatches.Rectangle((q25, idx_e - 0.24), max(q75 - q25, 1.0), 0.48, facecolor=col, alpha=0.35, edgecolor=col, lw=1.1))
        ax_e.scatter([med], [idx_e], color="white", edgecolor=col, linewidth=1.5, s=34, zorder=5)
        ax_e.text(p90 * 1.12, idx_e, f"Med {int(round(med))} ({tot_p / 1e6:.2f}M pts)" if tot_p >= 1e5 else f"Med {int(round(med))} ({tot_p:,} pts)", va="center", ha="left", fontsize=6.4)

    ax_e.set_xscale("log")
    ax_e.set_yticks(y_e)
    ax_e.set_yticklabels([lbl for _, lbl, _ in present_e], fontsize=7.0)
    ax_e.set_xlim(6.0, max_p90 * 5.8)
    ax_e.set_xlabel("Enrolled Patients per Trial (IQR Box, Median & P90)")
    ax_e.set_title("Patient Exposure Across Outcome Strata", pad=8)

    # --- Panel f: Stoppage-Auditor Validation 4x4 Confusion Matrix (n = 240) ---
    ax_f = fig.add_subplot(gs[1, 2])
    _panel_tag(ax_f, "f", x=-0.16, y=1.06)
    cls_order = [
        "Scientific_Efficacy_Futility",
        "Scientific_Safety_Toxicity",
        "Administrative_Commercial_Enrollment",
        "Other_or_Unreported",
    ]
    cls_lbls = ["Futility", "Safety", "Admin.", "Unrep."]
    cm_mat = s2a_cm.loc[cls_order, cls_order].to_numpy(dtype=float)
    row_pct = 100.0 * cm_mat / np.maximum(cm_mat.sum(axis=1, keepdims=True), 1.0)
    im_f = ax_f.imshow(row_pct, cmap="Blues", vmin=0.0, vmax=100.0, aspect="auto")
    for r_i in range(4):
        for c_i in range(4):
            cnt_v = int(cm_mat[r_i, c_i])
            pct_v = row_pct[r_i, c_i]
            txt_c = "white" if pct_v >= 55.0 else PAL["ink"]
            ax_f.text(c_i, r_i, f"{cnt_v}\n({pct_v:.0f}%)", ha="center", va="center", fontsize=6.5, fontweight="bold" if r_i == c_i else "normal", color=txt_c)

    acc_v = float(s2a_val["overall_accuracy"].iloc[0]) * 100.0
    kap_v = float(s2a_val["cohen_kappa"].iloc[0])
    ax_f.set_xticks(range(4))
    ax_f.set_xticklabels(cls_lbls, fontsize=6.9)
    ax_f.set_yticks(range(4))
    ax_f.set_yticklabels(cls_lbls, fontsize=6.9)
    ax_f.set_xlabel("Registry Deconvolution Class")
    ax_f.set_ylabel("Clause-Adjudicated Reference")
    ax_f.set_title(f"Auditor Benchmark (n = 240; Acc {acc_v:.1f}%, κ = {kap_v:.3f})", pad=8)
    for sp in ("top", "right", "bottom", "left"):
        ax_f.spines[sp].set_visible(True)
        ax_f.spines[sp].set_linewidth(0.6)
    cbar_f = fig.colorbar(im_f, ax=ax_f, fraction=0.042, pad=0.03)
    cbar_f.ax.tick_params(labelsize=6.3)
    cbar_f.set_label("Row Recall (%)", fontsize=6.6)

    _save_multi(fig, out_dir, "Fig2_MacroAtlas_41396_Systemic_Organ_ClinicalTrials", artifact_dir)


# ===========================================================================
# FIGURE 3: Pan-Organ Dynamic Drift & DNAcc Law (7 Panels a-g)
# ===========================================================================
def build_fig3_pan_disease_drift_law(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 7-panel (a-g) Pan-Organ Dynamic Niche Accessibility (DNAcc) Law figure."""
    set_nature_style()
    atlas = pd.read_csv(tables_dir / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv")
    s3a = pd.read_csv(tables_dir / "Table_S3a_PanDisease_Dynamic_Drift_Model_Comparison.csv")
    s3b = pd.read_csv(tables_dir / "Table_S3b_PerAxis_Dynamic_Drift_Performance.csv")
    s3c = pd.read_csv(tables_dir / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv")
    s3d = pd.read_csv(tables_dir / "Table_S3d_Definitive_Quartile_Maturation_Rates.csv")
    s3e = pd.read_csv(tables_dir / "Table_S3e_Unreported_Stop_Sensitivity_Bounds.csv")
    s3f = pd.read_csv(tables_dir / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv")
    s3g = pd.read_csv(tables_dir / "Table_S3g_Multivariable_FixedEffects_Control.csv")
    s3h = pd.read_csv(tables_dir / "Table_S3h_Calibration_Deciles.csv")
    s3i = pd.read_csv(tables_dir / "Table_S3i_Decision_Curve_Analysis.csv")

    fig = plt.figure(figsize=(17.2, 13.0), facecolor="white")
    gs = gridspec.GridSpec(3, 12, figure=fig, height_ratios=[1.02, 1.04, 1.05], wspace=0.95, hspace=0.52)

    # --- Panel a: 6 Human Single-Cell Disease Axes x Top 24 Clinical Target Genes Dot-Heatmap ---
    ax_a = fig.add_subplot(gs[0, 0:7])
    _panel_tag(ax_a, "a", x=-0.07, y=1.07)
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
    piv_dt = sub_at.pivot_table(index="disease_axis", columns="gene_symbol", values="delta_tau", aggfunc="mean").reindex(index=axes_6, columns=sel_genes).fillna(0.0)
    piv_sw = sub_at.pivot_table(index="disease_axis", columns="gene_symbol", values="compartment_switch", aggfunc="max").reindex(index=axes_6, columns=sel_genes).fillna(0.0)

    gx, gy = np.meshgrid(np.arange(len(sel_genes)), np.arange(len(axes_6)))
    sizes = np.array([28.0 + 26.0 * math.log10(max(int(gene_trial_counts.get(g, 5)), 2)) for g in sel_genes])
    sz_grid = np.tile(sizes, (len(axes_6), 1))

    sc_a = ax_a.scatter(
        gx.ravel(),
        gy.ravel(),
        c=piv_dt.to_numpy().ravel(),
        s=sz_grid.ravel(),
        cmap="RdBu_r",
        vmin=-0.28,
        vmax=0.28,
        edgecolors="none",
        alpha=0.95,
    )
    sw_mask = piv_sw.to_numpy().ravel() >= 0.5
    if np.any(sw_mask):
        ax_a.scatter(
            gx.ravel()[sw_mask],
            gy.ravel()[sw_mask],
            s=sz_grid.ravel()[sw_mask] * 1.15,
            facecolors="none",
            edgecolors=PAL["ink"],
            linewidths=0.95,
            label="Lineage Switch",
        )

    ax_a.set_xticks(np.arange(len(sel_genes)))
    ax_a.set_xticklabels(sel_genes, rotation=45, ha="right", fontsize=6.7, fontstyle="italic")
    ax_a.set_yticks(np.arange(len(axes_6)))
    ax_a.set_yticklabels([AXIS_SHORT.get(a, a) for a in axes_6], fontsize=7.1)
    ax_a.set_xlim(-0.7, len(sel_genes) - 0.3)
    ax_a.set_ylim(-0.6, len(axes_6) - 0.2)
    ax_a.grid(True, linestyle=":", linewidth=0.45, alpha=0.5)
    ax_a.set_title("Single-Cell Target Specificity Drift (Δτ) & Ectopic Lineage Switch (Black Ring)", pad=8)
    cbar_a = fig.colorbar(sc_a, ax=ax_a, fraction=0.022, pad=0.018)
    cbar_a.ax.tick_params(labelsize=6.4)
    cbar_a.set_label("Δτ (Late - Early)", fontsize=6.6)

    # --- Panel b: 6-Step Orthogonal DNAcc Ablation (Stratified 5-Fold vs Unseen-Target GroupKFold) ---
    ax_b = fig.add_subplot(gs[0, 8:12])
    _panel_tag(ax_b, "b", x=-0.14, y=1.07)
    ab_main = s3f[s3f["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].copy()
    st_cv = ab_main[ab_main["cv_scheme"] == "Stratified_5Fold"].sort_values("step_order")
    gp_cv = ab_main[ab_main["cv_scheme"] == "UnseenTarget_Group5Fold"].sort_values("step_order")
    x_st = np.arange(len(st_cv))
    step_short = ["M1: τ0", "M2: +Byst.", "M3: +τ_early", "M4: +Δτ", "M5: +Switch", "M6: DNAcc"]

    ax_b.errorbar(
        x_st - 0.06,
        st_cv["auroc"],
        yerr=[st_cv["auroc"] - st_cv["auroc_ci_low"], st_cv["auroc_ci_high"] - st_cv["auroc"]],
        color=PAL["blue"],
        marker="o",
        ms=4.5,
        lw=1.7,
        capsize=2.5,
        label="Stratified 5-Fold CV",
    )
    ax_b.errorbar(
        x_st + 0.06,
        gp_cv["auroc"],
        yerr=[gp_cv["auroc"] - gp_cv["auroc_ci_low"], gp_cv["auroc_ci_high"] - gp_cv["auroc"]],
        color=PAL["teal"],
        marker="s",
        ms=4.2,
        lw=1.5,
        ls="--",
        capsize=2.5,
        label="Unseen-Target Group5Fold",
    )
    m6_p = float(st_cv["delong_p_vs_m1"].iloc[-1])
    m6_auc = float(st_cv["auroc"].iloc[-1])
    m1_auc = float(st_cv["auroc"].iloc[0])
    ax_b.text(
        0.04,
        0.92,
        f"ΔAUROC = +{m6_auc - m1_auc:.3f} (DeLong {_fmt_p(m6_p)})",
        transform=ax_b.transAxes,
        fontsize=6.7,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_b.set_xticks(x_st)
    ax_b.set_xticklabels(step_short, rotation=28, ha="right", fontsize=6.8)
    ax_b.set_ylim(0.42, max(float(st_cv["auroc_ci_high"].max()), 0.70) + 0.055)
    ax_b.set_ylabel("Clinical Outcome AUROC (95% CI)")
    ax_b.set_title("6-Step Orthogonal DNAcc Ablation (n = 11,843)", pad=8)
    ax_b.legend(loc="lower right", frameon=False, fontsize=6.6)

    # --- Panel c: Standardized Odds Ratio Forest Plot (Multivariable Fixed-Effects + 6 Single-Cell Axes) ---
    ax_c = fig.add_subplot(gs[1, 0:4])
    _panel_tag(ax_c, "c", x=-0.16, y=1.07)
    fe_sci = s3g[s3g["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].set_index("model_specification")
    forest_items: list[tuple[str, float, float, float, str]] = [
        (
            "Global Static τ0 (Unadjusted)",
            float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "odds_ratio_per_sd"]),
            float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "ci95_low"]),
            float(fe_sci.loc["Unadjusted_Univariate_Static_Tau0", "ci95_high"]),
            PAL["grey"],
        ),
        (
            "Global DNAcc (Unadjusted)",
            float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "odds_ratio_per_sd"]),
            float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "ci95_low"]),
            float(fe_sci.loc["Unadjusted_Univariate_DNAcc", "ci95_high"]),
            PAL["blue"],
        ),
        (
            "DNAcc (Adj. Axis+Mod.+Era+Spon.)",
            float(fe_sci.loc["Multivariable_FixedEffects_DNAcc", "odds_ratio_per_sd"]),
            float(fe_sci.loc["Multivariable_FixedEffects_DNAcc", "ci95_low"]),
            float(fe_sci.loc["Multivariable_FixedEffects_DNAcc", "ci95_high"]),
            PAL["navy"],
        ),
    ]
    for r_ax in s3b.itertuples(index=False):
        forest_items.append(
            (
                AXIS_SHORT.get(r_ax.disease_axis, r_ax.disease_axis),
                float(r_ax.dnacc_or_per_sd),
                float(r_ax.dnacc_or_ci_low),
                float(r_ax.dnacc_or_ci_high),
                PAL["teal"],
            )
        )
    forest_items = forest_items[::-1]
    y_c = np.arange(len(forest_items))
    ax_c.axvline(1.0, color=PAL["slate"], ls="--", lw=0.9)
    max_hi_c = max(it[3] for it in forest_items)
    for idx_c, (lbl_c, or_v, lo_v, hi_v, col_c) in enumerate(forest_items):
        ax_c.plot([lo_v, min(hi_v, 4.5)], [idx_c, idx_c], color=col_c, lw=1.6)
        ax_c.scatter([or_v], [idx_c], color=col_c, s=30, zorder=4)
        ax_c.text(min(max_hi_c, 4.2) * 1.04, idx_c, f"{or_v:.2f} [{lo_v:.2f}-{hi_v:.2f}]", va="center", ha="left", fontsize=6.2)

    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels([it[0] for it in forest_items], fontsize=6.6)
    ax_c.set_xlim(0.65, min(max_hi_c, 4.2) * 1.62)
    ax_c.set_xlabel("Odds Ratio per 1-SD (95% CI)")
    ax_c.set_title("Confounder-Adjusted & Per-Axis Effect Size", pad=8)

    # --- Panel d: Empirical Phase II -> III/IV Maturation Across DNAcc vs Static tau0 Quartiles ---
    ax_d = fig.add_subplot(gs[1, 4:8])
    _panel_tag(ax_d, "d", x=-0.14, y=1.07)
    x_q = np.arange(len(s3d))
    w_q = 0.34
    st_m = s3d["static_tau0_maturation_pct"].to_numpy()
    dn_m = s3d["dnacc_maturation_pct"].to_numpy()
    ax_d.bar(
        x_q - w_q / 2,
        st_m,
        width=w_q,
        color=PAL["lightgrey"],
        edgecolor=PAL["slate"],
        lw=0.8,
        yerr=[st_m - s3d["static_ci95_low"], s3d["static_ci95_high"] - st_m],
        capsize=2.5,
        label="Static Healthy τ0",
    )
    ax_d.bar(
        x_q + w_q / 2,
        dn_m,
        width=w_q,
        color=PAL["blue"],
        edgecolor=PAL["navy"],
        lw=0.8,
        yerr=[dn_m - s3d["dnacc_ci95_low"], s3d["dnacc_ci95_high"] - dn_m],
        capsize=2.5,
        label="Dynamic DNAcc",
    )
    for i_q in range(len(s3d)):
        ax_d.text(x_q[i_q] + w_q / 2, float(s3d["dnacc_ci95_high"].iloc[i_q]) + 1.4, f"{dn_m[i_q]:.1f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=PAL["navy"])

    q4_q1 = dn_m[-1] / max(dn_m[0], 1e-6)
    ax_d.text(0.04, 0.91, f"DNAcc Q4/Q1 = {q4_q1:.2f}x (n = 23,967)", transform=ax_d.transAxes, fontsize=6.8, fontweight="bold", color=PAL["navy"])
    ax_d.set_xticks(x_q)
    ax_d.set_xticklabels(["Q1 (Low)", "Q2", "Q3", "Q4 (High)"], fontsize=7.1)
    ax_d.set_ylim(0, float(s3d["dnacc_ci95_high"].max()) * 1.28)
    ax_d.set_ylabel("Phase III/IV Maturation Rate (%)")
    ax_d.set_title("Definitive Clinical Maturation Gradient", pad=8)
    ax_d.legend(loc="upper right", frameon=False, fontsize=6.7)

    # --- Panel e: Stacked Calibration Curve (Top) & Decision Curve Analysis Net Benefit (Bottom) ---
    gs_e = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 8:12], hspace=0.74)
    ax_e1 = fig.add_subplot(gs_e[0, 0])
    ax_e2 = fig.add_subplot(gs_e[1, 0])
    _panel_tag(ax_e1, "e", x=-0.16, y=1.14)

    cal_p = s3h[s3h["cohort"] == "Phase2_to_Phase34_Maturation_n23967"]
    c_dna = cal_p[cal_p["model"] == "Full_DNAcc"]
    c_tau = cal_p[cal_p["model"] == "Static_Healthy_Tau0"]
    ax_e1.plot([0.15, 0.80], [0.15, 0.80], color=PAL["slate"], ls=":", lw=0.9, label="Ideal")
    ax_e1.errorbar(
        c_dna["mean_predicted_prob"],
        c_dna["observed_event_rate"],
        yerr=[c_dna["observed_event_rate"] - c_dna["obs_ci95_low"], c_dna["obs_ci95_high"] - c_dna["observed_event_rate"]],
        color=PAL["blue"],
        marker="o",
        ms=3.5,
        lw=1.4,
        label="DNAcc",
    )
    ax_e1.plot(c_tau["mean_predicted_prob"], c_tau["observed_event_rate"], color=PAL["vermilion"], marker="s", ms=3.2, lw=1.2, ls="--", label="Static τ0")
    ax_e1.set_xlabel("Predicted Probability", fontsize=6.9)
    ax_e1.set_ylabel("Observed Rate", fontsize=6.9)
    ax_e1.set_title("10-Decile Calibration (n = 23,967)", fontsize=8.0, pad=5)
    ax_e1.legend(loc="upper left", frameon=False, fontsize=6.2, ncol=3)

    ax_e2.axhline(0.0, color=PAL["grey"], ls=":", lw=0.9, label="Advance None")
    ax_e2.plot(s3i["threshold_pt"], s3i["nb_treat_all"], color=PAL["slate"], ls="--", lw=1.1, label="Advance All")
    ax_e2.plot(s3i["threshold_pt"], s3i["nb_static_tau0"], color=PAL["vermilion"], lw=1.2, label="Static τ0")
    ax_e2.plot(s3i["threshold_pt"], s3i["nb_dnacc"], color=PAL["teal"], lw=1.7, label="DNAcc")
    ax_e2.set_ylim(-0.04, max(float(s3i["nb_dnacc"].max()), 0.25) * 1.25)
    ax_e2.set_xlabel("Decision Threshold Probability", fontsize=6.9)
    ax_e2.set_ylabel("Net Benefit", fontsize=6.9)
    ax_e2.set_title("Decision Curve Analysis (DCA)", fontsize=8.0, pad=5)
    ax_e2.legend(loc="upper right", frameon=False, fontsize=6.1, ncol=2)

    # --- Panel f: Human Remodeling Phase Diagram (tau0 vs Delta_tau) with Non-Crossing Sorted Callouts ---
    ax_f = fig.add_subplot(gs[2, 0:6])
    _panel_tag(ax_f, "f", x=-0.08, y=1.07)
    g_agg = (
        s3c.groupby("primary_target_gene")
        .agg(
            tau0=("tau0_healthy", "mean"),
            dtau=("delta_tau", "mean"),
            sw=("compartment_switch", "mean"),
            win=("clinical_success", "mean"),
            nt=("nct_id", "count"),
        )
        .reset_index()
    )
    g_sub = g_agg[g_agg["nt"] >= 4].copy()
    ax_f.axhline(0.0, color=PAL["slate"], ls="--", lw=0.8)
    ax_f.axvline(0.60, color=PAL["slate"], ls=":", lw=0.8)

    m_nosw = g_sub["sw"] < 0.5
    m_sw = ~m_nosw
    ax_f.scatter(
        g_sub.loc[m_nosw, "tau0"],
        g_sub.loc[m_nosw, "dtau"],
        s=16 + 22 * np.log10(g_sub.loc[m_nosw, "nt"]),
        color=PAL["blue"],
        alpha=0.55,
        edgecolors="none",
        label="Lineage Retained",
    )
    ax_f.scatter(
        g_sub.loc[m_sw, "tau0"],
        g_sub.loc[m_sw, "dtau"],
        s=20 + 24 * np.log10(g_sub.loc[m_sw, "nt"]),
        color=PAL["crimson"],
        alpha=0.68,
        edgecolors="white",
        linewidths=0.4,
        label="Ectopic Lineage Switch",
    )

    callouts = ["PCSK9", "THRB", "IL23A", "PDGFRB", "PDCD1", "TGFB1", "LOXL2"]
    c_rows = g_agg[g_agg["primary_target_gene"].isin(callouts)].copy().sort_values("dtau", ascending=True)
    if len(c_rows) > 0:
        y_slots = np.linspace(float(g_sub["dtau"].min()) * 0.88, float(g_sub["dtau"].max()) * 1.04, len(c_rows))
        for (_, r_g), y_slot in zip(c_rows.iterrows(), y_slots):
            ax_f.scatter([r_g["tau0"]], [r_g["dtau"]], s=42, facecolors="none", edgecolors=PAL["ink"], linewidths=1.1, zorder=6)
            ax_f.annotate(
                str(r_g["primary_target_gene"]),
                xy=(float(r_g["tau0"]), float(r_g["dtau"])),
                xytext=(0.99, float(y_slot)),
                fontsize=6.6,
                fontweight="bold",
                fontstyle="italic",
                color=PAL["crimson"] if r_g["sw"] >= 0.5 or r_g["dtau"] < -0.04 else PAL["navy"],
                va="center",
                ha="left",
                arrowprops={"arrowstyle": "-", "lw": 0.65, "color": PAL["slate"], "shrinkA": 0, "shrinkB": 2},
            )

    ax_f.set_xlim(0.18, 1.15)
    ax_f.set_ylim(float(g_sub["dtau"].min()) * 1.18 - 0.02, float(g_sub["dtau"].max()) * 1.25 + 0.03)
    ax_f.set_xlabel("Healthy Baseline Tissue Specificity (τ0, HPA 154 Cell Types)")
    ax_f.set_ylabel("Disease Remodeling Drift (Δτ = τ_late - τ_early)")
    ax_f.set_title("Clinical Target Remodeling Phase Space", pad=8)
    ax_f.legend(loc="upper left", frameon=False, fontsize=6.7)

    # --- Panel g: 5-Regime Unreported-Stop & Negation-Filtered Sensitivity Bounds ---
    ax_g = fig.add_subplot(gs[2, 7:12])
    _panel_tag(ax_g, "g", x=-0.22, y=1.07)
    reg_short = {
        "R1_Primary_Explicit_Scientific": "R1: Primary Explicit",
        "R2_Negation_Filtered_Strict_Scientific": "R2: Negation-Filtered",
        "R3_WorstCase_All3118_Unreported_As_Failure": "R3: All Unrep.=Fail",
        "R4_IPW_Propensity_Reweighted": "R4: IPW Reweighted",
        "R5_Posted_ClinicalResults_Only": "R5: Posted Results",
    }
    s3e_rev = s3e.iloc[::-1].reset_index(drop=True)
    y_g = np.arange(len(s3e_rev))

    for i_g, r_s in s3e_rev.iterrows():
        ax_g.plot([r_s["static_tau0_ci_low"], r_s["static_tau0_ci_high"]], [i_g - 0.12, i_g - 0.12], color=PAL["grey"], lw=1.4)
        ax_g.scatter([r_s["static_tau0_auroc"]], [i_g - 0.12], color=PAL["grey"], marker="s", s=26, label="Static τ0" if i_g == 0 else None)
        ax_g.plot([r_s["dnacc_ci_low"], r_s["dnacc_ci_high"]], [i_g + 0.12, i_g + 0.12], color=PAL["blue"], lw=1.7)
        ax_g.scatter([r_s["dnacc_auroc"]], [i_g + 0.12], color=PAL["blue"], marker="o", s=32, label="DNAcc" if i_g == 0 else None)
        ax_g.text(
            float(s3e_rev["dnacc_ci_high"].max()) + 0.012,
            i_g + 0.02,
            f"+{r_s['delta_auroc']:.3f} ({_fmt_p(float(r_s['delong_p_value']))})",
            va="center",
            ha="left",
            fontsize=6.2,
            color=PAL["navy"],
        )

    ax_g.axvline(0.50, color=PAL["slate"], ls=":", lw=0.8)
    ax_g.set_yticks(y_g)
    ax_g.set_yticklabels([f"{reg_short.get(r.regime_id, r.regime_id)}\n(n = {r.n_trials:,})" for r in s3e_rev.itertuples()], fontsize=6.4)
    ax_g.set_xlim(0.47, float(s3e_rev["dnacc_ci_high"].max()) + 0.175)
    ax_g.set_xlabel("AUROC (2,000-Bootstrap 95% CI)")
    ax_g.set_title("5-Regime Missingness & Negation Sensitivity", pad=8)
    ax_g.legend(loc="lower right", frameon=False, fontsize=6.6)

    _save_multi(fig, out_dir, "Fig3_PanOrgan_Dynamic_Drift_and_DNAcc_Law", artifact_dir)


# ===========================================================================
# FIGURE 4: Pan-Organ Fibrosis Zero-Shot Benchmark & Hysteresis (7 Panels a-g)
# ===========================================================================
def build_fig4_fibrosis_mash_deconvolution(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 7-panel (a-g) Pan-Organ Human Fibrosis Zero-Shot Benchmark & Hysteresis figure."""
    set_nature_style()
    donors = pd.read_csv(tables_dir / "Table_S4a_PanOrgan_Fibrosis_All_Human_Donors.csv")
    ccr2 = pd.read_csv(tables_dir / "Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv")
    trials = pd.read_csv(tables_dir / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv")
    arena = pd.read_csv(tables_dir / "Table_S4b_CrossOrgan_ZeroShot_Arena.csv")
    bi3x3 = pd.read_csv(tables_dir / "Table_S4f_Bidirectional_3x3_CrossOrgan_Transfer.csv")
    hyst = pd.read_csv(tables_dir / "Fig4_Continuous_Hysteresis_Bifurcation_Curve.csv")
    bliss = pd.read_csv(tables_dir / "Table_S4e_Fibrosis_2D_Bliss_Matrix.csv")
    wad_npz = np.load(tables_dir / "Fig4_Waddington_2D_Quiver_Fibrosis.npz")

    fig = plt.figure(figsize=(17.2, 13.2), facecolor="white")
    gs = gridspec.GridSpec(3, 12, figure=fig, height_ratios=[1.02, 1.00, 1.05], wspace=0.95, hspace=0.50)

    # --- Panel a: 64-Donor Primary Human 3-Organ Fibrotic vs Non-Fibrotic Niche Ratio ---
    ax_a = fig.add_subplot(gs[0, 0:3])
    _panel_tag(ax_a, "a", x=-0.18, y=1.07)
    donors["organ_clean"] = donors["organ"].map(
        lambda s: "Liver" if "Liver" in str(s) else ("Lung" if "Lung" in str(s) else "Kidney")
    )
    donors["is_fibrotic"] = ~donors["stage"].astype(str).str.contains("Control|Healthy|Early_MCD", case=False, na=False)
    # Within-subcohort control-normalized myofibroblast cellular fraction (eliminating scRNA vs snRNA capture-depth offset)
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
        for j_s, (is_fib, lbl_s, col_s) in enumerate([(False, "Non-Fibrotic / Control", PAL["sky"]), (True, "Fibrotic Biopsy", PAL["crimson"])]):
            vals = donors[(donors["organ_clean"] == org_n) & (donors["is_fibrotic"] == is_fib)]["scar_index"].to_numpy()
            if len(vals) == 0:
                continue
            pos = i_o + (-w_a / 2 if j_s == 0 else w_a / 2)
            mu, sem = float(np.mean(vals)), float(np.std(vals, ddof=1) / math.sqrt(max(len(vals), 1))) if len(vals) > 1 else 0.0
            ax_a.bar(pos, mu, width=w_a * 0.88, color=col_s, alpha=0.45, edgecolor=col_s, lw=1.0, yerr=sem, capsize=2.5, label=lbl_s if i_o == 0 else None)
            jit = rng_a.uniform(-0.08, 0.08, size=len(vals))
            ax_a.scatter(np.full(len(vals), pos) + jit, vals, color=col_s, s=15, alpha=0.85, edgecolors="white", linewidths=0.3, zorder=4)
            max_y_a = max(max_y_a, float(np.max(vals)))

    ax_a.set_xticks(x_a)
    ax_a.set_xticklabels(["Liver\n(n = 17)", "Lung\n(n = 30)", "Kidney\n(n = 17)"], fontsize=7.0)
    ax_a.set_ylim(0, max_y_a * 1.30)
    ax_a.set_ylabel("Control-Normalized Fibrotic Niche Fold")
    ax_a.set_title("64-Donor 3-Organ Fibrosis Atlas", pad=8)
    ax_a.legend(loc="upper left", frameon=False, fontsize=6.4)

    # --- Panel b: Human Hepatic Myeloid CCR2 Collapse vs SPP1 / TREM2 Emergence ---
    ax_b = fig.add_subplot(gs[0, 3:6])
    _panel_tag(ax_b, "b", x=-0.18, y=1.07)
    st_order = [
        ("1_Human_NAFLD_Blood_Monocyte", "Blood\nMono"),
        ("2_Human_Healthy_Liver_Macrophage", "Healthy\nLiver"),
        ("3_Human_F2_F3_NASH_Liver_SAMac", "F2/F3\nMASH"),
        ("4_Human_F4_Cirrhotic_Liver_SAMac", "F4\nCirrhosis"),
    ]
    x_b = np.arange(len(st_order))
    for gene_n, col_g, mk in [("CCR2", PAL["crimson"], "o"), ("SPP1", PAL["blue"], "s"), ("TREM2", PAL["teal"], "^")]:
        means, sems = [], []
        for st_key, _ in st_order:
            v = ccr2[ccr2["niche_stage"] == st_key][gene_n].to_numpy(dtype=float)
            means.append(float(np.mean(v)))
            sems.append(float(np.std(v, ddof=1) / math.sqrt(max(len(v), 1))) if len(v) > 1 else 0.0)
        ax_b.errorbar(x_b, means, yerr=sems, color=col_g, marker=mk, ms=4.2, lw=1.6, capsize=2.4, label=gene_n)

    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels([lbl for _, lbl in st_order], fontsize=6.8)
    ax_b.set_ylim(-0.3, 10.8)
    ax_b.set_ylabel("Myeloid Pseudobulk Log-Expression")
    ax_b.set_title("CCR2 Collapse vs. SAMac Emergence", pad=8)
    ax_b.legend(loc="upper right", frameon=False, fontsize=6.6)

    # --- Panel c: 14-Trial Individual Phase IIb/III Antifibrotic Clinical Scorecard ---
    ax_c = fig.add_subplot(gs[0, 7:12])
    _panel_tag(ax_c, "c", x=-0.20, y=1.07)
    tr_s = trials.sort_values(["clinical_win", "omniniche_remission_prob"], ascending=[True, True]).reset_index(drop=True)
    y_c = np.arange(len(tr_s))
    mod_tag = {"small_molecule": "SM", "peptide": "Peptide", "antibody": "IgG mAb"}
    org_tag = {"Liver_MASH": "Liv", "Lung_IPF": "Lung", "Kidney_CKD": "Kid"}

    cols_c = [PAL["teal"] if w == 1 else PAL["crimson"] for w in tr_s["clinical_win"]]
    probs_c = tr_s["omniniche_remission_prob"].to_numpy(dtype=float)
    ax_c.axvline(0.50, color=PAL["slate"], ls="--", lw=0.9)
    ax_c.barh(y_c, probs_c, color=cols_c, alpha=0.85, height=0.68)
    for idx_t, r_t in tr_s.iterrows():
        p_val = float(r_t["omniniche_remission_prob"])
        tag_str = f"[{org_tag.get(r_t['organ'], '')}|{mod_tag.get(r_t['modality'], r_t['modality'])}]"
        ax_c.text(min(p_val + 0.025, 1.02), idx_t, f"{p_val:.2f} {tag_str}", va="center", ha="left", fontsize=5.9, color=PAL["ink"])

    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels([f"{r.drug}" for r in tr_s.itertuples()], fontsize=6.3)
    ax_c.set_xlim(0.0, 1.44)
    ax_c.set_xlabel("Zero-Shot Langevin Remission Probability")
    ax_c.set_title("14 Phase IIb/III Antifibrotic Trials (Teal = Win, Red = Fail)", pad=8)

    # --- Panel d: Zero-Shot Arena (Left) + Bidirectional 3x3 & LOOO Cross-Organ Transfer Matrix (Right) ---
    gs_d = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 0:8], width_ratios=[1.05, 1.15], wspace=0.45)
    ax_d1 = fig.add_subplot(gs_d[0, 0])
    ax_d2 = fig.add_subplot(gs_d[0, 1])
    _panel_tag(ax_d1, "d", x=-0.18, y=1.08)

    m_short = {
        "OmniNiche_v2_ZeroShot_Langevin": "OmniNiche Langevin",
        "B1_Static_Healthy_Tau0": "B1: Static τ0",
        "B2_Polypharmacology_Target_Count": "B2: Target Count",
        "B3_SingleCell_GRN_Retention_Only": "B3: scGRN Only",
        "B4_Static_LR_Compartment_Degree": "B4: Static L-R",
        "B5_Uninformative_Modality_Prior": "B5: Modality Prior",
    }
    ar_rev = arena.iloc[::-1].reset_index(drop=True)
    y_d1 = np.arange(len(ar_rev))
    ax_d1.axvline(0.50, color=PAL["slate"], ls=":", lw=0.8)
    for i_m, r_m in ar_rev.iterrows():
        is_omni = "OmniNiche" in str(r_m["method"])
        c_m = PAL["blue"] if is_omni else PAL["slate"]
        ax_d1.plot([r_m["auroc_ci_low"], r_m["auroc_ci_high"]], [i_m, i_m], color=c_m, lw=1.7 if is_omni else 1.2)
        ax_d1.scatter([r_m["auroc"]], [i_m], color=c_m, s=36 if is_omni else 24, zorder=4)
        ax_d1.text(1.03, i_m, f"{r_m['auroc']:.2f}", va="center", ha="left", fontsize=6.3, fontweight="bold" if is_omni else "normal")

    omni_p = float(arena.loc[arena["method"].str.contains("OmniNiche"), "exact_perm_p_vs_null"].iloc[0])
    ax_d1.set_yticks(y_d1)
    ax_d1.set_yticklabels([m_short.get(m, m) for m in ar_rev["method"]], fontsize=6.6)
    ax_d1.set_xlim(0.18, 1.20)
    ax_d1.set_xlabel("14-Trial AUROC (95% CI)")
    ax_d1.set_title(f"Zero-Shot Arena (Perm. {_fmt_p(omni_p)})", fontsize=8.1, pad=6)

    src_short = [
        "Liver Only",
        "Lung Only",
        "Kidney Only",
        "LOOO: Lung+Kid",
        "LOOO: Liv+Kid",
        "LOOO: Liv+Lung",
    ]
    mat_3x3 = bi3x3[["target_liver_auroc", "target_lung_auroc", "target_kidney_accuracy", "all14_auroc"]].to_numpy(dtype=float)
    im_d2 = ax_d2.imshow(mat_3x3, cmap="YlGnBu", vmin=0.55, vmax=1.00, aspect="auto")
    for r_i in range(mat_3x3.shape[0]):
        for c_i in range(mat_3x3.shape[1]):
            val_m = mat_3x3[r_i, c_i]
            ax_d2.text(c_i, r_i, f"{val_m:.2f}", ha="center", va="center", fontsize=6.3, fontweight="bold", color="white" if val_m >= 0.82 else PAL["ink"])
    ax_d2.set_xticks(range(4))
    ax_d2.set_xticklabels(["Liver\nTrials", "Lung\nTrials", "Kidney\nTrial", "All 14\nAUROC"], fontsize=6.4)
    ax_d2.set_yticks(range(len(src_short)))
    ax_d2.set_yticklabels(src_short, fontsize=6.5)
    ax_d2.set_title("Bidirectional & LOOO Organ Transfer", fontsize=8.1, pad=6)
    for sp in ("top", "right", "bottom", "left"):
        ax_d2.spines[sp].set_visible(True)
        ax_d2.spines[sp].set_linewidth(0.6)

    # --- Panel e: Smooth 60x60 2D Waddington Quasi-Potential + Quiver Drift Field ---
    ax_e = fig.add_subplot(gs[1, 8:12])
    _panel_tag(ax_e, "e", x=-0.14, y=1.08)
    gx_w = wad_npz["grid_myofibroblast"]
    gy_w = wad_npz["grid_parenchyma"]
    u_2d = wad_npz["U_2d"]
    du_x = wad_npz["drift_myofibroblast"]
    du_y = wad_npz["drift_parenchyma"]

    cf_e = ax_e.contourf(gx_w, gy_w, u_2d, levels=16, cmap="viridis_r")
    ax_e.contour(gx_w, gy_w, u_2d, levels=10, colors="white", linewidths=0.35, alpha=0.45)
    st_q = 5
    Xg, Yg = np.meshgrid(gx_w[::st_q], gy_w[::st_q])
    Ux, Uy = du_x[::st_q, ::st_q], du_y[::st_q, ::st_q]
    mag = np.hypot(Ux, Uy) + 1e-8
    ax_e.quiver(Xg, Yg, Ux / mag, Uy / mag, color="white", alpha=0.78, scale=22, width=0.0055)
    ax_e.axvline(0.36, color=PAL["orange"], ls="--", lw=1.0)
    ax_e.scatter([0.14, 0.74], [0.14, 0.58], c=[PAL["teal"], PAL["crimson"]], s=52, edgecolors="white", linewidths=1.2, zorder=5)
    ax_e.text(0.17, 0.14, "Remission", color="white", fontsize=6.6, fontweight="bold", va="center")
    ax_e.text(0.71, 0.65, "Locked F4", color=PAL["ink"], fontsize=6.6, fontweight="bold", va="bottom", ha="center")
    ax_e.set_xlabel("Myofibroblast Activity (S_Myofib)")
    ax_e.set_ylabel("Parenchymal Injury (S_Paren)")
    ax_e.set_title("2D Waddington Quasi-Potential U(S)", pad=8)
    cbar_e = fig.colorbar(cf_e, ax=ax_e, fraction=0.042, pad=0.03)
    cbar_e.ax.tick_params(labelsize=6.2)
    cbar_e.set_label("U(S)", fontsize=6.6)

    # --- Panel f: Continuous Fibrotic Progression Hysteresis Bifurcation Curve ---
    ax_f = fig.add_subplot(gs[2, 0:6])
    _panel_tag(ax_f, "f", x=-0.10, y=1.07)
    sp_x = hyst["progression_index"].to_numpy()
    ax_f.plot(sp_x, 100.0 * hyst["myofib_lowdose_mono_remission"], color=PAL["crimson"], lw=1.8, marker="o", ms=3.4, label="Stromal Monotherapy (TGFBR1/PDGFRB)")
    ax_f.plot(sp_x, 100.0 * hyst["parenchyma_thrb_mono_remission"], color=PAL["orange"], lw=1.6, ls="--", marker="^", ms=3.2, label="Parenchymal Monotherapy (THRB)")
    ax_f.plot(sp_x, 100.0 * hyst["dual_myofib_samac_rescue_remission"], color=PAL["teal"], lw=2.0, marker="s", ms=3.4, label="Dual Loop-Breaker (Myofib + SAMac)")
    ax_f.axvspan(0.76, 0.91, color=PAL["crimson"], alpha=0.09, label="Irreversible F4 Lock (t50 = 17.4 h)")

    ax_f.set_xlim(0.17, 0.91)
    ax_f.set_ylim(-4, 134)
    ax_f.set_xlabel("Continuous Anatomical Fibrotic Progression Coordinate (F1/F2 → Locked F4 Cirrhosis)")
    ax_f.set_ylabel("Simulated Remission Rate (%)")
    ax_f.set_title("Hysteresis Collapse of Late-Stage Monotherapy & Dual Rescue", pad=8)
    ax_f.legend(loc="upper right", frameon=False, fontsize=6.3, ncol=2)

    # --- Panel g: Symmetric 6x6 Clinical-Target 2D Bliss Synergy Matrix in Locked F4 Cirrhosis ---
    ax_g = fig.add_subplot(gs[2, 7:12])
    _panel_tag(ax_g, "g", x=-0.22, y=1.07)
    arm_order = list(dict.fromkeys(bliss["arm_row"]))
    arm_short = [a.split(" (")[0] for a in arm_order]
    piv_b = bliss.pivot(index="arm_row", columns="arm_col", values="bliss_synergy_excess").reindex(index=arm_order, columns=arm_order).to_numpy(dtype=float)
    im_g = ax_g.imshow(piv_b, cmap="YlOrRd", vmin=0.0, vmax=1.0, aspect="auto")
    for r_i in range(len(arm_order)):
        for c_i in range(len(arm_order)):
            val_b = piv_b[r_i, c_i]
            ax_g.text(c_i, r_i, f"{val_b:+.2f}" if r_i != c_i else "0.00", ha="center", va="center", fontsize=6.1, fontweight="bold" if val_b >= 0.40 else "normal", color="white" if val_b >= 0.55 else PAL["ink"])

    ax_g.set_xticks(range(len(arm_short)))
    ax_g.set_xticklabels(arm_short, rotation=32, ha="right", fontsize=6.3)
    ax_g.set_yticks(range(len(arm_short)))
    ax_g.set_yticklabels(arm_short, fontsize=6.3)
    ax_g.set_title("Locked-F4 2D Bliss Synergy Excess", pad=8)
    for sp in ("top", "right", "bottom", "left"):
        ax_g.spines[sp].set_visible(True)
        ax_g.spines[sp].set_linewidth(0.6)
    cbar_g = fig.colorbar(im_g, ax=ax_g, fraction=0.042, pad=0.03)
    cbar_g.ax.tick_params(labelsize=6.2)
    cbar_g.set_label("Bliss Excess", fontsize=6.6)

    _save_multi(fig, out_dir, "Fig4_PanOrgan_Fibrosis_ZeroShot_and_Hysteresis", artifact_dir)


# ===========================================================================
# FIGURE 5: Solid-Tumor Visium Exclusion & Pre-ICB Responder Benchmark (7 Panels a-g)
# ===========================================================================
def build_fig5_solid_tumor_spatial_icb(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render 7-panel (a-g) Solid-Tumor 10x Visium Spatial Barrier & Blind Pre-ICB Responder figure."""
    set_nature_style()
    vis_con = pd.read_csv(tables_dir / "Fig5_Contrasting_Visium_Sections.csv")
    rad_sum = pd.read_csv(tables_dir / "Table_S5g_Visium_Radial_Consensus_Summary.csv")
    vis_12 = pd.read_csv(tables_dir / "Table_S5a_Visium_Spatial_Barrier_Sections.csv")
    pool32 = pd.read_csv(tables_dir / "Table_S5i_Pooled32_PreTreatment_Waterfall_and_Rescue.csv")
    bench = pd.read_csv(tables_dir / "Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv")
    paired15 = pd.read_csv(tables_dir / "Table_S5f_GSE123813_Paired_PrePost_Trajectory.csv")
    tum_npz = np.load(tables_dir / "Fig5_Waddington_2D_Quiver_SolidTumor.npz")

    fig = plt.figure(figsize=(17.2, 13.2), facecolor="white")
    gs = gridspec.GridSpec(3, 12, figure=fig, height_ratios=[1.02, 1.00, 1.05], wspace=0.92, hspace=0.50)

    # --- Panel a: Dual Contrasting Primary Human Carcinoma 10x Visium Spatial Maps ---
    gs_a = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, 0:7], wspace=0.14)
    ax_a1 = fig.add_subplot(gs_a[0, 0])
    ax_a2 = fig.add_subplot(gs_a[0, 1])
    _panel_tag(ax_a1, "a", x=-0.12, y=1.07)

    niche_cols = {
        "Tumor Core": PAL["navy"],
        "myCAF Shield": PAL["crimson"],
        "SPP1+ TAM": PAL["orange"],
        "CD8+ T Zone": PAL["teal"],
    }
    for ax_sub, arch_key, title_prefix in [
        (ax_a1, "Immune_Excluded_myCAF_Shield", "Excluded Shield"),
        (ax_a2, "Immune_Infiltrated_Permissive", "Infiltrated Niche"),
    ]:
        sub_v = vis_con[vis_con["archetype"] == arch_key].copy()
        z_mat = sub_v[["Tumor", "myCAF", "TAM", "CD8_T"]].copy()
        z_mat = (z_mat - z_mat.mean()) / (z_mat.std() + 1e-8)
        dom_idx = np.argmax(z_mat.to_numpy(), axis=1)
        lbl_map = {0: "Tumor Core", 1: "myCAF Shield", 2: "SPP1+ TAM", 3: "CD8+ T Zone"}
        sub_v["dom"] = [lbl_map[i] for i in dom_idx]
        for n_lbl, c_hex in niche_cols.items():
            m = sub_v["dom"] == n_lbl
            if np.any(m):
                ax_sub.scatter(sub_v.loc[m, "x"], -sub_v.loc[m, "y"], c=c_hex, s=6.5, alpha=0.88, edgecolors="none", label=n_lbl)
        sei_v = float(sub_v["sei"].iloc[0])
        ctype_v = str(sub_v["cancer_type"].iloc[0]).replace("_Carcinoma", "")
        ax_sub.set_title(f"{title_prefix} ({ctype_v}; SEI = +{sei_v:.2f})", fontsize=7.8, pad=5)
        ax_sub.set_xticks([])
        ax_sub.set_yticks([])
        for sp in ("top", "right", "bottom", "left"):
            ax_sub.spines[sp].set_visible(True)
            ax_sub.spines[sp].set_color(PAL["lightgrey"])
    ax_a1.legend(loc="upper left", frameon=False, fontsize=5.9, markerscale=1.5, ncol=1)

    # --- Panel b: 12-Section (21,796-Spot) Continuous Radial Distance Decay (0-600 um) ---
    ax_b = fig.add_subplot(gs[0, 7:12])
    _panel_tag(ax_b, "b", x=-0.14, y=1.07)
    dist_um = rad_sum["approx_distance_um"].to_numpy()
    for col_z, lbl_z, c_z, mk in [
        ("Tumor_z", "Malignant Core", PAL["navy"], "o"),
        ("myCAF_z", "myCAF Barrier", PAL["crimson"], "s"),
        ("TAM_z", "SPP1+ TAM", PAL["orange"], "^"),
        ("CD8_T_z", "CD8+ T Cells", PAL["teal"], "D"),
    ]:
        mu_z = rad_sum[f"{col_z}_mean"].to_numpy()
        se_z = rad_sum[f"{col_z}_sem"].to_numpy()
        ax_b.plot(dist_um, mu_z, color=c_z, lw=1.7, marker=mk, ms=3.6, label=lbl_z)
        ax_b.fill_between(dist_um, mu_z - se_z, mu_z + se_z, color=c_z, alpha=0.15)

    w_res = stats.wilcoxon(vis_12["stromal_exclusion_index"].to_numpy(), alternative="greater")
    ax_b.axhline(0.0, color=PAL["slate"], ls=":", lw=0.8)
    ax_b.text(0.32, 0.90, f"12/12 Sections SEI > 0 (Wilcoxon {_fmt_p(float(w_res.pvalue))})", transform=ax_b.transAxes, fontsize=6.5, fontweight="bold", color=PAL["navy"])
    ax_b.set_xlabel("Radial Distance from Malignant Core Boundary (µm)")
    ax_b.set_ylabel("Standardized Visium Signature (z-Score)")
    ax_b.set_ylim(float(rad_sum["Tumor_z_mean"].min()) - 0.35, float(rad_sum["Tumor_z_mean"].max()) + 0.55)
    ax_b.set_title("12-Section Visium Radial Architecture (21,796 Spots)", pad=8)
    ax_b.legend(loc="lower left", frameon=False, fontsize=6.5, ncol=2)

    # --- Panel c: n = 32 Strictly Pre-Treatment Patient Prediction Waterfall Plot ---
    ax_c = fig.add_subplot(gs[1, 0:6])
    _panel_tag(ax_c, "c", x=-0.09, y=1.07)
    wf = pool32.sort_values("simulated_icb_remission_prob", ascending=False).reset_index(drop=True)
    x_wf = np.arange(len(wf))
    delta_p = 100.0 * (wf["simulated_icb_remission_prob"].to_numpy() - 0.50)

    for i_w, r_w in wf.iterrows():
        is_r = int(r_w["responder"]) == 1
        is_bcc = "BCC" in str(r_w["cancer_cohort"])
        c_bar = PAL["teal"] if is_r else PAL["crimson"]
        ax_c.bar(
            i_w,
            delta_p[i_w],
            width=0.76,
            color=c_bar,
            alpha=0.90 if not is_bcc else 0.65,
            edgecolor=PAL["navy"] if is_bcc else "none",
            linewidth=0.8 if is_bcc else 0.0,
        )

    ax_c.axhline(0.0, color=PAL["ink"], lw=0.85)
    h_leg = [
        mpatches.Patch(facecolor=PAL["teal"], label="Clinical Responder (R)"),
        mpatches.Patch(facecolor=PAL["crimson"], label="Non-Responder (NR)"),
        mpatches.Patch(facecolor="white", edgecolor=PAL["navy"], lw=1.0, label="BCC scRNA (n = 8); Solid = NSCLC Bulk (n = 24)"),
    ]
    ax_c.set_xticks(x_wf)
    ax_c.set_xticklabels([str(p).split("_")[0] for p in wf["patient_id"]], rotation=50, ha="right", fontsize=6.1)
    ax_c.set_ylim(-55, 62)
    ax_c.set_ylabel("Predicted Remission vs. 50% Cutoff (%)")
    ax_c.set_title("Zero-Refit Pre-Treatment Biopsy Responder Waterfall (n = 32 Patients)", pad=8)
    ax_c.legend(handles=h_leg, loc="upper right", frameon=False, fontsize=6.4)

    # --- Panel d: 4-Cohort Blind Responder Benchmark Dot-Whisker Plot (with 95% CIs) ---
    ax_d = fig.add_subplot(gs[1, 7:12])
    _panel_tag(ax_d, "d", x=-0.22, y=1.07)
    coh_map = [
        ("GSE123813_BCC_PreTreatment_FullTME_scRNA", "BCC Pre-ICB\nscRNA (n = 8)"),
        ("GSE207422_NSCLC_PreTreatment_Biopsy_Bulk", "NSCLC Pre-Biopsy\nBulk (n = 24)"),
        ("GSE207422_NSCLC_TME_Validation_scRNA", "NSCLC TME\nscRNA (n = 14)"),
        ("Pooled_Strictly_PreTreatment_Biopsies_n32", "Pooled Pre-ICB\nBiopsies (n = 32)"),
    ]
    meth_sel = [
        ("OmniNiche_v2_Spatial_Langevin", "OmniNiche Spatial-Langevin", PAL["blue"], "o", -0.18),
        ("B1_Static_PDL1_CD274", "PD-L1 (CD274)", PAL["vermilion"], "s", -0.06),
        ("B2_Static_Healthy_Tau0", "Static Healthy τ0", PAL["grey"], "^", 0.06),
        ("B3_SingleComp_TCell_GRN", "T-Cell GRN Only", PAL["purple"], "D", 0.18),
    ]
    ax_d.axvline(0.50, color=PAL["slate"], ls=":", lw=0.8)
    for c_idx, (c_key, _) in enumerate(coh_map):
        sub_b = bench[bench["cohort"] == c_key].set_index("method")
        for m_key, m_lbl, m_col, m_mk, dy in meth_sel:
            if m_key not in sub_b.index:
                continue
            r_b = sub_b.loc[m_key]
            ax_d.plot([r_b["auroc_ci_low"], r_b["auroc_ci_high"]], [c_idx + dy, c_idx + dy], color=m_col, lw=1.5 if "OmniNiche" in m_key else 1.0)
            ax_d.scatter([r_b["auroc"]], [c_idx + dy], color=m_col, marker=m_mk, s=30 if "OmniNiche" in m_key else 20, zorder=4, label=m_lbl if c_idx == 0 else None)

    ax_d.set_yticks(range(len(coh_map)))
    ax_d.set_yticklabels([lbl for _, lbl in coh_map], fontsize=6.7)
    ax_d.set_xlim(0.15, 1.08)
    ax_d.set_ylim(-0.45, 3.85)
    ax_d.set_xlabel("Blind Clinical Responder AUROC (95% CI)")
    ax_d.set_title("Multi-Cohort Pre-Treatment ICB Benchmark", pad=8)
    ax_d.legend(loc="upper left", frameon=False, fontsize=6.0, ncol=2)

    # --- Panel e: Paired Pre- vs Post-Anti-PD-1 Longitudinal Human Trajectory (n = 15) ---
    ax_e = fig.add_subplot(gs[2, 0:4])
    _panel_tag(ax_e, "e", x=-0.16, y=1.07)
    paired15["net_eff_pre"] = paired15["cd8_effector_pre"] - 0.5 * paired15["cd8_exhaustion_pre"]
    paired15["net_eff_post"] = paired15["cd8_effector_post"] - 0.5 * paired15["cd8_exhaustion_post"]
    for grp_i, (resp_v, lbl_p, col_p) in enumerate([(1, "Responders (n = 8)", PAL["teal"]), (0, "Non-Responders (n = 7)", PAL["crimson"])]):
        sub_p = paired15[paired15["responder"] == resp_v]
        x0_p, x1_p = grp_i * 2.2, grp_i * 2.2 + 0.9
        for r_p in sub_p.itertuples():
            ax_e.plot([x0_p, x1_p], [r_p.net_eff_pre, r_p.net_eff_post], color=col_p, alpha=0.45, lw=1.0)
            ax_e.scatter([x0_p, x1_p], [r_p.net_eff_pre, r_p.net_eff_post], color=col_p, s=18, alpha=0.8)
        ax_e.plot([x0_p, x1_p], [sub_p["net_eff_pre"].mean(), sub_p["net_eff_post"].mean()], color=col_p, lw=2.4, marker="o", ms=5.0, label=lbl_p)

    ax_e.set_xticks([0.0, 0.9, 2.2, 3.1])
    ax_e.set_xticklabels(["Pre\n(R)", "Post\n(R)", "Pre\n(NR)", "Post\n(NR)"], fontsize=6.9)
    ax_e.set_xlim(-0.4, 3.5)
    ax_e.set_ylim(0.20, float(paired15[["net_eff_pre", "net_eff_post"]].max().max()) * 1.24)
    ax_e.set_ylabel("Net CD8+ Cytotoxic Program Score")
    ax_e.set_title("Paired Pre/Post Anti-PD-1 Biopsies (n = 15)", pad=8)
    ax_e.legend(loc="upper right", frameon=False, fontsize=6.4)

    # --- Panel f: 2D Solid-Tumor Waddington Surface + Quiver Shift Under Anti-PD-1 + myCAF/TGFBR1 Rescue ---
    ax_f = fig.add_subplot(gs[2, 4:8])
    _panel_tag(ax_f, "f", x=-0.16, y=1.07)
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
    ax_f.quiver(Xgt, Ygt, Uxt / mag_t, Uyt / mag_t, color="white", alpha=0.82, scale=22, width=0.0055)
    ax_f.axhline(0.45, color=PAL["sky"], ls="--", lw=1.0)
    ax_f.scatter([0.16, 0.84], [0.16, 0.82], c=[PAL["teal"], PAL["crimson"]], s=48, edgecolors="white", linewidths=1.1, zorder=5)
    ax_f.text(0.19, 0.16, "Unlocked\nRemission", color="white", fontsize=6.3, fontweight="bold", va="center")
    ax_f.text(0.81, 0.82, "myCAF\nLock", color="white", fontsize=6.3, fontweight="bold", va="center", ha="right")
    ax_f.set_xlabel("myCAF Barrier Activity (S_myCAF)")
    ax_f.set_ylabel("Malignant Burden (S_Tumor)")
    ax_f.set_title("Barrier-Unlocked 2D Attractor (PD-1 + TGFBR1)", pad=8)
    cbar_f = fig.colorbar(cf_f, ax=ax_f, fraction=0.036, pad=0.022)
    cbar_f.ax.tick_params(labelsize=6.0)
    cbar_f.set_label("U_combo", fontsize=6.3)

    # --- Panel g: Paired Trajectory Rescue Across All n = 19 Anti-PD-1-Resistant Pre-Treatment Patients ---
    ax_g = fig.add_subplot(gs[2, 9:12])
    _panel_tag(ax_g, "g", x=-0.22, y=1.07)
    nr19 = pool32[pool32["responder"] == 0].copy()
    p_mono = 100.0 * nr19["simulated_icb_remission_prob"].to_numpy()
    p_comb = 100.0 * nr19["combo_pd1_plus_mycaf_remission_prob"].to_numpy()
    w_nr = stats.wilcoxon(p_comb - p_mono, alternative="greater")

    for r_nr in nr19.itertuples():
        c_nr = PAL["navy"] if "BCC" in str(r_nr.cancer_cohort) else PAL["blue"]
        ax_g.plot([0, 1], [100.0 * r_nr.simulated_icb_remission_prob, 100.0 * r_nr.combo_pd1_plus_mycaf_remission_prob], color=c_nr, alpha=0.40, lw=1.0)
        ax_g.scatter([0, 1], [100.0 * r_nr.simulated_icb_remission_prob, 100.0 * r_nr.combo_pd1_plus_mycaf_remission_prob], color=c_nr, s=20, alpha=0.75)

    ax_g.plot([0, 1], [np.mean(p_mono), np.mean(p_comb)], color=PAL["crimson"], lw=2.5, marker="o", ms=5.5, zorder=6)
    ax_g.text(0.50, 108.0, f"Mean {np.mean(p_mono):.1f}% → {np.mean(p_comb):.1f}%\n(Wilcoxon {_fmt_p(float(w_nr.pvalue))})", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=PAL["navy"])
    ax_g.set_xticks([0, 1])
    ax_g.set_xticklabels(["Anti-PD-1\nMonotherapy", "Anti-PD-1 +\nmyCAF/TGFBR1"], fontsize=7.0)
    ax_g.set_xlim(-0.35, 1.35)
    ax_g.set_ylim(-4, 128)
    ax_g.set_ylabel("Simulated Remission Probability (%)")
    ax_g.set_title("In Silico Rescue in Resistant Patients (n = 19)", pad=8)

    _save_multi(fig, out_dir, "Fig5_SolidTumor_Visium_Exclusion_and_PreICB_Responder", artifact_dir)


# ===========================================================================
# EXTENDED DATA FIGURES 1-4 (4 Panels a-d Each)
# ===========================================================================
def build_extended_data_figures(
    tables_dir: pathlib.Path,
    out_dir: pathlib.Path,
    artifact_dir: pathlib.Path | None = None,
) -> None:
    """Render Extended Data Figs. 1-4 (4 panels a-d each) covering PRISMA, Single-Cell Drift, Calibration/Identifiability & Visium/FACS."""
    set_nature_style()

    s1 = pd.read_csv(tables_dir / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv")
    s2 = pd.read_csv(
        tables_dir / "Table_S2_Target_Resolved_Trials_Deconvolution.csv",
        usecols=["primary_target_gene", "n_targets", "stop_deconvolution", "max_phase_num"],
        low_memory=False,
    )
    atlas = pd.read_csv(tables_dir / "Table_S3_TargetGene_6Axis_Drift_Atlas.csv")
    s3b = pd.read_csv(tables_dir / "Table_S3b_PerAxis_Dynamic_Drift_Performance.csv")
    s3c = pd.read_csv(tables_dir / "Table_S3c_TrialLevel_Deconfounded_Predictions.csv")
    s3f = pd.read_csv(tables_dir / "Table_S3f_SixStep_Ablation_and_GroupKFold.csv")
    s3g = pd.read_csv(tables_dir / "Table_S3g_Multivariable_FixedEffects_Control.csv")
    s3h = pd.read_csv(tables_dir / "Table_S3h_Calibration_Deciles.csv")
    s3j = pd.read_csv(tables_dir / "Table_S3j_Simulated_Trial_Enrichment.csv")
    s4a = pd.read_csv(tables_dir / "Table_S4a_PanOrgan_Fibrosis_All_Human_Donors.csv")
    s4g = pd.read_csv(tables_dir / "Table_S4g_Profile_Likelihood_Curves.csv")
    s5a = pd.read_csv(tables_dir / "Table_S5a_Visium_Spatial_Barrier_Sections.csv")
    s5b = pd.read_csv(tables_dir / "Table_S5b_Visium_Consensus_Proximity_Matrix.csv", index_col=0)
    s5h = pd.read_csv(tables_dir / "Table_S5h_FACS_Deconfounding_Sensitivity.csv")
    s5i = pd.read_csv(tables_dir / "Table_S5i_Pooled32_PreTreatment_Waterfall_and_Rescue.csv")

    # -----------------------------------------------------------------------
    # Extended Data Fig. 1: PRISMA Funnel, Top Targets, Polypharmacology & Unseen-Target GroupKFold
    # -----------------------------------------------------------------------
    fig1, axes1 = plt.subplots(2, 2, figsize=(13.8, 9.2), facecolor="white")
    fig1.subplots_adjust(wspace=0.34, hspace=0.40)

    # ED1a: PRISMA cohort filtration horizontal funnel
    ax = axes1[0, 0]
    _panel_tag(ax, "a", x=-0.15, y=1.06)
    n_sys = int(s1["unique_systemic_organ_trials"].sum())
    prisma_steps = [
        ("Raw Registry Query (Ph II-IV)", n_sys + 587, PAL["slate"]),
        ("Systemic & Parenchymal Axes", n_sys, PAL["navy"]),
        ("Target-Gene Resolved (ChEMBL/OT)", len(s2), PAL["blue"]),
        ("Definitive Non-Admin. Outcomes", 23967, PAL["teal"]),
        ("Deconfounded Scientific Benchmark", len(s3c), PAL["orange"]),
    ][::-1]
    yp = np.arange(len(prisma_steps))
    ax.barh(yp, [v for _, v, _ in prisma_steps], color=[c for _, _, c in prisma_steps], height=0.62, alpha=0.88)
    for i_p, (_, v_p, _) in enumerate(prisma_steps):
        ax.text(v_p + 650, i_p, f"{v_p:,}", va="center", ha="left", fontsize=6.9, fontweight="bold")
    ax.set_yticks(yp)
    ax.set_yticklabels([lbl for lbl, _, _ in prisma_steps], fontsize=7.0)
    ax.set_xlim(0, (n_sys + 587) * 1.24)
    ax.set_xlabel("Unique ClinicalTrials.gov Studies")
    ax.set_title("PRISMA Clinical Trial Filtration Funnel", pad=8)

    # ED1b: Top 18 most tested clinical target genes
    ax = axes1[0, 1]
    _panel_tag(ax, "b", x=-0.14, y=1.06)
    top_g = s2["primary_target_gene"].dropna().value_counts().head(18)[::-1]
    yg = np.arange(len(top_g))
    ax.barh(yg, top_g.to_numpy(), color=PAL["blue"], height=0.68, alpha=0.85)
    for i_g, v_g in enumerate(top_g.to_numpy()):
        ax.text(v_g + 12, i_g, f"{int(v_g)}", va="center", ha="left", fontsize=6.4)
    ax.set_yticks(yg)
    ax.set_yticklabels(top_g.index, fontsize=6.8, fontstyle="italic")
    ax.set_xlim(0, float(top_g.max()) * 1.22)
    ax.set_xlabel("Target-Resolved Phase II-IV Trials")
    ax.set_title("Top 18 Clinical Target Genes Across 12 Axes", pad=8)

    # ED1c: Polypharmacology target count spectrum & single- vs multi-target DNAcc AUROC
    ax = axes1[1, 0]
    _panel_tag(ax, "c", x=-0.15, y=1.06)
    sub_sg = s3f[
        (s3f["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes")
        & (s3f["cv_scheme"].str.startswith("Subgroup_"))
    ].copy()
    nt_bins = s2["n_targets"].clip(upper=5).value_counts().sort_index()
    xb = np.arange(len(nt_bins))
    ax.bar(xb, nt_bins.to_numpy(), color=PAL["sky"], edgecolor=PAL["navy"], lw=0.8, width=0.58)
    for i_b, v_b in enumerate(nt_bins.to_numpy()):
        ax.text(i_b, v_b + 300, f"{int(v_b):,}", ha="center", va="bottom", fontsize=6.5)
    if len(sub_sg) >= 2:
        r_single = sub_sg[sub_sg["cv_scheme"].str.contains("SingleTarget")].iloc[0]
        r_multi = sub_sg[sub_sg["cv_scheme"].str.contains("Polypharmacology")].iloc[0]
        ax.text(
            0.30,
            0.84,
            f"Single-Target DNAcc AUROC = {r_single['auroc']:.3f} [{r_single['auroc_ci_low']:.3f}-{r_single['auroc_ci_high']:.3f}]\n"
            f"Multi-Target DNAcc AUROC = {r_multi['auroc']:.3f} [{r_multi['auroc_ci_low']:.3f}-{r_multi['auroc_ci_high']:.3f}]",
            transform=ax.transAxes,
            fontsize=6.6,
            fontweight="bold",
            color=PAL["navy"],
        )
    ax.set_xticks(xb)
    ax.set_xticklabels(["1 Target", "2 Targets", "3 Targets", "4 Targets", "5+ Targets"], fontsize=7.0)
    ax.set_ylim(0, float(nt_bins.max()) * 1.28)
    ax.set_ylabel("Target-Resolved Trials")
    ax.set_title("Polypharmacology Distribution & Subgroup Robustness", pad=8)

    # ED1d: Core 6 Single-Cell Axes Stratified vs Unseen-Target GroupKFold Ablation
    ax = axes1[1, 1]
    _panel_tag(ax, "d", x=-0.14, y=1.06)
    ab_6 = s3f[s3f["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_Core6SingleCellAxes"]
    st6 = ab_6[ab_6["cv_scheme"] == "Stratified_5Fold"].sort_values("step_order")
    gp6 = ab_6[ab_6["cv_scheme"] == "UnseenTarget_Group5Fold"].sort_values("step_order")
    x6 = np.arange(len(st6))
    ax.errorbar(x6 - 0.05, st6["auroc"], yerr=[st6["auroc"] - st6["auroc_ci_low"], st6["auroc_ci_high"] - st6["auroc"]], color=PAL["blue"], marker="o", lw=1.6, capsize=2.4, label="Stratified 5-Fold (n = 5,815)")
    ax.errorbar(x6 + 0.05, gp6["auroc"], yerr=[gp6["auroc"] - gp6["auroc_ci_low"], gp6["auroc_ci_high"] - gp6["auroc"]], color=PAL["teal"], marker="s", ls="--", lw=1.5, capsize=2.4, label="Unseen-Target Group5Fold")
    ax.set_xticks(x6)
    ax.set_xticklabels(["M1: τ0", "M2: +Byst.", "M3: +τ_early", "M4: +Δτ", "M5: +Switch", "M6: DNAcc"], rotation=25, ha="right", fontsize=6.8)
    ax.set_ylim(0.41, float(st6["auroc_ci_high"].max()) + 0.05)
    ax.set_ylabel("Core 6-Axis AUROC (95% CI)")
    ax.set_title("Core 6 Single-Cell Axes Unseen-Target Generalization", pad=8)
    ax.legend(loc="lower right", frameon=False, fontsize=6.6)

    _save_multi(fig1, out_dir, "Extended_Data_Fig1_PRISMA_Target_Modality_Audit", artifact_dir)

    # -----------------------------------------------------------------------
    # Extended Data Fig. 2: 6-Axis Single-Cell Donor Architecture, Switch Frequency & Per-Axis AUROC
    # -----------------------------------------------------------------------
    fig2, axes2 = plt.subplots(2, 2, figsize=(13.8, 9.2), facecolor="white")
    fig2.subplots_adjust(wspace=0.36, hspace=0.42)

    # ED2a: Primary human single-cell / spatial sample & cell census
    ax = axes2[0, 0]
    _panel_tag(ax, "a", x=-0.16, y=1.06)
    coh_census = [
        ("Hepatic MASH/Cirrhosis (GSE136103+212837)", 17),
        ("Pulmonary IPF (GSE135893)", 30),
        ("Renal CKD/DKD (GSE271188)", 17),
        ("Solid Tumor NSCLC+BCC (GSE207422+123813)", 46),
        ("Ulcerative Colitis Biopsy (GSE250487)", 28),
        ("Longitudinal Sepsis (GSE167363)", 12),
        ("10x Visium Carcinoma Sections (12 Sections)", 12),
    ][::-1]
    yc2 = np.arange(len(coh_census))
    ax.barh(yc2, [v for _, v in coh_census], color=PAL["teal"], height=0.62, alpha=0.85)
    for i_c, (_, v_c) in enumerate(coh_census):
        ax.text(v_c + 0.8, i_c, f"n = {v_c}", va="center", ha="left", fontsize=6.6, fontweight="bold")
    ax.set_yticks(yc2)
    ax.set_yticklabels([lbl for lbl, _ in coh_census], fontsize=6.7)
    ax.set_xlim(0, 56)
    ax.set_xlabel("Primary Human Donors / Tissue Biopsies")
    ax.set_title("Primary Human Single-Cell & Spatial Cohort Census", pad=8)

    # ED2b: Ectopic lineage compartment-switch rate across 6 single-cell axes
    ax = axes2[0, 1]
    _panel_tag(ax, "b", x=-0.16, y=1.06)
    sw_rates = atlas.groupby("disease_axis")["compartment_switch"].mean().mul(100.0).sort_values()
    ysw = np.arange(len(sw_rates))
    ax.barh(ysw, sw_rates.to_numpy(), color=PAL["crimson"], alpha=0.82, height=0.62)
    for i_s, v_s in enumerate(sw_rates.to_numpy()):
        ax.text(v_s + 0.7, i_s, f"{v_s:.1f}%", va="center", ha="left", fontsize=6.6, fontweight="bold")
    ax.set_yticks(ysw)
    ax.set_yticklabels([AXIS_SHORT.get(a, a) for a in sw_rates.index], fontsize=6.9)
    ax.set_xlim(0, float(sw_rates.max()) * 1.25)
    ax.set_xlabel("Expressed Targets Switching Primary Lineage (%)")
    ax.set_title("Ectopic Compartment-Switch Frequency by Axis", pad=8)

    # ED2c: Delta_tau IQR distributions across 6 single-cell axes
    ax = axes2[1, 0]
    _panel_tag(ax, "c", x=-0.16, y=1.06)
    ax_list = list(sw_rates.index)
    data_dt = [atlas[atlas["disease_axis"] == a]["delta_tau"].dropna().to_numpy() for a in ax_list]
    bp = ax.boxplot(data_dt, vert=False, patch_artist=True, showfliers=False, widths=0.55)
    for patch in bp["boxes"]:
        patch.set_facecolor(PAL["sky"])
        patch.set_alpha(0.65)
        patch.set_edgecolor(PAL["navy"])
    for med_l in bp["medians"]:
        med_l.set_color(PAL["crimson"])
        med_l.set_linewidth(1.5)
    ax.axvline(0.0, color=PAL["slate"], ls="--", lw=0.8)
    ax.set_yticks(range(1, len(ax_list) + 1))
    ax.set_yticklabels([AXIS_SHORT.get(a, a) for a in ax_list], fontsize=6.9)
    ax.set_xlabel("Single-Cell Specificity Drift (Δτ = τ_late - τ_early)")
    ax.set_title("Target Specificity Drift Distributions Across 6 Axes", pad=8)

    # ED2d: Per-Axis DNAcc vs Static tau0 AUROC with 2,000-Bootstrap 95% CIs
    ax = axes2[1, 1]
    _panel_tag(ax, "d", x=-0.16, y=1.06)
    s3b_s = s3b.sort_values("dnacc_auroc", ascending=True).reset_index(drop=True)
    y_ax = np.arange(len(s3b_s))
    ax.axvline(0.50, color=PAL["slate"], ls=":", lw=0.8)
    for i_x, r_x in s3b_s.iterrows():
        ax.plot([r_x["static_tau0_ci_low"], r_x["static_tau0_ci_high"]], [i_x - 0.13, i_x - 0.13], color=PAL["grey"], lw=1.3)
        ax.scatter([r_x["static_tau0_auroc"]], [i_x - 0.13], color=PAL["grey"], marker="s", s=24, label="Static τ0" if i_x == 0 else None)
        ax.plot([r_x["dnacc_ci_low"], r_x["dnacc_ci_high"]], [i_x + 0.13, i_x + 0.13], color=PAL["blue"], lw=1.6)
        ax.scatter([r_x["dnacc_auroc"]], [i_x + 0.13], color=PAL["blue"], marker="o", s=30, label="DNAcc" if i_x == 0 else None)
        ax.text(float(s3b_s["dnacc_ci_high"].max()) + 0.012, i_x, f"+{r_x['delta_auroc']:.3f}", va="center", ha="left", fontsize=6.3, color=PAL["navy"])
    ax.set_yticks(y_ax)
    ax.set_yticklabels([AXIS_SHORT.get(a, a) for a in s3b_s["disease_axis"]], fontsize=6.9)
    ax.set_xlim(0.32, float(s3b_s["dnacc_ci_high"].max()) + 0.095)
    ax.set_xlabel("Per-Axis Clinical Maturation AUROC (95% CI)")
    ax.set_title("Per-Axis Predictive Discrimination (n = 6 Core Axes)", pad=8)
    ax.legend(loc="upper left", frameon=False, fontsize=6.6)

    _save_multi(fig2, out_dir, "Extended_Data_Fig2_SixAxis_SingleCell_Drift_and_Ablation", artifact_dir)

    # -----------------------------------------------------------------------
    # Extended Data Fig. 3: Profile Likelihood, Multivariable Control, Portfolio Enrichment & Scientific Calibration
    # -----------------------------------------------------------------------
    fig3, axes3 = plt.subplots(2, 2, figsize=(14.4, 9.2), facecolor="white")
    fig3.subplots_adjust(wspace=0.45, hspace=0.42)

    # ED3a: 25-Point Continuous Bilateral Profile-Likelihood Identifiability Curves
    ax = axes3[0, 0]
    _panel_tag(ax, "a", x=-0.15, y=1.06)
    edge_cols = [PAL["blue"], PAL["crimson"], PAL["teal"], PAL["orange"]]
    edge_short = {
        "v_Parenchyma_Injury_to_SPP1_SAMac": "Parenchyma → SAMac",
        "v_SPP1_SAMac_to_Myofibroblast": "SAMac → Myofibroblast",
        "v_Myofibroblast_to_SPP1_SAMac": "Myofibroblast → SAMac",
        "v_Myofibroblast_to_Parenchyma_Injury": "Myofibroblast → Parenchyma",
    }
    for (p_name, p_sub), col_e in zip(s4g.groupby("parameter", sort=False), edge_cols):
        p_sub_s = p_sub.sort_values("multiplier")
        ax.plot(p_sub_s["multiplier"], p_sub_s["delta_chi2"], color=col_e, lw=1.7, label=edge_short.get(p_name, p_name))
    ax.axhline(3.84, color=PAL["slate"], ls="--", lw=1.0, label="95% CI Cutoff (Δχ² = 3.84)")
    ax.set_xlabel("Parameter Fold Multiplier Relative to MLE (v / v_MLE)")
    ax.set_ylabel("Profile Likelihood Δχ² (Nuisance-Reoptimized)")
    ax.set_ylim(-1.5, 68.0)
    ax.set_title("25-Point Structural Parameter Identifiability", pad=8)
    ax.legend(loc="upper center", bbox_to_anchor=(0.51, 0.99), frameon=False, fontsize=5.7, ncol=1)

    # ED3b: Multivariable Fixed-Effects & Joint Head-to-Head Control
    ax = axes3[0, 1]
    _panel_tag(ax, "b", x=-0.16, y=1.06)
    fe_sub = s3g[s3g["cohort"] == "Deconfounded_Scientific_Ph34_vs_SciFail_All12Axes"].iloc[::-1].reset_index(drop=True)
    spec_short = {
        "Unadjusted_Univariate_DNAcc": "DNAcc (Unadjusted)",
        "Unadjusted_Univariate_Static_Tau0": "Static τ0 (Unadjusted)",
        "Multivariable_FixedEffects_DNAcc": "DNAcc (Adj. Axis+Mod+Era+Spon)",
        "Multivariable_FixedEffects_Static_Tau0": "Static τ0 (Adj. Axis+Mod+Era+Spon)",
        "Joint_FixedEffects_DNAcc_controlling_for_Tau0": "Joint Model: DNAcc (+τ0 +FE)",
        "Joint_FixedEffects_Tau0_controlling_for_DNAcc": "Joint Model: Static τ0 (+DNAcc +FE)",
    }
    yf3 = np.arange(len(fe_sub))
    ax.axvline(1.0, color=PAL["slate"], ls="--", lw=0.9)
    for i_f, r_f in fe_sub.iterrows():
        is_dna = "DNAcc" in str(r_f["model_specification"]) and not str(r_f["model_specification"]).startswith("Joint_FixedEffects_Tau0")
        c_f = PAL["blue"] if is_dna else PAL["grey"]
        ax.plot([r_f["ci95_low"], r_f["ci95_high"]], [i_f, i_f], color=c_f, lw=1.6)
        ax.scatter([r_f["odds_ratio_per_sd"]], [i_f], color=c_f, s=30, zorder=4)
        ax.text(float(fe_sub["ci95_high"].max()) * 1.03, i_f, f"OR {r_f['odds_ratio_per_sd']:.2f}", va="center", ha="left", fontsize=6.3)
    ax.set_yticks(yf3)
    ax.set_yticklabels([spec_short.get(s, s) for s in fe_sub["model_specification"]], fontsize=6.6)
    ax.set_xlim(0.75, float(fe_sub["ci95_high"].max()) * 1.28)
    ax.set_xlabel("Odds Ratio per 1-SD (95% CI; n = 11,843)")
    ax.set_title("Multivariable Fixed-Effects & Joint Control", pad=8)

    # ED3c: Simulated Phase II -> III Portfolio Enrichment & Spared Patients
    ax = axes3[1, 0]
    _panel_tag(ax, "c", x=-0.15, y=1.06)
    xg3 = s3j["bottom_percentile_gated_out"].to_numpy()
    ax.plot(xg3, s3j["dnacc_success_rate_pct"], color=PAL["blue"], lw=1.8, marker="o", ms=4.2, label="DNAcc-Gated Success Rate (%)")
    ax.plot(xg3, s3j["static_tau0_success_rate_pct"], color=PAL["grey"], lw=1.5, ls="--", marker="s", ms=3.8, label="Static τ0-Gated Success Rate (%)")
    r25 = s3j[s3j["bottom_percentile_gated_out"] == 25]
    if len(r25) > 0:
        sp_pts = int(r25["patients_spared_futile_exposure"].iloc[0])
        prev_f = int(r25["futile_trials_prevented_dnacc"].iloc[0])
        ax.text(0.04, 0.86, f"At 25% DNAcc Gate: {prev_f:,} Futile Trials Avoided\n(~{sp_pts:,} Patients Spared Exposure)", transform=ax.transAxes, fontsize=6.6, fontweight="bold", color=PAL["navy"])
    ax.set_xlabel("Bottom DNAcc Percentile Gated Out Before Phase III (%)")
    ax.set_ylabel("Phase III/IV Clinical Success Rate (%)")
    ax.set_ylim(float(s3j["static_tau0_success_rate_pct"].min()) * 0.90, float(s3j["dnacc_success_rate_pct"].max()) * 1.15)
    ax.set_title("Simulated Clinical Portfolio Enrichment Utility", pad=8)
    ax.legend(loc="lower right", frameon=False, fontsize=6.6)

    # ED3d: 10-Decile Calibration on Primary Deconfounded Scientific Cohort (n = 11,843)
    ax = axes3[1, 1]
    _panel_tag(ax, "d", x=-0.15, y=1.06)
    cal_s = s3h[s3h["cohort"] == "Deconfounded_Scientific_Outcomes_n11843"]
    cd_s = cal_s[cal_s["model"] == "Full_DNAcc"]
    ct_s = cal_s[cal_s["model"] == "Static_Healthy_Tau0"]
    ax.plot([0.55, 0.99], [0.55, 0.99], color=PAL["slate"], ls=":", lw=0.9, label="Ideal Diagonal")
    ax.errorbar(
        cd_s["mean_predicted_prob"],
        cd_s["observed_event_rate"],
        yerr=[cd_s["observed_event_rate"] - cd_s["obs_ci95_low"], cd_s["obs_ci95_high"] - cd_s["observed_event_rate"]],
        color=PAL["blue"],
        marker="o",
        ms=4.0,
        lw=1.6,
        capsize=2.2,
        label="DNAcc (10 Deciles)",
    )
    ax.plot(ct_s["mean_predicted_prob"], ct_s["observed_event_rate"], color=PAL["vermilion"], marker="s", ms=3.5, ls="--", lw=1.3, label="Static Healthy τ0")
    ax.set_xlabel("Mean Predicted Scientific Success Probability")
    ax.set_ylabel("Observed Scientific Success Rate (95% CI)")
    ax.set_title("Reliability Calibration (Scientific Cohort, n = 11,843)", pad=8)
    ax.legend(loc="lower right", frameon=False, fontsize=6.6)

    _save_multi(fig3, out_dir, "Extended_Data_Fig3_Calibration_DCA_FixedEffects_Sensitivity", artifact_dir)

    # -----------------------------------------------------------------------
    # Extended Data Fig. 4: All 12 Visium Sections, 4x4 Proximity Matrix, FACS Audit & Platform Scaling
    # -----------------------------------------------------------------------
    fig4, axes4 = plt.subplots(2, 2, figsize=(13.8, 9.2), facecolor="white")
    fig4.subplots_adjust(wspace=0.36, hspace=0.42)

    # ED4a: All 12 Human Solid-Tumor 10x Visium Sections Stromal Exclusion Index (SEI)
    ax = axes4[0, 0]
    _panel_tag(ax, "a", x=-0.16, y=1.06)
    v12_s = s5a.sort_values("stromal_exclusion_index", ascending=True).reset_index(drop=True)
    yv = np.arange(len(v12_s))
    ax.barh(yv, v12_s["stromal_exclusion_index"], color=PAL["crimson"], alpha=0.84, height=0.64)
    for i_v, r_v in v12_s.iterrows():
        ax.text(float(r_v["stromal_exclusion_index"]) + 0.004, i_v, f"+{r_v['stromal_exclusion_index']:.3f} ({r_v['n_spots']:,} spots)", va="center", ha="left", fontsize=6.2)
    ax.set_yticks(yv)
    ax.set_yticklabels([f"{r.section_id.replace('GSM6177', 'GSM').replace('GSE171351_', 'Bladder_')}" for r in v12_s.itertuples()], fontsize=6.4)
    ax.set_xlim(0, float(v12_s["stromal_exclusion_index"].max()) * 1.44)
    ax.set_xlabel("Stromal Exclusion Index (SEI = Prox[myCAF-TAM] - Prox[Tumor-CD8])")
    ax.set_title("Per-Section Stromal Exclusion Across 12 Human Visium Tumors", pad=8)

    # ED4b: Consensus 4x4 Bivariate Moran Spatial Proximity Matrix (21,796 Visium Spots)
    ax = axes4[0, 1]
    _panel_tag(ax, "b", x=-0.15, y=1.06)
    pmat = s5b.to_numpy(dtype=float)
    im4 = ax.imshow(pmat, cmap="PuBuGn", vmin=0.35, vmax=0.85, aspect="auto")
    for r_i in range(pmat.shape[0]):
        for c_i in range(pmat.shape[1]):
            ax.text(c_i, r_i, f"{pmat[r_i, c_i]:.3f}", ha="center", va="center", fontsize=6.8, fontweight="bold", color="white" if pmat[r_i, c_i] >= 0.65 else PAL["ink"])
    lbl4 = [c.replace("_", " ") for c in s5b.columns]
    ax.set_xticks(range(len(lbl4)))
    ax.set_xticklabels(lbl4, fontsize=6.9)
    ax.set_yticks(range(len(lbl4)))
    ax.set_yticklabels(lbl4, fontsize=6.9)
    ax.set_title("Consensus k = 6 Spatial Proximity Matrix (12 Sections)", pad=8)
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(True)
        ax.spines[sp].set_linewidth(0.6)
    cb4 = fig4.colorbar(im4, ax=ax, fraction=0.042, pad=0.03)
    cb4.ax.tick_params(labelsize=6.3)

    # ED4c: GSE123813 FACS CD3+ Sort Dilution vs Non-CD3+ Deconfounded Stromal Fraction
    ax = axes4[1, 0]
    _panel_tag(ax, "c", x=-0.15, y=1.06)
    xf = np.arange(len(s5h))
    wf_b = 0.35
    ax.bar(xf - wf_b / 2, 100.0 * s5h["naive_pooled_caf_frac"], width=wf_b, color=PAL["lightgrey"], edgecolor=PAL["slate"], lw=0.8, label="Naive All-Cell Pool (CD3+-Diluted)")
    ax.bar(xf + wf_b / 2, 100.0 * s5h["deconfounded_nonT_caf_frac"], width=wf_b, color=PAL["crimson"], alpha=0.82, label="FACS-Deconfounded Non-CD3+ TME")
    med_fold = float(s5h["caf_dilution_fold"].median())
    max_fold = float(s5h["caf_dilution_fold"].max())
    ax.text(0.03, 0.83, f"Median Recovery = {med_fold:.1f}x\n(Max {max_fold:.1f}x After CD3+ Deconfounding)", transform=ax.transAxes, fontsize=6.4, fontweight="bold", color=PAL["navy"])
    ax.set_xticks(xf)
    ax.set_xticklabels([f"{r.patient_id}\n({'R' if r.responder == 1 else 'NR'})" for r in s5h.itertuples()], fontsize=6.7)
    ax.set_ylim(0, float(100.0 * s5h["deconfounded_nonT_caf_frac"].max()) * 1.34)
    ax.set_ylabel("Pre-Treatment myCAF Fraction (%)")
    ax.set_title("GSE123813 FACS CD45+CD3+ Enrichment Deconfounding (n = 8)", pad=8)
    ax.legend(loc="upper right", frameon=False, fontsize=6.5)

    # ED4d: Platform-Invariant Baseline Niche State Alignment (scRNA BCC vs Bulk Pre-Biopsy NSCLC)
    ax = axes4[1, 1]
    _panel_tag(ax, "d", x=-0.15, y=1.06)
    grps_d = [
        ("BCC_scRNA_GSE123813", 1, "BCC scRNA\nResp.", PAL["teal"]),
        ("BCC_scRNA_GSE123813", 0, "BCC scRNA\nNon-Resp.", PAL["crimson"]),
        ("NSCLC_PreBiopsy_GSE207422", 1, "NSCLC Bulk\nResp.", PAL["blue"]),
        ("NSCLC_PreBiopsy_GSE207422", 0, "NSCLC Bulk\nNon-Resp.", PAL["vermilion"]),
    ]
    rng_d = np.random.default_rng(99)
    for idx_d, (ckey, rv, lbl_d, col_d) in enumerate(grps_d):
        vals_d = 100.0 * s5i[(s5i["cancer_cohort"] == ckey) & (s5i["responder"] == rv)]["simulated_icb_remission_prob"].to_numpy()
        mu_d = float(np.mean(vals_d))
        se_d = float(np.std(vals_d, ddof=1) / math.sqrt(max(len(vals_d), 1))) if len(vals_d) > 1 else 0.0
        ax.bar(idx_d, mu_d, width=0.56, color=col_d, alpha=0.45, edgecolor=col_d, lw=1.0, yerr=se_d, capsize=2.8)
        ax.scatter(np.full(len(vals_d), idx_d) + rng_d.uniform(-0.12, 0.12, size=len(vals_d)), vals_d, color=col_d, s=22, edgecolors="white", linewidths=0.4, zorder=4)
    ax.axhline(50.0, color=PAL["slate"], ls="--", lw=0.8)
    ax.set_xticks(range(4))
    ax.set_xticklabels([g[2] for g in grps_d], fontsize=6.8)
    ax.set_ylim(0, 118)
    ax.set_ylabel("Zero-Refit Simulated Remission Probability (%)")
    ax.set_title("Cross-Platform Pre-Treatment Responder Alignment (n = 32)", pad=8)

    _save_multi(fig4, out_dir, "Extended_Data_Fig4_Visium12_FACS_and_Platform_Scaling", artifact_dir)
