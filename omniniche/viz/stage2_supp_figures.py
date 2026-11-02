"""Extended Data / Supplementary Figure Ecosystem Builder for Nature Computational Science.

Organizes the clean Main vs Extended Data figure hierarchy:
- Main Figures (Fig. 1–5): outputs/stage2_main_figures/
- Extended Data Figures (Extended Data Fig. 1–6): outputs/stage2_supp_figures/
- Full Combined Archive (backwards-compatible with all Stage II pytest suites): outputs/stage2_figures/

Renders Extended Data Fig. 5 (Table_S6a–Table_S6d) and Extended Data Fig. 6 (Table_S6e:
wet-lab human fibroblast/TGFBR1i perturbation triad, non-ICB melanoma survival specificity
control & regime interaction, and 24-trial macro DNAcc -> meso Langevin consilience)
at 183-mm double-column width (600 DPI PNG + PDF + SVG).
"""
from __future__ import annotations

import math
import pathlib
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from omniniche.viz.stage2_flagship_figures import _save_multi, build_extended_data_figures
from omniniche.viz.visual_tokens import (
    AXIS_COMPACT,
    PAL,
    _fmt_p,
    _panel_tag,
    assert_zero_axes_frame_collision,
    set_nature_style,
)


def build_supp_ed_fig5_orthogonal_validations(
    tables_dir: pathlib.Path,
    supp_out_dir: pathlib.Path,
    legacy_out_dir: pathlib.Path | None = None,
) -> pathlib.Path:
    """Render Extended Data Fig. 5 (a-e): Shrinkage/LODO, ChEMBL/LOOO, Visium Toroidal Null, LOCO & Multivariable Cox."""
    set_nature_style()

    s6a = pd.read_csv(tables_dir / "Table_S6a_Core6_Regularized_Shrinkage_and_LODO_Axis.csv")
    s6b = pd.read_csv(tables_dir / "Table_S6b_Antifibrotic24_ChEMBL_Potency_and_LOOO_Calibration.csv")
    s6c = pd.read_csv(tables_dir / "Table_S6c_Visium12_Admixture_Residualization_and_Toroidal_Null.csv")
    s6d = pd.read_csv(tables_dir / "Table_S6d_Oncology107_LOCO_and_Survival74_Multivariable_Cox.csv")

    fig = plt.figure(figsize=(7.20, 5.95), facecolor="white")
    outer = gridspec.GridSpec(
        2,
        1,
        figure=fig,
        height_ratios=[1.0, 1.08],
        hspace=0.46,
        left=0.08,
        right=0.98,
        top=0.94,
        bottom=0.09,
    )
    top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[1.05, 1.15], wspace=0.34)
    bot_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1], width_ratios=[1.22, 0.88, 1.08], wspace=0.40)

    # -----------------------------------------------------------------------
    # Panel a: Core-6 Regularized Shrinkage & 6-Axis LODO Zero-Shot Transfer
    # -----------------------------------------------------------------------
    ax_a = fig.add_subplot(top_gs[0, 0])
    _panel_tag(ax_a, "a", x=-0.18, y=1.09)

    lodo_rows = s6a[s6a["analysis_block"] == "Leave_One_Disease_Axis_Out_ZeroShot"].copy()
    lodo_rows = lodo_rows.sort_values("auroc", ascending=False).reset_index(drop=True)
    x_ax = np.arange(len(lodo_rows))
    m1_lodo = (lodo_rows["auroc"] - lodo_rows["delta_auroc_vs_m1"]).to_numpy(dtype=float)
    m6_lodo = lodo_rows["auroc"].to_numpy(dtype=float)
    ci_l = np.maximum(m6_lodo - lodo_rows["auroc_ci_low"].to_numpy(dtype=float), 0.0)
    ci_h = np.maximum(lodo_rows["auroc_ci_high"].to_numpy(dtype=float) - m6_lodo, 0.0)

    ax_a.axhline(0.50, color=PAL["lightgrey"], ls="--", lw=0.65, zorder=1)
    ax_a.plot(x_ax, m1_lodo, marker="s", ms=3.4, lw=0.85, ls=":", color=PAL["midgrey"], label="Healthy $\\tau_0$ (LODO)")
    ax_a.fill_between(
        x_ax,
        lodo_rows["auroc_ci_low"].to_numpy(dtype=float),
        lodo_rows["auroc_ci_high"].to_numpy(dtype=float),
        color=PAL["sky"],
        alpha=0.22,
        zorder=2,
    )
    ax_a.errorbar(
        x_ax,
        m6_lodo,
        yerr=[ci_l, ci_h],
        fmt="o-",
        ms=4.0,
        lw=1.15,
        capsize=1.8,
        color=PAL["blue"],
        label="Shrunk DNAcc (LODO)",
        zorder=4,
    )
    pool_row = s6a[s6a["stratum_or_heldout_axis"].str.startswith("Pooled_6Axis_LODO")].iloc[0]
    ax_a.text(
        0.03,
        0.06,
        f"Pooled LODO AUROC = {pool_row['auroc']:.4f}\n({_fmt_p(float(pool_row['delong_p_vs_m1']))} vs $\\tau_0$)",
        transform=ax_a.transAxes,
        fontsize=4.9,
        color=PAL["navy"],
        fontweight="bold",
    )
    short_lbls = []
    for raw_s in lodo_rows["stratum_or_heldout_axis"]:
        key = str(raw_s).replace("HeldOut_", "")
        short_lbls.append(AXIS_COMPACT.get(key, key.replace("_", " ")).replace(" ", "\n", 1))
    ax_a.set_xticks(x_ax)
    ax_a.set_xticklabels(short_lbls, fontsize=4.7)
    ax_a.set_ylim(0.22, 0.88)
    ax_a.set_ylabel("Held-Out Organ Axis AUROC")
    ax_a.set_title("Core-6 Zero-Shot LODO Transfer ($n$=5,662)", pad=5)
    ax_a.legend(loc="upper right", frameon=False, fontsize=4.8)

    # -----------------------------------------------------------------------
    # Panel b: 24 Antifibrotic Trials — ChEMBL v33 pXC50 vs LOOO Recalibration
    # -----------------------------------------------------------------------
    ax_b = fig.add_subplot(top_gs[0, 1])
    _panel_tag(ax_b, "b", x=-0.16, y=1.09)

    tr24 = s6b[s6b["record_type"] == "Trial_Level_24"].copy()
    sum24 = s6b[s6b["record_type"] == "Model_Regime_Summary"].set_index("nct_id")

    wins = tr24[tr24["clinical_win"] == 1]
    fails = tr24[tr24["clinical_win"] == 0]

    ax_b.axhline(0.50, color=PAL["slate"], ls="--", lw=0.65, zorder=1)
    ax_b.axhspan(0.50, 1.03, color=PAL["pale"], alpha=0.55, zorder=0)
    ax_b.scatter(
        fails["chembl_pxc50"],
        fails["looo_platt_cal_prob"],
        s=24,
        marker="X",
        color=PAL["crimson"],
        alpha=0.86,
        edgecolors="white",
        linewidths=0.35,
        label="Phase IIb/III Failure ($n$=14)",
        zorder=4,
    )
    ax_b.scatter(
        wins["chembl_pxc50"],
        wins["looo_platt_cal_prob"],
        s=26,
        marker="o",
        color=PAL["teal"],
        alpha=0.90,
        edgecolors="white",
        linewidths=0.40,
        label="Phase IIb/III Success ($n$=10)",
        zorder=5,
    )
    # Annotate key pharmacological exemplars (high pXC50 failures vs moderate pXC50 wins)
    callouts = {
        "Lebrikizumab": (-0.18, -0.085),
        "Tropifexor": (-0.02, -0.075),
        "Pirfenidone": (0.24, -0.085),
        "Resmetirom": (-0.10, 0.06),
        "Dapagliflozin": (-0.20, 0.06),
    }
    for r_b in tr24.itertuples():
        dname = str(r_b.drug).split(" (")[0]
        if dname in callouts:
            dx, dy = callouts[dname]
            ax_b.text(
                float(r_b.chembl_pxc50) + dx,
                float(r_b.looo_platt_cal_prob) + dy,
                dname,
                fontsize=4.6,
                color=PAL["ink"],
                ha="center",
            )

    r_chem = sum24.loc["Summary_ChEMBL_pXC50_Alpha0.15_Raw"]
    r_looo = sum24.loc["Summary_LOOO_Platt_Recalibrated_Uniform"]
    ax_b.text(
        0.03,
        0.84,
        f"ChEMBL $pXC_{{50}}$ AUROC = {r_chem['summary_auroc']:.4f}\nLOOO Cal. Bal. Acc. = {100*r_looo['bal_acc_at_050']:.1f}%",
        transform=ax_b.transAxes,
        fontsize=4.8,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_b.set_xlim(3.0, 10.6)
    ax_b.set_ylim(-0.02, 1.05)
    ax_b.set_xlabel(r"Human Target Potency ($pXC_{50}$, $-\log_{10}$ M)")
    ax_b.set_ylabel("LOOO-Calibrated Remission")
    ax_b.set_title("ChEMBL v33 Potency & LOOO Calibration ($n$=24)", pad=5)
    ax_b.legend(loc="lower left", frameon=False, fontsize=4.7)

    # -----------------------------------------------------------------------
    # Panel c: 12-Section Visium (20,880 Spots) 2D Toroidal Shift Z-Score Forest/Bar Chart
    # -----------------------------------------------------------------------
    ax_c = fig.add_subplot(bot_gs[0, 0])
    _panel_tag(ax_c, "c", x=-0.26, y=1.09)

    sec12 = s6c[s6c["dataset"] != "All_12_Sections_Meta"].sort_values("toroidal_z_score", ascending=True).reset_index(drop=True)
    meta_c = s6c[s6c["dataset"] == "All_12_Sections_Meta"].iloc[0]

    sec_lbls = [
        str(s).replace("GSM6177599_", "").replace("GSM6177601_", "").replace("GSM6177603_", "")
        .replace("GSM6177607_", "").replace("GSM6177609_", "").replace("GSM6177612_", "")
        .replace("GSM6177614_", "").replace("GSM6177617_", "").replace("GSE171351_", "Blad-")
        .replace("NYU_", "").replace("_", "-")
        for s in sec12["section_id"]
    ] + ["Stouffer Meta\n(20,880 Spots)"]
    z_vals = np.append(sec12["toroidal_z_score"].to_numpy(dtype=float), float(meta_c["toroidal_z_score"]))
    yc = np.arange(len(z_vals))
    bar_c = []
    for idx_z, zv in enumerate(z_vals):
        if idx_z == len(z_vals) - 1:
            bar_c.append(PAL["navy"])
        elif zv >= 1.96:
            bar_c.append(PAL["crimson"])
        elif zv >= 0.0:
            bar_c.append(PAL["blue"])
        else:
            bar_c.append(PAL["midgrey"])

    ax_c.axhspan(len(z_vals) - 1.50, len(z_vals) - 0.35, color="#EFF6FF", zorder=0)
    ax_c.axvline(0.0, color=PAL["slate"], ls=":", lw=0.60, zorder=1)
    ax_c.axvline(1.96, color=PAL["crimson"], ls="--", lw=0.70, zorder=2, label="$Z = 1.96$ ($P < 0.05$)")
    ax_c.barh(yc, z_vals, height=0.66, color=bar_c, alpha=0.90, edgecolor="none", zorder=3)

    for idx_z, zv in enumerate(z_vals):
        x_txt = max(zv, 0.0) + 0.16 if zv >= 0.0 else 0.16
        use_pill = 0.25 <= zv < 1.96
        ax_c.text(
            x_txt,
            idx_z,
            f"{zv:+.2f}" if idx_z < len(z_vals) - 1 else f"+{zv:.2f} ({_fmt_p(float(meta_c['toroidal_empirical_p']))})",
            va="center",
            ha="left",
            fontsize=3.95 if idx_z < len(z_vals) - 1 else 4.15,
            fontweight="bold" if (zv >= 1.96 or idx_z == len(z_vals) - 1) else "normal",
            color=PAL["navy"] if idx_z == len(z_vals) - 1 else PAL["ink"],
            bbox={"boxstyle": "round,pad=0.04", "fc": "white", "ec": "none", "alpha": 0.92} if use_pill else None,
            zorder=5,
        )

    ax_c.text(
        0.98,
        0.04,
        f"12/12 Raw SEI > 0\n48.3% Admix.-Retained",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=4.15,
        fontweight="bold",
        color=PAL["teal"],
    )
    ax_c.set_yticks(yc)
    ax_c.set_yticklabels(sec_lbls, fontsize=4.15)
    ax_c.set_xlim(-1.45, 13.6)
    ax_c.set_ylim(-0.65, len(z_vals) - 0.15)
    ax_c.set_xlabel("2D Toroidal Shift $Z$-Score (24,000 Permutations)")
    ax_c.set_title("Visium Toroidal Shift Null ($n$=12)", pad=5)
    ax_c.legend(loc="lower right", bbox_to_anchor=(0.99, 0.20), frameon=False, fontsize=4.10)

    # -----------------------------------------------------------------------
    # Panel d: 4-Cohort Anti-PD-1 Leave-One-Cohort-Out (LOCO) Generalization
    # -----------------------------------------------------------------------
    ax_d = fig.add_subplot(bot_gs[0, 1])
    _panel_tag(ax_d, "d", x=-0.24, y=1.09)

    loco_df = s6d[s6d["analysis_section"] == "LOCO_4Cohort_PreTreatment_107"].copy().reset_index(drop=True)
    xd = np.arange(len(loco_df))
    w_d = 0.25

    ax_d.axhline(0.50, color=PAL["lightgrey"], ls="--", lw=0.65, zorder=1)
    ax_d.bar(xd - w_d, loco_df["pdl1_auroc"], width=w_d, color=PAL["lightgrey"], edgecolor=PAL["slate"], lw=0.45, label="Static PD-L1")
    ax_d.bar(xd, loco_df["gep18_auroc"], width=w_d, color=PAL["orange"], alpha=0.75, label="18-Gene GEP")
    ax_d.bar(xd + w_d, loco_df["auroc_or_cindex"], width=w_d, color=PAL["blue"], alpha=0.88, label="OmniNiche")

    coh_short = ["BCC\n(8)", "NSCLC\n(24)", "Mel78\n(26)", "Mel91\n(49)", "Pooled\n(107)"]
    ax_d.set_xticks(xd)
    ax_d.set_xticklabels(coh_short[: len(xd)], fontsize=4.5)
    ax_d.set_ylim(0.30, 1.02)
    ax_d.set_ylabel("Held-Out Cohort AUROC")
    ax_d.set_title("Pre-ICB 4-Cohort LOCO ($n$=107)", pad=5)
    ax_d.legend(loc="upper right", frameon=False, fontsize=4.5)

    # -----------------------------------------------------------------------
    # Panel e: Multivariable Cox Proportional-Hazards Forest & E-Value (n=74)
    # -----------------------------------------------------------------------
    ax_e = fig.add_subplot(bot_gs[0, 2])
    _panel_tag(ax_e, "e", x=-0.26, y=1.09)

    cox2 = s6d[s6d["model_or_cohort"] == "Cox_Model2_FullMultivariable_BinaryBarrier"].copy().reset_index(drop=True)
    # Display top-to-bottom with High Stromal Barrier Load at top
    cox2 = cox2.iloc[::-1].reset_index(drop=True)
    ye = np.arange(len(cox2))
    hrs = cox2["hazard_ratio"].to_numpy(dtype=float)
    cil = cox2["ci95_low"].to_numpy(dtype=float)
    cih = cox2["ci95_high"].to_numpy(dtype=float)

    ax_e.axvline(1.0, color=PAL["slate"], ls="--", lw=0.65, zorder=1)
    # Single horizontal forest errorbar call across the entire script (<= 1 per figure rule)
    ax_e.errorbar(
        hrs,
        ye,
        xerr=[np.maximum(hrs - cil, 0.0), np.maximum(cih - hrs, 0.0)],
        fmt="o",
        ms=4.0,
        lw=1.05,
        capsize=1.8,
        color=PAL["crimson"],
        ecolor=PAL["ink"],
        zorder=4,
    )
    pretty_cov = {
        "High_Stromal_Barrier_Load": "Stromal Barrier\n(Locked)",
        "Cohort_GSE91061_vs_GSE78220": "Cohort\n(Mel91/78)",
        "Tumor_Burden_per_SD": "Tumor Burden\n(+1 SD)",
        "Static_PDL1_per_SD": "Static PD-L1\n(+1 SD)",
        "Inflamed_GEP18_per_SD": "18-Gene GEP\n(+1 SD)",
    }
    ax_e.set_yticks(ye)
    ax_e.set_yticklabels([pretty_cov.get(str(c), str(c).replace("_", " ")) for c in cox2["covariate"]], fontsize=4.6)
    for idx_e, r_e in cox2.iterrows():
        ax_e.text(
            min(float(r_e["ci95_high"]) + 0.20, 7.15),
            idx_e,
            f"HR {r_e['hazard_ratio']:.2f}\n($P$={r_e['p_value']:.3f})",
            va="center",
            ha="left",
            fontsize=4.2,
            color=PAL["crimson"] if "Barrier" in str(r_e["covariate"]) else PAL["slate"],
            fontweight="bold" if "Barrier" in str(r_e["covariate"]) else "normal",
        )
    top_bar = cox2[cox2["covariate"] == "High_Stromal_Barrier_Load"].iloc[0]
    ax_e.text(
        0.03,
        0.03,
        f"E-value = {top_bar['e_value_point']:.2f} (CI: {top_bar['e_value_ci95_bound']:.2f})",
        transform=ax_e.transAxes,
        fontsize=4.7,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_e.set_xlim(0.25, 9.80)
    ax_e.set_ylim(-0.92, 4.55)
    ax_e.set_xlabel("Multivariable OS Hazard Ratio (95% CI)")
    ax_e.set_title("Multivariable Cox Survival ($n$=74)", pad=5)

    collisions = assert_zero_axes_frame_collision(fig)
    if collisions:
        raise RuntimeError(f"Extended Data Fig. 5 axes collision detected: {collisions}")

    stem = "Extended_Data_Fig5_Orthogonal_Shrinkage_Potency_Toroidal_Cox"
    png_path = _save_multi(fig, supp_out_dir, stem, legacy_out_dir)
    return png_path


def build_supp_ed_fig6_wetlab_perturbation_nonicb_consilience(
    tables_dir: pathlib.Path,
    supp_out_dir: pathlib.Path,
    legacy_out_dir: pathlib.Path | None = None,
) -> pathlib.Path:
    """Render Extended Data Fig. 6 (a-d): Wet-lab TGFBR1i/antifibrotic perturbation, non-ICB melanoma control & 24-trial consilience."""
    set_nature_style()

    s6e = pd.read_csv(tables_dir / "Table_S6e_WetLab_TGFBR1i_Perturbation_and_NonICB_Control.csv")
    s4 = pd.read_csv(tables_dir / "Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv")

    fig = plt.figure(figsize=(7.20, 5.65), facecolor="white")
    gs = gridspec.GridSpec(
        2,
        2,
        figure=fig,
        height_ratios=[1.0, 1.08],
        width_ratios=[1.04, 1.08],
        hspace=0.46,
        wspace=0.38,
        left=0.09,
        right=0.98,
        top=0.93,
        bottom=0.09,
    )

    # -----------------------------------------------------------------------
    # Panel a: Primary Human IPF & Lung Fibroblast Pirfenidone Perturbation (GSE226249, n=27)
    # -----------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    _panel_tag(ax_a, "a", x=-0.19, y=1.09)

    traj_226 = s6e[
        (s6e["validation_pillar"] == "WetLab_Sample_Trajectory")
        & (s6e["cohort_or_model"] == "GSE226249_Primary_Human_IPF_Fibroblasts")
    ].copy().reset_index(drop=True)
    sum_226 = s6e[
        (s6e["validation_pillar"] == "WetLab_Perturbation_Summary")
        & (s6e["subgroup_or_covariate"] == "Locked_21Gene_Myofibroblast_PFD_vs_Activated_FBS")
    ].iloc[0]

    donor_styles = {
        "N1": ("o", PAL["teal"], "Normal Lung (N1)"),
        "PF2": ("s", PAL["blue"], "IPF Donor 1 (PF2)"),
        "PF3": ("^", PAL["crimson"], "IPF Donor 2 (PF3)"),
    }
    state_x = {"Basal_SF": 0.0, "Activated_FBS": 1.0, "Pirfenidone_1mM": 2.0}
    rng_a = np.random.default_rng(226249)

    for donor_id, (mk, col, lbl) in donor_styles.items():
        b_vals = traj_226[traj_226["subgroup_or_covariate"] == f"{donor_id}_Basal_SF"]["metric_primary"].to_numpy(dtype=float)
        a_vals = traj_226[traj_226["subgroup_or_covariate"] == f"{donor_id}_Activated_FBS"]["metric_primary"].to_numpy(dtype=float)
        p_vals = traj_226[traj_226["subgroup_or_covariate"] == f"{donor_id}_Pirfenidone_1mM"]["metric_primary"].to_numpy(dtype=float)
        for rep_i in range(min(len(b_vals), len(a_vals), len(p_vals))):
            jit = float(rng_a.uniform(-0.075, 0.075))
            xs = np.array([0.0 + jit, 1.0 + jit, 2.0 + jit])
            ys = np.array([b_vals[rep_i], a_vals[rep_i], p_vals[rep_i]])
            ax_a.plot(xs, ys, lw=0.65, color=col, alpha=0.42, zorder=2)
            ax_a.scatter(
                xs,
                ys,
                s=17,
                marker=mk,
                color=col,
                alpha=0.86,
                edgecolors="white",
                linewidths=0.35,
                label=lbl if rep_i == 0 else None,
                zorder=4,
            )

    grp_means = [
        traj_226[traj_226["subgroup_or_covariate"].str.endswith(st)]["metric_primary"].mean()
        for st in ["Basal_SF", "Activated_FBS", "Pirfenidone_1mM"]
    ]
    grp_sems = [
        traj_226[traj_226["subgroup_or_covariate"].str.endswith(st)]["metric_primary"].sem()
        for st in ["Basal_SF", "Activated_FBS", "Pirfenidone_1mM"]
    ]
    ax_a.plot([0.0, 1.0, 2.0], grp_means, color=PAL["navy"], lw=1.35, marker="D", ms=4.2, zorder=5, label="Mean ± SEM")
    for xi, mu_i, se_i in zip([0.0, 1.0, 2.0], grp_means, grp_sems):
        ax_a.vlines(xi, mu_i - se_i, mu_i + se_i, color=PAL["navy"], lw=1.15, zorder=5)

    # Significance bracket between Activated (1) and + Pirfenidone (2)
    yb_a = 8.13
    ax_a.plot([1.0, 1.0, 2.0, 2.0], [yb_a - 0.03, yb_a, yb_a, yb_a - 0.03], color=PAL["ink"], lw=0.65, zorder=5)
    ax_a.text(
        1.50,
        yb_a + 0.018,
        f"$\\Delta = {float(sum_226['metric_primary']):.2f}$ ({_fmt_p(float(sum_226['p_value']))})",
        ha="center",
        va="bottom",
        fontsize=4.4,
        fontweight="bold",
        color=PAL["crimson"],
    )
    ax_a.text(
        0.03,
        0.04,
        "ACTA2 −76.7% | POSTN −67.2% CPM",
        transform=ax_a.transAxes,
        fontsize=4.35,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_a.set_xticks([0.0, 1.0, 2.0])
    ax_a.set_xticklabels(["Basal\n(SF, $n$=9)", "Activated\n(20% FBS, $n$=9)", "+ Pirfenidone\n(1 mM, $n$=9)"], fontsize=4.55)
    ax_a.set_xlim(-0.26, 2.26)
    ax_a.set_ylim(7.14, 8.32)
    ax_a.set_ylabel("Locked 21-Gene Myofibroblast ($\\log_2$ CPM+1)")
    ax_a.set_title("Human IPF Fibroblast Pirfenidone (GSE226249)", pad=5)
    ax_a.legend(loc="lower left", bbox_to_anchor=(0.01, 0.13), frameon=False, fontsize=4.15)

    # -----------------------------------------------------------------------
    # Panel b: Selective TGFBR1/ALK5 Blockade (RepSox, GSE154646, n=24) & Nintedanib (GSE210020, n=22)
    # -----------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    _panel_tag(ax_b, "b", x=-0.18, y=1.09)

    t154 = s6e[
        (s6e["validation_pillar"] == "WetLab_Sample_Trajectory")
        & (s6e["cohort_or_model"] == "GSE154646_Human_Carcinoma_RepSox_TGFBR1i")
    ].copy().reset_index(drop=True)
    sum154 = s6e[
        (s6e["validation_pillar"] == "WetLab_Perturbation_Summary")
        & (s6e["subgroup_or_covariate"] == "Core_14Gene_ECM_Myofibroblast_RepSox_vs_TGFB1")
    ].iloc[0]

    t210 = s6e[
        (s6e["validation_pillar"] == "WetLab_Sample_Trajectory")
        & (s6e["cohort_or_model"] == "GSE210020_Primary_Human_Fibroblasts_Nintedanib")
    ].copy().reset_index(drop=True)
    sum210 = s6e[
        (s6e["validation_pillar"] == "WetLab_Perturbation_Summary")
        & (s6e["subgroup_or_covariate"] == "Core_14Gene_ECM_Myofibroblast_Nintedanib_vs_DMSO")
    ].iloc[0]

    ax_b.axhline(0.0, color=PAL["slate"], ls="--", lw=0.60, zorder=1)
    ax_b.axvline(2.60, color=PAL["lightgrey"], ls=":", lw=0.75, zorder=1)

    # Left block: 8 paired triplets in GSE154646 centered on each triplet's TGF-B1 stimulated peak
    dmso_154 = t154[t154["subgroup_or_covariate"] == "DMSO"]["metric_secondary"].to_numpy(dtype=float)
    tgfb_154 = t154[t154["subgroup_or_covariate"] == "TGFB1"]["metric_secondary"].to_numpy(dtype=float)
    reps_154 = t154[t154["subgroup_or_covariate"] == "RepSox_TGFBR1i"]["metric_secondary"].to_numpy(dtype=float)

    d_dmso = dmso_154 - tgfb_154
    d_tgfb = np.zeros_like(tgfb_154)
    d_reps = reps_154 - tgfb_154

    rng_b = np.random.default_rng(154646)
    for idx_t in range(len(tgfb_154)):
        jit = float(rng_b.uniform(-0.07, 0.07))
        xb = np.array([0.0 + jit, 1.0 + jit, 2.0 + jit])
        yb = np.array([d_dmso[idx_t], d_tgfb[idx_t], d_reps[idx_t]])
        ax_b.plot(xb, yb, lw=0.60, color=PAL["slate"], alpha=0.38, zorder=2)
        ax_b.scatter(
            xb,
            yb,
            s=15,
            c=[PAL["midgrey"], PAL["orange"], PAL["blue"]],
            edgecolors="white",
            linewidths=0.35,
            zorder=4,
        )
    mu_154 = [float(np.mean(d_dmso)), 0.0, float(np.mean(d_reps))]
    ax_b.plot([0.0, 1.0, 2.0], mu_154, color=PAL["navy"], lw=1.30, marker="D", ms=3.9, zorder=5)

    # Right block: GSE210020 primary human fibroblasts centered on DMSO vehicle mean (n=22)
    dmso_210 = t210[t210["subgroup_or_covariate"] == "DMSO"]["metric_secondary"].to_numpy(dtype=float)
    nint_210 = t210[t210["subgroup_or_covariate"] == "Nintedanib_5uM"]["metric_secondary"].to_numpy(dtype=float)
    ref_210 = float(np.mean(dmso_210))
    d_c210 = dmso_210 - ref_210
    d_n210 = nint_210 - ref_210

    ax_b.bar(
        [3.25, 4.25],
        [float(np.mean(d_c210)), float(np.mean(d_n210))],
        width=0.52,
        color=[PAL["lightgrey"], PAL["teal"]],
        alpha=0.32,
        edgecolor=[PAL["slate"], PAL["teal"]],
        lw=0.65,
        zorder=2,
    )
    ax_b.scatter(
        3.25 + rng_b.uniform(-0.11, 0.11, size=len(d_c210)),
        d_c210,
        s=15,
        color=PAL["slate"],
        edgecolors="white",
        linewidths=0.35,
        zorder=4,
    )
    ax_b.scatter(
        4.25 + rng_b.uniform(-0.12, 0.12, size=len(d_n210)),
        d_n210,
        s=15,
        color=PAL["teal"],
        edgecolors="white",
        linewidths=0.35,
        zorder=4,
    )
    ax_b.plot([3.25, 4.25], [0.0, float(np.mean(d_n210))], color=PAL["navy"], lw=1.25, marker="D", ms=3.9, zorder=5)

    # Brackets & dataset headers
    ax_b.plot([1.0, 1.0, 2.0, 2.0], [0.22, 0.27, 0.27, 0.22], color=PAL["ink"], lw=0.62)
    ax_b.text(
        1.50,
        0.29,
        f"$\\Delta={float(sum154['metric_primary']):.2f}$\n({_fmt_p(float(sum154['p_value']))})",
        ha="center",
        va="bottom",
        fontsize=4.15,
        fontweight="bold",
        color=PAL["blue"],
    )
    ax_b.plot([3.25, 3.25, 4.25, 4.25], [0.52, 0.57, 0.57, 0.52], color=PAL["ink"], lw=0.62)
    ax_b.text(
        3.75,
        0.59,
        f"$\\Delta={float(sum210['metric_primary']):.2f}$\n({_fmt_p(float(sum210['p_value']))})",
        ha="center",
        va="bottom",
        fontsize=4.15,
        fontweight="bold",
        color=PAL["teal"],
    )
    ax_b.text(1.0, -1.49, "GSE154646 ($n$=24, ALK5i)", ha="center", va="bottom", fontsize=4.2, fontweight="bold", color=PAL["navy"])
    ax_b.text(3.75, -1.49, "GSE210020 ($n$=22)", ha="center", va="bottom", fontsize=4.2, fontweight="bold", color=PAL["teal"])

    ax_b.set_xticks([0.0, 1.0, 2.0, 3.25, 4.25])
    ax_b.set_xticklabels(["DMSO\n(8)", "+TGF-$\\beta$1\n(8)", "+RepSox\n(8)", "DMSO\n(8)", "+Ninted.\n(14)"], fontsize=4.35)
    ax_b.set_xlim(-0.48, 4.78)
    ax_b.set_ylim(-1.56, 0.92)
    ax_b.set_ylabel("$\\Delta$ Core ECM–Myofibroblast ($\\log_2$ CPM+1)")
    ax_b.set_title("Selective TGFBR1 Blockade & Nintedanib ($n$=46)", pad=5)

    # -----------------------------------------------------------------------
    # Panel c: Melanoma Anti-PD-1 vs Non-ICB Survival Specificity & Interaction (n=338, 192 deaths)
    # -----------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    _panel_tag(ax_c, "c", x=-0.19, y=1.09)

    mel = s6e[s6e["validation_pillar"] == "NonICB_Melanoma_Specificity_and_Interaction"].copy()
    row_specs = [
        ("GSE22153_NonICB_StageIV_Melanoma_OS", "Univariate_Locked_Barrier", "GSE22153 Non-ICB Univ.\n(Stage IV, $n$=54)", PAL["slate"], False),
        ("GSE65904_NonICB_Melanoma_DSS", "Univariate_Locked_Barrier", "GSE65904 Non-ICB Univ.\n(DSS, $n$=210)", PAL["slate"], False),
        ("Pooled_NonICB_Melanoma_Control_n264", "Multivariable_CTL_PDL1_Cohort_Adjusted_Barrier", "Pooled Non-ICB Multivar.\n(Adj. CTL/PD-L1, $n$=264)", PAL["teal"], True),
        ("AntiPD1_Monotherapy_Melanoma_n74", "Univariate_Locked_Barrier", "Anti-PD-1 Monother. Univ.\n($n$=74, 43 deaths)", PAL["orange"], True),
        ("AntiPD1_Monotherapy_Melanoma_n74", "Multivariable_GEP18_PDL1_Tumor_Cohort_Adjusted_Barrier", "Anti-PD-1 Multivariable\n(Adj. GEP/PD-L1, $n$=74)", PAL["crimson"], True),
    ]
    yc = np.arange(len(row_specs), dtype=float)
    hrs_c, cil_c, cih_c, pvs_c, lbls_c, cols_c, bolds_c = [], [], [], [], [], [], []
    for coh_k, sub_k, disp_lbl, col_k, is_b in row_specs:
        r_m = mel[(mel["cohort_or_model"] == coh_k) & (mel["subgroup_or_covariate"] == sub_k)].iloc[0]
        hrs_c.append(float(r_m["metric_primary"]))
        cil_c.append(float(r_m["ci95_low"]))
        cih_c.append(float(r_m["ci95_high"]))
        pvs_c.append(float(r_m["p_value"]))
        lbls_c.append(disp_lbl)
        cols_c.append(col_k)
        bolds_c.append(is_b)

    hrs_arr = np.asarray(hrs_c)
    cil_arr = np.asarray(cil_c)
    cih_arr = np.asarray(cih_c)

    ax_c.axhspan(2.52, 4.55, color="#FEF2F2", alpha=0.65, zorder=0)
    ax_c.axvline(1.0, color=PAL["slate"], ls="--", lw=0.65, zorder=1)
    # Single horizontal errorbar call across Extended Data Fig. 6 (<= 1 horizontal errorbar rule)
    ax_c.errorbar(
        hrs_arr,
        yc,
        xerr=[np.maximum(hrs_arr - cil_arr, 0.0), np.maximum(cih_arr - hrs_arr, 0.0)],
        fmt="none",
        lw=1.05,
        capsize=1.8,
        ecolor=PAL["ink"],
        zorder=3,
    )
    ax_c.scatter(hrs_arr, yc, s=26, c=cols_c, edgecolors="white", linewidths=0.45, zorder=4)

    for idx_r in range(len(row_specs)):
        pv_r = pvs_c[idx_r]
        p_txt = f"$P$={pv_r:.3f}" if pv_r < 0.05 else f"$P$={pv_r:.2f}, n.s."
        ax_c.text(
            min(cih_arr[idx_r] + 0.16, 6.95),
            yc[idx_r],
            f"HR {hrs_arr[idx_r]:.2f} ({p_txt})",
            va="center",
            ha="left",
            fontsize=4.15,
            fontweight="bold" if bolds_c[idx_r] else "normal",
            color=cols_c[idx_r],
        )

    inter_row = mel[mel["cohort_or_model"] == "TwoRegime_AntiPD1_vs_NonICB_Interaction_n338"].iloc[0]
    ax_c.text(
        0.10,
        0.03,
        f"Anti-PD-1 $\\times$ Barrier Interaction ($n$=338): HR = {float(inter_row['metric_primary']):.2f} ($P$ = {float(inter_row['p_value']):.4f})",
        transform=ax_c.transAxes,
        fontsize=4.15,
        fontweight="bold",
        color=PAL["navy"],
        bbox={"boxstyle": "round,pad=0.08", "fc": "white", "ec": "none", "alpha": 0.94},
        zorder=5,
    )
    ax_c.set_yticks(yc)
    ax_c.set_yticklabels(lbls_c, fontsize=4.30)
    ax_c.set_xlim(0.25, 10.60)
    ax_c.set_ylim(-0.88, 4.55)
    ax_c.set_xlabel("Locked Stromal Barrier Survival Hazard Ratio (95% CI)")
    ax_c.set_title("Melanoma ICB Specificity vs Non-ICB ($n$=338)", pad=5)

    # -----------------------------------------------------------------------
    # Panel d: 24-Trial Macro DNAcc Law -> Meso Zero-Shot Langevin Consilience
    # -----------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    _panel_tag(ax_d, "d", x=-0.18, y=1.09)

    cons_row = s6e[s6e["validation_pillar"] == "MultiScale_24Trial_DNAcc_Langevin_Consilience"].iloc[0]
    x_dnacc = 0.55 * s4["organ_dynamic_retention"].to_numpy(dtype=float) + 0.45 * s4["scar_barrier_access"].to_numpy(dtype=float)
    y_lang = s4["omniniche_remission_prob"].to_numpy(dtype=float)
    wins = s4["clinical_win"].to_numpy(dtype=int) == 1

    # Linear regression + 95% confidence ribbon
    xg = np.linspace(float(x_dnacc.min()) - 0.025, float(x_dnacc.max()) + 0.025, 120)
    slope_d, inter_d = np.polyfit(x_dnacc, y_lang, 1)
    yg = slope_d * xg + inter_d
    resid_d = y_lang - (slope_d * x_dnacc + inter_d)
    s_err = math.sqrt(float(np.sum(resid_d ** 2)) / max(len(x_dnacc) - 2, 1))
    sxx = float(np.sum((x_dnacc - np.mean(x_dnacc)) ** 2))
    ci_band = 2.0739 * s_err * np.sqrt(1.0 / len(x_dnacc) + ((xg - np.mean(x_dnacc)) ** 2) / max(sxx, 1e-9))

    ax_d.fill_between(xg, np.clip(yg - ci_band, -0.04, 1.05), np.clip(yg + ci_band, -0.04, 1.05), color=PAL["sky"], alpha=0.22, zorder=1)
    ax_d.plot(xg, yg, color=PAL["navy"], lw=1.15, ls="-", zorder=2)

    ax_d.scatter(
        x_dnacc[~wins],
        y_lang[~wins],
        s=25,
        marker="s",
        color=PAL["crimson"],
        alpha=0.86,
        edgecolors="white",
        linewidths=0.45,
        label="Phase IIb/III Failure ($n$=14)",
        zorder=4,
    )
    ax_d.scatter(
        x_dnacc[wins],
        y_lang[wins],
        s=28,
        marker="o",
        color=PAL["teal"],
        alpha=0.90,
        edgecolors="white",
        linewidths=0.45,
        label="Phase IIb/III Winner ($n$=10)",
        zorder=5,
    )

    callouts = {
        "Pamrevlumab": (0.015, 0.075, "left"),
        "Simtuzumab (F4)": (-0.042, 0.115, "center"),
        "Resmetirom": (-0.018, 0.032, "right"),
        "Nintedanib": (-0.010, -0.095, "right"),
        "Semaglutide": (-0.016, -0.075, "right"),
        "Pirfenidone": (0.010, -0.088, "center"),
    }
    for idx_t, r_t in s4.iterrows():
        d_name = str(r_t["drug"])
        if d_name in callouts:
            dx_c, dy_c, ha_c = callouts[d_name]
            short_d = d_name.replace(" (F4)", "")
            ax_d.annotate(
                short_d,
                xy=(x_dnacc[idx_t], y_lang[idx_t]),
                xytext=(x_dnacc[idx_t] + dx_c, y_lang[idx_t] + dy_c),
                fontsize=4.10,
                fontweight="bold" if bool(r_t["clinical_win"]) else "normal",
                color=PAL["teal"] if bool(r_t["clinical_win"]) else PAL["crimson"],
                ha=ha_c,
                va="center",
                arrowprops={"arrowstyle": "-", "lw": 0.45, "color": PAL["midgrey"]},
                zorder=6,
            )

    ax_d.text(
        0.03,
        0.95,
        f"Spearman $\\rho$ = {float(cons_row['metric_primary']):.2f} ({_fmt_p(float(cons_row['p_value']))})\n"
        f"Pearson $r$ = {float(cons_row['metric_secondary']):.2f} ($P$ = {float(cons_row['ci95_high']):.4f})",
        transform=ax_d.transAxes,
        va="top",
        ha="left",
        fontsize=4.35,
        fontweight="bold",
        color=PAL["navy"],
    )
    ax_d.set_xlim(0.28, 0.84)
    ax_d.set_ylim(-0.05, 1.08)
    ax_d.set_xlabel("Macro Structural DNAcc Index (Retention + Scar Access)")
    ax_d.set_ylabel("Zero-Shot Meso Langevin Remission Prob.")
    ax_d.set_title("Macro DNAcc $\\to$ Meso Langevin Consilience ($n$=24)", pad=5)
    ax_d.legend(loc="center left", bbox_to_anchor=(0.02, 0.68), frameon=False, fontsize=4.20)

    collisions = assert_zero_axes_frame_collision(fig)
    if collisions:
        raise RuntimeError(f"Extended Data Fig. 6 axes collision detected: {collisions}")

    stem = "Extended_Data_Fig6_WetLab_Perturbation_NonICB_Control_Consilience"
    png_path = _save_multi(fig, supp_out_dir, stem, legacy_out_dir)
    return png_path


def organize_and_render_r3_figure_ecosystem(root_dir: pathlib.Path | None = None) -> dict[str, list[str]]:
    """Sync Main Figs (1-5) into outputs/stage2_main_figures/ and render ED Figs (1-6) into outputs/stage2_supp_figures/."""
    root = root_dir or pathlib.Path("F:/OmniNiche_Universal")
    tables_dir = root / "outputs" / "stage2_tables"
    legacy_dir = root / "outputs" / "stage2_figures"
    main_dir = root / "outputs" / "stage2_main_figures"
    supp_dir = root / "outputs" / "stage2_supp_figures"

    main_dir.mkdir(parents=True, exist_ok=True)
    supp_dir.mkdir(parents=True, exist_ok=True)

    # Copy existing verified Main Figures (Fig1..Fig5) into dedicated outputs/stage2_main_figures/
    for p in sorted(legacy_dir.glob("Fig*.*")):
        shutil.copy2(p, main_dir / p.name)

    # Copy existing verified Extended Data Figs (1..4) into dedicated outputs/stage2_supp_figures/
    # (or rebuild if missing)
    if len(list(legacy_dir.glob("Extended_Data_Fig[1-4]*.png"))) < 4:
        build_extended_data_figures(tables_dir, legacy_dir)
    for p in sorted(legacy_dir.glob("Extended_Data_Fig[1-4]*.*")):
        shutil.copy2(p, supp_dir / p.name)

    # Render Extended Data Fig. 5 & Extended Data Fig. 6 into BOTH outputs/stage2_supp_figures/ and outputs/stage2_figures/
    ed5_png = build_supp_ed_fig5_orthogonal_validations(tables_dir, supp_dir, legacy_dir)
    print("[Module 4] Rendered Extended Data Fig. 5 ->", ed5_png)
    ed6_png = build_supp_ed_fig6_wetlab_perturbation_nonicb_consilience(tables_dir, supp_dir, legacy_dir)
    print("[Module 5] Rendered Extended Data Fig. 6 ->", ed6_png)

    return {
        "main_pngs": sorted(p.name for p in main_dir.glob("*.png")),
        "supp_pngs": sorted(p.name for p in supp_dir.glob("*.png")),
    }


if __name__ == "__main__":
    info = organize_and_render_r3_figure_ecosystem()
    print("Main PNGs:", len(info["main_pngs"]), info["main_pngs"])
    print("Supp PNGs:", len(info["supp_pngs"]), info["supp_pngs"])
