"""Round-3 (Station 4) Orthogonal Validations & 5-Piece Nature Medicine Submission Package Tests."""
from __future__ import annotations

import csv
import json
import pathlib
import zipfile

import openpyxl
from docx import Document

ROOT = pathlib.Path(__file__).resolve().parents[1]
TAB_DIR = ROOT / "outputs" / "stage2_tables"
MAIN_FIG_DIR = ROOT / "outputs" / "stage2_main_figures"
SUPP_FIG_DIR = ROOT / "outputs" / "stage2_supp_figures"
MS_DIR = ROOT / "manuscript"


def _read_csv_dicts(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def test_r3_orthogonal_core6_lodo_and_shrinkage_table() -> None:
    rows = _read_csv_dicts(TAB_DIR / "Table_S6a_Core6_Regularized_Shrinkage_and_LODO_Axis.csv")
    assert len(rows) == 17
    by_stratum = {r["stratum_or_heldout_axis"]: r for r in rows}
    nested_l1 = by_stratum["Core6_Group5Fold_M6_NestedCV_L1_Shrinkage"]
    nested_l2 = by_stratum["Core6_Group5Fold_M6_L2_Ridge_C0.10"]
    pooled_lodo = by_stratum["Pooled_6Axis_LODO_M6_L2_C0.25"]
    assert int(pooled_lodo["n_trials"]) == 5662
    assert int(pooled_lodo["n_failures"]) == 386
    assert abs(float(nested_l1["auroc"]) - 0.6325) < 1e-3
    assert abs(float(nested_l2["auroc"]) - 0.6363) < 1e-3
    assert abs(float(pooled_lodo["auroc"]) - 0.6306) < 1e-3
    assert float(pooled_lodo["delta_auroc_vs_m1"]) > 0.14
    assert float(pooled_lodo["delong_p_vs_m1"]) < 1e-14


def test_r3_orthogonal_antifibrotic24_chembl_and_looo_platt() -> None:
    rows = _read_csv_dicts(TAB_DIR / "Table_S6b_Antifibrotic24_ChEMBL_Potency_and_LOOO_Calibration.csv")
    trial_rows = [r for r in rows if r["record_type"] == "Trial_Level_24"]
    summary_rows = {r["nct_id"]: r for r in rows if r["record_type"] == "Model_Regime_Summary"}
    assert len(trial_rows) == 24
    pxc50_vals = [float(r["chembl_pxc50"]) for r in trial_rows]
    assert min(pxc50_vals) == 3.52
    assert max(pxc50_vals) == 10.00
    a015 = summary_rows["Summary_ChEMBL_pXC50_Alpha0.15_Raw"]
    looo = summary_rows["Summary_LOOO_Platt_Recalibrated_Uniform"]
    assert abs(float(a015["summary_auroc"]) - 0.8643) < 1e-3
    assert float(a015["summary_perm_p"]) < 0.005
    assert abs(float(looo["summary_brier"]) - 0.1735) < 1e-3
    assert float(looo["bal_acc_at_050"]) > 0.75


def test_r3_orthogonal_visium12_admixture_and_toroidal_null() -> None:
    rows = _read_csv_dicts(TAB_DIR / "Table_S6c_Visium12_Admixture_Residualization_and_Toroidal_Null.csv")
    sec_rows = [r for r in rows if r["dataset"] != "All_12_Sections_Meta"]
    meta = next(r for r in rows if r["dataset"] == "All_12_Sections_Meta")
    assert len(sec_rows) == 12
    assert int(meta["n_spots_qc"]) == 20880
    assert abs(float(meta["sei_retention_pct_after_admixture_control"]) - 48.3) < 0.5
    assert float(meta["toroidal_z_score"]) > 7.80
    assert float(meta["toroidal_empirical_p"]) < 1e-14


def test_r3_orthogonal_oncology_loco_and_multivariable_cox_evalue() -> None:
    rows = _read_csv_dicts(TAB_DIR / "Table_S6d_Oncology107_LOCO_and_Survival74_Multivariable_Cox.csv")
    loco_rows = [r for r in rows if r["analysis_section"] == "LOCO_4Cohort_PreTreatment_107"]
    cox_m2 = [
        r for r in rows
        if r["analysis_section"] == "Survival_74_Multivariable_Cox_and_EValue"
        and r["model_or_cohort"] == "Cox_Model2_FullMultivariable_BinaryBarrier"
    ]
    cox_m3 = [
        r for r in rows
        if r["analysis_section"] == "Survival_74_Multivariable_Cox_and_EValue"
        and r["model_or_cohort"] == "Cox_Model3_FullMultivariable_ContinuousBarrier"
    ]
    assert len(loco_rows) == 5
    barrier_row = next(r for r in cox_m2 if r["covariate"] == "High_Stromal_Barrier_Load")
    cont_row = next(r for r in cox_m3 if r["covariate"] == "Continuous_CAF_Barrier_per_SD")
    assert int(barrier_row["n_patients"]) == 74
    assert int(barrier_row["n_events_or_responders"]) == 43
    assert abs(float(barrier_row["hazard_ratio"]) - 3.02) < 0.05
    assert float(barrier_row["p_value"]) < 0.01
    assert abs(float(barrier_row["e_value_point"]) - 3.68) < 0.05
    assert float(barrier_row["e_value_ci95_bound"]) > 1.70
    assert float(cont_row["p_value"]) > 0.50


def test_r3_orthogonal_wetlab_perturbation_and_nonicb_control() -> None:
    rows = _read_csv_dicts(TAB_DIR / "Table_S6e_WetLab_TGFBR1i_Perturbation_and_NonICB_Control.csv")
    assert len(rows) == 87

    sample_rows = [r for r in rows if r["validation_pillar"] == "WetLab_Sample_Trajectory"]
    assert len(sample_rows) == 73

    pfd = next(
        r for r in rows
        if r["validation_pillar"] == "WetLab_Perturbation_Summary"
        and r["subgroup_or_covariate"] == "Locked_21Gene_Myofibroblast_PFD_vs_Activated_FBS"
    )
    repsox = next(
        r for r in rows
        if r["validation_pillar"] == "WetLab_Perturbation_Summary"
        and r["subgroup_or_covariate"] == "Core_14Gene_ECM_Myofibroblast_RepSox_vs_TGFB1"
    )
    nint = next(
        r for r in rows
        if r["validation_pillar"] == "WetLab_Perturbation_Summary"
        and r["subgroup_or_covariate"] == "Core_14Gene_ECM_Myofibroblast_Nintedanib_vs_DMSO"
    )
    assert abs(float(pfd["metric_primary"]) - (-0.4779)) < 0.01
    assert float(pfd["p_value"]) < 1e-5
    assert abs(float(repsox["metric_primary"]) - (-0.3021)) < 0.01
    assert float(repsox["p_value"]) < 0.002
    assert abs(float(nint["metric_primary"]) - (-0.5382)) < 0.01
    assert float(nint["p_value"]) < 0.02

    non_icb = next(r for r in rows if r["cohort_or_model"] == "Pooled_NonICB_Melanoma_Control_n264")
    inter = next(r for r in rows if r["cohort_or_model"] == "TwoRegime_AntiPD1_vs_NonICB_Interaction_n338")
    consil = next(r for r in rows if r["validation_pillar"] == "MultiScale_24Trial_DNAcc_Langevin_Consilience")
    assert int(non_icb["n_samples_or_patients"]) == 264
    assert abs(float(non_icb["metric_primary"]) - 0.8373) < 0.02
    assert float(non_icb["p_value"]) > 0.30
    assert int(inter["n_samples_or_patients"]) == 338
    assert abs(float(inter["metric_primary"]) - 3.6035) < 0.05
    assert float(inter["p_value"]) < 0.01
    assert abs(float(consil["metric_primary"]) - 0.7861) < 1e-3
    assert float(consil["p_value"]) < 1e-5


def test_r3_figure_ecosystem_and_69_verified_pubmed_refs() -> None:
    main_pngs = sorted(MAIN_FIG_DIR.glob("Fig*.png"))
    supp_pngs = sorted(SUPP_FIG_DIR.glob("Extended_Data_Fig*.png"))
    unified_figs = sorted((ROOT / "outputs" / "stage2_figures").glob("*.*"))
    assert len(main_pngs) == 5
    assert len(supp_pngs) == 6
    assert len(unified_figs) == 33
    assert all(p.stat().st_size > 100_000 for p in main_pngs + supp_pngs)

    ledger = _read_csv_dicts(MS_DIR / "numbers_ledger.csv")
    assert len(ledger) == 23
    metric_ids = {r["metric_id"] for r in ledger}
    assert {
        "core6_lodo_auroc",
        "antifibrotic24_looo_brier",
        "visium12_toroidal_z",
        "surv74_multivar_hr",
        "surv74_evalue",
        "wetlab_pirfenidone_delta",
        "nonicb264_multivar_hr",
        "icb_nonicb338_interaction_hr",
        "consilience24_spearman_rho",
    } <= metric_ids
    for row in ledger:
        assert (ROOT / row["source_csv"]).is_file(), f"Missing source_csv: {row['source_csv']}"
        scr_rel = row["script_line"].split(":L")[0]
        assert (ROOT / scr_rel).is_file(), f"Missing script_line target: {row['script_line']}"

    refs = json.loads((TAB_DIR / "verified_pubmed_references.json").read_text(encoding="utf-8"))
    assert len(refs) == 69
    pmids = {str(r.get("pmid")) for r in refs}
    assert {"42761496", "38374043", "28693043", "36585256", "37620362", "36055997", "25909218", "20460471"} <= pmids


def test_s5_nine_station_ledgers_and_track_contract() -> None:
    state = json.loads((MS_DIR / "MANUSCRIPT_STATE.json").read_text(encoding="utf-8"))
    assert state["research_track"] == "multiomics_spatial_mechanism"
    assert state["target_journal"] == "Nature Computational Science"
    assert state["story_contract"]["single_protagonist"]

    reg = json.loads((ROOT / "DATASETS_REGISTRY.json").read_text(encoding="utf-8"))
    assert reg["research_track"] == "multiomics_spatial_mechanism"
    prisma = reg["prisma_census_summary"]
    assert prisma["total_unique_records"] == prisma["included_count"] + prisma["excluded_count"]
    assert prisma["included_count"] == 20

    matrix = json.loads((ROOT / "results" / "CLAIM_EVIDENCE_MATRIX.json").read_text(encoding="utf-8"))
    assert matrix["research_track"] == "multiomics_spatial_mechanism"
    assert len(matrix["claims"]) == 6
    assert all(len(c["evidence_lines"]) >= 2 for c in matrix["claims"])


def test_r3_five_piece_submission_package_integrity() -> None:
    state = json.loads((MS_DIR / "MANUSCRIPT_STATE.json").read_text(encoding="utf-8"))
    locked_title = state["story_contract"]["locked_title"].lower()
    forbidden = [s.lower() for s in state["story_contract"]["zero_forbidden_jargon"]]

    r3_status = json.loads((ROOT / "R3_REVISION_STATUS.json").read_text(encoding="utf-8"))
    modules = r3_status.get("modules") or r3_status.get("gates")
    assert modules and all(bool(v) for v in modules.values()), f"Unclosed Round-3 status modules: {modules}"

    for fname, min_media in [
        ("OmniNiche_Flagship_Manuscript_EN.docx", 5),
        ("Supplementary_Information_EN.docx", 6),
        ("Cover_Letter_NatCompSci.docx", 0),
    ]:
        dpath = MS_DIR / fname
        assert dpath.exists()
        doc = Document(str(dpath))
        full_txt = "\n".join(
            [p.text for p in doc.paragraphs]
            + [c.text for t in doc.tables for row in t.rows for c in row.cells]
        )
        assert locked_title in " ".join(full_txt.split()).lower()
        assert "\u2014" not in full_txt
        assert "$" not in full_txt and "**" not in full_txt and "`" not in full_txt
        for fj in forbidden:
            assert fj not in full_txt.lower()
        with zipfile.ZipFile(str(dpath), "r") as zf:
            media_n = sum(1 for n in zf.namelist() if n.startswith("word/media/"))
        assert media_n == min_media

    wb = openpyxl.load_workbook(str(MS_DIR / "Supplementary_Tables.xlsx"), read_only=True)
    assert len(wb.sheetnames) == 13
    assert {
        "README_Index",
        "Supplementary_Table_1",
        "Table_6a_LODO",
        "Table_6b_ChEMBL",
        "Table_6c_Visium",
        "Table_6d_Cox",
        "Table_6e_WetLab_NonICB",
    } <= set(wb.sheetnames)
    wb.close()
