"""Execute Stage II Task 1: OpenTargets Drug-Target Dictionary + Exhaustive Non-CNS ClinicalTrials.gov v2 Census."""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(r"F:\OmniNiche_Universal")
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omniniche.stage2.macro_clinical_census import (
    build_opentargets_drug_target_dictionary,
    harvest_non_cns_clinical_trials_macro_atlas,
)


def main() -> None:
    vault_ct = ROOT / "data_vault" / "raw_clinical_trials"
    out_tables = ROOT / "outputs" / "stage2_tables"
    vault_ct.mkdir(parents=True, exist_ok=True)
    out_tables.mkdir(parents=True, exist_ok=True)

    ot_map_path = vault_ct / "opentargets_drug_target_map.json"
    if ot_map_path.exists() and ot_map_path.stat().st_size > 50_000:
        lookup = json.loads(ot_map_path.read_text(encoding="utf-8"))
        print(f"[OpenTargets] Loaded cached dictionary with {len(lookup)} entries.")
    else:
        lookup = build_opentargets_drug_target_dictionary(ot_map_path)

    raw_jsonl = vault_ct / "non_cns_phase234_trials.jsonl"
    # max_pages_per_axis=80 ensures 100% exhaustive full pagination without truncation (up to 80,000 per axis)
    df, prisma = harvest_non_cns_clinical_trials_macro_atlas(
        raw_jsonl_path=raw_jsonl,
        drug_target_lookup=lookup,
        max_pages_per_axis=80,
    )

    full_csv = vault_ct / "non_cns_phase234_trials_master.csv"
    df.to_csv(full_csv, index=False, encoding="utf-8")

    # Table S1: Per-Axis PRISMA & Stop-Reason Deconvolution Summary
    axis_rows = []
    for axis_id, sub in df.groupby("disease_axis", sort=False):
        vc = sub["stop_deconvolution"].value_counts().to_dict()
        axis_rows.append(
            {
                "disease_axis": axis_id,
                "api_total_reported": prisma["per_axis_counts"].get(axis_id, {}).get("api_total_available", len(sub)),
                "unique_systemic_organ_trials": int(len(sub)),
                "phase2_trials": int((sub["max_phase_num"] == 2).sum()),
                "phase3_trials": int((sub["max_phase_num"] == 3).sum()),
                "phase4_trials": int((sub["max_phase_num"] == 4).sum()),
                "completed_trials": int(vc.get("Completed", 0)),
                "scientific_efficacy_futility": int(vc.get("Scientific_Efficacy_Futility", 0)),
                "scientific_safety_toxicity": int(vc.get("Scientific_Safety_Toxicity", 0)),
                "administrative_commercial_enrollment": int(vc.get("Administrative_Commercial_Enrollment", 0)),
                "other_or_unreported_stop": int(
                    vc.get("Other_Unclassified_Stop", 0) + vc.get("Unreported_Stop_Reason", 0)
                ),
                "target_resolved_trials": int((sub["n_targets"] > 0).sum()),
                "unique_target_genes": int(
                    len({g for s in sub["target_genes"] if s for g in s.split(";") if g})
                ),
            }
        )
    import pandas as pd

    t_s1 = pd.DataFrame(axis_rows)
    t_s1.to_csv(out_tables / "Table_S1_Systemic_Organ_ClinicalTrials_PRISMA_Summary.csv", index=False, encoding="utf-8")

    # Table S2: Target-Resolved Clinical Trials Deconvolution Ledger
    df_resolved = df[df["n_targets"] > 0].copy()
    df_resolved.to_csv(out_tables / "Table_S2_Target_Resolved_Trials_Deconvolution.csv", index=False, encoding="utf-8")

    (vault_ct / "manifest_clinical_macro_atlas.json").write_text(
        json.dumps(prisma, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\n=== STAGE II TASK 1 COMPLETE ===")
    print(f"Total Unique Non-CNS Phase II-IV Drug/Biologic Trials: {len(df):,}")
    print(f"Target-Resolved Trials (ChEMBL/OpenTargets MoA):       {len(df_resolved):,}")
    print(f"CNS/Neuro Firewall Excluded Trials:                    {prisma['cns_neuro_firewall_excluded']:,}")
    print(t_s1.to_string(index=False))


if __name__ == "__main__":
    main()
