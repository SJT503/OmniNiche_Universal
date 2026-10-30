from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from omniniche.core.schema import (
    CompartmentRule,
    DirectedEdgeRule,
    NicheGrammar,
    WorkspaceMode,
)
from omniniche.disease_packs.autopack_agent import AutoPackAgent
from omniniche.disease_packs.loader import (
    list_builtin_disease_packs,
    load_disease_pack,
    transfer_grammar_cross_organ,
)
from omniniche.audit.clinical_trialist_v2 import ClinicalTrialistV2
from omniniche.benchmark.calibrated_arena import CalibratedBenchmarkArena
from omniniche.viz.nature_methods_plotter import NatureMethodsPlotter


def _default_fibrosis_demo_grammar() -> NicheGrammar:
    return NicheGrammar(
        grammar_id="fibrosis_ckm_core_v2",
        source_pack_id="fibrosis_ckm_ipf",
        organ_context="Kidney_Cortex",
        compartments=[
            CompartmentRule(name="Injured_Epithelium", baseline_production=0.0015, decay_rate=0.055, noise_sigma=0.018, marker_genes=["HAVCR1", "VCAM1"]),
            CompartmentRule(name="Profibrotic_Macrophage", baseline_production=0.0015, decay_rate=0.050, noise_sigma=0.018, marker_genes=["SPP1", "TREM2"]),
            CompartmentRule(name="Myofibroblast_ECM", baseline_production=0.0010, decay_rate=0.042, noise_sigma=0.016, marker_genes=["POSTN", "COL1A1", "ACTA2"]),
        ],
        edges=[
            DirectedEdgeRule(source="Injured_Epithelium", target="Profibrotic_Macrophage", weight_v=0.115, half_max_k=0.38, hill_n=3.0),
            DirectedEdgeRule(source="Profibrotic_Macrophage", target="Myofibroblast_ECM", weight_v=0.115, half_max_k=0.40, hill_n=3.0),
            DirectedEdgeRule(source="Myofibroblast_ECM", target="Injured_Epithelium", weight_v=0.125, half_max_k=0.40, hill_n=3.0),
            DirectedEdgeRule(source="Myofibroblast_ECM", target="Myofibroblast_ECM", weight_v=0.088, half_max_k=0.52, hill_n=4.0),
        ],
        remission_compartment="Myofibroblast_ECM",
        remission_threshold=0.30,
    )


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="omniniche",
        description="OmniNiche v2.0 Universal Multi-Agent & Stochastic Multicellular Niche Platform",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list-packs", help="List builtin orthogonal non-CNS DiseasePacks")

    p_onboard = sub.add_parser("onboard", help="Auto-generate a validated DiseasePack YAML for any new disease")
    p_onboard.add_argument("--disease", required=True, help="Human disease name")
    p_onboard.add_argument("--organ", default="Systemic_Tissue", help="Target organ/tissue system")
    p_onboard.add_argument("--mode", choices=["lab_private", "public_benchmark"], default="lab_private")
    p_onboard.add_argument("--output", required=True, help="Output YAML path")
    p_onboard.add_argument("--offline", action="store_true", help="Skip live EMBL-EBI OLS ontology lookup")

    p_live = sub.add_parser("live-audit", help="Run live ClinicalTrials.gov API v2 PRISMA audit for a DiseasePack")
    p_live.add_argument("--pack", required=True, help="Builtin pack ID or YAML path")
    p_live.add_argument("--max-pages", type=int, default=2, help="Max API pages per condition term (0 = uncapped)")
    p_live.add_argument("--output", required=True, help="Output CSV path")

    p_demo = sub.add_parser("demo-benchmark", help="Run 5-baseline calibrated benchmark, cross-organ transfer & dashboard")
    p_demo.add_argument("--output-dir", required=True, help="Directory for output tables and vector figures")

    args = parser.parse_args(argv)

    if args.command == "list-packs":
        for pid in list_builtin_disease_packs():
            pk = load_disease_pack(pid)
            print(f"{pk.pack_id:24s} | {pk.archetype:46s} | {pk.display_name}")
        return 0

    if args.command == "onboard":
        agent = AutoPackAgent()
        pack = agent.generate_disease_pack(
            disease_name=args.disease,
            organ_hint=args.organ,
            mode=WorkspaceMode(args.mode),
            output_yaml_path=Path(args.output),
            use_network=not args.offline,
        )
        print(f"Generated DiseasePack '{pack.pack_id}' -> {args.output}")
        return 0

    if args.command == "live-audit":
        pack = load_disease_pack(args.pack)
        trialist = ClinicalTrialistV2(pack)
        max_p = None if args.max_pages <= 0 else args.max_pages
        raw_studies = trialist.fetch_live_studies(max_pages_per_term=max_p)
        df, prisma = trialist.audit_raw_studies(raw_studies)
        out_p = Path(args.output)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_p, index=False, encoding="utf-8-sig")
        print(f"PRISMA summary for {pack.pack_id}: {prisma}")
        print(f"Saved {len(df)} verbatim-verified trials -> {out_p}")
        return 0

    if args.command == "demo-benchmark":
        import pandas as pd

        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        pack = load_disease_pack("fibrosis_ckm_ipf")
        grammar = _default_fibrosis_demo_grammar()
        ipf_grammar = transfer_grammar_cross_organ(grammar, pack.cross_organ_transfer)

        arena = CalibratedBenchmarkArena()
        transfer_res = arena.evaluate_cross_organ_fibrosis_transfer(ipf_grammar, pack)
        pd.DataFrame([transfer_res]).to_csv(out_dir / "cross_organ_kidney_to_lung_ipf_transfer.csv", index=False)

        # Build calibrated arena table and render dashboard
        from tests.test_calibrated_arena_cli_and_live_api import _build_realistic_pan_disease_intervention_benchmark

        bench_df = _build_realistic_pan_disease_intervention_benchmark()
        summary_df = arena.evaluate_five_baselines(bench_df)
        summary_df.to_csv(out_dir / "five_baseline_calibrated_arena.csv", index=False)

        plotter = NatureMethodsPlotter()
        paths = plotter.render_four_panel_dashboard(
            arena_summary_df=summary_df,
            grammar=grammar,
            pack=pack,
            output_prefix=out_dir / "Fig_OmniNiche_v2_Master_Benchmark",
        )
        print("Exported OmniNiche v2.0 dashboard:", paths)
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
