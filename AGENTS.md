# OmniNiche Universal — Workspace Rules (Always Active)

## 1. Strict Project Isolation & Writable Whitelist
- **TARGET / WRITABLE ROOT**: `F:/OmniNiche_Universal` is the **ONLY** writable project directory.
- **FROZEN / DO NOT TOUCH**: `D:/AI工具/NeuroImmune_Virtual_Biotech` and `C:/Users/sheng/Desktop/NeuroImmune_Virtual_Biotech_投稿包_v1.0` belong to an already-submitted *Nature Biomedical Engineering* manuscript (`nBME-26-6122`). **NEVER** read from, copy from, modify, or reference that project.
- **Silent CNS Exclusion Filter, Zero External `"Non-CNS"` Jargon**:
  - Keep the automated neurological/CNS lexical exclusion filter active inside Python data-ingestion code (`0%` CNS trial/cohort/target overlap).
  - **NEVER** use internal partitioning labels (`"Non-CNS"`, `"non-CNS"`, `"non-neurological"`, `"Firewall-2"`) in manuscript Title, Abstract, Introduction, Results, Discussion, Tables, Figure titles/axes, or Figure 1 schematic.
  - Frame the scientific scope positively as **"12 systemic and parenchymal organ disease axes"**, explaining in *Methods* (1 sentence) that the cohort focuses on peripheral systemic and solid-organ pathologies governed by interstitial extracellular matrix (ECM) remodeling rather than specialized blood-brain barrier neurovascular transport.

## 2. Mandatory Startup Sequence in New Window
1. Read `F:/OmniNiche_Universal/PROJECT_CONSTITUTION.md` and `F:/OmniNiche_Universal/STORY_FRAMEWORK.md` in full before taking any action.
2. Run `python -m pytest "F:/OmniNiche_Universal/tests" -v` to confirm the `27 passed` baseline.
3. **Execute Phase 0 Deep *Nature Medicine* / *Science* Benchmark Teardown FIRST**: Before modifying calculation scripts, figures, or prose, study 6–8 real 2024–2026 *Nature Medicine* and *Science* exemplar articles in translational AI, clinical-trial genetics/attrition, single-cell/spatial tissue remodeling, and immunotherapy response prediction, documenting concrete benchmarks for statistical rigor, clinical utility metrics, panel density (`6–8` panels/figure), and narrative structure in `manuscript/Phase0_NatMed_Science_Deep_Benchmark_Teardown.md`.

## 3. Hardware, Data Integrity & Quality Gates
- **RAM Hard Ceiling**: Peak Python process memory MUST stay `< 7 GB` at all times (use streaming / chunked `scipy.sparse` or Donor-level pseudobulk; never load dense multi-cohort cell matrices simultaneously).
- **100% Real Empirical Data & Live-Verified References**: Zero synthetic/mock numbers; every quantitative claim in prose/tables/figures must trace directly to deterministic CSV/NPZ outputs in `F:/OmniNiche_Universal/outputs/stage2_tables/`; every citation must be verified via NCBI PubMed E-utilities (`PMID`, `DOI`, `Journal`, `Year`, `Volume/Pages`).
- **Mandatory Visual & Automated Gates Before Completion**:
  - Visually inspect **100% of rendered PNG figures** via `read_file` to catch any text overlap, clipping, low-contrast heatmap, or visual artifact.
  - Pass `python "C:/Users/sheng/.gemini/skills/figure-workflow-discipline/scripts/figure_clean_aesthetic_gate.py" --dir "F:/OmniNiche_Universal/outputs/stage2_figures"`.
  - Pass `python "C:/Users/sheng/.gemini/skills/manuscript-three-lock-enforcer/scripts/manuscript_convention_and_parity_gate.py" --en "F:/OmniNiche_Universal/manuscript/OmniNiche_Flagship_Manuscript_EN.md" --zh "F:/OmniNiche_Universal/manuscript/OmniNiche_Flagship_Manuscript_ZH_Preview.md"`.
  - Run a full 5-reviewer simulated *Nature Medicine* review (`biomed-review-panel`) after v2.0 assembly and close all identified major concerns before final delivery.
