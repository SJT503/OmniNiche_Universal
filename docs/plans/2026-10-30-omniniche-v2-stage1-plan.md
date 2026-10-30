# OmniNiche v2.0 Stage I Universal Platform Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:test-driven-development and superpowers:verification-before-completion to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and unit-test `OmniNiche v2.0`, a disease-agnostic multi-agent clinical trial/target auditing, single-cell GRN perturbation, and deconfounded stochastic multicellular grammar auto-compilation platform in `D:/AI工具/OmniNiche_Universal/`, plus register the global `omniniche-universal` skill.

**Architecture:** Five-layer decoupled Python package (`omniniche`): (0) Runtime Integrity & Public/Private Isolation Firewall, (1) Five Orthogonal Non-CNS `DiseasePacks` YAML plugins + Cross-Organ Zero-Shot Grammar Transfer + `AutoPack` onboarding generator, (2) Two-tier Macro/Meso ClinicalTrials.gov v2 PRISMA deconvolution + $K$-stage Dynamic $\tau(t)$/$BC$/Pseudobulk + 6-Gate BTL/Elo Tournament, (3) Intra-compartment Ridge GRN Perturbation ($\eta$) + Deconfounded Causal BIC-Pruned `AutoGrammar` Compiler + Dimensionless Langevin SDE Quasi-Potential ($U(\mathbf{S})$) & Profile-Likelihood Identifiability Solver, and (4) 5-Baseline Calibrated Benchmark Arena + Publication Vector Plotter + CLI.

**Tech Stack:** Python 3.10+, `pydantic>=2.0`, `pyyaml`, `numpy`, `scipy`, `pandas`, `scikit-learn`, `matplotlib`, `pytest`, NCBI E-utilities / ClinicalTrials.gov REST API v2 / Open Targets GraphQL.

## Global Constraints

- Never touch or modify `D:/AI工具/NeuroImmune_Virtual_Biotech/` (`nBME-26-6122` frozen archive).
- Enforce `Firewall-1` (verbatim `NCT_ID` SHA-256 lock), `Firewall-2` (disease ontology & matrix SHA-256 provenance lock, blocking formula-interpolated columns), and `Firewall-3` (`public_benchmark` physically blocks reading `pr8_vgx_reanalysis` or private lab paths).
- Every production module MUST follow TDD (`RED -> Verify RED -> GREEN -> Verify GREEN`).

---

### Task 1: Package Scaffolding, Core Pydantic Schemas & Three-Gate Integrity Firewall

**Files:**
- Create: `D:/AI工具/OmniNiche_Universal/pyproject.toml`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/core/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/core/schema.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/core/integrity_firewall.py`
- Test: `D:/AI工具/OmniNiche_Universal/tests/test_core_schema_and_firewall.py`

- [ ] **Step 1: Write failing tests for Pydantic schemas (`DiseasePack`, `NicheGrammar`, `TrialRecord`, `MatrixProvenance`) and `IntegrityFirewall` (catching altered NCT title, synthetic/unhashed matrix metric, and private path access in `public_benchmark` mode).**
- [ ] **Step 2: Run `pytest tests/test_core_schema_and_firewall.py -v` and confirm RED failure.**
- [ ] **Step 3: Implement `schema.py` and `integrity_firewall.py`.**
- [ ] **Step 4: Run `pytest tests/test_core_schema_and_firewall.py -v` and confirm GREEN pass.**

---

### Task 2: Five Orthogonal Non-CNS `DiseasePacks`, Cross-Organ Grammar Transfer Adapter & `AutoPack` Generator

**Files:**
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/sepsis_ards.yaml`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/fibrosis_ckm_ipf.yaml`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/solidtumor_time.yaml`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/autoimmune_ibd_ra.yaml`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/cardiometabolic_mash.yaml`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/loader.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/disease_packs/autopack_agent.py`
- Test: `D:/AI工具/OmniNiche_Universal/tests/test_disease_packs_and_transfer.py`

- [ ] **Step 1: Write failing tests validating all 5 orthogonal non-CNS YAML packs, barrier accessibility curves $\mathcal{B}(m, t)$, cross-organ grammar zero-shot transfer (`kidney_ckm -> lung_ipf`), and `AutoPack` generation.**
- [ ] **Step 2: Run `pytest tests/test_disease_packs_and_transfer.py -v` and confirm RED failure.**
- [ ] **Step 3: Implement the 5 YAML packs, `loader.py` (`load_disease_pack`, `compute_barrier_access`, `transfer_grammar_cross_organ`), and `autopack_agent.py`.**
- [ ] **Step 4: Run `pytest tests/test_disease_packs_and_transfer.py -v` and confirm GREEN pass.**

---

### Task 3: `OmniBiotech-Audit` — Macro/Meso ClinicalTrials.gov v2 PRISMA Deconvolution, $K$-Stage `ChronoTau`, and 6-Gate BTL/Elo Tournament

**Files:**
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/audit/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/audit/clinical_trialist_v2.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/audit/chronotau_tracker.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/audit/sixgate_tournament.py`
- Test: `D:/AI工具/OmniNiche_Universal/tests/test_audit_pipeline.py`

- [ ] **Step 1: Write failing tests for PRISMA 4-step trial filtering & administrative-vs-scientific stop deconvolution, $K$-stage $\tau(t)$ & $\Delta\tau$ drift & Sarle's $BC$ & donor-level pseudobulk (rejecting cell-level pseudoreplication), and 6-gate ablation + BTL/Elo ranking.**
- [ ] **Step 2: Run `pytest tests/test_audit_pipeline.py -v` and confirm RED failure.**
- [ ] **Step 3: Implement `clinical_trialist_v2.py`, `chronotau_tracker.py`, and `sixgate_tournament.py`.**
- [ ] **Step 4: Run `pytest tests/test_audit_pipeline.py -v` and confirm GREEN pass.**

---

### Task 4: Phase 3 Core A — Intra-Compartment Ridge GRN Perturbation & Deconfounded BIC-Pruned `AutoGrammar` Compiler

**Files:**
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/engine/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/engine/grn_perturb.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/engine/autogrammar_compiler.py`
- Test: `D:/AI工具/OmniNiche_Universal/tests/test_grn_and_autogrammar.py`

- [ ] **Step 1: Write failing tests verifying: (a) `RidgeGRNPerturber` fits TF-masked regularized operator $\mathbf{W}^{(c)}$ and predicts bounded KO/drug efficacy $\eta \in [-1, 1]$; (b) `DeconfoundedAutoGrammarCompiler` eliminates spurious common-confounder edges ($Z \to X$, $Z \to Y$ with no $X \to Y$ causal edge!) while recovering true lagged positive-feedback loops, applying spatial proximity weights and BIC pruning, and exporting valid `NicheGrammar.yaml`.**
- [ ] **Step 2: Run `pytest tests/test_grn_and_autogrammar.py -v` and confirm RED failure.**
- [ ] **Step 3: Implement `grn_perturb.py` and `autogrammar_compiler.py`.**
- [ ] **Step 4: Run `pytest tests/test_grn_and_autogrammar.py -v` and confirm GREEN pass.**

---

### Task 5: Phase 3 Core B — Dimensionless Langevin SDE, Waddington Quasi-Potential Landscape $U(\mathbf{S})$, Profile-Likelihood Identifiability & Chrono-Synergy Solver

**Files:**
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/engine/langevin_solver.py`
- Test: `D:/AI工具/OmniNiche_Universal/tests/test_langevin_and_identifiability.py`

- [ ] **Step 1: Write failing tests verifying: (a) dimensionless $[0, 1]$ invariance and ensemble Langevin simulation outputting simulated clinical remission rate %; (b) Waddington quasi-potential $U(\mathbf{S})$ & saddle barrier $\Delta U_{\text{barrier}}$ drop under rescue combinations; (c) monotherapy window cliff $t_{50}$ & Bliss staged synergy matrix; (d) Profile-Likelihood identifiability scan flagging unidentifiable redundant parameters and LOCO-CV calibration.**
- [ ] **Step 2: Run `pytest tests/test_langevin_and_identifiability.py -v` and confirm RED failure.**
- [ ] **Step 3: Implement `langevin_solver.py`.**
- [ ] **Step 4: Run `pytest tests/test_langevin_and_identifiability.py -v` and confirm GREEN pass.**

---

### Task 6: Five-Baseline Calibrated Benchmark Arena, Publication Vector Plotter, CLI & Live End-to-End Verification

**Files:**
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/benchmark/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/benchmark/calibrated_arena.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/viz/__init__.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/viz/nature_methods_plotter.py`
- Create: `D:/AI工具/OmniNiche_Universal/omniniche/cli.py`
- Create: `C:/Users/sheng/.gemini/config/skills/omniniche-universal/SKILL.md`
- Test: `D:/AI工具/OmniNiche_Universal/tests/test_calibrated_arena_cli_and_live_api.py`

- [ ] **Step 1: Write failing tests for the 5-Baseline Calibrated Arena (`B1_Static_Healthy_Tau`, `B2_Uninformative_Linear_Mean`, `B3_SingleCell_GRN_Only`, `B4_Static_LR_Centrality`, `OmniNiche_v2_Full`), cross-organ fibrosis zero-shot transfer evaluation, 4-panel *Nature Methods* vector figure generation (`.pdf`/`.svg`/`.png`), CLI subcommands (`list-packs`, `onboard`, `demo-benchmark`), and live ClinicalTrials.gov API v2 `Firewall-1` smoke test.**
- [ ] **Step 2: Run `pytest tests/test_calibrated_arena_cli_and_live_api.py -v` and confirm RED failure.**
- [ ] **Step 3: Implement `calibrated_arena.py`, `nature_methods_plotter.py`, `cli.py`, and global `SKILL.md`.**
- [ ] **Step 4: Run full `pytest -v` across the entire repository + live CLI verification and confirm 100% GREEN.**
