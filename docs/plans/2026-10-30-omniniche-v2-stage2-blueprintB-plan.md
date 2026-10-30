# OmniNiche v2.0 Stage II (Blueprint B + Solid Tumor Anti-PD-1 Primary Battlefield) Implementation Plan

> **For agentic workers:** Execute this plan task-by-task using TDD, strict `lab-constitution` authenticity/census gates, and `< 8 GB` peak RAM streaming on `F:\OmniNiche_Universal\`. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute the full real-data empirical validation and 5-figure flagship discovery pipeline of `OmniNiche v2.0` under Blueprint B (*Science* `Virtual Biotech` / *Nature Medicine* 5-figure translational paradigm: `35,000+` non-CNS interventional clinical trials macro-atlas + 6-axis human single-cell dynamic drift law + Pan-organ MASH/Fibrosis Phase III failure deconvolution & zero-shot cross-organ transfer + Solid Tumor Anti-PD-1 paired patient single-cell & 10x Visium spatial barrier blind-test).

**Architecture:** All raw clinical trial payloads, OpenTargets drug-target mappings, and GEO/HPA single-cell & Visium matrices reside in `F:\OmniNiche_Universal\data_vault\` with verbatim SHA-256 provenance manifests (`Firewall-1` & `Firewall-2`) and zero overlap with CNS/neurological disorders (`Firewall-3` / `nBME-26-6122` isolation). Large single-cell matrices are streamed one cohort at a time via sparse CSR / `anndata` backed mode into donor-level pseudobulk and compartment-stage state trajectories (`< 7 GB` peak RAM on a `31.6 GB` machine), feeding `OmniNiche`'s `RidgeGRNPerturber`, `DeconfoundedAutoGrammarCompiler`, `LangevinNicheSolver`, and `5-Baseline Calibrated Arena`.

**Tech Stack:** Python 3.12, `numpy 2.4`, `scipy 1.16`, `pandas 2.3`, `h5py 3.16`, `anndata 0.12`, `scanpy 1.12`, `scikit-learn 1.8`, `matplotlib 3.10`, `uv 0.11.8` (OpenTargets Platform GraphQL skill wrapper), `pytest`.

## Global Constraints

- **Workspace & Storage:** `F:\OmniNiche_Universal\` (`388.5 GB` free); total Stage II raw + processed footprint `< 12 GB`; peak RAM `< 8 GB`.
- **Zero CNS Overlap (`Firewall-2`):** Strictly exclude all CNS/neurological/psychiatric trials and datasets; never touch `D:\AI工具\NeuroImmune_Virtual_Biotech\flagship_paper\submission_package_v1.0\`.
- **100% Data Authenticity (`lab-constitution` Rules 1, 2, 13):** Zero synthetic, simulated, or interpolated data points in any figure, table, or statistic; every `NCT` and `GSE` verified live against official APIs with SHA-256 manifests.
- **Clean Figure Gate (`lab-constitution` Rules 8, 9):** Zero explanatory prose boxes inside data panels, `frameon=False` legends, `figure_clean_aesthetic_gate.py` PASS + `read_file` visual inspection on every rendered PNG.

---

### Task 1: Macro-Atlas Non-CNS ClinicalTrials.gov v2 Census & OpenTargets Drug-Target Deconvolution (`Fig. 2`)

**Files:**
- Create: `F:\OmniNiche_Universal\.licenses\opentargets_database_LICENSE.txt`
- Create: `F:\OmniNiche_Universal\omniniche\stage2\macro_clinical_census.py`
- Create: `F:\OmniNiche_Universal\scripts\run_stage2_task1_clinical_atlas.py`
- Test: `F:\OmniNiche_Universal\tests\test_stage2_clinical_atlas.py`

**Interfaces:**
- Produces:
  - `F:\OmniNiche_Universal\data_vault\raw_clinical_trials\non_cns_phase234_trials.jsonl` (verbatim API records + SHA-256)
  - `F:\OmniNiche_Universal\data_vault\raw_clinical_trials\opentargets_drug_target_map.json`
  - `F:\OmniNiche_Universal\outputs\stage2_tables\Table_S1_NonCNS_ClinicalTrials_PRISMA_Summary.csv`
  - `F:\OmniNiche_Universal\outputs\stage2_tables\Table_S2_Target_Resolved_Trials_Deconvolution.csv`

- [ ] **Step 1: Write failing unit tests in `tests/test_stage2_clinical_atlas.py`** verifying CNS exclusion, administrative vs scientific failure classification (`Scientific_Efficacy_Futility`, `Scientific_Safety_Toxicity`, `Administrative_Commercial_Enrollment`), and drug-to-target gene mapping normalization.
- [ ] **Step 2: Run pytest to verify RED state.**
- [ ] **Step 3: Create `.licenses/opentargets_database_LICENSE.txt` and query OpenTargets drug-target associations across the 12 non-CNS EFO/MONDO disease roots via official `query_opentargets.py` (`uv run`).**
- [ ] **Step 4: Implement `omniniche/stage2/macro_clinical_census.py` & execute full paginated ClinicalTrials.gov API v2 harvest across 12 non-CNS disease axes (`>= 30,000` unique Phase II–IV drug/biologic trials).**
- [ ] **Step 5: Run pytest to verify GREEN state and commit Task 1.**

---

### Task 2: Pan-Organ Normal Baseline `τ₀` + 6-Axis Human Single-Cell Dynamic Drift `Δτ(t) × Barrier` Statistical Law (`Fig. 3`)

**Files:**
- Create: `F:\OmniNiche_Universal\omniniche\stage2\streaming_sc_loader.py`
- Create: `F:\OmniNiche_Universal\omniniche\stage2\pan_disease_drift_law.py`
- Create: `F:\OmniNiche_Universal\scripts\download_stage2_geo_cohorts.py`
- Create: `F:\OmniNiche_Universal\scripts\run_stage2_task2_drift_law.py`
- Test: `F:\OmniNiche_Universal\tests\test_stage2_drift_law.py`

**Verified Cohorts (`100%` Primary Human Non-CNS):**
1. **Normal Pan-Organ Baseline (`τ₀`)**: Human Protein Atlas (`20,151` genes × `154` single-cell types, SHA-256 `2ddd3d90fb050a8e...`).
2. **Axis 1 — Hepatic Fibrosis / MASH**: `GSE136103` (Ramachandran et al., *Nature* 2019 human healthy vs cirrhotic liver `Homo sapiens` GSMs) + `GSE269926` (human MASLD/MASH liver snRNA-seq).
3. **Axis 2 — Pulmonary Fibrosis (IPF)**: `GSE135893` (Habermann et al., *Sci. Adv.* 2020 human IPF vs control lung scRNA-seq).
4. **Axis 3 — Renal CKD / Glomerular Fibrosis**: `GSE271188` (human glomerular CKD vs control kidney snRNA-seq).
5. **Axis 4 — Solid Tumor ICB**: `GSE207422` (Hu et al., *Cell Res.* 2023 human NSCLC anti-PD-1) + `GSE123813` (Yost et al., *Nat. Med.* 2019 human BCC/SCC pre/post anti-PD-1).
6. **Axis 5 — Autoimmune Mucosal IBD**: `GSE250487` (human ulcerative colitis mucosal biopsy `.h5ad.gz`).
7. **Axis 6 — Systemic Bacterial Sepsis**: `GSE167363` (human healthy vs surviving vs fatal longitudinal sepsis scRNA-seq).

- [ ] **Step 1: Write failing tests in `tests/test_stage2_drift_law.py`** checking donor-level pseudobulk aggregation, HPA `τ₀` computation, `Δτ(t)` drift directionality, and logistic/AUROC comparison against static `τ₀`.
- [ ] **Step 2: Download & SHA-256-verify all 6 human GEO disease-axis cohorts into `data_vault/raw_single_cell/` via `scripts/download_stage2_geo_cohorts.py`.**
- [ ] **Step 3: Stream each cohort (< 7 GB RAM) into donor-level compartment pseudobulk matrices in `data_vault/processed_pseudobulk/` and compute `τ₀`, `τ(t)`, `Δτ(t)`, `BC`, and `Barrier(m,t)`.**
- [ ] **Step 4: Link target-resolved Phase II–IV clinical trials (`Task 1`) with `τ₀` vs `OmniNiche Dynamic Drift × Barrier` scores, compute calibrated Odds Ratios, DeLong AUROC P-values, and save `Table_S3_PanDisease_Dynamic_Drift_Law.csv`.**
- [ ] **Step 5: Run pytest to verify GREEN state and commit Task 2.**

---

### Task 3: Deep Case #1 — MASH & Pan-Organ Fibrosis Phase III Failure Deconvolution + Cross-Organ Zero-Shot Grammar Transfer (`Fig. 4`)

**Files:**
- Create: `F:\OmniNiche_Universal\omniniche\stage2\fibrosis_mash_deconvolution.py`
- Create: `F:\OmniNiche_Universal\scripts\run_stage2_task3_fibrosis_mash.py`
- Test: `F:\OmniNiche_Universal\tests\test_stage2_fibrosis_mash.py`

- [ ] **Step 1: Pull verbatim ClinicalTrials.gov API v2 records (`Firewall-1` SHA-256) for real Phase IIb/III antifibrotic/MASH benchmark trials** (Winners: `Resmetirom` `NCT03900429`, `Finerenone` `NCT02540993`, `Nintedanib` `NCT01335464`, `Pirfenidone` `NCT01366209`, `Lanifibranor` `NCT03008070`; Scientific Failures: `Simtuzumab` `NCT01672866`/`NCT01672879`, `Selonsertib` `NCT03053050`/`NCT03053063`, `Cenicriviroc` `NCT03028740`, `Elafibranor` `NCT02704403`, `Pamrevlumab` `NCT03955146`, `Ziritaxestat` `NCT03711162`).
- [ ] **Step 2: Compile deconfounded 3-compartment fibrotic lock-in grammar on Human Liver (`GSE136103` + `GSE269926`), compute `RidgeGRNPerturber` `η` from real single-cell matrices, fit dimensionless `LangevinNicheSolver` `U(S)`, verify Profile-Likelihood identifiability, and quantify the `F2 -> F3` monotherapy window cliff `t₅₀`.**
- [ ] **Step 3: Perform Zero-Shot Cross-Organ Grammar Transfer (`Human Liver Fibrosis -> Human Lung IPF [GSE135893]` & `Human Kidney CKD [GSE271188]`) without refitting network topology, evaluate 5-Baseline Calibrated Arena on Phase III Winners vs Failures, and rank synergistic barrier-breaking dual-target regimens (`Table_S4_Fibrosis_Phase3_Benchmark_and_Synergy.csv`).**
- [ ] **Step 4: Run pytest to verify GREEN state and commit Task 3.**

---

### Task 4: Deep Case #2 — Solid Tumor Anti-PD-1 Paired Patient Single-Cell + 10x Visium Spatial Collagen Barrier & Blind Responder Prediction (`Fig. 5`)

**Files:**
- Create: `F:\OmniNiche_Universal\omniniche\stage2\solid_tumor_spatial_icb.py`
- Create: `F:\OmniNiche_Universal\scripts\run_stage2_task4_solid_tumor_icb.py`
- Test: `F:\OmniNiche_Universal\tests\test_stage2_solid_tumor_icb.py`

- [ ] **Step 1: Download & process human solid-tumor 10x Visium spatial transcriptomics (`GSE171351` `47 MB` Visium `.h5ad.gz` + `GSE203612` *Nat. Genet.* human primary carcinoma Visium `.h5` sections) to quantify physical neighborhood contact & `FAP+/COL1A1+/MMP11+` myCAF collagen barrier excluding `CD8+` T cells ($\mathbf{\Omega}_{\text{spatial}}$).**
- [ ] **Step 2: Couple $\mathbf{\Omega}_{\text{spatial}}$ with paired human Anti-PD-1 clinical cohorts `GSE207422` (NSCLC neoadjuvant Anti-PD-1 scRNA-seq + independent Bulk RNA-seq cohort) and `GSE123813` (BCC/SCC paired Pre/Post Anti-PD-1 scRNA-seq) to compile the spatial-deconfounded 4-compartment TME grammar (`Tumor` – `SPP1+ TAM` – `FAP+/MMP11+ myCAF` – `CXCL13+/TOX+ CD8+ T`).**
- [ ] **Step 3: Blind-predict individual patient clinical response (Responder vs Non-Responder AUROC, AUPRC, Brier Score) from baseline TME projection onto `U(S)` across `GSE207422` (scRNA + bulk) and `GSE123813`, benchmarking against `Static PD-L1 (CD274)`, `Static τ₀`, `T-cell Exhaustion/TIDE`, `LIANA+`, and `CellOracle GRN-only`, and derive optimal barrier-unlocking combination therapies (`Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv`).**
- [ ] **Step 4: Run pytest to verify GREEN state and commit Task 4.**

---

### Task 5: Publication-Grade 5 Main Figures + 3 Extended Data Figures + Clean Aesthetic Audit + Flagship Manuscript Package

**Files:**
- Create: `F:\OmniNiche_Universal\omniniche\viz\stage2_flagship_figures.py`
- Create: `F:\OmniNiche_Universal\scripts\generate_all_stage2_figures.py`
- Create: `F:\OmniNiche_Universal\manuscript\OmniNiche_Flagship_Manuscript_EN.md`
- Create: `F:\OmniNiche_Universal\manuscript\OmniNiche_Flagship_Manuscript_ZH_Preview.md`

- [ ] **Step 1: Render `Fig1`–`Fig5` + `Extended_Data_Fig1`–`Extended_Data_Fig3` (`PDF` + `SVG` + `600 DPI PNG`) strictly from verified Stage II CSV/JSON outputs.**
- [ ] **Step 2: Run `figure_clean_aesthetic_gate.py` + `accession_metadata_authenticity_gate.py` + `read_file` visual inspection on every single figure, fixing any overlap or clutter until `0` issues remain.**
- [ ] **Step 3: Write the complete, verified flagship manuscript (`EN` + `ZH` preview) with 100% live-verified five-field references and Three-Lock figure/table numbering check.**
