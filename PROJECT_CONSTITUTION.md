# PROJECT_CONSTITUTION — OmniNiche Universal (`F:/OmniNiche_Universal`)

- **Project Name**: `OmniNiche v2.0` (*Dynamic Niche Accessibility & Spatial-Stromal Hysteresis Across Systemic and Parenchymal Organ Diseases*)
- **Canonical Root**: `F:/OmniNiche_Universal`
- **Primary Target Journal**: ***Nature Medicine*** (benchmarked simultaneously against *Nature Medicine* clinical-translational standards and *Science* mechanistic/dynamical rigor)
- **Stage II v1 Git Baseline**: `commit 52d08ec` (`27/27` unit tests passing; Stage II v1 baseline outputs in `outputs/stage2_tables/` and `outputs/stage2_figures/`)

---

## 1. Writable Whitelist & Cross-Project Firewall (`cross-project-handoff-checker` Step 6)

| Role | Path | Permission |
| :--- | :--- | :--- |
| **TARGET / WRITABLE** | `F:/OmniNiche_Universal/` | **Sole writable directory** for all code, data tables, figures, manuscripts, and git commits |
| **FROZEN / FORBIDDEN** | `D:/AI工具/NeuroImmune_Virtual_Biotech/` & `C:/Users/sheng/Desktop/NeuroImmune_Virtual_Biotech_投稿包_v1.0/` | **STRICTLY DO NOT TOUCH** (submitted *Nat. Biomed. Eng.* manuscript `nBME-26-6122`; `0%` data/code/prose touch) |

---

## 2. Verified Real-Data Inventory (`F:/OmniNiche_Universal/data_vault/`)

1. **Macro-Clinical Trial & Target Knowledge Base (`data_vault/raw_clinical_trials/`)**:
   - ClinicalTrials.gov API v2 SHA-256-locked interventional Phase II–IV census across 12 systemic and parenchymal organ disease axes: **`41,396` unique trials** (`23,536` Phase II, `11,936` Phase III, `5,924` Phase IV; `587` records with secondary neurological terms excluded silently at ingestion).
   - Open Targets Platform (`24.09` Parquet schemas: `molecule`, `mechanismOfAction`, `targets`) + ChEMBL target mapping: **`29,791` target-resolved trials (`72.0%`)** spanning **`1,155` unique human target genes** (`586` primary target genes).
   - Multi-agent stoppage audit across **`6,595` prematurely stopped trials**: `2,706` (`41.0%`) Administrative/Commercial/Enrollment, `3,118` (`47.3%`) Unreported/Other, `466` (`7.1%`) Safety/Toxicity, `305` (`4.6%`) Efficacy/Futility; **`23,967` definitive non-administrative trials** (`11,072` Completed Phase III/IV, `12,124` Completed Phase II, `771` explicit biological failures).
2. **Healthy Human Single-Cell Reference**:
   - Human Protein Atlas (HPA) consensus single-cell transcriptomics: **`20,151` human protein-coding genes × `154` healthy cell types** ($\tau_0$ and vital-organ bystander expression share).
3. **Six Primary Human Single-Cell Tissue-Remodeling Axes (`data_vault/raw_single_cell/` & `processed_pseudobulk/`)**:
   - **Hepatic MASH & Cirrhosis**: `GSE136103` (`10` human donors) + `GSE212837` (`7` human MASH donors + blood monocyte trajectories; `17` liver tissue donors / `19` myeloid profiles).
   - **Pulmonary IPF / ILD**: `GSE135893` (`30` human donors, `114,396` cells).
   - **Renal CKD / DKD**: `GSE271188` (`17` human kidney biopsies).
   - **Solid-Tumor NSCLC**: `GSE207422` (`~92,000` TME single cells across `14` patients + `24` strictly pre-treatment diagnostic bulk RNA-seq biopsies).
   - **Mucosal Autoimmunity (Ulcerative Colitis)**: `GSE250487` (`93,900` human colon biopsy cells).
   - **Systemic Bacterial Sepsis**: `GSE167363` (`12` longitudinal human PBMC single-cell samples).
4. **Solid-Tumor Spatial 10x Visium & Pre-Treatment Anti-PD-1 Clinical Cohorts**:
   - **12 Primary Human Carcinoma 10x Visium Sections (`21,796` spatial spots)**: `GSE171351` (`4` human bladder carcinoma sections) + `GSE203612` (`8` human BRCA, OVCA, LIHC, and GIST sections).
   - **Blind Strictly Pre-Treatment Anti-PD-1 Biopsy Cohorts**: `GSE123813` (`n = 8` strictly pre-ICB full-TME BCC scRNA-seq biopsies + `n = 15` paired pre/post clonal trajectories) and `GSE207422` (`n = 14` NSCLC scRNA-seq + `n = 24` strictly pre-treatment baseline NSCLC biopsies; **`n = 32` pooled strictly pre-treatment biopsies** with `0%` post-surgical leakage: `13` Responders, `19` Non-Responders).

---

## 3. Canonical Stage II v1 Output Tables (`F:/OmniNiche_Universal/outputs/stage2_tables/`)

- `Table_S1_NonCNS_ClinicalTrials_PRISMA_Summary.csv` *(rename external display headers in v2.0 to remove `"NonCNS"`, keeping backwards-compatible alias or updating pipeline cleanly)*
- `Table_S2_Target_Resolved_Trials_Deconvolution.csv` (`29,791` rows)
- `Table_S3_TargetGene_6Axis_Drift_Atlas.csv`, `Table_S3a..S3d`
- `Table_S4_Fibrosis_Phase3_Benchmark_Trials.csv` (`14` landmark Phase IIb/III antifibrotic trials), `Table_S4a..S4e`, `Fig4_*.csv/.npz`
- `Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv`, `Table_S5a..S5f`, `Fig5_*.csv/.npz`
- `verified_pubmed_references.json` (`30` live-verified PubMed citations; expand and live-verify any newly added *Nature Medicine* / *Science* citations in v2.0)

---

## 4. Mandatory Entry Verification Commands

```powershell
python -m pytest "F:\OmniNiche_Universal\tests" -v
```
All `27` baseline tests (plus new v2.0 statistical rigor & ablation unit tests) must be `100% GREEN` before and after Stage II v2.0 upgrades.
