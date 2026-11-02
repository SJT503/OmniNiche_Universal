# OmniNiche: Modeling dynamic tissue accessibility predicts clinical trial outcomes across human organs

**Jiangtao Sheng<sup>1,\*</sup>, Xueer Liu<sup>1</sup>, Runmin Lin<sup>1</sup>, Qiaochu Sun<sup>1</sup>, Xiaoxuan Chen<sup>1</sup>, Kangsheng Li<sup>1,\*</sup> & Weiqiang Chen<sup>2,\*</sup>**

<sup>1</sup>Department of Microbiology and Immunology, Shantou University Medical College, 22 Xinling Road, Shantou 515041, Guangdong, China.  
<sup>2</sup>Department of Neurosurgery, The First Affiliated Hospital of Shantou University Medical College, 57 Changping Road, Shantou 515041, Guangdong, China.  
<sup>\*</sup>Correspondence: `jtsheng@stu.edu.cn` (J.S.; ORCID: `0000-0003-2438-7345`); `ksli2013@yeah.net` (K.L.); `wqchen@stu.edu.cn` (W.C.; ORCID: `0000-0002-2612-0890`)

---

## Overview

Human genetics and healthy single-cell atlases guide drug target selection, yet over 85% of clinical drugs fail because diseased tissues remodel cellular target expression and erect physical extracellular matrix (ECM) barriers against macromolecules. **OmniNiche** is a deterministic, multi-scale computational and biophysical modeling framework that bridges germline target validity and diseased-tissue accessibility across human organs:

1. **Macro-Registry Deconvolution (`41,396` Phase II–IV Trials Across `12` Systemic and Parenchymal Organ Disease Axes)**:
   - Deterministic two-stage clause-aware stoppage grammar (`omniniche/stage2/clinical_trials_macro_atlas.py`) validated against `n = 240` clinically adjudicated trial records (`95.83%` accuracy, Cohen's $\kappa = 0.9413$).
   - Stripping operational and commercial closures (`51.4%` of all trial stoppages) exposes a `1.57`-fold excess in biological efficacy futility for ~150-kDa monoclonal antibodies relative to small molecules ($P = 2.26 \times 10^{-4}$).
2. **Dynamic Niche Accessibility Law (`DNAcc`)**:
   - Parsimonious three-parameter white-box formulation coupling single-cell target specificity drift ($\Delta\tau = \tau_{\mathrm{late}} - \tau_{\mathrm{early}}$) with drug-modality steric permeability ($B_{\mathrm{acc}}$) across `29,791` target-resolved clinical trials and `1,155` human target genes.
   - Outperforms static healthy specificity (`AUROC = 0.6710` versus `0.5387`, $P = 1.11 \times 10^{-21}$), generalizes zero-shot to unseen target genes (`Group5Fold AUROC = 0.6541`) and held-out organ axes (`LODO AUROC = 0.6306`, $P = 5.19 \times 10^{-16}$), and combines additively with Open Targets human genetics (`AUROC = 0.7023` scientific; `0.7210` Phase II-to-III/IV maturation).
3. **Zero-Label Multicellular Stochastic Itô–Langevin Waddington Engine**:
   - Compiled strictly from primary human single-cell and 10x Visium spatial biopsies (`20` cohorts across `17` NCBI GEO accessions; `>530,000` single cells and `20,880` quality-filtered Visium spots) with **zero clinical trial outcome supervision**.
   - Ranks `24` SHA-256-locked Phase IIb/III antifibrotic trials across liver, lung, and kidney (`AUROC = 0.8714`, exact permutation $P = 0.0010$; Core-14 `AUROC = 0.9583`), mirrors wet-lab pharmacological reversal in human fibroblasts (`n = 73` across pirfenidone, RepSox, and nintedanib), stratifies `107` strictly pre-treatment anti-PD-1 tumor biopsies (`AUROC = 0.6367`, $P = 0.0187$; single-cell BCC `AUROC = 0.8750`), and isolates a `3.60`-fold immunotherapy-specific stromal survival hazard across `338` melanoma patients ($P = 0.0059$).

---

## Repository Structure

```text
OmniNiche_Universal/
├── omniniche/                         # Core OmniNiche Python package
│   ├── stage2/                        # Multi-scale computational & biophysical engines
│   │   ├── clinical_trials_macro_atlas.py   # Clause-aware stoppage parser & 41,396-trial deconvolution
│   │   ├── pan_disease_drift_law.py         # Single-cell specificity drift & white-box DNAcc law
│   │   ├── fibrosis_mash_deconvolution.py   # 64-donor closure-free coupling & 24-trial Langevin arena
│   │   ├── solid_tumor_spatial_icb.py       # 12-section 10x Visium barrier & 107-patient pre-ICB engine
│   │   ├── r3_orthogonal_validations.py     # LODO, ChEMBL potency, 24k toroidal shifts, Cox E-value, wet-lab & non-ICB controls
│   │   └── run_stage2_all.py                # End-to-end deterministic pipeline runner
│   └── viz/
│       └── stage2_flagship_figures.py       # 100% pure-code Matplotlib vector figure renderers (Figs. 1–5 & Extended Data Figs. 1–6)
├── outputs/
│   ├── stage2_tables/                 # 65 deterministic empirical CSV/JSON/MD archive files
│   └── stage2_figures/                # High-resolution 600-DPI PNG, vector PDF, and SVG figures
├── manuscript/                        # Full manuscript, Supplementary Information, and 13-sheet Excel workbook
├── scripts/                           # Data ingestion, GEO streaming, and audit utilities
├── tests/                             # 55-test automated pytest verification suite
├── DATASETS_REGISTRY.json             # 1,038 screened records & 20 included cohorts across 17 GEO accessions
├── PROJECT_CONSTITUTION.md            # Governance & data-integrity rules
└── STORY_FRAMEWORK.md                 # Multi-scale biophysical formulation & narrative architecture
```

---

## System Requirements & Installation

- **Operating System**: Windows, macOS, or Linux (64-bit; peak RAM `< 7 GB`)
- **Python**: `>= 3.11` (tested on `Python 3.14.4`)
- **Core Dependencies**:
  - `numpy >= 2.0` (`2.4.4`)
  - `scipy >= 1.13` (`1.16.3`)
  - `pandas >= 2.2` (`2.3.3`)
  - `scikit-learn >= 1.5` (`1.8.0`)
  - `matplotlib >= 3.9` (`3.10.9`)
  - `pytest >= 8.0`

```bash
git clone https://github.com/SJT503/OmniNiche_Universal.git
cd OmniNiche_Universal
pip install -e .
```

---

## Reproducing All Tables, Figures, and Verification Gates

### 1. Run the 55-Test Automated Verification Suite
```bash
python -m pytest tests/ -v
```

### 2. Re-execute the End-to-End Empirical & Biophysical Pipeline
```bash
python -m omniniche.stage2.run_stage2_all
python -m omniniche.stage2.r3_orthogonal_validations
```

### 3. Re-render All Main Figures (`Figs. 1–5`) and Extended Data Figures (`Extended Data Figs. 1–6`)
```bash
python -c "from omniniche.viz.stage2_flagship_figures import build_all_stage2_figures, build_supplementary_figures; build_all_stage2_figures(); build_supplementary_figures()"
```

---

## Public Data Accessions

All analyzed clinical trial and human transcriptomic datasets are publicly available without restriction:
- **Clinical Trials**: ClinicalTrials.gov API v2 (`https://clinicaltrials.gov/api/v2/studies`; `41,396` unique Phase II–IV studies across 12 organ axes)
- **Drug Mechanisms & Potency**: ChEMBL v33 (`https://www.ebi.ac.uk/chembl/`)
- **Human Genetics**: Open Targets Platform v24.09 (`https://platform.opentargets.org/`, CC0 1.0)
- **Healthy Single-Cell Reference**: Human Protein Atlas v23 (`https://www.proteinatlas.org/`)
- **NCBI GEO Transcriptomic Cohorts (`17` Accessions Across `20` Datasets)**:
  - *Hepatic, Pulmonary & Renal Fibrosis Single-Cell Atlases*: `GSE136103`, `GSE212837`, `GSE135893`, `GSE271188`
  - *10x Visium Spatial Transcriptomics (`12` Carcinoma Sections, `20,880` QC Spots)*: `GSE171351`, `GSE203612`
  - *Strictly Pre-Treatment Anti-PD-1 Cohorts (`n = 107` Biopsies; `n = 74` OS)*: `GSE123813`, `GSE91061`, `GSE78220`, `GSE207422`
  - *Autoimmune IBD & Intravascular Sepsis Boundary*: `GSE250487`, `GSE167363`
  - *Human Wet-Lab Fibroblast / TGFBR1i Perturbation Triad (`n = 73`)*: `GSE226249`, `GSE154646`, `GSE210020`
  - *Non-ICB Metastatic Melanoma Survival Specificity Controls (`n = 264`)*: `GSE65904`, `GSE22153`

---

## License

This project is released under the MIT License for academic and research use.
