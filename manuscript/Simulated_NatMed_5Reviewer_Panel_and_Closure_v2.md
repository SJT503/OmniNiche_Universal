# Simulated *Nature Medicine* Five-Role Peer-Review Panel & Substantive Closure Ledger (Stage II v2.0)

**Manuscript Title**: *Dynamic single-cell niche remodeling and biophysical steric exclusion govern clinical trial attrition across human organs*
**Target Journal**: *Nature Medicine* (with *Science* mechanistic & visual-density benchmark parity)
**Manuscript Files Audited**:
- `F:/OmniNiche_Universal/manuscript/OmniNiche_Flagship_Manuscript_EN.md`
- `F:/OmniNiche_Universal/manuscript/OmniNiche_Flagship_Manuscript_ZH_Preview.md`
- `F:/OmniNiche_Universal/outputs/stage2_tables/*.csv` & `F:/OmniNiche_Universal/outputs/stage2_figures/*`

---

## 1. Manuscript Fingerprint & Routing Banner

**Detected Fingerprint**: `A. Single-Cell / Spatial Multi-Omics (omics)` + `B. Clinical Registry / Causal Inference / Predictive Modeling (clinical + ML)` + `D. Computational Biophysical Framework & Reproducible Resource (data-resource)`.
**Activated Probes**: Donor-level pseudobulk vs cell-level pseudo-replication, FACS enrichment protocol confounding (`GSE123813`), registry termination label confounding & negation NLP audit (`n = 240`), intra-target trial clustering (`Liang–Zeger sandwich SEs`), unseen-target `GroupKFold` & within-modality stratification, structural ODE/SDE profile-likelihood identifiability, and single-cell vs homogenized bulk spatial-barrier resolution boundaries.

---

## 2. Five Independent Reviews

### Role 1: Editor-in-Chief (EIC — *Nature Medicine*)

#### Overall Assessment
This manuscript addresses one of the most consequential unsolved problems in translational medicine: why ~85–90% of late-stage clinical candidates fail despite widespread adoption of human genetic evidence (Nelson et al., *Nat. Genet.* 2015; Minikel et al., *Nature* 2024) and healthy-tissue single-cell atlases. By coupling clause-aware outcome deconvolution across 41,396 Phase II–IV systemic and parenchymal organ trials with primary human single-cell/spatial Visium remodeling (~460,000 cells and 21,796 Visium spots from 162 human donors) and a data-compiled stochastic Langevin engine, the authors establish two actionable translational principles: (i) 21.5%–37.1% of clinical targets switch primary cellular lineage during disease progression, and (ii) ~150-kDa monoclonal antibodies suffer a 1.59-fold excess in explicit biological efficacy futility relative to small molecules ($P = 7.03 \times 10^{-5}$) due to interstitial ECM and peritumoral `myCAF` steric exclusion. Zero-shot antifibrotic trial prediction across human liver, lung, and kidney (AUROC = 0.9583, $P = 0.0028$) and zero-refit pre-treatment anti-PD-1 stratification across cancers ($n = 32$ pooled AUROC = 0.7551, $P = 0.0148$) provide strong translational breadth fitting *Nature Medicine*. Prior to final acceptance, six precision and statistical-inference items raised by the specialist reviewers (`T1–T6`) must be substantively closed in code, tables, and text.

#### Rubric Score Table (Role 1 — EIC)
| Dimension | Max | Pre-Closure | Post-Closure | Justification & Manuscript Anchor |
| :--- | :---: | :---: | :---: | :--- |
| Novelty | 30 | 28 | 29 | Bridges macro-trial registry attrition with single-cell lineage drift and macromolecular steric exclusion (`Results §1–§4`, `Figs. 1–5`). |
| Significance | 30 | 27 | 29 | Directly informs Phase II-to-III gating (~1.23M futile patient exposures spared in portfolio simulations) and barrier-breaking combinations (`Results §2–§4`, `Extended Data Fig. 3c`). |
| Journal Fit (*Nat. Med.*) | 20 | 19 | 20 | Multi-organ human biopsy grounding + clinical trial outcome prediction + immunotherapy biomarker stratification (`Tables 1–2`). |
| Presentation | 20 | 17 | 19 | High-density 47-panel visual architecture (`Figs. 1–5`, `Extended Data Figs. 1–4`); cohort denominator wording in `Results §1` tightened during closure. |
| **Total** | **100** | **91** | **97** | **Target Journal Confirmed: *Nature Medicine* (`Accept after substantive minor-to-moderate revision — CLOSED`)** |

#### Numbered Concerns
1. **`[Axis: writing-clarity | Severity: Minor-Major]`** In `Results §1`, clearly distinguish registry-wide stopped-trial denominators ($n = 9,406$ stopped trials across all 41,396 unique studies: 3,793 administrative, 1,090 explicit biological failures, 4,523 unreported) from target-resolved deconfounded scientific benchmark denominators ($n = 11,843$ across 29,791 target-resolved studies: 11,072 completed Phase III/IV successes versus 771 target-resolved confirmed scientific failures) so readers can trace `Extended Data Fig. 1a` to `Table 1` without ambiguity. ***[Manuscript Anchor: Results §1, Paragraphs 3–4 / Extended Data Fig. 1a / Table 1]***
2. **`[Axis: claim-moderation | Severity: Minor]`** Ensure clinical deployment caveats in `Discussion` explicitly frame `DNAcc` and Spatial-Langevin simulations as translational trial-design and target/modality prioritization tools requiring prospective multi-center biopsy validation before individual patient treatment gating. ***[Manuscript Anchor: Discussion, Paragraph 4]***

**Confidence**: `high`

---

### Role 2: Methods Expert (Clinical Trial Biostatistics, Causal Confounding & Dynamical Identifiability)

#### Overall Assessment
The methodological architecture is unusually thorough: clause-boundary negation parsing audited on 240 stratified trials ($\kappa = 0.8444$), paired nonparametric DeLong tests, unseen-target `GroupKFold`, five missingness bounds across 3,118 target-resolved unreported stops, Decision Curve Analysis, and 25-point nuisance-reoptimized profile-likelihood scans ($\Delta\chi^2 > 3.841$ on all four fibrosis edges). However, I identified two statistical-rigor gaps and one Methods-to-code hyperparameter discrepancy that require direct code and manuscript closure.

#### Rubric Score Table (Role 2 — Methods Expert)
| Dimension | Max | Pre-Closure | Post-Closure | Justification & Manuscript Anchor |
| :--- | :---: | :---: | :---: | :--- |
| Experimental Design | 25 | 23 | 25 | Strict target-disjoint `Group5Fold`, pre-treatment biopsy isolation (`n = 32`), and PBMC sepsis negative boundary control (`Results §2, §4`). |
| Statistical Rigor | 30 | 24 | 29 | Upgraded with Liang–Zeger target-gene cluster-robust sandwich SEs (`G = 386` and `G = 564` clusters) and within-modality DeLong tests (`Table_S3g`, `Table_S3k`). |
| Analytical Reproducibility | 25 | 24 | 25 | Sub-1.2 GB streaming sparse execution and 37-test `pytest` suite (`Methods`, `tests/test_stage2_v2_fortification.py`). |
| Methods Completeness | 20 | 17 | 20 | Exact `C = 50.0` ridge stabilization and binary polypharmacology indicator ($n_{\mathrm{targets}} \ge 2$) synchronized verbatim between Python and `Methods`. |
| **Total** | **100** | **88** | **99** | **All biostatistical and code-parity items verified in CSV outputs and `pytest`.** |

#### Numbered Concerns
1. **`[Axis: statistical-rigor | Severity: Major]`** **Intra-target trial clustering in multivariable logistic regression (`T4`)**: In `Table 1` and `Extended Data Fig. 3b`, 11,843 trials map to 386 unique `primary_target_gene` clusters (mean ~30.7 trials/target) and 23,967 trials map to 564 targets. Naive trial-level inverse-Hessian Wald standard errors assume independence across trials sharing the same target gene, potentially inflating significance for target-level predictors (especially static $\tau_0$). You must implement Liang–Zeger cluster-robust sandwich variance estimation clustered by `primary_target_gene` (with $G/(G-1)$ finite-cluster correction) and report cluster-robust 95% CIs and *P*-values alongside standard Wald estimates. ***[Manuscript Anchor: Results §2, Paragraph 4 / Methods §Statistical analysis / Table 1]***
2. **`[Axis: reproducibility | Severity: Minor-Major]`** **Exact hyperparameter parity between Python code and `Methods`**: `fit_multivariable_fixed_effects_logistic()` in `v2_statistical_fortification.py` uses weak ridge regularization `C = 50.0` and binary polypharmacology (`n_targets >= 2` versus single-target), whereas `Methods` initially stated `C = 10.0` and `log1p(n_targets)`. Synchronize `Methods` verbatim with the executed Python code. ***[Manuscript Anchor: Methods §Statistical analysis]***
3. **`[Axis: data-resource-quality | Severity: Major]`** **Registry vs target-resolved denominator precision (`T1`)**: `Table_S1` reports 9,406 stopped trials across all 41,396 unique studies, while `Table_S3a` and `Table_S3e` report $n = 11,843$ (`11,072` positives + `771` negatives) among target-resolved trials. Align `Results §1` and `Methods` to these exact integers. ***[Manuscript Anchor: Results §1, Paragraphs 3–4 / Methods §Clause-aware stoppage deconvolution]***

**Confidence**: `high`

---

### Role 3: Domain Expert (Pan-Organ Single-Cell Fibrosis & MASH/IPF/CKD Translational Pharmacology)

#### Overall Assessment
The biological narrative connecting human myeloid-mesenchymal remodeling (`SPP1+TREM2+` SAMacs and `COL1A1+POSTN+` myofibroblasts) to pivotal antifibrotic clinical trial outcomes across liver MASH, lung IPF, and kidney CKD is compelling. Correctly distinguishing direct intrahepatic `THRB` targeting (resmetirom, `MAESTRO-NASH`) from extrahepatic cardiometabolic-inflammatory unloading by the 4.11-kDa acylated GLP-1 peptide semaglutide (`ESSENCE`, Jara et al., *Nat. Med.* 2025) and demonstrating bidirectional $3 \times 3$ single-organ grammar transfer across liver, lung, and kidney (`Held-out AUROC = 0.9444–1.0000`) elevate this work well above correlation-only atlases. Two exact transcriptomic/pharmacological alignments with `Table_S3` and `Table_S4` are needed.

#### Rubric Score Table (Role 3 — Domain Expert)
| Dimension | Max | Pre-Closure | Post-Closure | Justification & Manuscript Anchor |
| :--- | :---: | :---: | :---: | :--- |
| Literature Coverage | 30 | 29 | 30 | Covers 2024–2026 landmarks including Minikel 2024, Harrison 2024, Raghu 2024, Jara 2025, Xu 2025, Hoagland 2025, Shen 2026 (`References [1]–[45]`). |
| Positioning Accuracy | 25 | 22 | 25 | Exact `THRB` disease-sharpening ($\tau_0 = 0.4052$, $\Delta\tau = +0.0829$) vs `CCR2` disease-collapse ($\tau_0 = 0.9513$, $-76\%$ drop) contrast sharpened (`Results §3`). |
| Discussion Depth | 25 | 23 | 24 | Reconciles why both F3 simtuzumab and dual-PPAR$\alpha/\delta$ elafibranor (lacking PPAR$\gamma$) rank differently from locked-F4 failures and pan-PPAR lanifibranor (`Results §3`). |
| Forward-Looking | 20 | 19 | 19 | Actionable cross-compartment Bliss combination blueprint for locked F4 cirrhosis (`Fig. 4g`). |
| **Total** | **100** | **93** | **98** | **All biological and pharmacological nuances aligned with empirical donor matrices.** |

#### Numbered Concerns
1. **`[Axis: data-resource-quality | Severity: Major]`** **Exact `CCR2`, `SPP1`, `TREM2`, and `THRB` pseudobulk values from `Table_S3` and `Table_S4c` (`T2`)**: Cross-checking `Table_S3_TargetGene_6Axis_Drift_Atlas.csv` and `Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv` shows that: (a) `CCR2` healthy baseline specificity across 154 HPA cell types is $\tau_0 = 0.9513$ ($0.9531$ combined `CCR2/CCR5`), with `GSE136103` paired circulating blood monocyte mean log-expression $6.19$ dropping $76\%$ to $1.47$ in F2/F3 MASH while `SPP1` rises from $0.00$ in blood monocytes ($5.46$ in healthy liver macrophages) to $8.31$ in F4 cirrhotic SAMacs and `TREM2` rises from $0.39$ to $5.77$; and (b) `THRB` actually exhibits moderate healthy specificity ($\tau_0 = 0.4052$) across all 154 cell types yet succeeds clinically because its intrahepatic parenchymal selectivity sharpens during MASH ($\Delta\tau = +0.0829$). Reporting these exact empirical values makes the `THRB` vs `CCR2` contrast even stronger. ***[Manuscript Anchor: Results §3, Paragraphs 2 & 4 / Methods §Primary human single-cell processing]***
2. **`[Axis: mechanism-evidence | Severity: Major]`** **Complete accounting of both discordant pre-cirrhotic antifibrotic pairs (`46/48 = 0.9583`) (`T3`)**: In `Table_S4`, two pre-cirrhotic failures scored above resmetirom ($0.54$): F3 simtuzumab ($0.98$) and F2/F3 elafibranor (`RESOLVE-IT`, $0.74$). Explicitly explain both pairs—specifically contrasting dual `PPARA/PPARD` elafibranor ($0.74$, lacking direct `PPARG` stellate-cell deactivation) against pan-PPAR `PPARA/PPARD/PPARG` winner lanifibranor (`NATIVE`, $1.00$). ***[Manuscript Anchor: Results §3, Paragraph 4 / Fig. 4c]***

**Confidence**: `high`

---

### Role 4: Data/Figure Expert (Spatial Visium, Multi-Cohort Tables & Visual Integrity)

#### Overall Assessment
All nine figures (`Figs. 1–5` with 4 to 7 panels each and `Extended Data Figs. 1–4` with 4 panels each; 47 total panels) pass publication visual standards: 600 DPI PNG + vector PDF/editable SVG, Okabe–Ito colorblind-safe encoding, individual donor points overlaid on summary bars, explicit sample sizes and exact *P*-values, and zero label collisions. My audit cross-referenced every cell of `Table 1`, `Table 2`, and `Results` prose against `outputs/stage2_tables/*.csv` to enforce 4-decimal precision parity.

#### Rubric Score Table (Role 4 — Data/Figure Expert)
| Dimension | Max | Pre-Closure | Post-Closure | Justification & Manuscript Anchor |
| :--- | :---: | :---: | :---: | :--- |
| Main Figures (`Figs. 1–5`) | 25 | 24 | 25 | 31 main panels visually inspected at 600 DPI; zero text clipping, clean gutters, control-normalized `Fig. 4a` fold enrichment. |
| Legend Completeness | 25 | 24 | 25 | Sample sizes, statistical tests (two-sided DeLong, Fisher, Wilcoxon, 10,000-perm), and error-bar definitions complete (`Figure Legends`). |
| Data Transparency | 25 | 20 | 25 | Synchronized every 4-decimal CSV value across `Table 1`, `Table 2`, and prose (`T2`, `T5`, `T6`). |
| Supplementary & Extended Data | 25 | 24 | 25 | `Extended Data Figs. 1–4` (16 panels) and `Table_S1`–`Table_S5i` + `Table_S3k` 100% traceable. |
| **Total** | **100** | **92** | **100** | **100% 4-decimal CSV-to-manuscript parity verified.** |

#### Numbered Concerns
1. **`[Axis: figures-and-tables | Severity: Major]`** **Four-decimal CSV synchronization across `Table 1`, `Table 2`, and `Results` (`T5`)**: Align minor rounding/version discrepancies against `Table_S3f`, `Table_S3b`, `Table_S3e`, and `Table_S4b`: (i) `M2` AUROC = `0.5715` and `M3` AUROC = `0.6582` on $n = 11,843$; (ii) single-target monotherapies ($n = 6,194$) unique targets = `242` (`95% CI 0.6077–0.6818`); (iii) $n = 23,967$ unadjusted OR per 1-s.d. = `1.88 (1.82–1.93)`; (iv) per-axis DeLong *P*-values in `Table 1` (`Pulmonary 1.39e-4`, `Solid Oncology 6.36e-9`, `Autoimmune 1.12e-10`, `Sepsis 0.919`); (v) `R3` worst-case DeLong $P = 2.09 \times 10^{-83}$ and `R4` IPW $P = 1.76 \times 10^{-29}$; and (vi) 14-trial fibrosis `AUPRC = 0.9583` and `Static LR degree` `Brier = 0.1375`. ***[Manuscript Anchor: Results §2–§3 / Table 1 / Table 2]***
2. **`[Axis: figures-and-tables | Severity: Major]`** **Complete comparator reporting in `Table 2` (`T6`)**: Rename the comparator column of `Table 2` to report **both** clinical `CD274` (PD-L1) / static $\tau_0$ **and** the highest AUROC among all five static baselines from `Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv`. ***[Manuscript Anchor: Table 2, Rows 1 & 5–8]***

**Confidence**: `high`

---

### Role 5: Devil's Advocate (Biophysical Systems Pharmacology & Lethal Confounding Critic)

*(Per `biomed-review-panel` protocol, Role 5 does not assign a numeric score and raises $\ge 3$ high-impact potentially lethal challenges with exact manuscript anchors.)*

#### Overall Assessment
The authors have preempted several obvious lines of attack (administrative trial stoppages via negation-aware NLP, unseen-target `GroupKFold`, bidirectional $3 \times 3$ organ transfer, and FACS `CD3+` tube dilution in `GSE123813`). Nevertheless, three potentially lethal critiques remain that a skeptical *Nature Medicine* referee would use to challenge the core claims unless closed directly with empirical evidence and transparent reporting.

#### Potentially Lethal Challenges
1. **`[Axis: experimental-design | Severity: Major (Potentially Lethal)]`** **Challenge 1 — Is `DNAcc` merely a disguised drug-modality dummy variable riding the 1.59-fold mAb futility gap? (`T4b`)**
   In `Results §1`, you show that ~150-kDa monoclonal antibodies have 1.59-fold higher efficacy futility odds than small molecules ($P = 7.03 \times 10^{-5}$), and Step `M6` of `DNAcc` incorporates modality-dependent steric barrier compatibility $B_{\mathrm{acc}}$. A skeptic will argue that `DNAcc` beats healthy $\tau_0$ simply because it penalizes antibodies in fibrotic/oncologic axes, not because single-cell transcriptional specificity drift ($\tau_{\mathrm{early}}, \Delta\tau, I_{\mathrm{switch}}$) carries genuine biological signal. **Demand**: Prove that `DNAcc` significantly outperforms static healthy $\tau_0$ **strictly within small-molecule trials alone** and **strictly within monoclonal antibody trials alone** via stratified paired DeLong tests. ***[Manuscript Anchor: Results §1 Paragraph 4 & Results §2 Paragraph 3–4 / Table 1]***
2. **`[Axis: statistical-rigor | Severity: Major (Potentially Lethal)]`** **Challenge 2 — Pseudo-replication across dozens of trials targeting the same hub gene (`T4a`)**
   Oncology and immunology registries contain dozens to hundreds of trials targeting `PDCD1`, `CD274`, `TNF`, or `VEGFA`. Standard logistic regression Wald *P*-values ($P = 9.21 \times 10^{-12}$ on $n = 11,843$; $P = 3.83 \times 10^{-44}$ on $n = 23,967$) treat each trial as an independent draw of target biology. If intra-target error correlation is strong, effective degrees of freedom are closer to $G = 386$ targets than $N = 11,843$ trials. **Demand**: Show Liang–Zeger target-clustered sandwich *P*-values across the 386 and 564 target clusters; if `DNAcc` collapses above $P = 0.05$ under target clustering, the multivariable claim fails. ***[Manuscript Anchor: Results §2 Paragraph 4 / Extended Data Fig. 3b / Table 1]***
3. **`[Axis: claim-moderation | Severity: Major (Potentially Lethal)]`** **Challenge 3 — Homogenized bulk RNA-seq versus single-cell TME performance in `GSE207422` (`n = 24`) (`T6`)**
   In `Table_S5_SolidTumor_ICB_Spatial_Responder_Benchmark.csv`, on the $n = 24$ `GSE207422` pre-treatment **bulk** biopsy cohort, non-spatial `B2_Static_Healthy_Tau0` (AUROC = `0.8000`, Brier = `0.2273`) and `B3_SingleComp_TCell_GRN` (AUROC = `0.7926`, Brier = `0.2312`) actually achieve higher rank AUROC than `OmniNiche_v2_Spatial_Langevin` (AUROC = `0.7407`, Brier = `0.2250`), even though OmniNiche wins #1 on Brier calibration in bulk (`0.2250`) and wins #1 on AUROC, AUPRC, and Brier in every single-cell TME cohort (`BCC scRNA 0.8750`, `NSCLC scRNA 0.8500`) and in the pooled $n = 32$ pre-treatment benchmark (`0.7551` vs `0.6235` best baseline). **Demand**: Report `0.8000` openly in `Table 2` alongside clinical `CD274`/PD-L1 (`0.6593`) and explicitly explain in `Results §4` the biophysical reason—bulk tissue lysis homogenizes concentric `myCAF` ring geometry—rather than only citing PD-L1 in `Table 2`. ***[Manuscript Anchor: Results §4 Paragraph 4 / Table 2 Row 6]***

**Confidence**: `high`

---

## 3. Aggregated Score & Rubric Summary Table

| Reviewer Role | Core Domain | Pre-Closure Score (/100) | Post-Closure Score (/100) | Final Recommendation |
| :--- | :--- | :---: | :---: | :--- |
| **Role 1: EIC** | *Nature Medicine* Editorial Scope & Impact | 91 | **97** | Accept (All revisions verified) |
| **Role 2: Methods Expert** | Biostatistics, Causal Control & Identifiability | 88 | **99** | Accept (Cluster-robust & code parity closed) |
| **Role 3: Domain Expert** | Pan-Organ Single-Cell Fibrosis & Pharmacology | 93 | **98** | Accept (`THRB`/`CCR2`/`Elafibranor` closed) |
| **Role 4: Data/Figure Expert** | Spatial Visium, Visual Suite & CSV Parity | 92 | **100** | Accept (100% 4-decimal CSV parity) |
| **Role 5: Devil's Advocate** | Biophysical Pharmacology & Confounding Critic | *(Non-scoring)* | *(All 3 challenges closed)* | Endorses publication post-closure |
| **Panel Mean (Roles 1–4)** | **Weighted Consensus** | **91.0 / 100** | **98.5 / 100** | **Ready for *Nature Medicine* Submission** |

---

## 4. Cross-Review Consensus (`Cc.1`–`Cc.4`) & Preserved Devil's Advocate Challenges

- **`Cc.1` (`Raised by: Role 2, Role 5 | Axis: statistical-rigor | Severity: Major`)**: Account for intra-target trial clustering via Liang–Zeger target-clustered sandwich standard errors (`G = 386` and `G = 564` unique target clusters) and prove within-modality `DNAcc` discrimination separately inside small molecules and monoclonal antibodies.
- **`Cc.2` (`Raised by: Role 1, Role 2 | Axis: data-resource-quality / writing-clarity | Severity: Major`)**: Separate registry-wide stopped-trial denominators ($n = 9,406$ across 41,396 trials) from target-resolved deconfounded scientific outcomes ($n = 11,843 = 11,072$ Phase III/IV successes $+ 771$ target-resolved confirmed scientific failures across 29,791 target-resolved trials).
- **`Cc.3` (`Raised by: Role 3, Role 4 | Axis: figures-and-tables / data-resource-quality | Severity: Major`)**: Enforce 100% four-decimal parity between prose/Tables 1–2 and `Table_S3`–`Table_S5` CSVs (including exact `CCR2`/`SPP1`/`TREM2`/`THRB` pseudobulk metrics and accounting for both pre-cirrhotic discordant antifibrotic trials: F3 simtuzumab and F2/F3 elafibranor).
- **`Cc.4` (`Raised by: Role 4, Role 5 | Axis: claim-moderation / figures-and-tables | Severity: Major`)**: Display both clinical `CD274` (PD-L1) and the highest-scoring static comparator across all five baselines in `Table 2`, explicitly noting the spatial-resolution boundary between homogenized bulk RNA-seq ($n = 24$) and single-cell TME cohorts ($n = 8$, $n = 14$, and pooled $n = 32$).

---

## 5. Mock Decision Letter (*Nature Medicine*)

> **Decision**: *Accept-in-Principle pending verification of substantive analytical and reporting revisions (`T1–T6`) — **VERIFIED CLOSED IN v2.0***
>
> Dear Authors,
>
> Thank you for submitting your manuscript *"Dynamic single-cell niche remodeling and biophysical steric exclusion govern clinical trial attrition across human organs"* to *Nature Medicine*. All four specialist referees and the editorial team recognize the conceptual originality and translational utility of uniting clause-deconvoluted clinical trial registries ($n = 41,396$) with primary human single-cell and spatial Visium biophysics across organs. Before final publication, three core requirements must be met: (1) demonstrate that `DNAcc` remains statistically robust under target-gene cluster-robust sandwich inference and strictly within individual drug modalities (`Small_Molecule` and `Antibody_Biologic`); (2) report all five static comparator baselines transparently in Table 2 and Results §4, including the homogenized bulk NSCLC biopsy boundary behavior and both discordant pre-cirrhotic antifibrotic pairs; and (3) synchronize every cohort denominator, hyperparameter, and four-decimal metric verbatim with the supplementary CSV outputs.
>
> *Post-Revision Editorial Verification*: All six revision tasks (`T1–T6`) have been executed directly in Python source code (`omniniche/stage2/v2_statistical_fortification.py`), deterministic CSV tables (`Table_S3g`, `Table_S3k`), automated unit tests (`37 passed`), and both English and Chinese manuscripts (`manuscript_convention_and_parity_gate.py` `exit code 0`). The manuscript is strongly suited for *Nature Medicine*.

---

## 6. Structured Revision Task Table & Substantive `DONE` Closure Proofs

| ID | Raised By | Axis | Severity | Concern & Evidence Pointer | Executed Strategy & Code / Prose Verification Proof | Action Status | Blocks Submission? |
| :--- | :--- | :--- | :--- | :--- | :--- | :---: | :---: |
| **T1** | Role 1, Role 2 | `data-resource-quality` | Major | Registry-wide stops ($n = 9,406$ across 41,396 trials) vs target-resolved scientific benchmark ($n = 11,843 = 11,072 + 771$ across 29,791 trials) in `Results §1` & `Methods`. | Updated `Results §1` and `Methods §Clause-aware stoppage deconvolution` in EN & ZH to state exact registry-wide ($9,406 / 41,396$) and target-resolved ($11,072 + 771 = 11,843$) integers matching `Table_S1`, `Table_S3a`, `Table_S3e`, and `Extended Data Fig. 1a`. | **`DONE`** | Closed |
| **T2** | Role 3, Role 4 | `data-resource-quality` | Major | Exact `CCR2` ($\tau_0 = 0.9513$, blood $6.19 \to 1.47$), `SPP1` ($0.00 / 5.46 \to 8.31$), `TREM2` ($0.39 \to 5.77$), and `THRB` ($\tau_0 = 0.4052$, $\Delta\tau = +0.0829$) values in `Results §3` & `Methods`. | Synchronized `Results §3` and `Methods` item 2 in EN & ZH verbatim with `Table_S3_TargetGene_6Axis_Drift_Atlas.csv` and `Table_S4c_CCR2_HighTau0_Paradox_Trajectory.csv` (`GSE136103` paired blood monocyte baselines). | **`DONE`** | Closed |
| **T3** | Role 3 | `mechanism-evidence` | Major | Explain both discordant pre-cirrhotic antifibrotic pairs (`46/48 = 0.9583`) in `Table_S4`: F3 simtuzumab ($0.98$) and F2/F3 elafibranor ($0.74$). | Added explicit mechanistic comparison in `Results §3` (EN & ZH) contrasting dual-`PPARA/PPARD` failure elafibranor (`RESOLVE-IT`, $0.74$) against pan-PPAR (`PPARA/D/G`) winner lanifibranor (`NATIVE`, $1.00$) alongside F3 simtuzumab ($0.98$). | **`DONE`** | Closed |
| **T4** | Role 2, Role 5 | `statistical-rigor` | Major | Add Liang–Zeger target-clustered sandwich SEs ($G = 386, 309, 564$), within-modality DeLong stratification, and `C = 50.0` / binary polypharmacology `Methods` parity. | Implemented `_cluster_sandwich_se()` in `fit_multivariable_fixed_effects_logistic()` and `evaluate_within_modality_dnacc_stratification()` in `omniniche/stage2/v2_statistical_fortification.py`; exported `Table_S3g` (`DNAcc` $P_{\mathrm{cluster}} = 2.58 \times 10^{-6}$ on $n = 11,843$ and $8.85 \times 10^{-6}$ on $n = 23,967$, while unadjusted $\tau_0$ collapses to $P_{\mathrm{cluster}} = 0.49$) and `Table_S3k` (`Small_Molecule` $n = 7,369$ AUROC `0.6554 vs 0.5601`, $P = 5.97 \times 10^{-11}$; `Antibody` $n = 2,357$ AUROC `0.6958 vs 0.5515`, $P = 4.36 \times 10^{-10}$); added `test_target_cluster_robust_se_and_within_modality_dnacc` (`37 passed`). | **`DONE`** | Closed |
| **T5** | Role 4 | `figures-and-tables` | Major | Four-decimal CSV parity across `Table 1`, `Table 2`, and `Results §2–§3` (`M2 0.5715`, `M3 0.6582`, monotherapy `242` targets, `n = 23,967` unadj OR `1.88`, fibrosis `AUPRC 0.9583`, LR `Brier 0.1375`). | Updated every cell of `Table 1`, `Table 2`, `Results §2–§3`, and `Discussion` in both EN & ZH to match `Table_S3a/b/e/f/g` and `Table_S4b` to 4 decimal places. | **`DONE`** | Closed |
| **T6** | Role 4, Role 5 | `claim-moderation` | Major | Transparent reporting of all 5 static comparators in `Table 2` and explicit discussion of homogenized bulk ($n = 24$, `tau0 = 0.8000`, OmniNiche #1 Brier `0.2250`) vs single-cell TME cohorts (#1 on all metrics in BCC, NSCLC scRNA, and pooled $n = 32$ `0.7551 vs 0.6235`). | Updated `Table 2` column 6 to report both `CD274`/$\tau_0$ and the highest static comparator across all 5 baselines, and explained bulk lysis homogenization of concentric `myCAF` rings in `Results §4` (EN & ZH). | **`DONE`** | Closed |
