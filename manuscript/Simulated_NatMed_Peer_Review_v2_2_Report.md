# Simulated *Nature Medicine* Five-Role Post-Revision Peer Review Report (Stage II v2.2)

**Manuscript Title:** *Dynamic single-cell niche remodeling and biophysical steric exclusion govern clinical trial attrition across human organs*  
**Stage / Version:** Stage II v2.2 (Post-`T1`–`T10` Full Code, Empirical Data, Figure, and Bilingual Manuscript Revision)  
**Step 0 Pre-Review Gate (`manuscript/REVISION_STATUS.json`):** `✅ PASSED — 10/10 revision tasks (T1–T10) completed, 0 pending, 0 failed`

---

## 1. Manuscript Fingerprint & Routing Banner

**Detected Fingerprint:** `A. Single-Cell / Spatial / Multi-Omics (omics)` + `B. Clinical Cohort / Deconfounded Registry Inference / Translational ML (clinical + ML)` + `D. Mechanistic Biophysics & Reproducible Translational Resource (data-resource)` — activating Role 2 probes for donor-level pseudobulk vs cell-level closure-free statistics, unseen-target `Group5Fold` leakage prevention, Liang–Zeger target-clustered sandwich inference, and missingness sensitivity bounds; Role 3 & Role 5 probes for human-genetics incremental value, response-label circularity, bulk-vs-spatial platform boundaries, and stromal-unlocking clinical toxicology; and Role 4 probes for native 183-mm vector/300-DPI rendering, individual donor/patient points, physical Visium µm calibration, and text-table-figure parity.

---

## 2. Five Independent Reviews

### Role 1: Editor-in-Chief (EIC — *Nature Medicine*)

#### Overall Assessment
In our v2.1 evaluation (`Score: 71.0 / 100`, `Major Revision`), this manuscript presented an original physical-biological thesis—that disease-stage single-cell lineage drift and interstitial macromolecular steric exclusion gate clinical trial realization—but suffered from two fatal reproducibility gaps (`T1` synthetic confusion-matrix counts; `T2` hand-coded drug-compartment routing in the fibrosis arena) and six major evidentiary gaps (`T3` 14-trial fibrosis sample size, `T4` absence of an Open Targets human-genetics benchmark, `T5` hand-assigned ECM weights and single-target `Group5Fold` gaps, `T6` oncology drift staging circularity, `T7` `n = 32` pre-ICB sample size and `umi_stride = 3` subsampling, and `T8` split NLP stoppage parsers).

In Stage II v2.2, the authors have executed every requested code-level, empirical-data, statistical, and visual upgrade without shortcuts. Most decisively for *Nature Medicine* readers, `DNAcc` is now benchmarked head-to-head against Open Targets v24.09 human genetic evidence across 11,637 deconfounded scientific trials and 23,761 maturation trials (`Joint Human Genetics + DNAcc` AUROC = `0.6977` vs `0.6362` for genetics alone, DeLong $P = 7.13 \times 10^{-10}$; maturation `0.7190` vs `0.6631`, $P = 2.19 \times 10^{-103}$), while static healthy tissue specificity ($\tau_0$) adds zero gain to genetics ($P = 0.451$). The translational validation cohorts now reach *Nature Medicine* evidentiary scale: 24 live ClinicalTrials.gov SHA-256-locked Phase IIb/III antifibrotic trials (`AUROC = 0.8714`, exact $P = 0.0010$; Core-14 `0.9583`) compiled via objective single-cell expression routing, and 107 strictly pre-treatment anti-PD-1 solid-tumor biopsies (`AUROC = 0.6367`, exact $P = 0.0187$; single-cell TME cohorts `AUROC = 0.8750` in both BCC and NSCLC). Scope boundaries (peripheral interstitial ECM remodeling vs blood-brain barrier neurovascular transport, intravascular PBMC sepsis null boundary, and bulk lysis attenuation in melanoma) are framed with exemplary scientific transparency.

#### Scoring Rubric (Role 1: 93 / 100)
| Criterion | Max | v2.1 | v2.2 | Justification |
| :--- | :---: | :---: | :---: | :--- |
| Novelty | 30 | 24 | **28** | Establishes orthogonal synergy between germline human genetics (causal target validity) and dynamic single-cell + steric tissue accessibility (`DNAcc`). |
| Significance | 30 | 22 | **28** | Deconfounds 41,396 Phase II–IV trials ($\kappa = 0.941$), resolves macromolecular futility ($\mathrm{OR} = 1.57$), and transfers zero-label to 24 fibrosis trials and 107 pre-ICB biopsies. |
| Journal Fit (*Nature Medicine*) | 20 | 14 | **18** | Directly informs target selection, modality matching (small molecule/peptide vs IgG), and biomarker-enriched Phase II/III trial design. |
| Presentation | 20 | 13 | **19** | Crisp 305-line 1:1 EN/ZH prose, explicit Block A/B Table 1 structure, and 9 publication-grade 183-mm multi-panel figures (44 panels). |

#### Numbered Comments
1. `[Axis: novelty-significance | Severity: Minor]` The conceptual distinction between germline human genetics ("does perturbing the target modify disease?") and `DNAcc` ("does the target remain in a physically accessible compartment when the chosen drug modality arrives?") is compelling and well supported by the `0.7012` AUROC strictly within genetically supported targets ($g_{\max} \ge 0.10$, $n = 7,174$). Keeping the responsible-AI guardrail that outputs are hypothesis-generating and must never be used alone to withhold standard-of-care therapy is appropriate for clinical readership. ***[原文锚点：Results §3 "Within already genetically supported targets (Open Targets genetic score >= 0.10, n = 7,174 trials), DNAcc reached AUROC = 0.7012"; Discussion §5 "OmniNiche outputs are intended as mechanistic hypothesis-generating instruments"]***

**Confidence:** `high`

---

### Role 2: Methods & Biostatistics Expert

#### Overall Assessment
I audited the Python analytical codebase (`macro_clinical_census.py`, `pan_disease_drift_law.py`, `fibrosis_mash_deconvolution.py`, `solid_tumor_spatial_icb.py`, `v2_statistical_fortification.py`), the 56 empirical files in `outputs/stage2_tables/`, and the 42-test `pytest` suite. All methodological flaws flagged in v2.1 (`T1, T2, T4, T5, T6, T7, T8`) have been rigorously eliminated at the code and data levels:
- **Zero synthetic PRNG validation (`T1, T8`)**: The synthetic confusion-matrix generator is completely removed. A single unified clause-aware negation-stripping parser operates across all 41,396 trials and is benchmarked head-to-head against a Naive Keyword Baseline on 240 manually adjudicated verbatim `whyStopped` narratives (`95.83%` vs `70.00%` accuracy; $\kappa = 0.9413$ vs $0.6105$).
- **Leakage-free single-cell & cross-validation architecture (`T5, T6`)**: `AXIS_FIBROTIC_BARRIER` hand constants are replaced by empirical 10-gene collagen-myofibroblast ECM densities $D_{\mathrm{ECM}} \in [0.070, 0.350]$ (`Table_S3m`). Solid-oncology $\Delta\tau$ is now computed strictly from response-label-free early primary BCC vs advanced resected NSCLC tumors. Unseen-target `Group5Fold` CV is now reported inside strictly single-target monotherapies ($n = 6,115$ scientific `Group5Fold` $P = 1.35 \times 10^{-4}$; $n = 10,507$ maturation `Group5Fold` AUROC = `0.6375`, $P = 7.39 \times 10^{-123}$), and within-modality DeLong tests are significant across all three modalities (small molecules $P = 2.04 \times 10^{-13}$, IgG antibodies $P = 0.0133$, proteins/peptides $P = 0.0045$).
- **Cluster-robust inference & sensitivity regimes (`T4, T5`)**: Multivariable fixed-effects logistic regression uses Liang–Zeger target-clustered sandwich covariance across $G = 374$ and $G = 561$ primary target clusters, including a simultaneous three-pillar specification (`DNAcc` + out-of-fold Open Targets genetics + static $\tau_0$ + 5 confounder families: `DNAcc` $\mathrm{OR} = 1.38$, $P_{\mathrm{cluster}} = 2.39 \times 10^{-4}$; $\tau_0$ $\mathrm{OR} = 1.00$, $P_{\mathrm{cluster}} = 0.98$). All seven sensitivity regimes (`R1`–`R7`, including worst-case failure assignment of all 2,618 target-resolved unreported stops in `R3`, $n = 14,255$, AUROC = `0.6541`, $P = 1.50 \times 10^{-70}$) are exact and reproducible.
- **Closure-free & LOPO zero-label biophysics (`T2, T7`)**: Donor fibrosis coupling uses cell-level $\log_2(\mathrm{CPM}+1)$ module scores rather than simplex-constrained proportions, `GSE123813` is streamed at exact single-UMI fidelity (`umi_stride = 1`), and multi-cohort oncology pooling ($n = 83$ and $n = 107$) uses zero-label leave-one-patient-out (`_lopo_robust_01`) rank scaling alongside raw physical Langevin probabilities (`13.7% -> 78.1%` / `74.6%`).

#### Scoring Rubric (Role 2: 92 / 100)
| Criterion | Max | v2.1 | v2.2 | Justification |
| :--- | :---: | :---: | :---: | :--- |
| Experimental Design | 25 | 17 | **23** | Clean separation of discovery/compilation vs zero-label external clinical trial (24 RCTs) and pre-treatment biopsy ($n = 107$) transfer. |
| Statistical Rigor | 30 | 21 | **28** | Paired DeLong, 2,000-bootstrap CIs, 10,000-permutation exact *P*-values, Liang–Zeger target-clustered sandwich SEs, and 7 sensitivity regimes (`R1`–`R7`). |
| Analytical Reproducibility | 25 | 16 | **23** | 100% deterministic empirical CSV/NPZ outputs (56 archive files), `< 2.3 GB` peak RAM, and 42/42 automated `pytest` checks passing. |
| Methods Completeness | 20 | 15 | **18** | Explicit equations, initial conditions, stage-blind `FixedS0` baseline, LOPO normalization, and denominator definitions documented. |

#### Numbered Comments
1. `[Axis: statistical-rigor | Severity: Minor]` In `Pooled Pure Anti-PD-1 Monotherapy Pre-Treatment Biopsies` ($n = 83$, 28 responders), OmniNiche achieves higher AUROC (`0.6016`) than `No-Visium` (`0.5815`) and all clinical biomarker signatures (`PD-L1` `0.5338`, `18-gene GEP` `0.5442`, `CAF/TIDE` `0.5571`, `CD8:myCAF` `0.5873`), though the 10,000-permutation *P*-value on the $n = 83$ subset alone is `0.1302` due to heavy weighting by homogenized melanoma bulk biopsies (`GSE91061` $n = 49$), reaching statistical significance ($P = 0.0187$, AUROC = `0.6367`) across the full $n = 107$ strictly pre-treatment cohort. Both numbers are transparently reported in Table 2 and Discussion §4. ***[原文锚点：Table 2 "Pooled Pure Anti-PD-1 Monotherapy Pre-Treatment Biopsies (LOPO Zero-Refit) | 83 patients | 0.6016 | P = 0.1302" & "Pooled All Strictly Pre-Treatment Human Tumor Biopsies | 107 patients | 0.6367 | P = 0.01870"]***

**Confidence:** `high`

---

### Role 3: Translational Fibrosis & Immuno-Oncology Domain Expert

#### Overall Assessment
From a disease-biology and clinical-pharmacology standpoint, v2.2 resolves every domain concern I raised in v2.1 (`T2, T3, T7, T10`). Expanding the fibrosis benchmark to 24 Phase IIb/III trials across liver MASH ($n = 13$), lung IPF ($n = 8$), and kidney CKD ($n = 3$)—including recent landmark readouts such as `SYNERGY-NASH` (tirzepatide), `HARMONY` (efruxifermin), `FIBRONEER-IPF` (nerandomilast), `DAPA-CKD` (dapagliflozin), `ZEPHYRUS-1` (pamrevlumab), `STELLAR-3/4` and `MOSAIC` (selonsertib), and `ISABELA 1` (ziritaxestat)—provides a balanced, clinically realistic stress test (`AUROC = 0.8714`, $P = 0.0010$). Crucially, replacing manual compartment rules with objective single-cell expression routing plus upstream metabolic/hemodynamic injury unloading (`GLP1R`, `GIPR`, `SLC5A2`) accurately aligns with Jara et al. (*Nat. Med.* 2025) showing absent intrahepatic `GLP1R` alongside systemic suppression of `TREM2` and fibrotic remodeling. Moreover, the explicit analysis of the `ISABELA 1` ziritaxestat false positive (`0.9667`, halted for systemic benefit-risk/off-niche toxicity despite pulmonary access) adds high clinical credibility.

In immuno-oncology, scaling to 107 strictly pre-treatment biopsies across cutaneous BCC (`GSE123813`), melanoma (`GSE91061`, `GSE78220`), and NSCLC (`GSE207422`), benchmarking directly against `18-gene GEP`, `CAF/TIDE`, and `CD8:myCAF`, and adding `Discussion` §3 on the historical cardiac valvulopathy of systemic pan-`TGFBR1`/ALK5 inhibitors versus inhaled ALK5 inhibitors, $LRRC15^+$ myofibroblast targeting, and `TREM2/SPP1` or `CXCL9:SPP1` myeloid reprogramming (`74.6%` simulated rescue without pan-ALK5 blockade) elevates the translational sophistication to *Nature Medicine* standards.

#### Scoring Rubric (Role 3: 91 / 100)
| Criterion | Max | v2.1 | v2.2 | Justification |
| :--- | :---: | :---: | :---: | :--- |
| Literature Coverage | 30 | 22 | **28** | 56 live-verified PubMed citations covering 2024–2026 genetics, MASH/IPF/CKD Phase III RCTs, spatial oncology, and clinical AI benchmarks. |
| Positioning Accuracy | 25 | 18 | **23** | Accurately positioned as a zero-label mechanistic complement to germline genetics (Open Targets) and supervised clinical AI (SCORPIO, COMPASS, HEX-MICA). |
| Discussion Depth | 25 | 17 | **22** | Nuanced treatment of ALK5 cardiotoxicity vs lineage/organ-selective stromal modulation, `ISABELA 1` toxicity boundary, and bulk-vs-scRNA spatial attenuation. |
| Forward-Looking Utility | 20 | 14 | **18** | Concrete implications for Phase II/III modality selection, histological stage gating, and prospective spatial-biopsy combination trials. |

#### Numbered Comments
1. `[Axis: clinical-validity | Severity: Minor]` In `GSE207422` ($n = 24$ pre-treatment NSCLC bulk biopsies), patients received neoadjuvant anti-PD-1 concurrently with platinum-doublet chemotherapy, whereas `GSE123813` ($n = 8$), `GSE91061` ($n = 49$), and `GSE78220` ($n = 26$) represent pure anti-PD-1 monotherapy ($n = 83$). Separating `Pooled Pure Anti-PD-1 Monotherapy` ($n = 83$) from `Pooled All Strictly Pre-Treatment Biopsies` ($n = 107$) in Results §6, Table 2, and Fig. 5c,d resolves regimen heterogeneity cleanly. ***[原文锚点：Results §6 "Pooled Pure Anti-PD-1 Monotherapy Pre-Treatment biopsies (n = 83, 28 responders: AUROC = 0.6016...) and Pooled All Strictly Pre-Treatment Biopsies (n = 107, 37 responders: AUROC = 0.6367...)"]***

**Confidence:** `high`

---

### Role 4: Data, Visualization & Reproducibility Expert

#### Overall Assessment
I inspected all 9 rendered publication figures (`Fig. 1`–`Fig. 5` and `Extended Data Fig. 1`–`Extended Data Fig. 4`; 44 panels total across PNG/PDF/SVG triplets at 183-mm width), passed them through `figure_clean_aesthetic_gate.py` (`ALL PASSED`), and cross-checked every displayed number against `OmniNiche_Flagship_Manuscript_EN.md`, `OmniNiche_Flagship_Manuscript_ZH_Preview.md`, and `Supplementary_Tables_Master_Index.md` (`56` archive files):
- **`Fig. 2f`**: Now faithfully displays the head-to-head `Manual_Clinical_Adjudication_240` confusion matrix (`95.8%`, $\kappa = 0.941$) and per-class comparison against the Naive Keyword Baseline (`70.0%`, $\kappa = 0.611$).
- **`Fig. 3b,c,d`**: Cleanly integrates the `Human Genetics` (`0.6362`) and `Genetics + DNAcc` (`0.6977`) benchmark bars, the 3-pillar Liang–Zeger target-clustered forest plot (`OR = 1.38`, $P_{\mathrm{cluster}} = 2.39 \times 10^{-4}$), and significant within-modality gains across all 3 drug classes.
- **`Fig. 4c,d,e`**: Displays all 24 Phase IIb/III antifibrotic trials with aligned `Org / Mod / Acc` tracks, the 24-trial + Stage-Blind classification arena (`0.8714` and `0.7643`), and the updated 24-trial bidirectional $3 \times 3$ + `LOOO` transfer matrix (`0.8357–0.8821`).
- **`Fig. 5c,d,g`**: Scales cleanly to the 107-patient pre-treatment waterfall + 4-row OncoPrint, the multi-cohort benchmark including `No-Visium`, `18-gene GEP`, `CAF/TIDE`, and `CD8:myCAF`, and the 3-arm physical Langevin rescue trajectory across $n = 70$ clinical non-responders (`13.7% -> 78.1%` and `74.6%`).
- **`Table 1` (`T9`)**: Explicit **Block A** ($n = 11,637$, $G = 374$) and **Block B** ($n = 23,761$, $G = 561$) headers completely eliminate any ambiguity between the Deconfounded Scientific Cohort and the per-axis Definitive Maturation rows.

#### Scoring Rubric (Role 4: 92 / 100)
| Criterion | Max | v2.1 | v2.2 | Justification |
| :--- | :---: | :---: | :---: | :--- |
| Main Figures (`Fig. 1`–`Fig. 5`) | 25 | 19 | **23** | High-density 4-to-8 panel *Nature Medicine* / *Science* layouts, zero text overlap or clipping, colorblind-safe Okabe–Ito/Viridis palettes, physical 500-µm Visium scale bars. |
| Legend Completeness | 25 | 18 | **23** | Every panel specifies exact sample size ($n$), test side/type, error bar definition, and exact *P*-values matching the CSVs. |
| Data Transparency | 25 | 17 | **23** | Individual donor/patient points overlaid on summary bars (`Fig. 4a,b`, `Fig. 5e,g`, `ED Fig. 4d`); 56 deterministic CSV/NPZ/JSON files archived. |
| Supplementary & Parity | 25 | 18 | **23** | 100% 1:1 EN/ZH line parity (`305 == 305` lines), 56/56 live-verified PubMed citations in strict `[1]..[56]` order, and complete Master Index (`Table S1`–`Table S5k`). |

#### Numbered Comments
1. `[Axis: figures-and-tables | Severity: Minor]` In `Fig. 5c`, with $n = 107$ individual patient bars and 4 OncoPrint rows across 183-mm width, individual patient IDs are appropriately omitted from the x-axis to prevent label crowding while full per-patient identifiers and predictions are archived in `Table_S5i_Pooled_PreTreatment_107Patients_Master.csv`. ***[原文锚点：Fig. 5c & Supplementary Table S5i]***

**Confidence:** `high`

---

### Role 5: Devil's Advocate (Adversarial Stress-Test — Non-Scoring)

#### Overall Assessment
My mandate is to attack the revised v2.2 manuscript at its hardest remaining joints. First, I confirm that all five potentially lethal vulnerabilities I raised against v2.1—synthetic PRNG validation (`T1`), hand-coded fibrosis compartment maps (`T2`), omitting human genetics (`T4`), oncology response-label drift circularity (`T6`), and small $n = 32$ / subsampled pre-ICB cohorts (`T7`)—have been genuinely fixed in code and data rather than papered over rhetorically. I now raise three sharp post-revision stress challenges to test whether any residual overclaiming remains in v2.2:

#### Adversarial Challenges
1. `[Axis: clinical-validity | Severity: Minor-Major (Properly Bounded in v2.2)]` **Bulk RNA-seq homogenization blunts spatial-Langevin advantage over algebraic linear summaries in cutaneous/metastatic melanoma.** While OmniNiche clearly beats every baseline including `Uninformative Linear Mean` in single-cell TME cohorts (`GSE123813` BCC scRNA: `0.8750` vs `0.7500`; `GSE207422` NSCLC scRNA: `0.8750` vs `0.8250`, ranking 1st of 10 methods) and beats `PD-L1`, `18-gene GEP`, `CAF/TIDE`, `CD8:myCAF`, and `No-Visium` in pooled pre-treatment biopsies ($n = 107$: `0.6367`), its AUROC in homogenized melanoma bulk RNA-seq (`GSE91061`: `0.5385`; `GSE78220`: `0.5982`) is modest and slightly below a simple uninformative linear mean on the pooled bulk-heavy $n = 107$ set (`0.6367` vs `0.6442` in `Table_S5`). *Verdict on v2.2 handling:* The authors do **not** hide this; `Table_S5` archives all 10 rows per cohort, `Results` §6 explicitly states *"Spatial-Langevin modeling excelled wherever single-cell microenvironment resolution was preserved"* and reports `0.5385` and `0.5982` verbatim, and `Discussion` §4 explicitly identifies bulk lysis homogenization of concentric `myCAF` geometry as Limitation 3. This honest boundary reporting prevents overclaiming. ***[原文锚点：Results §6 "0.5982 in GSE78220 melanoma... and 0.5385 in homogenized GSE91061 melanoma biopsies"; Discussion §4 "bulk RNA-seq homogenization in cutaneous and metastatic melanoma biopsies (GSE91061 AUROC = 0.5385; GSE78220 0.5982) attenuates concentric myCAF geometry"]***
2. `[Axis: mechanism-evidence | Severity: Minor (Properly Bounded in v2.2)]` **Interstitial accessibility cannot predict off-niche chemical or systemic benefit-risk toxicity (`ISABELA 1` ziritaxestat).** In the expanded 24-trial fibrosis arena, the small-molecule autotaxin (`ENPP2`) inhibitor ziritaxestat (`ISABELA 1`) receives a high simulated remission score (`0.9667`) yet failed in Phase III due to benefit-risk/mortality imbalance. *Verdict on v2.2 handling:* Highlighting `ISABELA 1` explicitly in both `Results` §4 and `Discussion` §4 as the sole high-scoring false positive—and explaining mechanistically why an interstitial-transport and target-engagement model cannot anticipate idiosyncratic systemic chemical toxicity—turns a potential critique into a clear domain-of-applicability boundary. ***[原文锚点：Results §4 "The sole high-scoring failure among all 24 trials was the small-molecule autotaxin (ENPP2) inhibitor ziritaxestat (ISABELA 1, 0.9667)"; Discussion §4 "OmniNiche models local target-niche engagement and interstitial transport rather than idiosyncratic off-target or systemic chemical toxicity"]***
3. `[Axis: claim-moderation | Severity: Minor (Closed via T10 in v2.2)]` **Translational feasibility of `myCAF/TGFBR1` blockade.** Systemic pan-ALK5 inhibitors caused hemorrhagic cardiac valvulopathy in preclinical and early clinical studies. *Verdict on v2.2 handling:* `Discussion` §3 now explicitly warns against systemic pan-`TGFBR1` cardiotoxicity and pairs this caveat with both literature-backed organ/lineage-selective alternatives (inhaled ALK5 inhibitors [52], $LRRC15^+$ myofibroblast targeting [53]) and an empirical computation showing that `anti-PD-1 + SPP1/TREM2` myeloid reprogramming achieves `74.6%` simulated rescue ($P = 1.78 \times 10^{-13}$) without pan-`TGFBR1` inhibition (`Fig. 5g`). ***[原文锚点：Discussion §3 "first-generation systemic pan-TGFBR1 (ALK5) small-molecule inhibitors historically encountered dose-limiting cardiac valvulopathy... myeloid TREM2/SPP1 or CXCL9:SPP1 polarity reprogramming [51], which achieved nearly comparable simulated responder rescue (74.6% versus 78.1%)"]***

**Confidence:** `high`

---

## 3. Aggregated Score & Rubric Summary Table

| Reviewer Role | Weight | v2.1 Score | v2.2 Score | Status of Prior Blocking Items (`T1`–`T8`) |
| :--- | :---: | :---: | :---: | :--- |
| **Role 1: Editor-in-Chief (*Nature Medicine*)** | 25% | 73 / 100 | **93 / 100** | All editorial scope, genetics-benchmark, and scale requirements met |
| **Role 2: Methods & Biostatistics Expert** | 25% | 69 / 100 | **92 / 100** | `T1, T4, T5, T6, T8` 100% resolved in Python code, CSVs & 42/42 `pytest` |
| **Role 3: Fibrosis & Immuno-Oncology Expert** | 25% | 71 / 100 | **91 / 100** | `T2, T3, T7, T10` 100% resolved (24 fibrosis RCTs, 107 pre-ICB biopsies, ALK5 safety) |
| **Role 4: Data, Visualization & Parity Expert** | 25% | 72 / 100 | **92 / 100** | `T9` + all 9 figures (44 panels) + 56 archive files + 5-lock manuscript gate green |
| **Role 5: Devil's Advocate** | Non-scoring | 5 Fatal/Major | **0 Fatal / 0 Major** | Residual bulk-homogenization & toxicity boundaries explicitly bounded in text |
| **Weighted Panel Total** | **100%** | **71.25 / 100** | **92.00 / 100** | **Decision Upgraded: `Major Revision` $\to$ `Accept in Principle`** |

---

## 4. Cross-Review Consensus (`Cc.1`–`Cc.3`)

- **`Cc.1` [Raised by: Role 1, Role 2, Role 3, Role 5 | Axis: `novelty-significance` & `statistical-rigor` | Status: `RESOLVED IN v2.2`]**: Head-to-head evaluation against Open Targets v24.09 human genetic evidence confirms that `DNAcc` provides orthogonal, statistically robust predictive value beyond germline genetics (`0.6977` vs `0.6362` scientific, $P = 7.13 \times 10^{-10}$; `0.7190` vs `0.6631` maturation, $P = 2.19 \times 10^{-103}$; 3-pillar target-clustered $\mathrm{OR} = 1.38$, $P_{\mathrm{cluster}} = 2.39 \times 10^{-4}$), whereas static healthy $\tau_0$ adds zero value ($P = 0.451$).
- **`Cc.2` [Raised by: Role 2, Role 3, Role 5 | Axis: `reproducibility` & `clinical-validity` | Status: `RESOLVED IN v2.2`]**: Eliminating synthetic validation counts (`T1`: independent 240-trial manual adjudication, $\kappa = 0.9413$) and hand-coded drug-compartment maps (`T2`: objective single-cell pseudobulk routing + closure-free donor module coupling) while expanding external validation to 24 Phase IIb/III antifibrotic trials (`T3`: `AUROC = 0.8714`, $P = 0.0010$) and 107 pre-treatment anti-PD-1 biopsies (`T7`: `AUROC = 0.6367`, $P = 0.0187$; scRNA TME `0.8750`) places all translational claims on auditable empirical footing.
- **`Cc.3` [Raised by: Role 3, Role 5 | Axis: `claim-moderation` | Status: `PROPERLY BOUNDED IN v2.2`]**: Domain-of-applicability limits—specifically bulk RNA-seq homogenization in melanoma (`0.5385–0.5982` vs `0.8750` in single-cell TME cohorts), off-niche systemic chemical toxicity (`ISABELA 1` ziritaxestat), intravascular PBMC sepsis (`0.3841`), and systemic pan-`TGFBR1` cardiotoxicity—are transparently quantified and bounded in `Results` and `Discussion`.

---

## 5. Mock Decision Letter (*Nature Medicine*)

**Decision:** **Accept in Principle (Zero Blocking Concerns Remaining)**

Dear Authors,

Thank you for submitting your revised Stage II v2.2 manuscript, *"Dynamic single-cell niche remodeling and biophysical steric exclusion govern clinical trial attrition across human organs"*, together with the updated source code, 56-file empirical data archive, 42-test verification suite, and revised display items.

All five members of our review panel have re-evaluated the manuscript and its underlying computational artifacts. The reviewers agree that your team has thoroughly resolved every methodological, statistical, and translational concern raised in the previous round (`T1`–`T10`). In particular, demonstrating orthogonal synergy between Dynamic Niche Accessibility (`DNAcc`) and Open Targets human genetic evidence across 11,637 deconfounded scientific trials and 23,761 maturation trials, replacing manual fibrosis compartment rules with objective single-cell expression routing across 24 Phase IIb/III trials, and expanding zero-label pre-treatment immunotherapy validation to 107 human tumor biopsies across single-cell and bulk cohorts establish a rigorous and broadly useful translational framework for *Nature Medicine* readers. All reviewers also commend the transparent framing of anatomical and platform boundary conditions (intravascular sepsis, off-niche chemical toxicity, bulk melanoma homogenization, and lineage-selective alternatives to systemic ALK5 inhibition).

No blocking analytical, visual, or textual revisions remain (`0` Fatal, `0` Major). Your Stage II v2.2 manuscript package is ready for final submission formatting.

Sincerely,  
**Chief Editor, *Nature Medicine* (Simulated Five-Role Review Panel)**

---

## 6. Structured Revision Task Audit Ledger (`T1`–`T10` Closure Verification)

| ID | Raised By (v2.1) | Axis | v2.1 Severity | Resolution Verified in Stage II v2.2 | Action Status | Blocks Submission? |
| :--- | :--- | :--- | :---: | :--- | :---: | :---: |
| **T1** | Role 2, Role 5 | `reproducibility` | Fatal | Removed synthetic PRNG confusion generator; audited unified clause-aware parser (`95.83%`, $\kappa = 0.9413$) vs Naive Keyword Baseline (`70.00%`, $\kappa = 0.6105$) on 240 manually adjudicated trials (`Table_S2a`, `Fig. 2f`). | **`DONE`** | No |
| **T2** | Role 2, Role 3, Role 5 | `mechanism-evidence` | Fatal | Replaced `DRUG_COMPARTMENT_MAP` & manual `pIC50` with objective organ scRNA pseudobulk routing ($\ge 58\%$ peak), $\tau_0$ sink dilution, unsigned ridge-GRN leverage (`pIC50 = 7.5`), closure-free module Pearson coupling, and `StageBlind_FixedS0` (`Table_S4b`, `Fig. 4a,d`). | **`DONE`** | No |
| **T3** | Role 1, Role 3 | `clinical-validity` | Major | Expanded fibrosis benchmark from 14 to 24 live ClinicalTrials.gov SHA-256-locked Phase IIb/III RCTs across liver ($n=13$), lung ($n=8$), and kidney ($n=3$) (`AUROC = 0.8714`, $P = 0.0010$; Core-14 `0.9583`; `Table_S4`, `Fig. 4c–e`). | **`DONE`** | No |
| **T4** | Role 1, Role 2, Role 5 | `novelty-significance` | Major | Benchmarked `DNAcc` head-to-head against Open Targets v24.09 human genetic evidence (`Genetics + DNAcc` `0.6977` vs `0.6362`, $P = 7.13 \times 10^{-10}$; maturation `0.7190` vs `0.6631`, $P = 2.19 \times 10^{-103}$; 3-pillar clustered $\mathrm{OR} = 1.38$; `Table_S3l`, `Table 1`, `Fig. 3b,c`). | **`DONE`** | No |
| **T5** | Role 2, Role 5 | `statistical-rigor` | Major | Derived empirical $D_{\mathrm{ECM}} \in [0.070, 0.350]$ (`Table_S3m`); added single-target `Group5Fold` (`0.5872` / `0.6375`), significant within-modality `DNAcc` gains across all 3 modalities (`Table_S3k`, `Fig. 3d`), and `R6`/`R7` sensitivity bounds (`Table_S3e`, `ED Fig. 3c`). | **`DONE`** | No |
| **T6** | Role 2, Role 5 | `experimental-design` | Major | Computed solid-oncology $\Delta\tau$ strictly from response-label-free early primary BCC vs advanced resected NSCLC tumors; reported ectopic switch across all profiled (`17.8%–37.1%`) and expressed (`21.0%–41.3%`) targets (`Table_S3`, `ED Fig. 2b`). | **`DONE`** | No |
| **T7** | Role 1, Role 2, Role 3, Role 5 | `clinical-validity` | Major | Set `umi_stride = 1`; added `GSE91061` ($n=49$) & `GSE78220` ($n=26$) melanoma monotherapy cohorts ($n = 107$ pre-ICB biopsies, $n = 83$ pure monotherapy); added `No-Visium`, `18-gene GEP`, `CAF/TIDE`, `CD8:myCAF`, and `SPP1/TREM2` TAM rescue (`74.6%`; `Table_S5`, `Table_S5i`, `Fig. 5c,d,g`). | **`DONE`** | No |
| **T8** | Role 2, Role 4 | `reproducibility` | Major | Unified clause-aware parser across `macro_clinical_census.py` & `v2_statistical_fortification.py` (`4,838` admin stops = `51.4%`; 2020–2025 `64.0%`, 2021 peak `66.3%`) and reconciled all `29,791` modality subcategories (`Results` §1, `Methods`). | **`DONE`** | No |
| **T9** | Role 4 | `figures-and-tables` | Minor-Major | Added explicit **Block A** ($n = 11,637$, $G = 374$) and **Block B** ($n = 23,761$, $G = 561$) headers in `Table 1` and synchronized all 9 figures & legends. | **`DONE`** | No |
| **T10** | Role 3, Role 5 | `claim-moderation` | Minor | Added explicit `Discussion` §3 caveat on systemic pan-ALK5 cardiac valvulopathy and live-verified citations for inhaled ALK5 [52], $LRRC15^+$ myofibroblast [53], and `TREM2/SPP1` / `CXCL9:SPP1` [51] alternatives. | **`DONE`** | No |
