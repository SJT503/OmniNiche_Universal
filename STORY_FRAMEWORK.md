# STORY_FRAMEWORK & STAGE II v2.0 *NATURE MEDICINE* UPGRADE BLUEPRINT

## 1. Strategic Anchor
- **First-Choice Target Journal**: ***Nature Medicine*** (Article;同時吸纳 *Science* 顶刊的机制动力学严密性与高密度视觉叙事标准).
- **Core Scientific Paradigm Shift**:
  - Human germline genetics (`Nelson 2015`; `Minikel et al., Nature 2024`) and healthy single-cell atlases (`HPA` $\tau_0$) answer **"Who"** (whether a gene participates in disease etiology or shows baseline cell-type selectivity), yet late-stage Phase II–III clinical attrition still exceeds 85%.
  - Using `41,396` interventional Phase II–IV clinical trials across **12 systemic and parenchymal organ disease axes**, `OmniNiche` proves *why*: (1) `41.0%` of premature trial closures are administrative/accrual noise that corrupts unstratified registry benchmarks; (2) once deconfounded (`n = 11,843` explicit biological outcomes; `n = 23,967` definitive maturation trials), healthy $\tau_0$ performs near chance (`AUROC = 0.528–0.541`) because active human tissue remodeling switches the dominant expressing lineage in `21.5%–37.1%` of clinical targets while assembling macromolecular stromal barriers (`myCAF` / cross-linked ECM).
  - We formulate **Dynamic Niche Accessibility (`DNAcc`)** — coupling single-cell remodeling drift with modality-dependent stromal barrier compatibility — and embed it into a profile-likelihood-audited multicellular **Itô-Langevin stochastic attractor engine**, achieving pan-disease clinical maturation stratification, **bidirectional zero-shot cross-organ antifibrotic Phase IIb/III trial prediction**, and **blind pre-treatment solid-tumor anti-PD-1 patient response & combination rescue prediction**.

---

## 2. Execution Roadmap (Strictly Ordered Phases for the New Window)

### Phase 0 — Deep *Nature Medicine* & *Science* Same-Genre Benchmark Teardown (MANDATORY FIRST STEP)
Before touching any script or figure, activate `biomed-literature`, `journal-benchmark-review`, `manuscript-lift-playbook`, and `top-journal-storycraft` to conduct a deep structural, statistical, and visual teardown of **6–8 real 2024–2026 *Nature Medicine* and *Science* exemplar articles** across:
1. Clinical-trial success / genetics / target-prioritization macro-analyses (e.g., `Minikel et al., Nature 2024` + recent *Nat. Med.* clinical trial attrition / AI trial enrichment studies);
2. Pan-organ single-cell fibrosis / macrophage-stromal niche atlases (*Science* / *Nat. Med.*);
3. Spatial TME / multi-modal immunotherapy responder predictors in *Nature Medicine* (examining how *Nat. Med.* reports calibration, Decision Curve Analysis net clinical benefit, trial sample-size enrichment simulations, multi-cohort forest plots, waterfall plots, CIs, and confounding control).
- **Deliverable**: Save `F:/OmniNiche_Universal/manuscript/Phase0_NatMed_Science_Deep_Benchmark_Teardown.md` and extract a concrete **"Gap-Zero Checklist"** (including any overlooked blind spots found during deep benchmarking) before coding Phase 1–3!

---

### Phase 1 — Complete De-Scarring of `"Non-CNS"` Framing & Unified Concept Handle (`DNAcc`)
1. Remove `"Non-CNS"`, `"non-CNS"`, `"non-neurological"`, and `"Firewall-2"` **100%** from external scientific framing (Title, Abstract, Introduction, Results, Discussion, Tables 1–2, Figure 1–5 titles/axes/legends, Extended Data Figs).
2. Keep the neurological exclusion filter silently enforced in code (`0%` CNS overlap), and state the cohort scope naturally in **one biophysical sentence in *Methods*** (focusing on 12 systemic and parenchymal organ disease axes where interstitial ECM architecture, rather than neurovascular blood-brain barrier transport, governs drug and immune access).
3. Introduce and consistently use **Dynamic Niche Accessibility (`DNAcc`)** as the unifying quantitative index across text, equations, tables, and figures.

---

### Phase 2 — Statistical & Methodological Fortification + Proactive Gap-Hunting ("查漏补缺")
Implement deterministic, unit-tested quantitative upgrades (`< 7 GB` RAM) outputting new/updated CSVs in `F:/OmniNiche_Universal/outputs/stage2_tables/`:
1. **Full Uncertainty & Exact Small-Sample Inference**:
   - 2,000-iteration stratified bootstrap **95% CIs** for every AUROC, AUPRC, and Brier score across Tables 1–2 and Figs. 3–5.
   - Paired **DeLong test *P*-values** for large-scale trial comparisons (`n = 11,843`, `n = 5,815`, `n = 23,967`).
   - **10,000-permutation exact two-sided *P*-values** for small-`n` translational cohorts (`n = 14` Phase IIb/III antifibrotic trials, `n = 8` BCC, `n = 14` NSCLC scRNA, `n = 24` NSCLC pre-biopsy bulk, `n = 32` pooled pre-treatment biopsies).
2. ** Multi-Agent Stoppage Classifier Audit & `3,118` Unreported-Stop Sensitivity Bounds**:
   - Stratified manual/rule-independent benchmark validation (`n >= 200` audited stopped trials with explicit ground-truth rationales) reporting **precision, recall, macro-F1, confusion matrix, and Cohen's $\kappa$**.
   - **Missingness sensitivity bounds** across the `3,118` (`47.3%`) `Unreported / Other` stopped trials: test (i) strict deconfounded baseline, (ii) worst-case all-unreported = failure, (iii) inverse-probability / propensity-weighted re-estimation, and (iv) results-posted-only subset — proving `DNAcc` significantly outperforms $\tau_0$ across every boundary assumption.
3. **6-Step Orthogonal Component Ablation, Unseen-Target `GroupKFold`, Multivariable Fixed-Effects Control & Clinical Decision Utility**:
   - Stepwise orthogonal ablation: `Healthy tau0` $\to$ `+Vital-Organ Bystander` $\to$ `+Disease Early tau` $\to$ `+Dynamic Drift (Delta tau)` $\to$ `+Lineage Compartment Switch` $\to$ `+Biophysical Barrier (Full DNAcc)`.
   - **Strict `GroupKFold` by primary target gene** (evaluating generalization to **unseen target genes** never observed in training folds!) + **multivariable logistic regression adjusting for disease-axis, clinical phase, modality, and calendar-era fixed effects** (proving `DNAcc` is not a proxy for indication baseline rates).
   - **10-decile Probability Calibration Curves** + ***Nature Medicine*-standard Decision Curve Analysis (DCA net clinical benefit)** + **Simulated Phase II-to-III Clinical Trial Portfolio Enrichment** (quantifying futile Phase III patient-exposure & trial-cost reduction when gating by `DNAcc`).
4. **Bidirectional / Leave-One-Organ-Out (`3×3`) Cross-Organ Fibrosis Transfer**:
   - Expand `Liver -> Lung/Kidney` to the full **bidirectional / Leave-One-Organ-Out cross-organ matrix** across all 64 primary human donors (`Liver -> Lung/Kidney`, `Lung -> Liver/Kidney`, `Kidney -> Liver/Lung`), proving true organ-invariant universality of the `Parenchyma – SPP1+ SAMac – Myofibroblast` Langevin attractor.
5. **Continuous Fibrosis Progression Hysteresis & Spatial Visium Radial Distance Decay**:
   - Compute continuous disease-progression bifurcation/hysteresis curves smoothly linking early reversible injury, the `t50 = 17.43 h` monotherapy cliff, and locked F4 cirrhosis (`0.9%`), plus the full `12-pair` 2D Bliss synergy matrix.
   - Across all `12` human carcinoma 10x Visium sections (`21,796` spots), compute **continuous radial distance profiles from malignant core boundaries (`0` to `5+` spot pitches / `0–500 µm`)** showing peritumoral `myCAF + SPP1+ TAM` peak rings physically intercepting `CD8+` T-cell infiltration, plus **2D Waddington quasi-potential surfaces $U(S_{\text{myCAF}}, S_{\text{Tumor}})$ with phase-velocity quiver arrows** (eliminating the v1 1D clipped-plateau artifact).

---

### Phase 3 — *Nature Medicine* / *Science* High-Density Figure Craftsmanship Rebuild (`6–7` Panels per Main Figure)
Eliminate uniform `2×2` bar-chart monotony; rebuild all figures via asymmetric multi-column `GridSpec` (`600 DPI` PNG + vector PDF + editable-text SVG, `svg.fonttype = 'none'`):
- **`Fig. 1` (`a–d` Architecture & Translational Paradigm)**: Regenerate clean high-resolution schematic **without `"Non-CNS"`**, depicting **(a)** 12 systemic & parenchymal organ axes + multi-agent clinical trial stoppage deconvolution (`41,396` trials), **(b)** 6-axis primary human single-cell remodeling drift & `DNAcc` biophysical barrier law, **(c)** coupled Itô-Langevin SDE, profile-likelihood audit & Waddington hysteresis saddle, and **(d)** zero-shot cross-organ antifibrotic trial scorecard & Visium-guided pre-treatment anti-PD-1 responder / combination rescue workflow.
- **`Fig. 2` (`6 panels a–f`, Macro-Clinical Deconvolution Atlas)**:
  - **a**, Multi-stage Alluvial/Sankey flow (`41,396` trials $\to$ Target Resolution $\to$ Modality $\to$ Multi-Agent Stop Class $\to$ Definitive Maturation);
  - **b**, 12 systemic/parenchymal organ axes stacked Phase II/III/IV census paired with Open Targets gene resolution rate & diversity;
  - **c**, 20-year (`2005–2025`) longitudinal dynamics of biologic modality share and deconfounded stoppage causes (with 95% confidence ribbons);
  - **d**, Drug modality $\times$ termination mechanism (Safety/Toxicity vs Efficacy/Futility vs Administrative) log-odds forest plot with 95% CIs;
  - **e**, Trial enrollment sample-size distributions across 12 axes / stoppage classes (half-violin raincloud + IQR boxplots);
  - **f**, Multi-agent stoppage auditor validation confusion matrix & Cohen's $\kappa$ concordance against audited ground-truth trials.
- **`Fig. 3` (`7 panels a–g`, Pan-Organ Single-Cell Drift & `DNAcc` Law)**:
  - **a**, **6 Human Single-Cell Axes $\times$ Top 36 Clinical Target Genes Clustered Dot-Heatmap** (encoding $\Delta\tau$ color + compartment-switch ring + trial count size across Liver MASH, Lung IPF, Kidney CKD, NSCLC, Ulcerative Colitis, Sepsis);
  - **b**, 6-step orthogonal ablation trajectory (`AUROC` + 2,000-bootstrap 95% CI whiskers + paired DeLong *P*-values, comparing standard CV vs unseen-target `GroupKFold`);
  - **c**, Full forest plot with **95% CI error bars** comparing Healthy $\tau_0$ vs `DNAcc` across global cohorts, 6 single-cell axes, and fixed-effects multivariable adjustment;
  - **d**, Empirical `Q1 -> Q4` Phase II-to-III/IV maturation curve (`2.56×` enrichment) with 95% binomial CIs;
  - **e**, 10-decile probability calibration curve + Decision Curve Analysis (DCA) net benefit inset;
  - **f**, Exemplar human remodeling quadrant phase diagram ($\tau_0$ vs $\Delta\tau$ with non-crossing radial fan-out target callouts);
  - **g**, Sensitivity forest plot across the `3,118` unreported-stop boundary regimes.
- **`Fig. 4` (`7 panels a–g`, Pan-Organ Fibrosis Zero-Shot Transfer & F4 Hysteresis Rescue)**:
  - **a**, 64-donor primary human 3-organ (`Liver, Lung, Kidney`) pseudobulk state distribution across `Parenchyma / SPP1+ SAMac / Myofibroblast` progression;
  - **b**, Human myeloid `CCR2` 4.2-fold collapse vs `SPP1 / TREM2` de novo emergence (`mean ± s.e.m.` overlaid with individual donor jittered points);
  - **c**, **14-Trial Individual Phase IIb/III Antifibrotic Clinical Scorecard Matrix** (`MAESTRO-NASH`, `ESSENCE`, `NATIVE`, `INPULSIS-1`, `ASCEND`, `FIDELIO-DKD`, `AURORA`, `STELLAR-3/4`, `RAINIER`, `RESOLVE-IT`, `ZEPHYRUS-1`, etc. resolved row-by-row with organ, modality, $\tau_0$, `DNAcc`, SDE probability, and actual clinical outcome!);
  - **d**, Zero-shot model arena (with 95% bootstrap CIs & exact permutation *P*) paired with **`3×3` Leave-One-Organ-Out bidirectional cross-organ transfer AUROC matrix**;
  - **e**, 2D Waddington quasi-potential surface $U(\text{Myofibroblast}, \text{Parenchyma})$ overlaid with **deterministic drift quiver vector field** and saddle-escape trajectory;
  - **f**, Continuous fibrotic progression hysteresis bifurcation curve connecting early reversible states, `t50 = 17.43 h` cliff, and locked F4 monotherapy failure (`0.9%`) vs dual-compartment rescue (`100%`);
  - **g**, **2D Bliss Synergy Matrix Heatmap** across systematic two-target antifibrotic combinations in locked F4 cirrhosis.
- **`Fig. 5` (`7 panels a–g`, 10x Visium `myCAF` Exclusion & Blind Pre-Treatment Anti-PD-1 Responder Prediction)**:
  - **a**, **Dual contrasting primary human carcinoma 10x Visium spatial maps** (`myCAF`-encapsulated immune-excluded tumor vs barrier-low `CD8+` T-infiltrated tumor);
  - **b**, **12-Section continuous spatial radial distance profile** (`myCAF`, `SPP1+ TAM`, and `CD8+ T` density as a function of distance from malignant core boundary) + paired proximity Wilcoxon comparison (`P = 0.0002`);
  - **c**, **`n = 32` Strictly Pre-Treatment Patient-Level Prediction Waterfall Plot** (bars colored by true clinical Responder vs Non-Responder, annotated by cancer cohort BCC vs NSCLC);
  - **d**, Multi-cohort blind responder performance (`BCC n=8`, `NSCLC scRNA n=14`, `NSCLC Pre-Biopsy Bulk n=24`, `Pooled Pre-Treatment n=32`) dot-whisker plot with 95% bootstrap CIs & exact permutation *P*-values;
  - **e**, Longitudinal paired pre- vs post-anti-PD-1 `CD8+` cytotoxic expansion slope/raincloud plot (`GSE123813`, `n = 15`);
  - **f**, **2D Solid-Tumor Waddington Quasi-Potential & Vector-Field Shift** ($U(S_{\text{myCAF}}, S_{\text{Tumor}})$ under anti-PD-1 monotherapy vs `anti-PD-1 + myCAF/TGFBR1` barrier unlocking — smooth 2D potential with zero 1D clipping plateau);
  - **g**, Individual patient trajectory rescue paired slope plot across all `n = 19` anti-PD-1-resistant pre-treatment biopsies (`16.3%` $\to$ `99.2%`).
- **`Extended Data Figs. 1–4` (`4–6` panels each)**:
  - `ED Fig. 1`: PRISMA audit, Open Targets modality/polypharmacology distributions, per-axis target resolution, and `GroupKFold` unseen-target generalization curves;
  - `ED Fig. 2`: 6-axis single-cell donor/cell QC summary, compartment-switch distributions, $\Delta\tau$ distributions, and vital-organ bystander hazard analysis;
  - `ED Fig. 3`: Bilateral profile-likelihood $\Delta\chi^2$ curves (`4` individual parameter profile scans!), bidirectional `3×3` organ calibration residuals, and trial enrichment / Decision Curve Analysis curves;
  - `ED Fig. 4`: All `12` human Visium section spatial barrier scores, FACS-deconfounding sensitivity on `GSE123813`, and platform-invariant MAD scaling calibration across scRNA vs bulk pre-treatment biopsies.

---

### Phase 4 — *Nature Medicine* Manuscript Rewrite, 5-Reviewer Simulated Panel & Dual Gates
1. Rewrite `manuscript/OmniNiche_Flagship_Manuscript_EN.md` and `manuscript/OmniNiche_Flagship_Manuscript_ZH_Preview.md` strictly to ***Nature Medicine* Article** structure and prose standards (`nature-polishing` + `humanizer_academic` + `manuscript-lift-playbook`), live-verifying 100% of old and new PubMed references.
2. Run `figure_clean_aesthetic_gate.py` (`0` errors, `0` warnings) + `read_file` visual inspection on every PNG + `manuscript_convention_and_parity_gate.py` (`exit code 0`).
3. Execute `biomed-review-panel` (5 independent *Nature Medicine* reviewer personas: Chief Editor, Biostatistics/Clinical Trials Reviewer, Single-Cell/Fibrosis Expert, Spatial Oncology/ICB Expert, Devil's Advocate) $\to$ fix every substantive critique in code/figures/text until consensus score reaches *Nature Medicine* external-review/accept tier $\to$ `pytest` 100% green $\to$ `git commit`.
