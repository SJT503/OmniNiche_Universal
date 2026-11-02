# Round-3 Four-Pass Flagship Benchmark Teardown (`Nature Medicine` / `Nature` / `Nature Communications`)

**Target Journal**: *Nature Medicine* (`US` spelling, `Extended Data Fig. 1–5`, `Extended Data Table 1`, `Supplementary Tables 1–6`, unstructured `<= 220`-word Abstract, `Declarations` inside main manuscript)  
**Round-2 Immutable Baseline Archive**: `F:/OmniNiche_Universal/v_round2_archive/` (`9` PNGs, `55` CSVs, Stage II v2.6 manuscript & rebuttal frozen prior to Round-3 edits)  
**Benchmark Corpus Directory**: `F:/OmniNiche_Universal/benchmarks_r3/` (`6` Open-Access Main PDFs, `5` Supplementary/ESM PDFs, `3` Open Peer Review File PDFs totaling `174` pages of verbatim referee critiques & author rebuttals, and `10` rendered 200-DPI figure pages visually inspected by the PI agent)

---

## Downloaded Benchmark Corpus Ledger

1. **Bench 1 — Clinical Trial Genetics & Attrition Anchor (*Nature*, 2024)**: Minikel, E. V., Painter, J. L., Dong, C. C. & Nelson, M. R. Refining the impact of genetic evidence on clinical success. *Nature* **629**, 624–629 (2024). `PMID: 38632401` | `DOI: 10.1038/s41586-024-07316-0`  
   - Local files: `Bench1_Minikel2024_Nature_GeneticEvidenceClinicalTrials.pdf` (`7.58 MB`, 20 pp), `ESM_Supp_Bench1_Nature_2024.pdf` (`3.54 MB`), `Peer_Review_File_Bench1_Nature_2024.pdf` (`378 KB`, 28 pp).
2. **Bench 2 — Clinical AI Immunotherapy Efficacy & Survival Anchor (*Nature Medicine*, 2025)**: Yoo, S.-K. et al. Prediction of checkpoint inhibitor immunotherapy efficacy for cancer using routine blood tests and clinical data (`SCORPIO`). *Nat. Med.* **31**, 869–880 (2025). `PMID: 39762425` | `DOI: 10.1038/s41591-024-03398-5`  
   - Local files: `Bench2_Yoo2025_NatMed_SCORPIO_ImmunotherapyPrediction.pdf` (`3.54 MB`, 25 pp), `ESM_Supp_Bench2_NatMed_2025.pdf` (`29.24 MB`).
3. **Bench 3 — AI Antifibrotic Translational Trial Anchor (*Nature Medicine*, 2025)**: Xu, F. et al. A generative AI-discovered TNIK inhibitor for idiopathic pulmonary fibrosis: a randomized phase 2a trial. *Nat. Med.* **31**, 2602–2610 (2025). `PMID: 40461817` | `DOI: 10.1038/s41591-025-03743-2`  
   - Local files: `Bench3_Xu2025_NatMed_TNIK_PulmonaryFibrosisPhase2a.pdf` (`14.12 MB`, 25 pp).
4. **Bench 3b — Zero-Shot Biomedical Foundation Model Anchor (*Nature Medicine*, 2024)**: Huang, K. et al. A foundation model for clinician-centered drug repurposing (`TxGNN`). *Nat. Med.* **30**, 3601–3613 (2024). `PMID: 39322717` | `DOI: 10.1038/s41591-024-03233-x`  
   - Local files: `Bench3b_Huang2024_NatMed_TxGNN_DrugRepurposing.pdf` (`4.48 MB`, 20 pp), `ESM_Supp_Bench3b_NatMed_2024.pdf` (`5.84 MB`).
5. **Bench 4 — Spatial Tumor–Stroma Boundary & Anti-PD-1 Barrier Anchor (*Nature Communications*, 2024)**: Feng, Y. et al. Spatially organized tumor-stroma boundary determines the efficacy of immunotherapy in colorectal cancer patients. *Nat. Commun.* **15**, 10259 (2024). `PMID: 39592630` | `DOI: 10.1038/s41467-024-54710-3`  
   - Local files: `Bench4_TumorStromaBoundary2024_NatCommun_Immunotherapy.pdf` (`6.30 MB`, 18 pp), `ESM_Supp_Bench4_NatCommun_2024.pdf` (`3.72 MB`), `Peer_Review_File_Bench4_NatCommun_2024.pdf` (`13.75 MB`, 79 pp).
6. **Bench 5 — Spatial Tumor–Immune Heterogeneity & Survival Anchor (*Nature Communications*, 2025)**: Klein, L. et al. Spatial tumor immune heterogeneity facilitates subtype co-existence and therapy response in pancreatic cancer. *Nat. Commun.* **16**, 335 (2025). `PMID: 39762215` | `DOI: 10.1038/s41467-024-55330-7`  
   - Local files: `Bench5_SpatialImmuneHeterogeneity2025_NatCommun_Therapy.pdf` (`10.27 MB`, 19 pp), `ESM_Supp_Bench5_NatCommun_2025.pdf` (`1.96 MB`), `Peer_Review_File_Bench5_NatCommun_2025.pdf` (`16.00 MB`, 67 pp).

---

## Pass 1 — Macro Story Arc, Structural Cadence & *Nature Medicine* House Conventions

1. **Title & Abstract Architecture**:
   - Across `Bench2` (*Nat. Med.* 2025 `SCORPIO`), `Bench3` (*Nat. Med.* 2025 TNIK IPF), and `Bench3b` (*Nat. Med.* 2024 `TxGNN`), titles are declarative or noun-phrase clinical mechanisms (`11–17` words) free of internal pipeline acronyms or negative scope exclusions. Unstructured single-paragraph Abstracts (`185–220` words) follow a strict 4-beat cadence: **(i) Clinical bottleneck in 1 sentence (`<= 28` words)** -> **(ii) Cohort & multi-modal architecture (`n` trials, `n` patients, `n` organs)** -> **(iii) Primary out-of-sample discrimination & orthogonal lift over standard-of-care (`AUROC`, `OR`, `P`)** -> **(iv) External clinical transportability & patient survival (`HR`, `P`)**.
2. **Introduction Funnel (`3–4` paragraphs, `550–720` words)**:
   - `Bench1` (Minikel et al.) and `Bench3b` (`TxGNN`) open directly on clinical attrition cost (`~90%` failure) and immediately pinpoint the structural blind spot of existing paradigms (`TxGNN` p. 3602: *"92% of 17,080 diseases... have no indications... ability of existing machine learning models... drops drastically"*). Neither paper uses defensive scope-exclusion labels; scope is defined positively by biological mechanism.
3. **Results Subheading & Citation House Style**:
   - *Nature Medicine* Results contain `4–6` declarative, mechanism-first subheadings (`<= 14` words each).
   - **Critical House-Style Finding (`Bench3` p. 2604 & `Bench2` p. 871)**: *Nature Medicine* pairs **`Extended Data Fig. 1–8`** (embedded in main online article / SI) with **`Supplementary Table 1–9`** (NO `"S"` prefix in `Supplementary Table 1`, matching `MANUSCRIPT_STATE.json -> supp_prefix: "Extended_Data"`!). Updating OmniNiche's internal `Supplementary Table S1`..`S6` citations to `Supplementary Table 1`..`6` aligns 100% with *Nature Medicine* house style and `manuscript_convention_and_parity_gate.py`.
4. **Discussion & Limitations Positioning**:
   - Exemplars close with `5–6` tight paragraphs (`900–1,200` words): synthesis of multi-scale advance -> mechanistic resolution of clinical paradoxes -> clinical trial design actionable implications -> **bias-directional Limitations paragraph** (stating whether each boundary condition biases effect estimates conservatively toward the null or requires prospective calibration, rather than apologizing for missing assays).

---

## Pass 2 — Visual Figure Autopsy & Chart-Type Grammar (From 7 Visually Inspected 200-DPI PNGs)

Direct visual inspection of `Bench1_Nature_p2_FigPage.png`, `Bench2_NatMed_p3_FigPage.png`, `Bench2_NatMed_p4_FigPage.png`, `Bench3_NatMed_p3_FigPage.png`, `Bench3b_NatMed_p2_FigPage.png`, `Bench4_NatCommun_p3_FigPage.png`, and `Bench5_NatCommun_p5_FigPage.png` established six concrete visual design laws:

1. **Asymmetric Multi-Column GridSpec & Right-Margin Denominator Columns (`Bench1` Fig. 1a–e)**:
   - `Bench1` avoids uniform `2x2` grids, combining stacked horizontal point-interval plots (`a, b`), an L2G threshold ribbon curve (`c`), a trait-count histogram (`e`), and a tall 24-stratum subgroup forest column (`d`) with right-aligned `numerator/denominator` (`189/1,519`, `103/412`) text columns outside the plot frame.
2. **80/20 Neutral-to-Protagonist Chromatic Restraint (`Bench1` Fig. 1, `Bench3b` Fig. 1, `Bench5` Fig. 2)**:
   - Over 80% of ink in `Bench1` and `Bench5` is charcoal (`#182026`), slate grey (`#64707D`), and porcelain fill (`#F3F6F9`), reserving saturated cobalt/teal/coral strictly for the protagonist predictor, confidence ribbon, or high-risk survival stratum.
3. **Cross-Cohort / Cross-Stratum Generalization Matrices (`Bench2` Fig. 2)**:
   - `SCORPIO` dedicates a full visual matrix to showing out-of-sample `AUC` across every individual cohort and clinical trial arm rather than hiding behind a single pooled `AUROC` — motivating a dedicated Leave-One-Cohort-Out (LOCO) & Leave-One-Disease-Axis-Out (LODO) generalization visual in our Extended Data ecosystem.
4. **Spatial Transcriptomics Boundary Proof (`Bench4` Fig. 1b–f)**:
   - Spatial barrier papers pair representative tissue coordinate maps (`d`) with continuous physical distance-to-boundary profiles (`e`) and regression confidence ribbons (`f`), annotated with exact spot denominators (`17,639` boundary spots) and two-sided corrected `P`-values.
5. **Translational Survival + Multivariable Cox Coupling (`Bench5` Fig. 2k,l)**:
   - Survival panels pair Kaplan–Meier curves and aligned `Number at risk` tables directly with multivariable/quartile Cox proportional-hazards forest intervals (`HR`, `95% CI`, `P`-value), proving prognostic independence beyond baseline tumor/clinical covariates.
6. **Zero Internal Code Tokens or Narrative Text Boxes**:
   - Across all 7 inspected pages, zero axis labels contain raw underscores (`snake_case`), zero panels embed multi-word prose sentences inside plot boxes, and panel tags are crisp `9–10 pt` bold lowercase (`a–l`).

---

## Pass 3 — Open Peer Review Attacks -> Author Defenses -> OmniNiche Round-3 Orthogonal Code Action Plan

Reading the `174` pages of downloaded Open Peer Review Files (`Peer_Review_File_Bench1_Nature_2024.pdf`, `Peer_Review_File_Bench4_NatCommun_2024.pdf`, `Peer_Review_File_Bench5_NatCommun_2025.pdf`) surfaced four recurring referee attacks that directly target translational multi-omics & clinical-trial papers:

1. **Attack 1 (From `Bench1` *Nature* 2024 Peer Review File, Referee #1 Comments 4–7 — Therapeutic-Area Heterogeneity, Finite-Sample Overparameterization & Imputation Sensitivity)**:
   - *Referee critique*: Aggregate trial-success enrichment can be confounded by therapeutic-area baseline differences, target-family clustering, or proxy imputation across indications lacking direct data.
   - *OmniNiche Round-3 Code Defense (`Table_S6a_Core6_Regularized_Shrinkage_and_LODO_Axis.csv`)*: In Stage II v2.6, unpenalized 9-feature `M6` (`0.6693`) trailed 3-feature `M3` (`0.6787`) on the smaller Core-6 direct-data subset (`n = 5,662`, `386` failures) due to finite-sample parameter variance (~42 events/parameter), and cross-axis transportability without consensus imputation had not been stress-tested via Leave-One-Disease-Axis-Out (LODO). We will execute an **L2-regularized Ridge/Elastic-Net shrinkage sweep (`C` grid) + strict 6-axis LODO cross-validation** on Core-6 unseen targets (`Group5Fold`) to test whether variance shrinkage eliminates the finite-sample penalty and whether `DNAcc` transfers zero-shot to held-out organ axes without target or axis leakage.
2. **Attack 2 (From `Bench1` Referee #2 & `Bench2` *Nat. Med.* `SCORPIO` Calibration Standards — Empirical Drug Potency vs Uniform Reference Potency & Cross-Organ Probability Calibration)**:
   - *Referee critique*: Uniform drug assumptions across heterogeneous molecules and miscalibrated raw risk scores (`0.05–0.23` compressed probabilities vs a `0.50` decision threshold) undermine clinical decision utility.
   - *OmniNiche Round-3 Code Defense (`Table_S6b_Antifibrotic24_ChEMBL_Potency_and_LOOO_Calibration.csv`)*: Replace the uniform reference potency (`pIC50 = 7.5`) across all `24` Phase IIb/III antifibrotic trials (`Liver`, `Lung`, `Kidney`) with **curated ChEMBL v33 / clinical human target `pXC50` bioactivity + a 2D potency-by-Hill-slope perturbation grid**, AND implement **Leave-One-Organ-Out (LOOO) Platt & isotonic monotonic probability recalibration** (fitting calibration strictly on the other 2 organs and predicting the held-out organ) to test whether out-of-organ recalibration corrects probability compression at the `0.50` threshold while preserving `AUROC = 0.8714` ranking.
3. **Attack 3 (From `Bench4` *Nat. Commun.* 2024 Peer Review File, Reviewer #1 — Spatial Transcriptomics "Barrier" vs Compositional Spot Admixture & Spatial Autocorrelation)**:
   - *Referee critique*: Reviewer #1 challenged whether co-localization of `Mac_SPP1` and `myCAF` ECM genes at the tumor–stroma boundary reflects a true physical exclusion architecture rather than 55-µm spot cellular admixture, sequencing depth gradients, or spatial autocorrelation.
   - *OmniNiche Round-3 Code Defense (`Table_S6c_Visium12_Admixture_Residualization_and_Toroidal_Null.csv`)*: Across all `20,880` quality-filtered 10x Visium spots in the `12` primary human carcinoma sections (`GSE203612` + `GSE171351`), compute **OLS admixture- and library-depth-residualized Stromal Exclusion Index (`SEI_resid`)** (regressing out total spot signature density and malignant/stromal admixture entropy) AND run **`2,000` 2D toroidal cyclic coordinate-shift permutations** per section (preserving exact within-lineage spatial autocorrelation variograms while breaking cross-lineage boundary alignment!).
4. **Attack 4 (From `Bench5` *Nat. Commun.* 2025 Peer Review File, Reviewers #1–#3 & `Bench2` `SCORPIO` — Leave-One-Cohort-Out Transportability & Multivariable Cox Confounding + E-Value)**:
   - *Referee critique*: Pooled immunotherapy survival associations (`HR = 1.95` across `n = 74` melanoma patients) can be driven by cohort batch effects, baseline tumor burden, `PD-L1` expression, or T-cell inflammation (`18-gene GEP`).
   - *OmniNiche Round-3 Code Defense (`Table_S6d_Oncology107_LOCO_and_Survival74_Multivariable_Cox.csv`)*: Execute **Leave-One-Cohort-Out (LOCO) validation across the 4 strictly pre-treatment anti-PD-1 cohorts (`n = 107`)** AND fit **multivariable Cox proportional-hazards models on the `n = 74` monotherapy OS cohort** adjusting simultaneously for cohort (`GSE91061` vs `GSE78220`), baseline `Tumor` proliferative signature, `PD-L1` (`CD274`), and `18-gene T-cell inflamed GEP`, accompanied by **VanderWeele E-value unmeasured-confounding tipping-point quantification**.

---

## Pass 4 — Micro-Sentence Rhythm, Three-Reader Cognitive Ergonomics & 10-Dimension Gap Action Table

Auditing our Round-2 (`v2.6`) prose against `Bench1`, `Bench2`, `Bench3`, and `three-reader-ergonomics.md`:
- **Strengths to preserve**: `200`-word single-paragraph Abstract, zero `"Non-CNS"` / `"Firewall-2"` jargon, 100% `<= 75`-word sentences in `.md` and `.docx`, `61` live-verified PubMed citations in strict `[1]..[61]` order.
- **Ergonomic upgrades for Round 3**:
  1. Convert all `Supplementary Table S1`..`S6` citations in EN/ZH prose and SI into *Nature Medicine* house-standard `Supplementary Table 1`..`6` (`non_S` / `Extended_Data` contract).
  2. Ensure every Results paragraph opens with a `<= 25`-word conclusion-first topic sentence and introduces specialized quantitative metrics (`DNAcc`, `SEI`, toroidal shift null, E-value) with a `6–12` word plain-scientific physical anchor on first mention.
  3. Upgrade the Discussion Limitations paragraph so every boundary condition states its **exact direction of statistical bias** (e.g., clinical registry `whyStopped` under-reporting biases unadjusted modality contrasts toward the null; 55-µm Visium spot diameter averages sub-spot contacts and conservatively attenuates measured `SEI`).
  4. Assemble the complete **5-Piece *Nature Medicine* Submission Package** (`OmniNiche_Flagship_Manuscript_EN.docx` with all 5 Main Figures embedded inline, `Supplementary_Information_EN.docx` with verbatim title mirror + all 5 Extended Data Figures embedded inline + Supplementary Notes 1–3, multi-sheet `Supplementary_Tables.xlsx` with `README_Index` + `Supplementary Tables 1–6`, and `Cover_Letter_NatMed.docx` with verbatim title mirror).

### 10-Dimension Round-2 vs Benchmark Gap & Round-3 Execution Matrix

| Dim # | Audit Dimension | *Nat. Med.* / *Nature* Benchmark Standard | OmniNiche Round-2 (`v2.6`) Baseline | Concrete Gap Identified | Round-3 Code / Figure / Prose Action |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | **House Citation & Supp Naming** | `Extended Data Fig. 1–8` + `Supplementary Table 1–9` (zero `"S"` prefix; `Bench3` p. 2604) | Mixed `Extended Data Fig. 1–4` with `Supplementary Table S1–S5` | `"S"` prefix in `Supplementary Table S1` violates *Nat. Med.* house style & `MANUSCRIPT_STATE.json` | Normalize 100% of EN/ZH/SI citations to `Supplementary Table 1–6` and build 7-sheet `Supplementary_Tables.xlsx` |
| **2** | **Core-6 Finite-Sample Variance & Cross-Axis LODO** | Penalized shrinkage + out-of-domain generalization (`Bench1` & `Bench2`) | Unpenalized 9-param `M6` (`0.6693`) vs 3-param `M3` (`0.6787`) on Core-6 (`n = 5,662`) | No regularized shrinkage sweep or Leave-One-Disease-Axis-Out (LODO) test on Core-6 | Compute L2 shrinkage sweep + 6-axis LODO in `Table_S6a` & plot in `Extended Data Fig. 5a` |
| **3** | **Antifibrotic Drug Potency & Probability Calibration** | Empirical target bioactivity + out-of-sample probability calibration (`Bench1`, `Bench2`) | Uniform `pIC50 = 7.5`; compressed Lung/Kidney probabilities (`0.0535–0.2330`) at `0.50` cutoff | Reviewers question uniform `pIC50` assumption and threshold calibration across organs | Compute ChEMBL v33 `pXC50` sensitivity + Leave-One-Organ-Out (LOOO) calibration in `Table_S6b` & `Extended Data Fig. 5b` |
| **4** | **Visium Spatial Barrier vs Admixture & Autocorrelation** | Spatial null preserving autocorrelation + admixture control (`Bench4` PRF Reviewer #1) | Spot-level `k = 6` Mann–Whitney `P`; no admixture residualization or toroidal shift null | Need proof that `SEI > 0` across `20,880` spots is true geometry, not spot density/autocorrelation | Run OLS admixture residualization + `2,000` 2D toroidal cyclic shifts/section in `Table_S6c` & `Extended Data Fig. 5c` |
| **5** | **ICB Multi-Cohort LOCO & Multivariable Cox + E-Value** | Cross-cohort holdout matrix (`Bench2` Fig. 2) + multivariable Cox forest (`Bench5` Fig. 2k,l) | Stratified log-rank `HR = 1.95` (`n = 74`); no simultaneous `PD-L1` + `GEP18` + `Tumor` Cox + E-value | Referees demand multivariable Cox independence over `PD-L1`/`GEP18` + E-value tipping point | Compute 4-cohort LOCO (`n = 107`) + multivariable Cox & E-value (`n = 74`) in `Table_S6d` & `Extended Data Fig. 5d` |
| **6** | **Main + Extended Data Figure Ecosystem Split** | Dedicated Main Figs (`5–8`) + Extended Data Figs (`4–8`) in clean submission hierarchy | All 9 PNGs stored together in `outputs/stage2_figures/` (`4` ED Figs) | Gate 2 requires separate `stage2_main_figures/` (`5` Figs) & `stage2_supp_figures/` (`5` ED Figs) | Create `visual_tokens.py` + `stage2_supp_figures.py` (`ED Figs 1–5`), preserving `outputs/stage2_figures/` for 47 existing tests |
| **7** | **Inline Figure Embedding in Submission `.docx`** | Initial peer-review `.docx` embeds high-res figures inline alongside legends (`r3_gate_checker` Gate 3) | `v2.6` `.docx` had `[INSERT FIGURE X ABOUT HERE]` text placeholders with `0` embedded PNGs | Reviewers must scroll without visual context if `word/media/` has `0` images | Embed all `5` 600-DPI Main Fig PNGs inside `OmniNiche_Flagship_Manuscript_EN.docx` & `5` ED Fig PNGs inside `Supplementary_Information_EN.docx` |
| **8** | **Standalone `Supplementary_Information_EN.docx`** | Complete SI Word document with verbatim title mirror, ED Figs + legends, and Supp Notes | Only `Supplementary_Tables_Master_Index.md` existed; no standalone `Supplementary_Information_EN.docx` | Missing piece #2 of the mandatory 5-piece submission package | Author and compile `Supplementary_Information_EN.md` + `Supplementary_Information_EN.docx` (`zero garble`, verbatim title) |
| **9** | **Multi-Sheet `Supplementary_Tables.xlsx`** | Consolidated multi-tab Excel workbook (`README_Index` + `Table 1..6`) | `55` individual CSV files in `outputs/stage2_tables/` without unified `.xlsx` workbook | Missing piece #3 of the mandatory 5-piece submission package | Compile formatted `Supplementary_Tables.xlsx` (`README_Index` + `Supplementary Table 1`..`6` + key sub-tables) |
| **10** | **Editorial `Cover_Letter_NatMed.docx` & Bias-Directional Limitations** | Non-abstract-repeating *Nat. Med.* Cover Letter with exact title mirror + directional Limitations | No Cover Letter `.docx` yet; Limitations qualitative in 2 clauses | Missing piece #4 of the 5-piece submission package + bias-direction sharpness | Author `Cover_Letter_NatMed.md`/`.docx` with verbatim title mirror and sharpen Limitations with explicit bias directions |
