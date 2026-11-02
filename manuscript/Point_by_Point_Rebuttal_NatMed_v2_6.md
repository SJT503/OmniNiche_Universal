# Point-by-Point Response to Reviewers (*Nature Medicine* Revision v2.6)

**Manuscript Title:** Dynamic single-cell niche remodeling and biophysical steric exclusion govern clinical trial attrition across human organs  
**Revision Version:** Stage II v2.6 (Response to Simulated *Nature Medicine* 5-Reviewer Panel v2.5, Aggregate Score 86.5/100)  
**Revised Deliverables:** `OmniNiche_Flagship_Manuscript_EN.md` (`335` lines), `OmniNiche_Flagship_Manuscript_ZH_Preview.md` (`335` lines, 1:1 parity), `OmniNiche_Flagship_Manuscript_EN.docx` (`67,031` bytes), `Supplementary_Tables_Master_Index.md` (`60` deterministic empirical files, `61` PubMed-verified references), and `9` re-rendered 183-mm vector/600-DPI figures (`Fig. 1`–`Fig. 5`, `Extended Data Fig. 1`–`4`).

---

## 1. Editor Cover Letter

Dear Chief Editor and Handling Editorial Team,

Thank you for the constructive **Major Revision** evaluation of our manuscript and for synthesizing the five reviewers' cross-disciplinary critiques into five consensus priorities (`Cc.1`–`Cc.5`) and ten actionable revision tasks (`T1`–`T10`). In Stage II v2.6, we executed four foundational empirical and methodological upgrades across the Python codebase, deterministic output tables, publication figures, and manuscript text, paired with explicit scope and claim moderations:

1. **Cross-Platform $\log_2(\mathrm{CPM}+1)$ Harmonization of Interstitial ECM Barrier Density (`T1` / `Cc.2`)**: We eliminated cross-platform scale distortion between Human Protein Atlas v23 normalized consensus transcripts per million (`nTPM`) and primary single-cell UMI pseudobulk profiles by harmonizing both inputs onto a common $\log_2(\mathrm{CPM}+1)$ stromal scale across the 10-gene collagen-myofibroblast signature prior to min-max normalization (`Table_S3m_Empirical_Axis_ECM_Density.csv`). On this unified biophysical scale, desmoplastic solid oncology (`0.3500`) and pulmonary fibrosis (`0.3265`) rank at the top while bloodstream sepsis (`0.0700`) anchors the minimum, lifting out-of-fold `DNAcc` scientific discrimination to **AUROC = 0.6710** ($P = 1.11 \times 10^{-21}$ versus static $\tau_0 = 0.5387$; unseen-target `Group5Fold` **0.6541**, $P = 3.91 \times 10^{-25}$) and `Human Genetics + DNAcc` synergy to **AUROC = 0.7023** ($P = 1.58 \times 10^{-10}$ versus genetics alone `0.6362`; `Group5Fold` **0.6844**; maturation **0.7210**, $P = 1.72 \times 10^{-111}$).
2. **Formal Pairwise `M6 vs M3` Step-Attribution DeLong Audit & Six-Axis Consensus-Drift Disclosure (`T2`, `T6` / `Cc.1`, `Cc.5`)**: We added exact two-sided paired DeLong tests comparing full `DNAcc` (`M6`) against early diseased-niche specificity (`M3`) and sequential steps (`M_{k}` versus `M_{k-1}`) across all 48 rows of `Table_S3f_SixStep_Ablation_and_GroupKFold.csv`, documented the exact six-axis consensus-drift assignment rule (`51.3%` of scientific trials) in *Methods*, and moderated our framing across *Abstract*, *Results*, and *Discussion* to state explicitly that `M3` ($\tau_{\mathrm{early}}$) captures `92.7%` of the global scientific lift while progression drift and steric compatibility provide statistically significant incremental refinement across single-target monotherapies ($+0.0394$, $P = 0.0353$) and the 23,761-trial maturation cohort ($+0.0137$, $P = 8.93 \times 10^{-13}$; `Group5Fold` $P = 1.70 \times 10^{-4}$), alongside transparent reporting of the Core-6 `Group5Fold` finite-sample variance trade-off (`M3` `0.6575` versus `M6` `0.6415`, $P = 0.0188$ across 386 failure events).
3. **Dual Rank vs Absolute Probability Calibration Audit in 24 Antifibrotic Trials (`T4`, `T7` / `Cc.3`)**: In *Results §4*, *Table 2*, *Fig. 4c*, and *Discussion (Limitations)*, we now contrast the empirical Youden-optimal rank cutoff ($p \ge 0.0535$: **10/10 [100.0%] sensitivity**, **10/14 [71.4%] specificity**) against the raw physical attractor-basin threshold ($p \ge 0.50$: **5/10 [50.0%] sensitivity** [4/5 hepatic, 1/3 pulmonary, 0/2 renal], **13/14 [92.9%] specificity**), explaining mechanistically why liver-compiled three-compartment Hill kinetics compress pulmonary and renal hemodynamic/pleiotropic winners into `0.0535–0.2330` while verifying renal tubular `SLC5A2` specificity ($\tau_0 = 0.8959$ in `GSE271188` epithelium) and citing the negative Phase IIb semaglutide RCT in compensated F4 cirrhosis (Loomba et al., *Lancet Gastroenterol Hepatol* 2023 `[61]`, `PMID: 36934740`).
4. **Symmetric `No-Visium`, `Uninformative Linear Mean`, and `CD8:myCAF` Benchmarks Across Solid-Tumor Anti-PD-1 Cohorts (`T3`, `T5`, `T8`–`T10` / `Cc.4`, `DA-U1`, `DA-U2`)**: We separated strictly pre-treatment pure monotherapy biopsies ($n = 83$), pre-treatment chemo-anti-PD-1 biopsies ($n = 24$; $n = 107$ total pre-treatment), and post-neoadjuvant surgical resection specimens ($n = 14$), restored `CD8:myCAF` as the sixth plotted series in `Fig. 5d`, added symmetric `No-Visium` and `Uninformative Linear Mean` comparisons across every oncology row of *Table 2*, moderated our bulk RNA-seq interpretation to reflect spatial-gradient erasure under tissue homogenization (`0.6367` versus `Linear Mean` `0.6442` on $n = 107$, whereas intact single-cell TME reaches `0.8750` versus `Linear Mean` `0.8250` in NSCLC and `0.7500` in BCC), reported exact non-significant *P*-values for the exploratory $n = 15$ paired longitudinal scRNA-seq trajectory (Mann-Whitney $P = 0.779$), and clarified deterministic in silico operator rescue statistics (`Op. P = 1.8e-13`) and `20,880` quality-filtered Visium spots (of `21,796` profiled).

All 47 automated unit tests (`pytest`), the 10-gate figure aesthetic audit, the 3-Layer Zero-Garble `.docx` audit, and the bilingual manuscript convention and parity gate pass without warning. Below we provide a cross-review consensus summary, atomic point-by-point responses with verified line numbers and tracked excerpts, and the completed Revision Task Table (`T1`–`T10`).

Sincerely,  
The OmniNiche Universal Authors

---

## 2. Cross-Review Consensus Summary (`Cc.1`–`Cc.5` & Red-Team Vectors `DA-U1`–`DA-U2`)

| Consensus / Vector ID | Raising Reviewers | Core Methodological or Translational Concern | Resolution in Revised v2.6 Manuscript, Code & Figures | Status |
| :--- | :--- | :--- | :--- | :---: |
| **`Cc.1`** | Role 1 (`#1`), Role 2 (`#1`), Role 5 (`#1`) | `M3` ($\tau_{\mathrm{early}}$) drives >90% of `DNAcc` scientific lift (`0.6614 -> 0.6710`), and `M3` (`0.6575`) exceeds `M6` (`0.6415`) under Core-6 `Group5Fold`. | Added exact pairwise `M6 vs M3` and sequential DeLong tests to `Table_S3f`; reported full step-attribution and Core-6 finite-sample variance (`386` failures, `3` vs `9` features) in *Abstract* (`Line 7`), *Results §3* (`Lines 49–51`), *Table 1* (`Lines 197–198, 202, 204, 206`), and *Discussion* (`Lines 105, 113`). | **DONE** |
| **`Cc.2`** | Role 2 (`#2`), Role 4 (`#1`) | Cross-platform `nTPM` vs scRNA UMI scale mismatch in `Table_S3m` depressed `Pulmonary_IPF_ARDS` (`0.1235`) below `Ophthalmic_Retinal` (`0.2286`). | Harmonized both HPA v23 stromal types and primary scRNA stromal pseudobulks onto common $\log_2(\mathrm{CPM}+1)$ before min-max scaling (`Solid_Oncology = 0.3500`, `Pulmonary_IPF_ARDS = 0.3265`, `Ophthalmic_Retinal = 0.0810`, `Infectious_Sepsis = 0.0700`); regenerated all downstream tables and `Fig. 1`/`Fig. 3`/`ED Figs. 1–3` (`Lines 41, 156`). | **DONE** |
| **`Cc.3`** | Role 1 (`#2`), Role 3 (`#1`), Role 5 (`#3`) | High 24-trial antifibrotic rank `AUROC = 0.8714` coexists with probability compression (`0.0535–0.2330`, `5/10` winners $\ge 0.50$) on pulmonary/renal drugs. | Added dual Youden-optimal ($p \ge 0.0535$: `10/10` sens., `10/14` spec.) vs physical `0.50` (`5/10` sens., `13/14` spec.) cutoffs to `Fig. 4c`, *Results §4* (`Lines 75–77`), *Table 2* (`Line 219`), *Methods* (`Line 162`), and *Limitations* (`Line 115`), plus verified tubular `SLC5A2` ($\tau_0 = 0.8959$) and cited Loomba et al. 2023 `[61]`. | **DONE** |
| **`Cc.4`** | Role 1 (`#3`), Role 3 (`#2`), Role 5 (`#2`) | `Uninformative Linear Mean` (`0.6442` on $n = 107$) matches bulk Langevin (`0.6367`); `No-Visium` (`1.0000`) exceeds full model (`0.8750`) on BCC; post-neoadjuvant $n = 14$ timing. | Separated $n = 83$ pre-treatment monotherapy, $n = 24$ pre-treatment chemo-anti-PD-1 ($n = 107$ total), and $n = 14$ post-neoadjuvant resections; tabulated symmetric `No-Visium`, `Linear Mean`, and `CD8:myCAF` across all *Table 2* rows (`Lines 224–229`) and `Fig. 5d`; explained bulk lysis gradient erasure and visceral-to-cutaneous Visium domain shift (`Lines 91–95, 115`). | **DONE** |
| **`Cc.5`** | Role 2 (`#3`) | `51.3%` of scientific trials span 6 registry axes lacking stage-resolved single-cell progression datasets; consensus imputation rule omitted from *Methods*. | Documented exact inheritance (`Dermatologic_Immune` $\leftarrow$ `Autoimmune_IBD_Rheum`) and 5-axis gene-level `consensus_drift` averaging rule in *Methods §4* (`Line 156`) and *Discussion (Limitations)* (`Line 113`). | **DONE** |
| **`DA-U1`** | Role 5 (`#4`) | `P = 1.78e-13` on $n = 70$ non-responders (`Fig. 5g`) tests deterministic mathematical operator displacement rather than observed clinical combination response. | Explicitly qualified `P = 1.78e-13` as deterministic in silico operator displacement nominating prospective combination hypotheses in *Results §6* (`Line 97`), *Limitations* (`Line 115`), *Fig. 5g* panel (`Op. P = 1.8e-13`), and *Fig. 5g* legend (`Line 261`). | **DONE** |
| **`DA-U2`** | Role 5 (`#5`) | Uniform reference potency $\mathrm{pIC}_{50} = 7.5$ (`31.6 nM`) and extrahepatic reference leverage $\eta = 0.5400$ (`0.7637` dual) in 24 antifibrotic trials. | Disclosed exact reference leverage values ($\eta = 0.5400$; $0.7637$ for `GLP1R;GIPR`) alongside $\mathrm{pIC}_{50} = 7.5$ in *Methods §5* (`Line 162`), clarifying zero-PK-tuning isolation of compartment topology and steric permeability. | **DONE** |

---

## 3. Point-by-Point Response to Reviewers

### 3.1 Role 1 — Senior Translational / Clinical-AI Editor (Score: 87/100)

#### Comment R1.1a [Major — Consensus `Cc.1`: Formal step-attribution of `M3` ($\tau_{\mathrm{early}}$) vs `M4`–`M6` progression & barrier terms]
> *"Report formal pairwise DeLong tests comparing `M6` directly against `M3` (`0.6682` vs `0.6614` on scientific; `0.6745` vs `0.6607` on maturation; `0.6385` vs `0.6575` on Core-6 `Group5Fold`)."*

**Strategy:** `add_statistical_test` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added exact two-sided paired nonparametric DeLong tests comparing every model step against `M3` (`delta_auroc_vs_m3`, `delong_p_vs_m3`) and against the immediately preceding step (`delta_auroc_vs_prev`, `delong_p_vs_prev`) across all 48 cohort-by-CV strata in `outputs/stage2_tables/Table_S3f_SixStep_Ablation_and_GroupKFold.csv` (`omniniche/stage2/v2_statistical_fortification.py`), and incorporated these statistics directly into *Results §3* (`Lines 49–51`) and *Table 1* (`Lines 197–198, 202, 204, 206`). Following cross-platform $\log_2(\mathrm{CPM}+1)$ harmonization of $D_{\mathrm{ECM}}$ (`T1`), `M3` ($\tau_{\mathrm{early}}$, stratified AUROC = `0.6614`) captures `92.7%` of the global `+0.1323` scientific lift over healthy $\tau_0$ (`0.5387`, $P = 7.38 \times 10^{-18}$), while `M6` (`0.6710`) provides incremental gains of `+0.0096` ($P = 0.205$) overall, `+0.0151` ($P = 0.0563$) under unseen-target `Group5Fold` (`0.6541` vs `0.6389`), `+0.0394` ($P = 0.0353$) within strictly single-target monotherapies ($n = 6,115$), `+0.0137` ($P = 8.93 \times 10^{-13}$) across the 23,761-trial Definitive Maturation Cohort (`0.6744` vs `0.6607`), and `+0.0076` ($P = 1.70 \times 10^{-4}$) under Maturation `Group5Fold` (`0.6509` vs `0.6432`; single-target maturation `+0.0295`, $P = 1.19 \times 10^{-12}$). Conversely, under Core-6 unseen-target `Group5Fold` ($n = 5,662$, `386` failures), three-feature `M3` (`0.6575`) exceeds nine-feature `M6` (`0.6415`, $\Delta = -0.0160$, $P = 0.0188$), which we report transparently alongside its finite-sample variance explanation.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §3, Lines 49–51):**
> "Pairwise step attribution showed that **replacing healthy $\tau_0$ with diseased-lesion specificity (`M3`) captured 92.7% of the global scientific gain ($+0.1227$, $P = 7.38 \times 10^{-18}$), while progression drift and steric exclusion (`M6` versus `M3`) added finer incremental refinement ($\\Delta = +0.0096$, $P = 0.205$ overall; $+0.0151$, $P = 0.0563$ under unseen-target `Group5Fold`; $+0.0394$, $P = 0.0353$ within single-target monotherapies).**
>
> Across the higher-powered 23,761-trial Definitive Maturation Cohort, out-of-fold discrimination rose monotonically at every ablation step (`0.5412 -> 0.5777 -> 0.6607 -> 0.6623 -> 0.6716 -> 0.6744`; $\Delta\mathrm{AUROC} = +0.1332$ versus `M1`, DeLong $P = 4.12 \times 10^{-180}$; Brier error reduced from 0.2477 to 0.2273). **Notably, `M6` significantly exceeded `M3` under both stratified cross-validation ($\Delta = +0.0137$, $P = 8.93 \times 10^{-13}$) and unseen-target `Group5Fold` ($\Delta = +0.0076$, $P = 1.70 \times 10^{-4}$; single-target maturation $+0.0295$, $P = 1.19 \times 10^{-12}$)... Under Core-6 `Group5Fold`, the three-feature `M3` model (0.6575) exceeded nine-feature `M6` ($\Delta = -0.0160$, $P = 0.0188$), reflecting higher parameter-estimation variance when fitting six progression-barrier terms across 386 failure events without cross-axis consensus regularization (Extended Data Fig. 1e and Supplementary Table S3f).**"

---

#### Comment R1.1b [Major — Consensus `Cc.1`: Moderate primary claim in Abstract & Discussion to reflect $\tau_{\mathrm{early}}$ dominance]
> *"Moderate the narrative in Abstract, Results §3, and Discussion §1 to frame the finding accurately: diseased-tissue target specificity ($\tau_{\mathrm{early}}$) is the primary biological driver overcoming healthy-atlas failure, whereas progression drift ($\Delta\tau$) and steric exclusion ($B_{\mathrm{acc}}$) provide a secondary refinement that reaches statistical significance across the larger 23,761-trial maturation cohort and in late-stage macromolecular trials."*

**Strategy:** `moderate_claim` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We revised the *Abstract* (`Line 7`), *Discussion §1* (`Line 105`), and *Discussion (Limitations)* (`Line 113`) to state explicitly that clinical translation is anchored primarily by diseased-niche target specificity ($\tau_{\mathrm{early}} = 0.6614$) and refined across clinical maturation and single-target regimens by progression drift and macromolecular steric exclusion ($P = 8.93 \times 10^{-13}$).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Abstract Line 7 & Discussion Line 105):**
> *Abstract (`Line 7`):* "Dynamic Niche Accessibility (`DNAcc`) lifts scientific discrimination over healthy specificity from AUROC 0.5387 to 0.6710 ($P = 1.11 \times 10^{-21}$), **anchored by early diseased-niche specificity (0.6614) with progression-barrier maturation refinement ($P = 8.93 \times 10^{-13}$)**, and synergizes with human genetics (AUROC 0.7023 versus 0.6362, $P = 1.58 \times 10^{-10}$; adjusted odds ratio = 1.42)."  
> *Discussion (`Line 105`):* "...we show that clinical translation across peripheral solid and mucosal organs **is gated primarily by diseased-niche target specificity ($\tau_{\mathrm{early}}$) and refined across clinical maturation by stage progression drift and macromolecular steric exclusion.**"

---

#### Comment R1.2a [Major — Consensus `Cc.3`: Dual absolute `0.50` vs empirical Youden-optimal threshold metrics in 24 antifibrotic trials]
> *"Report sensitivity, specificity, PPV, and NPV at both the physical `0.50` basin threshold and the Youden-optimal threshold in Results §4 and Table 2."*

**Strategy:** `recalibrate_threshold` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We computed exact operating characteristics at both thresholds across the 24 SHA-256-locked Phase IIb/III antifibrotic trials (`Table_S4d_Clinical_24Trials_ZeroShot_Predictions.csv`) and added them to *Results §4* (`Lines 75–77`), *Table 2* (`Line 219`), *Methods §5* (`Line 162`), and *Fig. 4c* (`Line 259`), plotting both vertical reference lines (`Youden p >= 0.054` and `Basin p >= 0.50`) directly on `Fig. 4c`:
- **Raw physical `0.50` attractor-basin threshold**: Sensitivity = `50.0%` (`5/10` winners: `4/5` hepatic, `1/3` pulmonary, `0/2` renal), Specificity = `92.9%` (`13/14` failures), PPV = `83.3%` (`5/6`), NPV = `72.2%` (`13/18`).
- **Empirical Youden-optimal rank threshold ($p \ge 0.0535$, $J = 0.7143$)**: Sensitivity = `100.0%` (`10/10` winners), Specificity = `71.4%` (`10/14` failures), PPV = `71.4%` (`10/14`), NPV = `100.0%` (`10/10`), because sterically excluded antibodies (`Simtuzumab`, `Pamrevlumab`, `Lebrikizumab`: `0.0366–0.0442`) and lineage-eroded inhibitors (`Cenicriviroc`, `Selonsertib`: `0.0409–0.0423`) cluster below `0.048`.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §4, Line 77):**
> "**Consequently, the raw physical 0.50 basin cutoff yielded 50.0% sensitivity (5/10 winners: 4/5 hepatic, 1/3 pulmonary, 0/2 renal) and 92.9% specificity (13/14 failures, with permeable autotaxin inhibitor ziritaxestat [`ISABELA 1`, 0.9667] as the sole high-scoring off-niche toxicity failure). By contrast, the empirical Youden-optimal rank threshold ($p \ge 0.0535$, $J = 0.7143$) captured 10 of 10 winners (100.0% sensitivity) at 71.4% specificity (10/14 failures), because sterically excluded antibodies (`Simtuzumab`, `Pamrevlumab`, `Lebrikizumab`, 0.0366-0.0442) and lineage-eroded `Cenicriviroc` or `Selonsertib` (0.0409-0.0423) clustered below 0.048.**"

---

#### Comment R1.2b [Major — Consensus `Cc.3`: Explain probability compression of pulmonary/renal winners under liver-compiled Hill parameters]
> *"Explicitly discuss why liver-compiled three-compartment Hill parameters compress absolute remission probabilities for non-hepatic hemodynamic/pleiotropic mechanisms (`pirfenidone`, `nerandomilast`, `dapagliflozin`, `finerenone`, `efruxifermin`) into `0.05–0.23`, making zero-shot scores an ordinal ranking tool rather than an out-of-the-box calibrated probability across organs."*

**Strategy:** `clarify_boundary_condition` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added a dedicated paragraph in *Results §4* (`Lines 75–77`), updated the *Abstract* (`Line 7`: *"Zero-label Langevin transfer ranks 24 Phase IIb/III antifibrotic trials"*) and *Introduction* (`Line 19`), and expanded *Discussion (Limitations)* (`Line 115`) to state explicitly that three-compartment hepatic Hill kinetics parameterize local parenchymal-SAMac-myofibroblast feedback and therefore compress absolute probabilities (`0.0535–0.2330`) for pulmonary and renal winners acting via pleiotropic, cAMP (`PDE4B` `nerandomilast`), systemic endocrine (`FGF21` `efruxifermin`), or tubuloglomerular/mineralocorticoid hemodynamic pathways (`dapagliflozin`, `finerenone`), rendering zero-label cross-organ scores a continuous rank prioritizer rather than a pre-calibrated absolute probability without organ-specific threshold calibration.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Discussion Limitations, Line 115):**
> "**Second, in zero-label cross-organ antifibrotic transfer, high continuous rank discrimination (`AUROC = 0.8714`; Youden cutoff $p \ge 0.0535$: 10/10 sensitivity, 10/14 specificity) coexists with compressed absolute probabilities (0.0535-0.2330, 50% sensitivity at the raw 0.50 basin threshold) for pulmonary and renal drugs acting through hemodynamic, endocrine, or cAMP axes weakly parameterized by three-compartment hepatic Hill kinetics (`pirfenidone`, `efruxifermin`, `dapagliflozin`, `nerandomilast`, and `finerenone`).**"

---

#### Comment R1.3a [Major — Consensus `Cc.4`: Surface `Uninformative Linear Mean` & `No-Visium` across all oncology cohorts and moderate bulk claims]
> *"Surface the `Uninformative_Linear_Combination` and `No-Visium` ablation numbers across all cohorts in Results §6 and Table 2, explicitly framing nonlinear spatial-Langevin dynamics as providing clear gains in single-cell TME profiling (`0.8750` vs `0.7500` linear in BCC; `0.8750` vs `0.8250` linear and `0.8500` No-Visium in NSCLC scRNA) while converging toward linear-combination parity in homogenized bulk biopsies where spatial gradients are lost."*

**Strategy:** `add_baseline_comparison` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We updated every oncology row of *Table 2* (`Lines 224–229`) and restructured *Results §6* into two contrasting paragraphs (`Lines 93–95`) plus *Discussion (Limitations)* (`Line 115`). We now report `OmniNiche (No-Visium)` and `Uninformative Linear Mean` side-by-side across all single-cell, individual bulk, and pooled cohorts: in intact single-cell TME cohorts, OmniNiche (`0.8750`) outperforms `Linear Mean` in both NSCLC (`0.8250`; `No-Visium` `0.8500`) and BCC (`0.7500`; with `No-Visium` `1.0000` reflecting visceral-to-cutaneous domain shift), whereas in homogenized bulk biopsies tissue lysis erases concentric `myCAF` distance gradients and compresses OmniNiche (`0.6016` on $n = 83$; `0.6367`, $P = 0.0187$ on $n = 107$) to parity with `Uninformative Linear Mean` (`0.6065` on $n = 83$; `0.6442`, $P = 0.0155$ on $n = 107$) while consistently exceeding `No-Visium` (`+0.0201` and `+0.0157`) and published clinical biomarkers (`PD-L1`, `18-gene GEP`, `CAF/TIDE`, and `CD8:myCAF`).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §6, Line 95):**
> "**In homogenized pre-treatment bulk RNA-seq biopsies, tissue lysis erased concentric `myCAF` distance gradients and compressed nonlinear Langevin trajectories toward linear-combination parity... Across Pooled All Strictly Pre-Treatment Biopsies ($n = 107$, 37 responders), OmniNiche reached AUROC = 0.6367 (95% CI 0.5253-0.7425, $P = 0.0187$), exceeding `No-Visium` (0.6210) and `CD8:myCAF` (0.6263) while matching the `Uninformative Linear Mean` (`0.6065` on $n = 83$; `0.6442`, $P = 0.0155$ on $n = 107$) after bulk homogenization (Table 2, Fig. 5c,d, and Supplementary Tables S5, S5i).**"

---

#### Comment R1.3b [Major — Consensus `Cc.4`: Strictly separate pre-treatment biopsies ($n = 107$) from post-neoadjuvant surgical resection specimens ($n = 14$)]
> *"Keep the post-neoadjuvant surgical resection cohort (`GSE207422` scRNA, $n = 14$) clearly separated from strictly pre-treatment predictive claims."*

**Strategy:** `stratify_by_biopsy_timing` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We removed the ambiguous `single-cell AUROC = 0.8750` parenthetical from the *Abstract* pre-treatment sentence (`Line 7`), and explicitly partitioned the `121` oncology specimens at the opening of *Results §6* (`Line 91`), *Table 2* (`Line 225`), and *Fig. 5d* legend (`Line 261`) into: (i) $n = 83$ strictly pre-treatment pure anti-PD-1 monotherapy biopsies, (ii) $n = 24$ strictly pre-treatment neoadjuvant chemo-anti-PD-1 biopsies ($n = 107$ total strictly pre-treatment biopsies), and (iii) $n = 14$ post-neoadjuvant surgical resection specimens reflecting on-treatment microenvironment remodeling.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §6, Line 91):**
> "...we evaluated five independent human anti-PD-1 cohorts ($n = 121$ profiled tumor samples) **separated strictly by biopsy timing and companion regimen: $n = 83$ strictly pre-treatment pure anti-PD-1 monotherapy biopsies, $n = 24$ strictly pre-treatment neoadjuvant chemo-anti-PD-1 biopsies ($n = 107$ total strictly pre-treatment biopsies), and $n = 14$ post-neoadjuvant surgical resection specimens reflecting on-treatment microenvironment remodeling.**"

---

### 3.2 Role 2 — Methods & Biostatistical Rigor Reviewer (Score: 86/100)

#### Comment R2.1 [Major — Consensus `Cc.1`: Sequential `M_{k}` vs `M_{k-1}` and `M6` vs `M3` DeLong tests + Core-6 `Group5Fold` non-monotonicity]
> *"Add pairwise sequential DeLong tests (`M_{k}` vs `M_{k-1}` and `M6` vs `M3`) to `Table_S3f_SixStep_Ablation_and_GroupKFold.csv` and Table 1, and explicitly address the Core-6 `Group5Fold` non-monotonicity (`M3 = 0.6575` vs `M6 = 0.6385`, likely variance inflation from 6 extra progression/barrier coefficients across 386 failure events when holding out target groups)."*

**Strategy:** `add_statistical_test` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added four exact paired two-sided DeLong columns (`delta_auroc_vs_m3`, `delong_p_vs_m3`, `delta_auroc_vs_prev`, and `delong_p_vs_prev`) to `Table_S3f_SixStep_Ablation_and_GroupKFold.csv`, documented the cumulative-and-incremental DeLong protocol in *Methods §6* (`Line 174`), reported the `M3` and `M6 vs M3` statistics across *Table 1* (`Lines 197–198, 202, 204, 206`), and analyzed the Core-6 `Group5Fold` variance trade-off in *Results §3* (`Line 51`) and *Discussion (Limitations)* (`Line 113`). Specifically, on the Core-6 scientific subset ($n = 5,662$ trials across 293 targets, `386` biological failure events), stratified CV rises from `M3 = 0.6766` to `M6 = 0.6792` ($\Delta = +0.0025$, $P = 0.696$), whereas holding out entire target-gene families under `Group5Fold` yields `M3 = 0.6575` versus `M6 = 0.6415` ($\Delta = -0.0160$, $P = 0.0188$) due to higher parameter-estimation variance when fitting six additional progression/barrier coefficients (`9` total features in `M6` versus `3` in `M3`) across 386 failure events without cross-axis consensus regularization. Once event power increases on the 23,761-trial Definitive Maturation Cohort, `M6` significantly outperforms `M3` under both Stratified 5-Fold (`0.6744` vs `0.6607`, $\Delta = +0.0137$, $P = 8.93 \times 10^{-13}$) and unseen-target `Group5Fold` (`0.6509` vs `0.6432`, $\Delta = +0.0076$, $P = 1.70 \times 10^{-4}$).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Methods §6, Line 174 & Table 1, Line 206):**
> *Methods §6 (`Line 174`):* "Correlated ROC curves on shared cohorts were compared using two-sided paired nonparametric DeLong tests [22], **evaluating both cumulative gains over `M1` (static $\tau_0$) and incremental step attributions (`M6` versus `M3` early-disease specificity and sequential `M_{k}` versus `M_{k-1}`).**"  
> *Table 1 (`Line 206`):* "| All 12 Axes Maturation: `DNAcc` vs Static $\tau_0$ (Stratified / `Group5Fold` CV) | 23,761 | 561 | $\tau_0$: 0.5412 / 0.4879 | **`DNAcc`: 0.6744 / 0.6509 (`M3`: 0.6607 / 0.6432)** | **+0.1332 / +0.1630** | **$4.12 \times 10^{-180}$ / $2.33 \times 10^{-244}$** | **`M6 vs M3`: $P = 8.93 \times 10^{-13}$ / $1.70 \times 10^{-4}$**; FE OR: 1.32 ($1.84 \times 10^{-52}$ / $3.09 \times 10^{-7}$) |"

---

#### Comment R2.2 [Major — Consensus `Cc.2`: Harmonize cross-platform HPA `nTPM` vs scRNA UMI scale in 10-gene `D_ECM` (`Table_S3m`)]
> *"Harmonize the 10-gene ECM signature onto a common scale across all 12 axes (e.g., transforming HPA consensus `nTPM` via $\log_2(\mathrm{nTPM}+1)$ on stromal/organ-resident types or computing standardized z-scores within a unified reference) and verify that `DNAcc` ablation and within-modality performance remain invariant."*

**Strategy:** `recompute_harmonized_metric` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** In `omniniche/stage2/pan_disease_drift_law.py`, we harmonized both primary advanced-disease single-cell stromal pseudobulk profiles (`scRNA_Late_Disease_10Gene_ECM`) and Human Protein Atlas v23 organ stromal cell types (`HPA_v23_Organ_Stromal_10Gene_ECM`) onto an identical $\log_2(\mathrm{CPM}+1)$ expression scale across the 10-gene collagen-myofibroblast signature (`COL1A1`, `COL1A2`, `COL3A1`, `FN1`, `POSTN`, `ACTA2`, `FAP`, `LUM`, `DCN`, `THBS2`) prior to min-max normalization onto $[0.070, 0.350]$, and regenerated `Table_S3m` and all downstream `DNAcc` tables and figures. As shown in revised `Table_S3m_Empirical_Axis_ECM_Density.csv`, desmoplastic and fibrotic organs now rank at the biological apex while non-fibrotic retinal and intravascular axes fall to the bottom:
- `Solid_Oncology`: `0.3500` (raw $\log_2(\mathrm{CPM}+1) = 9.048$)
- `Pulmonary_IPF_ARDS`: `0.3265` (`8.431`, corrected from `0.1235`)
- `Autoimmune_IBD_Rheum`: `0.2921` (`7.530`)
- `Cardiovascular_HF_Athero`: `0.2535` (`6.517`)
- `Dermatologic_Immune` & `Musculoskeletal_Bone`: `0.2417` (`6.208`)
- `Renal_CKD_AKI`: `0.2240` (`5.744`)
- `Hepatic_MASH_Fibrosis`: `0.2107` (`5.394`)
- `Metabolic_Endocrine`: `0.1325` (`3.344`)
- `Hematologic_Malignancy`: `0.1237` (`3.112`)
- `Ophthalmic_Retinal`: `0.0810` (`1.992`, corrected from `0.2286`)
- `Infectious_Sepsis`: `0.0700` (`1.705`)

Harmonizing `D_ECM` improved out-of-fold `DNAcc` discrimination across every primary benchmark (`0.6682 -> 0.6710` scientific stratified, $P = 1.11 \times 10^{-21}$; `0.6456 -> 0.6541` `Group5Fold`, $P = 3.91 \times 10^{-25}$; `0.6977 -> 0.7023` joint `Human Genetics + DNAcc`, $P = 1.58 \times 10^{-10}$; `1.38 -> 1.42` 3-pillar adjusted OR; `0.6663` small molecules, `0.6690` IgG antibodies, `0.6190` proteins/peptides).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Methods §4, Line 156):**
> "Axis-specific interstitial ECM barrier density $D_{\mathrm{ECM}}(\mathrm{axis}) \in [0.070, 0.350]$ was derived objectively from a 10-gene collagen-myofibroblast signature (`COL1A1`, `COL1A2`, `COL3A1`, `FN1`, `POSTN`, `ACTA2`, `FAP`, `LUM`, `DCN`, `THBS2`) **after harmonizing both advanced human single-cell stromal pseudobulk profiles and Human Protein Atlas v23 organ stromal cell types onto a common $\log_2(\mathrm{CPM}+1)$ expression scale prior to min-max normalization (`Solid_Oncology` = 0.3500, `Pulmonary_IPF_ARDS` = 0.3265, `Autoimmune_IBD_Rheum` = 0.2921, down to `Infectious_Sepsis` = 0.0700; `Supplementary Table S3m`)**, yielding modality-by-matrix barrier compatibility $B_{\mathrm{acc}}(m, d)$ across small molecules, peptides, and ~150-kDa IgG monoclonal antibodies."

---

#### Comment R2.3 [Major — Consensus `Cc.5`: Explicit *Methods* disclosure of six-axis consensus-drift assignment across `51.3%` of scientific trials]
> *"Document the exact cross-axis drift imputation/consensus mapping rule for the six non-single-cell axes (`Dermatologic_Immune` inheriting `Autoimmune_IBD_Rheum`; 5 axes using `consensus_drift`) directly in Methods ('Dynamic specificity drift...') and mention in Limitations that `51.3%` of All-12-Axes scientific trials rely on pan-organ consensus progression features alongside axis-specific `D_ECM`, which is why Core-6 validation (`AUROC = 0.6752`) is reported in parallel."*

**Strategy:** `disclose_methodological_protocol` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added full, verbatim disclosure of the six-axis mapping rule in *Methods §4* (`Line 156`) and explicitly noted this boundary condition alongside the parallel Core-6 validation (`AUROC = 0.6792`) in *Discussion (Limitations)* (`Line 113`).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Methods §4, Line 156):**
> "**For the six registry disease axes without dedicated stage-resolved single-cell progression datasets (`51.3%` of scientific trials), `Dermatologic_Immune` inherited mucosal epithelial-immune barrier drift features from `Autoimmune_IBD_Rheum`, whereas `Cardiovascular_HF_Athero`, `Metabolic_Endocrine`, `Hematologic_Malignancy`, `Musculoskeletal_Bone`, and `Ophthalmic_Retinal` assigned gene-specific pan-organ mean progression features (`consensus_drift`) averaged across the six directly profiled single-cell axes while preserving exact gene-level healthy $\tau_0$, vital-organ bystander hazard $H_{\mathrm{vital}}(g)$, drug modality, and axis-specific ECM density.**"

---

#### Comment R2.4 [Minor: State exact LOPO logistic mapping formula and `[0.08, 0.92]` clipping in *Methods*]
> *"State the exact algebraic formula of the LOPO rank-to-probability transformation and upper-bound clipping (`0.9200` tie across 26 patients in `Table_S5i`) in Methods."*

**Strategy:** `disclose_exact_formula` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added the exact LOPO median/MAD robust z-score logistic projection $(1 + \exp(-0.90\,z_{\mathrm{LOPO}}))^{-1}$ and $[0.08, 0.92]$ clipping interval (noting the observed range `0.2714–0.9200`, the 26 high-scoring patients [12 responders and 14 non-responders] tied at `0.9200` in `Fig. 5c`, and the unclipped `raw_icb_remission_prob` column archived in `Supplementary Table S5i`) directly in *Methods §5* (`Line 164`) and *Fig. 5c* legend (`Line 261`).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Methods §5, Line 164):**
> "Multi-cohort pre-treatment pooling ($n = 83$ pure monotherapy and $n = 107$ total pre-treatment biopsies) used **zero-label leave-one-patient-out (LOPO) robust median/MAD z-scoring mapped onto a logistic probability scale $(1 + \exp(-0.90\,z_{\mathrm{LOPO}}))^{-1}$ clipped to $[0.08, 0.92]$ (yielding observed aligned probabilities of $0.2714\text{-}0.9200$, with 26 high-scoring patients [12 responders and 14 non-responders] tied at the $0.9200$ upper bound in Fig. 5c while preserving unclipped `raw_icb_remission_prob` in Supplementary Table S5i).**"

---

#### Comment R2.5 [Minor: Explain Open Targets Human Genetics attenuation (`OR = 1.70 -> 1.12 -> 1.06`) under disease-axis fixed effects]
> *"Report and briefly interpret why `Open Targets Human Genetics` attenuates from unadjusted `OR = 1.70` ($P_{\mathrm{cluster}} = 4.47 \times 10^{-17}$) to `OR = 1.12` ($P_{\mathrm{cluster}} = 0.0477$) with 5 confounder fixed effects and `OR = 1.06` ($P_{\mathrm{cluster}} = 0.328$) in the 3-pillar fixed-effects model."*

**Strategy:** `interpret_statistical_attenuation` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added the exact unadjusted (`OR = 1.70`, $P_{\mathrm{cluster}} = 4.47 \times 10^{-17}$), five-confounder fixed-effects (`OR = 1.12`, 95% CI `1.00–1.25`, $P_{\mathrm{cluster}} = 0.0477$), and three-pillar fixed-effects (`OR = 1.06`, 95% CI `0.95–1.18`, $P_{\mathrm{cluster}} = 0.328$) odds ratios for Open Targets Human Genetics to *Results §3* (`Line 59`), *Table 1* (`Line 201`), and *Fig. 3c* legend (`Line 257`), explaining that the 11 disease-axis indicator variables absorb substantial between-indication variation in genetic and somatic mutation score prevalence (notably oncology somatic enrichment versus non-malignant axes), whereas `DNAcc` (`OR = 1.42`, $P_{\mathrm{cluster}} = 1.04 \times 10^{-4}$) discriminates target accessibility strictly within individual disease axes and drug modalities.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §3, Line 59):**
> "...whereas **Open Targets human genetic support attenuated from unadjusted OR = 1.70 ($P_{\mathrm{cluster}} = 4.47 \times 10^{-17}$) to adjusted OR = 1.12 (target-clustered 95% CI 1.00-1.25, $P_{\mathrm{cluster}} = 0.0477$)... Conversely, static $\tau_0$ collapsed to unity (OR = 0.97, $P_{\mathrm{cluster}} = 0.719$) and human genetics attenuated to OR = 1.06 (target-clustered 95% CI 0.95-1.18, $P_{\mathrm{cluster}} = 0.328$), because the 11 disease-axis indicators absorb between-axis genetic and somatic score prevalence while `DNAcc` operates within individual axes and modalities (Fig. 3c and Extended Data Fig. 3b).**"

---

### 3.3 Role 3 — Fibrosis, Immunology & Tumor Microenvironment Domain Expert (Score: 88/100)

#### Comment R3.1 [Major — Consensus `Cc.3`: Cross-organ probability compression & accurate biological routing of renal tubular `SLC5A2`]
> *"Clarify in Results §4 and Discussion why pulmonary and renal winners (`pirfenidone`, `nerandomilast`, `dapagliflozin`, `finerenone`) receive low absolute remission probabilities (`0.0535–0.2330`) under liver-compiled parameters, and refine the biological description of `SLC5A2` (which IS expressed in renal proximal tubular epithelium in `GSE271188`, $\tau_0 = 0.8959$, rather than being 'absent' from the kidney — its indirect routing reflects tubular hemodynamic unloading rather than absent organ expression)."*

**Strategy:** `refine_biological_mechanism` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** In addition to detailing the cross-organ probability compression (`0.0535–0.2330`) and Youden vs `0.50` calibration (`R1.2a`–`R1.2b`), we corrected the wording for `SLC5A2` in *Results §4* (`Line 69`) and *Methods §5* (`Line 162`) to distinguish extrahepatic incretin receptors (`GLP1R` and `GIPR`, which lack intrahepatic receptor expression) from renal proximal-tubule `SLC5A2` ($\tau_0 = 0.8959$), which is strongly expressed in renal tubular epithelium in `GSE271188` rather than in `Myofibroblast` or `SAMac` scar compartments and unloads `Parenchyma_Injury` via tubuloglomerular feedback hemodynamic reduction.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §4, Line 69):**
> "**Extrahepatic incretin peptides lacking intrahepatic receptor expression (`GLP1R` and `GIPR`) and renal proximal-tubule `SLC5A2` ($\tau_0 = 0.8959$, expressed in tubular epithelium rather than `Myofibroblast` or `SAMac` scar compartments in `GSE271188` and acting via tubuloglomerular hemodynamic unloading) were routed objectively to upstream `Parenchyma_Injury` unloading.**"

---

#### Comment R3.2 [Major — Consensus `Cc.4`: Cutaneous BCC vs visceral Visium domain shift (`0.8750` vs `1.0000`) and NSCLC chemo-debulking $\tau_0$ (`0.8000`)]
> *"Discuss the histological domain mismatch between visceral pan-carcinoma Visium priors (`bladder, breast, ovarian, liver, GI`) and cutaneous BCC (`0.8750` with Visium vs `1.0000` without Visium, $P = 0.0266$), and acknowledge that in pre-treatment NSCLC chemo-immunotherapy biopsies (`GSE207422`, $n = 24$), concurrent platinum-doublet chemotherapy alters the predictive ranking (where `Static_Healthy_Tau0` reaches `0.8000` vs `OmniNiche = 0.7704`)."*

**Strategy:** `clarify_domain_boundary` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added both biological boundary observations explicitly to *Results §6* (`Lines 93–95`), *Table 2* (`Lines 224, 226`), and *Discussion (Limitations)* (`Line 115`): (i) in cutaneous BCC (`GSE123813`, $n = 8$), `OmniNiche (No-Visium)` reaches `1.0000` ($P = 0.0266$) versus `0.8750` with pan-carcinoma Visium priors because the 12 Visium sections derive from visceral carcinomas (`bladder, breast, ovarian, hepatocellular, gastrointestinal`) rather than cutaneous palisading basaloid nests; and (ii) in pre-treatment NSCLC neoadjuvant chemo-anti-PD-1 bulk biopsies (`GSE207422`, $n = 24$), concurrent platinum-doublet cytotoxic chemotherapy directly debulks proliferating tumor-stroma compartments, elevating static $\tau_0$ to `0.8000` alongside `OmniNiche` (`0.7704`, $P = 0.0288$) and `No-Visium` (`0.7556`).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §6, Lines 93–95):**
> *Line 93:* "In pre-treatment BCC monotherapy scRNA-seq ($n = 8$), OmniNiche achieved AUROC = 0.8750 (95% CI 0.5000-1.0000, AUPRC = 0.9167, Brier = 0.1423, exact $P = 0.1093$), exceeding `Linear Mean` (0.7500), `18-gene GEP` (0.6875), `CD8:myCAF` (0.5625), `PD-L1` (0.4375), and `CAF/TIDE` (0.3750), **whereas `OmniNiche (No-Visium)` reached AUROC = 1.0000 ($P = 0.0266$), reflecting histological domain mismatch between visceral pan-carcinoma Visium priors (bladder, breast, ovarian, liver, gastrointestinal) and cutaneous BCC architecture.**"  
> *Line 95:* "In pre-treatment NSCLC chemo-anti-PD-1 biopsies, OmniNiche achieved AUROC = 0.7704 (95% CI 0.5592-0.9407, exact $P = 0.0288$) versus `No-Visium` 0.7556, `Linear Mean` 0.7630, `18-gene GEP` 0.7704, `CAF/TIDE` 0.7630, `CD8:myCAF` 0.7481, `PD-L1` 0.6593, **and static $\tau_0$ 0.8000 under concurrent cytotoxic chemotherapy debulking.**"

---

#### Comment R3.3 [Minor-Major: Report exact paired *P*-values & T-cell sorting composition for $n = 15$ longitudinal `GSE123813` (`Fig. 5e`) and cite negative Phase IIb semaglutide F4 cirrhosis RCT]
> *"1) In `Fig. 5e` and Results §6 (`GSE123813` longitudinal $n = 15$), report the paired Wilcoxon / group-comparison P-value, define `Net CD8+ Cytotoxic Score` ($S_{\mathrm{CD8,eff}} - 0.5\,S_{\mathrm{CD8,exh}}$) in Methods, and note that 7/15 longitudinal patients come from T-cell-sorted fractions lacking paired `myCAF` counts (`pre_mycaf = 0.0` in `Table_S5f`). 2) Cite the published negative Phase IIb RCT of semaglutide in compensated NASH cirrhosis (Loomba et al., Lancet Gastroenterol. Hepatol. 2023) to clinically validate the predicted pre-saddle F2/F3 vs post-saddle F4 cirrhosis boundary for incretin therapy."*

**Strategy:** `add_statistical_test_and_clinical_citation` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** First, we computed exact nonparametric longitudinal statistics on `Table_S5f_BCC_GSE123813_Longitudinal_15Patients.csv` (`8` responders: mean $\Delta = +0.0345$, two-sided paired Wilcoxon $P = 0.945$; `7` non-responders: mean $\Delta = -0.0625$, paired Wilcoxon $P = 0.219$; responder vs non-responder $\Delta$ two-sided Mann-Whitney $P = 0.779$), annotated `Mann-Whitney P = 0.78 (n.s.)` directly onto `Fig. 5e`, moderated *Results §6* (`Line 97`), *Discussion (Limitations)* (`Line 115`), and *Fig. 5e* legend (`Line 261`) from *"confirmed"* to *"exploratory... non-significant directional divergence"*, and documented both the exact formula `Net CD8+ Cytotoxic Score` $= S_{\mathrm{CD8,eff}} - 0.5\,S_{\mathrm{CD8,exh}}$ and the `7/15` T-cell-sorted fraction caveat in *Methods §3 item 5* (`Line 147`). Second, we live-verified Loomba et al. (*Lancet Gastroenterol Hepatol* 2023;8(6):511–522, `PMID: 36934740`, `DOI: 10.1016/S2468-1253(23)00068-7`) via NCBI PubMed E-utilities, added it to `verified_pubmed_references.json` (`61` entries) and *References* (`[61]`, `Line 335`), and cited it in *Discussion (Limitations)* (`Line 115`) as direct clinical validation of the pre-cirrhotic stage boundary for incretin unloading.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §6 Line 97 & Discussion Line 115):**
> *Results §6 (`Line 97`):* "**We also examined exploratory paired pre- versus post-anti-PD-1 longitudinal scRNA-seq across 15 cutaneous carcinoma patients (`GSE123813`, 11 BCC and 4 cutaneous SCC, 7 of whom were profiled from T-cell-sorted fractions lacking paired non-`CD3+` `myCAF` counts). Longitudinal `Net CD8+ Cytotoxic Score` ($S_{\mathrm{CD8,eff}} - 0.5\,S_{\mathrm{CD8,exh}}$) showed a non-significant directional divergence between responders (mean $\Delta = +0.0345$, paired Wilcoxon $P = 0.945$) and non-responders ($\Delta = -0.0625$, paired Wilcoxon $P = 0.219$; responder versus non-responder Mann-Whitney $P = 0.779$) (Fig. 5e and Supplementary Table S5f).**"  
> *Discussion (`Line 115`):* "...and incretin unloading exhibits strict pre-cirrhotic stage dependence **consistent with negative Phase IIb semaglutide histology in compensated F4 cirrhosis [61].**"

---

### 3.4 Role 4 — Data Integrity, Reproducibility & Visual Forensics Reviewer (Score: 89/100)

#### Comment R4.1 [Major — Consensus `Cc.2`: Propagate harmonized `Table_S3m` ECM densities across all deterministic CSVs, tables & figures]
> *"Reconcile `Table_S3m_Empirical_Axis_ECM_Density.csv` by harmonizing HPA v23 and single-cell pseudobulk onto a shared log-transformed stromal scale so cross-organ ECM ordering is biologically coherent (`IPF` and `MASH` not artificially depressed relative to `Ophthalmic_Retinal`), and update any downstream `DNAcc` table/figure digits if they shift."*

**Strategy:** `recompute_and_synchronize_archive` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** Following $\log_2(\mathrm{CPM}+1)$ harmonization in `omniniche/stage2/pan_disease_drift_law.py` (`R2.2`), we re-executed the entire deterministic Stage II Task 2 statistical fortification pipeline (`Peak RAM = 2,306.0 MB < 7,000 MB`), re-rendered all 9 figures (`Fig. 1`–`Fig. 5`, `Extended Data Figs. 1–4`, passing all 10 aesthetic gates and 100% PNG visual inspection), and synchronized every affected digit across `OmniNiche_Flagship_Manuscript_EN.md`, `OmniNiche_Flagship_Manuscript_ZH_Preview.md`, `OmniNiche_Flagship_Manuscript_EN.docx`, and `Supplementary_Tables_Master_Index.md` (verified zero stale v2.5 numbers via automated sweep and regression test `test_v26_harmonized_ecm_m6_vs_m3_and_manuscript_sync`).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §2, Line 41):**
> "...while empirical interstitial ECM barrier density $D_{\mathrm{ECM}}(\mathrm{axis}) \in [0.070, 0.350]$ was **derived on a harmonized $\log_2(\mathrm{CPM}+1)$ 10-gene structural collagen-myofibroblast scale across all 12 axes (highest in solid oncology at 0.3500 and pulmonary fibrosis at 0.3265, lowest in bloodstream sepsis at 0.0700; Supplementary Table S3m).**"

---

#### Comment R4.2 [Minor-Major: Restore `CD8:myCAF` 6th series in `Fig. 5d` & make *Table 2* baseline columns symmetric across oncology rows]
> *"1) In `Fig. 5d` and its legend (Line 253), restore strict 1:1 parity: either plot `CD8:myCAF` as a 6th bar series in `Fig. 5d` or remove `CD8:myCAF` from the `Fig. 5d` legend description so legend and graphic match 100%. 2) Make Table 2 baseline columns symmetric across all oncology rows (reporting `No-Visium`, `Uninformative Linear Mean`, and `CD8:myCAF` consistently across `BCC`, `NSCLC`, `Melanoma`, `Pooled n = 83`, and `Pooled n = 107`)."*

**Strategy:** `fix_figure_and_table_parity` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** In `omniniche/viz/stage2_flagship_figures.py`, we added `CD8:myCAF` (`CD8_to_myCAF_Algebraic_Ratio`: `0.6263` on $n = 107$, `0.5873` on $n = 83$, `0.7481` on NSCLC pre-biopsy, `0.7750` on NSCLC scRNA, `0.5625` on BCC scRNA) as the sixth plotted horizontal series in `Fig. 5d` alongside its dedicated legend entry, verified visually on the rendered 600-DPI PNG and asserted in `pytest` (`assert "CD8:myCAF" in fig5_svg`). In parallel, we standardized *Table 2* (`Lines 224–229`) so every oncology row reports `No-Visium`, `Linear`, `PD-L1`, `GEP18`, `CAF/TIDE`, and `CD8:myCAF` in identical order.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Table 2, Lines 224 & 229):**
> "| **Pre-Treatment BCC Anti-PD-1 Monotherapy Full-TME scRNA-seq (`GSE123813`)** | 8 patients | 4 responders | 0.8750 (0.5000-1.0000) | 0.9167 / 0.1423 | **`No-Visium`: 1.0000 ($P = 0.0266$); Linear: 0.7500; `PD-L1`: 0.4375; `GEP18`: 0.6875; `CAF/TIDE`: 0.3750; `CD8:myCAF`: 0.5625** | $P = 0.1093$ (`umi_stride = 1` scRNA) |  
> ...  
> | **Pooled All Strictly Pre-Treatment Human Tumor Biopsies (LOPO Zero-Refit)** | 107 patients | 37 responders | 0.6367 (0.5253-0.7425) | 0.4489 / 0.2876 | **`No-Visium`: 0.6210; Linear: 0.6442 ($P = 0.0155$); `PD-L1`: 0.5610; `GEP18`: 0.5948; `CAF/TIDE`: 0.6023; `CD8:myCAF`: 0.6263** | $P = 0.01870$ (+0.0157 vs `No-Visium`) |"

---

#### Comment R4.3 [Minor: Harmonize `20,880` quality-filtered Visium spots (of `21,796` profiled) & clarify directed $k = 6$ Moran boundary asymmetry]
> *"1) Harmonize the 10x Visium spot count across Results §6 (`Line 83`), Extended Data Table 1 (`Line 236`), Figure Legends (`Lines 253, 261`), and `Fig. 5b` / `ED Fig. 4a` annotations so `20,880` quality-filtered spots (of `21,796` total profiled spots) is stated consistently everywhere. 2) Clarify in Methods (`Line 156`) that $k = 6$ nearest-neighbor Moran lag is directed (`0.533` vs `0.535` off-diagonal boundary asymmetry in `ED Fig. 4b`) rather than strictly 'symmetric'."*

**Strategy:** `harmonize_counts_and_terminology` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** First, we updated *Results §2* (`Line 41`), *Results §6* (`Line 87`), *Methods §3 item 5* (`Line 147`), *Extended Data Table 1* (`Line 244`), *Fig. 1d* (`20,880` QC spots), *Fig. 5b* (`20,880` QC spots), *Extended Data Fig. 4a* (`20,880` QC spots), and *Figure Legends* (`Lines 253, 261, 269`) so `20,880` quality-filtered spots (of `21,796` total profiled spots) is stated uniformly across text, tables, and graphics. Second, we replaced `"symmetric bivariate Moran spatial lag"` with `"directed k = 6 nearest-neighbor bivariate Moran spatial lag"` in *Methods §5* (`Line 164`) and *Extended Data Fig. 4b* legend (`Line 269`), noting the minor boundary-induced off-diagonal directional asymmetry ($\le 0.002$) on irregular tissue margins.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Methods §5, Line 164):**
> "For each 10x Visium section ($21,796$ total profiled spots across 12 carcinoma sections; $20,880$ quality-filtered spots in Supplementary Tables S5a,g), standardized spot-level compartment scores $z_c(u)$ were linked via physical $k = 6$ Euclidean nearest-neighbor graphs, mapping **directed $k = 6$ nearest-neighbor bivariate Moran spatial lag $I_{ij}$ (exhibiting minor boundary-induced off-diagonal directional asymmetry $\le 0.002$ on irregular tissue margins in Extended Data Fig. 4b)** to proximity $P_{ij} = (1 + \exp(-4.5\,I_{ij}))^{-1}$..."

---

#### Comment R4.4 [Minor: Fix `cancer_cohort` -> `cohort` column name in `Supplementary_Tables_Master_Index.md`]
> *"Update `Supplementary_Tables_Master_Index.md` entry for `Table_S5i` from `cancer_cohort` to `cohort` to match the exact CSV header."*

**Strategy:** `fix_documentation_schema` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We corrected `cancer_cohort` to `cohort` in `manuscript/Supplementary_Tables_Master_Index.md` (entry `54`, `Table_S5i_Pooled_PreTreatment_107Patients_Master.csv`) and synchronized all harmonized v2.6 metrics, the four new DeLong step-attribution columns in `Table_S3f`, the LOPO `[0.08, 0.92]` clipping metadata in `Table_S5i`, and the 61-reference archive count (`verified_pubmed_references.json`).

**Tracked Changes Excerpt (`Supplementary_Tables_Master_Index.md`, Section 5, Item 54):**
> "54. **`Table_S5i_Pooled_PreTreatment_107Patients_Master.csv`** — Patient-level zero-refit LOPO master prediction matrix across all `107` strictly pre-treatment anti-PD-1 solid-tumor biopsies (`patient_id`, **`cohort`**, `platform`, `clinical_responder`, `pdl1_status`, `mycaf_barrier_score`, `lopo_aligned_prob`, `raw_icb_remission_prob`... **(LOPO logistic projection clipped to `[0.08, 0.92]`, observed `0.2714-0.9200` alongside unclipped `raw_icb_remission_prob`)**)."

---

### 3.5 Role 5 — Devil's Advocate / Red-Team Challenger (Score: 81/100)

#### Comment R5.1 [Major — Consensus `Cc.1`: Is `DNAcc` mostly static early-disease specificity (`M3`) dressed in dynamical terminology?]
> *"Report `M6 vs M3` DeLong P-values across all cohorts and explicitly acknowledge in Abstract, Results §3, and Discussion that replacing healthy $\tau_0$ with diseased-lesion specificity (`M3`) accounts for >90% of the global registry gain, while dynamic drift and steric compatibility act as a finer-grained late-stage and maturation refinement."*

**Strategy:** `add_statistical_test_and_moderate_claim` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** As detailed in `R1.1a`–`R1.1b` and `R2.1`, we executed the exact `M6 vs M3` DeLong audit across all cohorts in `Table_S3f` and incorporated the results into *Abstract* (`Line 7`), *Results §3* (`Lines 49–51`), *Table 1* (`Lines 197–198, 202, 204, 206`), and *Discussion* (`Lines 105, 113`). We state upfront that `M3` ($\tau_{\mathrm{early}} = 0.6614$) captures `92.7%` of the global scientific gain over healthy $\tau_0$ (`0.5387`), while `M4`–`M6` progression drift and steric exclusion contribute statistically significant incremental discrimination in single-target scientific monotherapies (`+0.0394`, $P = 0.0353$) and across the 23,761-trial Definitive Maturation Cohort (`+0.0137`, $P = 8.93 \times 10^{-13}$ stratified; `+0.0076`, $P = 1.70 \times 10^{-4}$ `Group5Fold`), alongside candid disclosure of the Core-6 `Group5Fold` finite-sample inversion (`M3` `0.6575` vs `M6` `0.6415`, $P = 0.0188$).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Discussion Limitations, Line 113):**
> "**Moreover, while three-feature `M3` ($\tau_{\mathrm{early}}$) captures >90% of `DNAcc`'s scientific lift over healthy $\tau_0$ and slightly exceeds nine-feature `M6` under Core-6 unseen-target `Group5Fold` (0.6575 versus 0.6415) due to finite-sample variance across 386 failure events, progression drift and steric compatibility provide highly significant incremental gains across the 23,761-trial maturation cohort ($P = 8.93 \times 10^{-13}$) and single-target monotherapies ($P = 0.0353$).**"

---

#### Comment R5.2 [Major — Consensus `Cc.4`: Does spatial-Langevin complexity beat `Uninformative Linear Mean` (`0.6442` vs `0.6367`) or `No-Visium` (`1.0000` on BCC)?]
> *"Report `Uninformative_Linear_Combination` and `OmniNiche_NoVisium` side-by-side in Table 2 for every single-cell, bulk, and pooled cohort, and explicitly restrict the claim of nonlinear spatial-Langevin superiority to spatially/single-cell resolved TME profiles while acknowledging linear-combination equivalence in homogenized bulk RNA-seq."*

**Strategy:** `add_baseline_comparison_and_restrict_scope` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** As detailed in `R1.3a` and `R3.2`, every oncology row of *Table 2* (`Lines 224–229`) and *Results §6* (`Lines 93–95`) now presents `OmniNiche`, `OmniNiche (No-Visium)`, and `Uninformative Linear Mean` side-by-side. We explicitly restrict nonlinear spatial-Langevin superiority to single-cell TME cohorts where cellular compartment geometry and histological domain alignment are preserved (`NSCLC` scRNA $n = 14$: `0.8750` vs `No-Visium` `0.8500` and `Linear Mean` `0.8250`; `BCC` scRNA $n = 8$: `0.8750` vs `Linear Mean` `0.7500`, noting `No-Visium = 1.0000` due to visceral-to-cutaneous Visium domain shift), and state directly in *Results §6* (`Line 95`) and *Discussion (Limitations)* (`Line 115`) that bulk tissue homogenization collapses nonlinear spatial trajectories to parity with `Uninformative Linear Mean` (`0.6367` vs `0.6442` on $n = 107$; `0.6016` vs `0.6065` on $n = 83$).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Discussion Limitations, Line 115):**
> "**Third, in solid-tumor checkpoint blockade, bulk tissue lysis collapses concentric `myCAF` geometry to linear-combination parity (`OmniNiche` 0.6367 versus `Linear Mean` 0.6442 on $n = 107$), while pan-visceral Visium priors slightly attenuate cutaneous BCC discrimination relative to `No-Visium` (0.8750 versus 1.0000).**"

---

#### Comment R5.3 [Major — Consensus `Cc.3`: `0/2` renal and `1/3` pulmonary winners crossing `p = 0.50` in 24 antifibrotic trials]
> *"Report both rank-based discrimination (`AUROC = 0.8714`) and threshold-level calibration (`5/10` sensitivity at `p = 0.50` vs `10/10` sensitivity, `10/14` specificity at Youden cutoff `p >= 0.0535`) in Results §4, Table 2, and Limitations."*

**Strategy:** `recalibrate_threshold_and_disclose_limit` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** Addressed in full across *Results §4* (`Lines 75–77`), *Table 2* (`Line 219`), *Fig. 4c* (`Line 259`, with dual `Youden p >= 0.054` and `Basin p >= 0.50` vertical reference lines), *Methods §5* (`Line 162`), and *Discussion (Limitations)* (`Line 115`) (see `R1.2a`–`R1.2b`). We explicitly break down the `5/10` raw `0.50` basin sensitivity by organ (`4/5` hepatic, `1/3` pulmonary, `0/2` renal) alongside the `10/10` sensitivity and `10/14` specificity at the empirical Youden cutoff ($p \ge 0.0535$).

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Table 2, Line 219):**
> "| **Pan-Organ Human Fibrosis Phase IIb/III Arena (All 24 Trials; Liver $\to$ 3 Organs)** | 24 trials | 10 winners | 0.8714 (0.7071-0.9929) | 0.7769 / 0.2072 | Stage-blind: 0.7643; Stage prior: 0.7714; GRN: 0.8357; $\tau_0$: 0.4643; **Youden ($p \ge 0.0535$): 10/10 sens., 10/14 spec. (`p >= 0.50`: 5/10 sens., 13/14 spec.)** | $P = 0.00100$ ($P = 0.0311$ stage-blind) |"

---

#### Comment R5.4 [Minor-Major — Red-Team Vector `DA-U1`: `P = 1.78e-13` on $n = 70$ non-responders (`Fig. 5g`) reflects deterministic operator displacement]
> *"Clarify in Results §6 and `Fig. 5g` legend that $P = 1.78 \times 10^{-13}$ tests consistency of the simulated operator shift across empirical patient initial states rather than clinical efficacy evidence, framing the combination rescue predictions strictly as prospective hypotheses."*

**Strategy:** `qualify_in_silico_inference` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We revised *Results §6* (`Line 97`), *Discussion (Limitations)* (`Line 115`), *Fig. 5g* panel annotation (`Op. P = 1.8e-13` instead of bare `P = 1.8e-13`), and *Fig. 5g* legend (`Line 261`) to state explicitly that $P = 1.78 \times 10^{-13}$ measures the paired consistency of deterministic in silico operator displacement across empirical patient baseline states, nominating prospective combination hypotheses rather than observed clinical trial efficacy.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Results §6, Line 97):**
> "**This paired Wilcoxon shift ($P = 1.78 \times 10^{-13}$) measures the consistency of in silico operator displacement across empirical patient baseline states, nominating prospective combination hypotheses rather than observed clinical trial efficacy.**"

---

#### Comment R5.5 [Minor — Red-Team Vector `DA-U2`: Disclose uniform reference potency $\mathrm{pIC}_{50} = 7.5$ and extrahepatic reference leverage $\eta = 0.5400$ (`0.7637` dual)]
> *"State the exact reference leverage constant ($\eta = 0.5400$; $0.7637$ for dual `GLP1R;GIPR`) and uniform reference potency ($\mathrm{pIC}_{50} = 7.5$) in Methods, noting that drug-specific human free-exposure / receptor-occupancy data could further calibrate cross-organ absolute probabilities."*

**Strategy:** `disclose_exact_parameter` | **Status:** `DONE` | **Follow-up Risk:** `low`

**Response:** We added exact disclosure of reference leverage $\eta = 0.5400$ ($0.7637 = 0.5400 \times \sqrt{2}$ for dual `GLP1R;GIPR` agonism) alongside uniform reference potency $\mathrm{pIC}_{50} = 7.5$ (`31.6 nM`) directly in *Methods §5* (`Line 162`), clarifying that holding potency fixed prevents post-hoc pharmacokinetic parameter tuning and isolates compartmental target topology, receptor-sink dilution, and macromolecular steric permeability.

**Tracked Changes Excerpt (`OmniNiche_Flagship_Manuscript_EN.md`, Methods §5, Line 162):**
> "...were routed to upstream `Parenchyma_Injury` unloading **(reference leverage $\eta = 0.5400$; $0.7637$ for dual `GLP1R;GIPR`)**."

---

## 4. Revision Task Table (`T1`–`T10`, 100% Verified `DONE`)

| Task ID | Priority | Consensus / Reviewer Source | Concrete Verification Target in Revised Code, Tables, Figures & Manuscript | Status |
| :---: | :---: | :--- | :--- | :---: |
| **`T1`** | **Major** | `Cc.2` (Roles 2 `#2`, 4 `#1`) | Harmonized 10-gene ECM signature onto common $\log_2(\mathrm{CPM}+1)$ stromal scale in `pan_disease_drift_law.py`; regenerated `Table_S3m` (`Solid_Oncology = 0.3500`, `Pulmonary_IPF_ARDS = 0.3265`, `Ophthalmic_Retinal = 0.0810`, `Infectious_Sepsis = 0.0700`) and all downstream `DNAcc` CSVs, `Fig. 1d`, `Fig. 3`, `ED Figs. 1–3`, and *Table 1* (`Lines 41, 156, 197–213`). | **DONE** |
| **`T2`** | **Major** | `Cc.1` (Roles 1 `#1`, 2 `#1`, 5 `#1`) | Added `delta_auroc_vs_m3`, `delong_p_vs_m3`, `delta_auroc_vs_prev`, and `delong_p_vs_prev` to `Table_S3f`; reported exact step attribution (`M3` `92.7%` of scientific lift; `M6 vs M3` $P = 0.0353$ single-target, $P = 8.93 \times 10^{-13}$ maturation; Core-6 `Group5Fold` variance trade-off $P = 0.0188$) and moderated claims (`Lines 7, 49–51, 105, 113, 174, 197–206`). | **DONE** |
| **`T3`** | **Major** | `Cc.4` (Roles 1 `#3`, 3 `#2`, 5 `#2`) | Separated $n = 83$ pre-treatment monotherapy, $n = 24$ pre-treatment chemo-anti-PD-1 ($n = 107$ total), and $n = 14$ post-neoadjuvant resections; tabulated symmetric `No-Visium` and `Uninformative Linear Mean` across all *Table 2* rows (`Lines 224–229`) and moderated bulk spatial-erasure & BCC visceral-to-cutaneous domain shift (`Lines 91–95, 115`). | **DONE** |
| **`T4`** | **Major** | `Cc.3` (Roles 1 `#2`, 3 `#1`, 5 `#3`) | Added dual Youden-optimal ($p \ge 0.0535$: `10/10` sens., `10/14` spec.) vs physical `0.50` (`5/10` sens., `13/14` spec.) vertical threshold lines to `Fig. 4c` and reported full calibration + hepatic Hill compression mechanics in *Results §4* (`Lines 75–77`), *Table 2* (`Line 219`), *Methods* (`Line 162`), *Fig. 4c* legend (`Line 259`), and *Limitations* (`Line 115`). | **DONE** |
| **`T5`** | **Major** | Role 4 (`#2`) | Plotted `CD8:myCAF` as the 6th bar series in `Fig. 5d` (`stage2_flagship_figures.py`), matching the 6-series legend and *Table 2* (`Lines 224–229, 261`); verified visually via PNG inspection and `pytest` SVG check. | **DONE** |
| **`T6`** | **Minor-Major** | `Cc.5` (Role 2 `#3`) | Documented exact 6-axis cross-axis drift assignment (`Dermatologic_Immune` mucosal inheritance + 5-axis `consensus_drift` across `51.3%` of scientific trials) in *Methods §4* (`Line 156`) and *Limitations* (`Line 113`). | **DONE** |
| **`T7`** | **Minor-Major** | Role 3 (`#1`, `#3`) | Clarified renal proximal-tubule `SLC5A2` ($\tau_0 = 0.8959$ in `GSE271188` epithelium, hemodynamic unloading) (`Lines 69, 162`); live-verified and cited Loomba et al. *Lancet Gastroenterol Hepatol* 2023 `[61]` (`PMID: 36934740`) (`Lines 115, 335`). | **DONE** |
| **`T8`** | **Minor-Major** | Role 3 (`#3`), `DA-U1` (Role 5 `#4`) | Reported exact paired Wilcoxon ($P = 0.945$, $P = 0.219$) and Mann-Whitney ($P = 0.779$, n.s.) statistics and `7/15` T-cell-sorted caveat for `GSE123813` longitudinal scRNA-seq (`Fig. 5e`, `Lines 97, 115, 147, 261`); qualified `Fig. 5g` $P = 1.78 \times 10^{-13}$ as deterministic operator displacement (`Lines 97, 115, 261`). | **DONE** |
| **`T9`** | **Minor** | Role 2 (`#4`, `#5`), `DA-U2` (Role 5 `#5`) | Added exact LOPO logistic mapping $(1 + \exp(-0.90\,z_{\mathrm{LOPO}}))^{-1}$ and `[0.08, 0.92]` clipping (`Line 164`), Open Targets genetics FE attenuation (`OR = 1.70 -> 1.12 -> 1.06`) (`Lines 59, 201`), and reference potency/leverage ($\mathrm{pIC}_{50} = 7.5$, $\eta = 0.5400$ / $0.7637$) (`Line 162`). | **DONE** |
| **`T10`** | **Minor** | Role 4 (`#3`, `#4`) | Harmonized `20,880` quality-filtered Visium spots (of `21,796` profiled) across text, *Extended Data Table 1*, legends, and `Fig. 1d`/`Fig. 5b`/`ED Fig. 4a`; clarified directed $k = 6$ Moran boundary asymmetry ($\le 0.002$) (`Lines 164, 269`); fixed `cancer_cohort` $\to$ `cohort` in `Supplementary_Tables_Master_Index.md`. | **DONE** |
