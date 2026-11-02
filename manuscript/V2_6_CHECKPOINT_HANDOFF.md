# Stage II v2.6 Revision Checkpoint & Instant Resume Handoff

**Saved At**: `2026-11-01 13:18 +08:00`  
**Active Review Reference**: `Simulated_NatMed_5Reviewer_Panel_v2_5.md` (`Aggregate = 86.5/100`, `Major Revision`, Tasks `T1`–`T10`)

---

## 1. Completed & Gate-Verified Work (Do NOT Re-Run Heavy Pipelines)

1. **Code & All 60 Deterministic Empirical Tables (`outputs/stage2_tables/`) — DONE**:
   - `omniniche/stage2/pan_disease_drift_law.py`: Harmonized $\log_2(\mathrm{CPM}+1)$ 10-gene `D_ECM` scaling (`Solid_Oncology = 0.3500`, `Pulmonary_IPF_ARDS = 0.3265`, `Autoimmune_IBD_Rheum = 0.2921`, ..., `Infectious_Sepsis = 0.0700` in `Table_S3m_Empirical_Axis_ECM_Density.csv`).
   - `omniniche/stage2/v2_statistical_fortification.py`: Added exact pairwise `M6 vs M3` (`delta_auroc_vs_m3`, `delong_p_vs_m3`) and sequential (`delta_auroc_vs_prev`, `delong_p_vs_prev`) two-sided DeLong tests to `Table_S3f_SixStep_Ablation_and_GroupKFold.csv`.
   - All downstream Task 2 statistical tables regenerated (`Peak RAM = 2306.0 MB < 7,000 MB`).
2. **All 9 Publication Figures (`.png`, `.pdf`, `.svg` in `outputs/stage2_figures/`) — DONE & GATE-PASSED**:
   - `Fig1`–`Fig5` and `Extended_Data_Fig1`–`4` re-rendered from v2.6 tables, 100% visually inspected via `read_file`, polished (`Fig. 4c` borderless `Youden p >= 0.053` vs `Basin p >= 0.50` legend; `Fig. 5b` `20,880 QC Spots`; `Fig. 5d` 6-series plot including `CD8:myCAF` with non-overlapping 2-col legend; `Fig. 5e` exact $\Delta$ and Mann–Whitney $P = 0.78$, n.s.).
   - Passed `python "C:/Users/sheng/.claude/scripts/figure_clean_aesthetic_gate.py" --allow-schematic "F:/OmniNiche_Universal/omniniche/viz/stage2_flagship_figures.py" "F:/OmniNiche_Universal/outputs/stage2_figures"` -> **`ALL 10 GATES PASSED (脚本: 1 个 | PNG: 9 张)`**.
3. **Live NCBI PubMed E-Utilities Verification for Compensation Reference `[61]` — DONE**:
   - `Loomba R, Abdelmalek MF, Armstrong MJ, et al. Semaglutide 2.4 mg once weekly in patients with non-alcoholic steatohepatitis-related cirrhosis: a randomised, placebo-controlled phase 2 trial. Lancet Gastroenterol Hepatol 2023;8(6):511-522. PMID: 36934740. DOI: 10.1016/S2468-1253(23)00068-7.` (Cite as `[61]` in `Discussion` Limitations after `[60]` to preserve strict monotonic order `[1]`..`[61]`).

---

## 2. Exact Verified v2.6 Metrics for Manuscript & Rebuttal Synchronization

- **All-12 Scientific (`n = 11,637`, `G = 374`)**:
  - Stratified 5-Fold `M1..M6`: `0.5387 (0.5152–0.5610) -> 0.5618 -> 0.6614 (0.6393–0.6840) -> 0.6589 -> 0.6650 -> 0.6710 (0.6480–0.6932)`; `M6 vs M1`: $\Delta = +0.1323$, DeLong $P = 1.11 \times 10^{-21}$; `M3 vs M1`: $\Delta = +0.1227$, $P = 7.38 \times 10^{-18}$ (92.7% of total lift); `M6 vs M3`: $\Delta = +0.0096$, $P = 0.205$.
  - `Group5Fold` `M1..M6`: `0.5030 (0.4805–0.5256) -> 0.5206 -> 0.6389 -> 0.6335 -> 0.6466 -> 0.6541 (0.6300–0.6763)`; `M6 vs M1`: $\Delta = +0.1510$, $P = 3.91 \times 10^{-25}$; `M6 vs M3`: $\Delta = +0.0151$, $P = 0.0563$.
  - Single-target scientific (`n = 6,115`, `G = 235`): Stratified `0.6404` vs `0.5073` ($P = 1.18 \times 10^{-5}$; `M6 vs M3` $\Delta = +0.0394$, $P = 0.0353$); `Group5Fold` `0.5989` vs `0.4697` ($P = 7.63 \times 10^{-5}$; `M6 vs M3` $+0.0367$, $P = 0.0540$). Multi-target (`n = 5,522`, `G = 290`): Stratified `0.6835` ($P = 4.59 \times 10^{-13}$), `Group5Fold` `0.6645` ($P = 2.63 \times 10^{-20}$).
- **Core-6 Single-Cell Scientific (`n = 5,662`, `G = 293`, `386` failures)**:
  - Stratified 5-Fold: `M1 = 0.5447`, `M3 = 0.6766`, `M6 = 0.6792 (0.6521–0.7054)` (`M6 vs M1` $\Delta = +0.1345$, $P = 2.11 \times 10^{-15}$; `M6 vs M3` $\Delta = +0.0025$, $P = 0.696$).
  - `Group5Fold`: `M1 = 0.5049`, `M3 = 0.6575`, `M6 = 0.6415 (0.6118–0.6699)` (`M6 vs M1` $\Delta = +0.1366$, $P = 4.97 \times 10^{-15}$; `M6 vs M3` $\Delta = -0.0160$, $P = 0.0188$).
- **Definitive Maturation (`n = 23,761`, `G = 561`)**:
  - Stratified 5-Fold `M1..M6`: `0.5412 -> 0.5777 -> 0.6607 -> 0.6623 -> 0.6716 -> 0.6744 (0.6674–0.6810)` ($\Delta = +0.1332$ vs `M1`, $P = 4.12 \times 10^{-180}$; `M6 vs M3` $\Delta = +0.0137$, $P = 8.93 \times 10^{-13}$; Brier `0.2477 -> 0.2273`).
  - `Group5Fold`: `M6 = 0.6509` vs `M1 = 0.4879` ($\Delta = +0.1630$, $P = 2.33 \times 10^{-244}$; `M6 vs M3` [`0.6432`] $\Delta = +0.0076$, $P = 1.70 \times 10^{-4}$). Single-target maturation (`n = 10,507`, `G = 374`): `Group5Fold` `0.6451` vs `0.4310` ($P = 3.81 \times 10^{-130}$; `M6 vs M3` $\Delta = +0.0295$, $P = 1.19 \times 10^{-12}$).
- **Human Genetics Orthogonality (`Table_S3l`)**:
  - Scientific Stratified: `G1 Genetics = 0.6362 (0.6085–0.6625)`, `G3 Genetics+tau0 = 0.6344` ($\Delta = -0.0018$, $P = 0.451$), `G4 DNAcc = 0.6710` ($\Delta = +0.0348$ vs `G1`, $P = 0.0224$), `G5 Joint = 0.7023 (0.6791–0.7244)` ($\Delta = +0.0661$ vs `G1`, $P = 1.58 \times 10^{-10}$).
  - Scientific `Group5Fold`: `G1 = 0.6175`, `G4 = 0.6541` ($P = 0.0215$), `G5 = 0.6844 (0.6614–0.7077)` ($\Delta = +0.0669$ vs `G1`, $P = 1.34 \times 10^{-9}$). Within `OT >= 0.10` (`n = 7,174`, `G = 291`): `0.6916 (0.6672–0.7159)` vs `0.5267` ($\Delta = +0.1650$, $P = 3.59 \times 10^{-25}$). Core-6 `G5 = 0.6974` ($P = 8.21 \times 10^{-5}$ vs `G1 = 0.6496`).
  - Maturation: Stratified `G5 = 0.7210 (0.7145–0.7273)` vs `G1 = 0.6631` ($\Delta = +0.0579$, $P = 1.72 \times 10^{-111}$); `Group5Fold` `G5 = 0.6998` vs `G1 = 0.6413` ($\Delta = +0.0584$, $P = 3.86 \times 10^{-115}$).
- **Within-Modality (`Table_S3k`)**: Small molecules `0.6663` vs `0.5451` ($\Delta = +0.1211$, $P = 1.58 \times 10^{-14}$); IgG antibodies `0.6690` vs `0.5959` ($\Delta = +0.0730$, $P = 0.0093$); Proteins/peptides `0.6190` vs `0.4870` ($\Delta = +0.1320$, $P = 0.0030$).
- **Multivariable Fixed-Effects (`Table_S3g`)**:
  - Scientific: Unadj. `DNAcc OR = 1.88` (Wald `1.73–2.05`, $P = 3.71 \times 10^{-46}$ / cluster $P = 5.68 \times 10^{-16}$); `+5 FE DNAcc OR = 1.42` (Wald `1.28–1.59`, $P = 1.69 \times 10^{-10}$; cluster `1.21–1.68`, $P_{\mathrm{cluster}} = 2.81 \times 10^{-5}$); `+5 FE Genetics OR = 1.12` (cluster `1.00–1.25`, $P_{\mathrm{cluster}} = 0.0477$, down from unadj. `1.70`, $P_{\mathrm{cluster}} = 4.47 \times 10^{-17}$); 3-Pillar FE: `DNAcc OR = 1.42` (Wald `1.27–1.60`, $P = 2.86 \times 10^{-9}$; cluster `1.19–1.70`, $P_{\mathrm{cluster}} = 1.04 \times 10^{-4}$), `Genetics OR = 1.06` (cluster `0.95–1.18`, $P_{\mathrm{cluster}} = 0.328$), `Static tau_0 OR = 0.97` (Wald `0.87–1.08`, $P = 0.529$; cluster $P_{\mathrm{cluster}} = 0.719$).
  - Maturation: `+5 FE DNAcc OR = 1.32` (Wald `1.28–1.37`, $P = 1.84 \times 10^{-52}$; cluster $P = 3.09 \times 10^{-7}$); 3-Pillar FE `DNAcc OR = 1.35` (Wald `1.30–1.40`, $P = 2.30 \times 10^{-53}$; cluster $P = 7.96 \times 10^{-7}$).
- **Per-Axis Maturation (`Table_S3b`)**: Hepatic `0.6166` ($\Delta +0.1238$, $P = 2.51 \times 10^{-3}$, `OR = 1.49`, $P = 2.01 \times 10^{-4}$); Pulmonary `0.6264` ($\Delta +0.0399$, $P = 6.27 \times 10^{-5}$, `OR = 1.54`, $P = 8.24 \times 10^{-19}$); Renal `0.5661` ($\Delta +0.0944$, $P = 7.04 \times 10^{-4}$, `OR = 1.25`, $P = 1.29 \times 10^{-3}$); Solid Oncology `0.5886` ($\Delta +0.0569$, $P = 2.05 \times 10^{-12}$, `OR = 1.41`, $P = 1.83 \times 10^{-45}$); Autoimmune IBD `0.7013` ($\Delta +0.0943$, $P = 1.04 \times 10^{-10}$, `OR = 2.12`, $P = 2.37 \times 10^{-35}$); Sepsis boundary `0.3822` vs `0.3616` ($\Delta +0.0206$, $P = 0.689$, `OR = 0.57`, $P = 3.79 \times 10^{-3}$).
- **Quartiles, `R1–R7` & Enrichment (`Tables S3d, S3e, S3j`)**:
  - Quartiles: `2.45-fold` from `27.10%` (`Q1`) to `66.28%` (`Q4`) vs `1.20-fold` (`44.96%` to `53.87%`).
  - `R1 = 0.6710`, `R2 = 0.6751`, `R3 = 0.6533` ($P = 3.89 \times 10^{-70}$), `R4 = 0.6713`, `R5 = 0.6742`, `R6 = 0.6767` vs `0.5177` ($P = 6.27 \times 10^{-28}$), `R7 = 0.6884` vs `0.5304` ($P = 1.21 \times 10^{-25}$).
  - Bottom-25% gate: `46.60% -> 53.10%`, preventing `4,326` failed trials (`~1.19 million` patient exposures [`1,193,976`]), retaining `85.51%` of winners.
- **Fibrosis Thresholds & Solid-Tumor Baselines (`Tables S4, S5, S5f, S5i`)**:
  - Fibrosis 24 trials: continuous `AUROC = 0.8714`; Youden-optimal $p \ge 0.0535$: `10/10` (`100.0%`) sensitivity, `10/14` (`71.4%`) specificity; raw basin $p \ge 0.50$: `5/10` (`50.0%`) sensitivity (`4/5` liver, `1/3` lung, `0/2` kidney), `13/14` (`92.9%`) specificity.
  - Solid-tumor `Table_S5` symmetric baselines: `BCC n=8`: `OmniNiche 0.8750 (P=0.1093)`, `No-Visium 1.0000 (P=0.0266)`, `Linear 0.7500`; `NSCLC scRNA n=14` (post-neoadjuvant resection): `OmniNiche 0.8750 (P=0.0324)`, `No-Visium 0.8500`, `Linear 0.8250`; `NSCLC bulk n=24` (chemo-anti-PD-1): `OmniNiche 0.7704 (P=0.0288)`, `No-Visium 0.7556`, `Linear 0.7630`, `tau_0 0.8000`; `Melanoma bulk GSE78220 / GSE91061`: `OmniNiche 0.5982 / 0.5385`, `No-Visium 0.5893 / 0.5500`, `Linear 0.6548 / 0.5385`; `Pooled monotherapy n=83`: `OmniNiche 0.6016 (P=0.1302)`, `No-Visium 0.5815`, `Linear 0.6065`; `Pooled all pre-treatment n=107`: `OmniNiche 0.6367 (P=0.0187)`, `No-Visium 0.6210`, `Linear 0.6442 (P=0.0155)`.
  - `Table_S5i` LOPO logistic projection $(1 + \exp(-0.90\,z_{\mathrm{LOPO}}))^{-1}$ clipped to `[0.08, 0.92]`: observed `0.2714–0.9200` (`0.3387–0.9200` in bulk), `26` patients tied at `0.9200` (`12` responders, `14` non-responders).

---

## 3. Immediate 4-Step Execution Sequence on Resume

1. Run `v26_update_manuscripts.py` updating `OmniNiche_Flagship_Manuscript_EN.md`, `OmniNiche_Flagship_Manuscript_ZH_Preview.md` (enforcing EN Abstract `<= 200` words, every narrative paragraph `<= 270` EN words / `<= 600` ZH chars, every EN sentence `<= 72` words, `0` `\u2014` em-dashes, `0` forbidden jargon, 1:1 line-parity), `Supplementary_Tables_Master_Index.md`, and `verified_pubmed_references.json` (`[61]` Loomba 2023).
2. Run `manuscript_convention_and_parity_gate.py` on `EN.md` + `--zh ZH_Preview.md` AND rebuild `OmniNiche_Flagship_Manuscript_EN.docx` (+ copy to `OmniNiche_Flagship_Manuscript_EN_v2_6.docx`) via `C:/Users/sheng/.gemini/antigravity/brain/9697612f-91b2-4dc2-916d-234e0cd73347/scratch/mod08b_build_submission_docx.py`.
3. Add `test_v26_harmonized_ecm_m6_vs_m3_and_manuscript_sync` to `tests/test_stage2_v2_fortification.py` and run `python -m pytest "F:/OmniNiche_Universal/tests" -v` (`47 passed`).
4. Write the 4-part `Point_by_Point_Rebuttal_NatMed_v2_6.md` (`biomed-rebuttal`, all atomic sub-questions `DONE` with verified line numbers and Tracked Changes Excerpts) to both the artifact directory and `F:/OmniNiche_Universal/manuscript/`, then mark `T1`–`T10` all `completed` in `REVISION_STATUS.json`.
