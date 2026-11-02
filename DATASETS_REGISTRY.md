# OmniNiche Universal — Station M1 Exhaustive Public Dataset Census & PRISMA Registry

- **Research Track**: `multiomics_spatial_mechanism` (hybrid multiscale biophysical coupling with `biomedical_digital_twin`)
- **Databases Queried (5)**: NCBI GEO, ClinicalTrials.gov API v2, EMBL-EBI ChEMBL v33, Open Targets Platform v24.09, Human Protein Atlas v23
- **PRISMA Arithmetic Closure**: `1038` unique candidate cohorts/registries screened = **`20` included** + **`1018` excluded** (`0` unscreened).
- **Data Integrity & RAM Locks**: 100% public deidentified human records; SHA-256 clinical-trial truth lock active; zero synthetic interpolation; exact single-UMI streaming (`umi_stride = 1`) with peak Python RAM `< 2.3 GB`.

## 1. Included Cohorts & Knowledgebases (`n = 20`)

| ID | Repository | Organ / Scope | Modality | Role | Official Samples (`n`) | Replicates / Group |
| :--- | :--- | :--- | :--- | :--- | :---: | :---: |
| `ClinicalTrials_gov_API_v2_Macro41396` | ClinicalTrials.gov API v2 | 12 systemic and parenchymal organ disease axes | `clinical_trial_registry` | `trial_truth_benchmark` | 41,396 | 565 |
| `ChEMBL_v33_OpenTargets_v24_09` | EMBL-EBI ChEMBL v33 + Open Targets v24.09 | 12 systemic and parenchymal organ disease axes | `pharmacogenomic_knowledgebase` | `orthogonal_genetics_and_potency_reference` | 29,791 | 374 |
| `Human_Protein_Atlas_v23` | Human Protein Atlas v23 | Multi-organ healthy human baseline | `scRNA-seq_consensus_reference` | `healthy_baseline_reference` | 154 | 15 |
| `GSE136103` | NCBI GEO | Hepatic_MASH_Cirrhosis | `scRNA-seq` | `discovery` | 26 | 5 |
| `GSE212837` | NCBI GEO | Hepatic_MASH_Cirrhosis | `single_nuclear_RNA_sequencing_HTS` | `discovery` | 20 | 3 |
| `GSE135893` | NCBI GEO | Pulmonary_IPF_ARDS | `scRNA-seq` | `primary_external_validation` | 34 | 10 |
| `GSE271188` | NCBI GEO | Renal_CKD_AKI | `snRNA-seq` | `primary_external_validation` | 37 | 5 |
| `GSE203612` | NCBI GEO | Solid_Oncology | `10x_Visium_spatial_transcriptomics` | `discovery` | 31 | 8 |
| `GSE171351` | NCBI GEO | Solid_Oncology | `10x_Visium_spatial_transcriptomics` | `spatial_prior_calibration` | 4 | 2 |
| `GSE123813` | NCBI GEO | Solid_Oncology | `scRNA-seq` | `external_validation` | 86 | 4 |
| `GSE207422` | NCBI GEO | Solid_Oncology | `scRNA-seq_and_bulk_RNA-seq` | `external_validation` | 39 | 9 |
| `GSE91061` | NCBI GEO | Solid_Oncology | `bulk_RNA-seq` | `external_validation` | 109 | 10 |
| `GSE78220` | NCBI GEO | Solid_Oncology | `bulk_RNA-seq` | `external_validation` | 28 | 12 |
| `GSE250487` | NCBI GEO | Autoimmune_IBD_Rheum | `scRNA-seq` | `external_validation` | 8 | 4 |
| `GSE167363` | NCBI GEO | Infectious_Sepsis | `scRNA-seq` | `external_validation` | 12 | 4 |
| `GSE226249` | NCBI GEO | Pulmonary_IPF_ARDS | `bulk_RNA-seq_perturbation` | `external_validation` | 27 | 9 |
| `GSE154646` | NCBI GEO | Solid_Oncology | `bulk_RNA-seq_perturbation` | `external_validation` | 24 | 8 |
| `GSE210020` | NCBI GEO | Pulmonary_IPF_ARDS | `bulk_RNA-seq_perturbation` | `external_validation` | 40 | 8 |
| `GSE65904` | NCBI GEO | Solid_Oncology | `bulk_microarray_clinical_survival` | `external_validation` | 214 | 102 |
| `GSE22153` | NCBI GEO | Solid_Oncology | `bulk_microarray_clinical_survival` | `external_validation` | 57 | 27 |

## 2. Summary of Excluded NCBI GEO Candidate Cohorts (`n = 1018`)

All `1018` excluded NCBI GEO series retrieved across the exhaustive multi-organ single-cell, spatial 10x Visium, and anti-PD-1 search queries are itemized with deterministic exclusion reasons in `DATASETS_REGISTRY.json` (first 25 representative exclusions shown below):

| GEO Accession | Official Samples (`n`) | Census Bucket | Exclusion Reason |
| :--- | :---: | :--- | :--- |
| `GSE341136` | 6 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: insufficient donor biological replication (n_samples=6 < 8) for stage-resolved pseudobulk drift or pre-treatment responder stratification. |
| `GSE350364` | 136 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE331320` | 11 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: non-primary human biopsy model (in vitro cell line, organoid, or xenograft/murine perturbation) lacking native human interstitial ECM architecture. |
| `GSE331214` | 38 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE349540` | 10 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE346465` | 14 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: non-primary human biopsy model (in vitro cell line, organoid, or xenograft/murine perturbation) lacking native human interstitial ECM architecture. |
| `GSE244745` | 355 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: lacks downloadable cell-level count matrix (.h5/.mtx/.h5ad) or spatial spot coordinates required for exact single-UMI streaming. |
| `GSE347204` | 6 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: non-primary human biopsy model (in vitro cell line, organoid, or xenograft/murine perturbation) lacking native human interstitial ECM architecture. |
| `GSE339010` | 14 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: central nervous system / blood-brain barrier neurovascular cohort outside the 12 peripheral systemic and parenchymal organ disease axes. |
| `GSE339006` | 10 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: central nervous system / blood-brain barrier neurovascular cohort outside the 12 peripheral systemic and parenchymal organ disease axes. |
| `GSE338999` | 4 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: central nervous system / blood-brain barrier neurovascular cohort outside the 12 peripheral systemic and parenchymal organ disease axes. |
| `GSE338996` | 8 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: central nervous system / blood-brain barrier neurovascular cohort outside the 12 peripheral systemic and parenchymal organ disease axes. |
| `GSE308435` | 42 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE308434` | 42 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE314616` | 3 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: insufficient donor biological replication (n_samples=3 < 8) for stage-resolved pseudobulk drift or pre-treatment responder stratification. |
| `GSE307398` | 12 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE337136` | 5 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: lacks downloadable cell-level count matrix (.h5/.mtx/.h5ad) or spatial spot coordinates required for exact single-UMI streaming. |
| `GSE337135` | 12 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE293914` | 1 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: insufficient donor biological replication (n_samples=1 < 8) for stage-resolved pseudobulk drift or pre-treatment responder stratification. |
| `GSE314906` | 81 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE301351` | 6 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: insufficient donor biological replication (n_samples=6 < 8) for stage-resolved pseudobulk drift or pre-treatment responder stratification. |
| `GSE337519` | 1 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: insufficient donor biological replication (n_samples=1 < 8) for stage-resolved pseudobulk drift or pre-treatment responder stratification. |
| `GSE327404` | 4 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded: insufficient donor biological replication (n_samples=4 < 8) for stage-resolved pseudobulk drift or pre-treatment responder stratification. |
| `GSE309595` | 270 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
| `GSE309593` | 91 | `solid_tumor_anti_pd1_scrna` | [solid_tumor_anti_pd1_scrna] Excluded during full-text metadata screening: lacks paired stage-resolved progression annotations or strictly pre-treatment clinical response/survival endpoints without regimen confounding. |
