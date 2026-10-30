import hashlib
import json
from pathlib import Path
import pytest

from omniniche.core.schema import (
    BarrierDeliverySpec,
    ClinicalQuerySpec,
    CompartmentRule,
    CrossOrganTransferMap,
    DirectedEdgeRule,
    DiseasePack,
    MatrixProvenance,
    NicheGrammar,
    OrganPuritySpec,
    TrialRecord,
    WorkspaceMode,
)
from omniniche.core.integrity_firewall import (
    FirewallViolationError,
    IntegrityFirewall,
)


def _sample_raw_ctg_study(nct_id: str = "NCT04567890", title: str = "Phase 3 Randomized Trial of DrugX in Sepsis") -> dict:
    return {
        "protocolSection": {
            "identificationModule": {
                "nctId": nct_id,
                "officialTitle": title,
                "briefTitle": "DrugX in Sepsis",
            },
            "statusModule": {
                "overallStatus": "TERMINATED",
                "whyStopped": "Lack of efficacy at interim futility analysis",
            },
            "designModule": {
                "studyType": "INTERVENTIONAL",
                "phases": ["PHASE3"],
                "enrollmentInfo": {"count": 420},
            },
            "armsInterventionsModule": {
                "interventions": [{"type": "DRUG", "name": "DrugX"}]
            },
            "conditionsModule": {
                "conditions": ["Septic Shock", "ARDS"]
            },
        }
    }


def test_firewall1_accepts_verbatim_trial_and_rejects_tampered_title():
    raw_study = _sample_raw_ctg_study()
    record = TrialRecord.from_ctg_api_v2(raw_study)
    assert record.nct_id == "NCT04567890"
    assert record.raw_payload_sha256 == hashlib.sha256(
        json.dumps(raw_study, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()

    # Verbatim check passes
    IntegrityFirewall.assert_trial_verbatim_integrity(record, raw_study)

    # Tampering with official_title (simulating regex title rewriting) MUST raise FirewallViolationError
    tampered = record.model_copy(update={"official_title": "Rewritten Fake Sepsis Title"})
    with pytest.raises(FirewallViolationError, match="Firewall-1"):
        IntegrityFirewall.assert_trial_verbatim_integrity(tampered, raw_study)


def test_firewall2_enforces_real_matrix_sha256_and_disease_ontology_match(tmp_path: Path):
    matrix_file = tmp_path / "GSE999999_sepsis_counts.tsv"
    matrix_file.write_text("gene\tcell1\tcell2\nTLR4\t5\t8\n", encoding="utf-8")

    prov = MatrixProvenance.from_file(
        accession="GSE999999",
        file_path=matrix_file,
        organism="Homo sapiens",
        tissue="Blood / BALF",
        condition_summary="Human septic shock longitudinal PBMC single-cell atlas",
        n_donors=24,
        n_cells_or_samples=48500,
        derivation_method="raw_matrix_direct_computation",
    )
    assert len(prov.source_matrix_sha256) == 64

    pack = DiseasePack(
        pack_id="sepsis_ards",
        display_name="Sepsis and ARDS",
        archetype="hyperacute_biphasic_flip",
        organ_system="Systemic_Vascular_Lung",
        time_scale="hours_to_days",
        clinical_query_spec=ClinicalQuerySpec(
            condition_terms=["Sepsis", "Septic Shock", "ARDS"],
            exclusion_terms=["Traumatic Brain Injury", "Glioblastoma"],
            stages=["t0_baseline", "t1_hyperacute_storm", "t2_vascular_leak", "t3_immunoparalysis"],
        ),
        organ_purity_spec=OrganPuritySpec(
            ambient_contamination_genes=["HBA1", "HBB", "PPBP"],
            dissociation_stress_genes=["HSPA1A", "FOS", "JUN"],
            vital_organ_toxicity_genes=["KCNH2", "MYH6"],
        ),
        barrier_delivery_spec=BarrierDeliverySpec(
            barrier_type="endothelial_glycocalyx_microthrombus_biphasic",
            stage_permeability={
                "small_molecule": [0.95, 0.92, 0.85, 0.72],
                "mAb": [0.45, 0.82, 0.60, 0.32],
                "ADC": [0.40, 0.75, 0.55, 0.28],
            },
        ),
        cell_ontology_map={
            "Initiator_Myeloid": ["CD14 Monocyte", "MS1 Monocyte", "Macrophage"],
            "Amplifier_Neutrophil": ["Neutrophil", "NET_Neutrophil"],
            "Barrier_Endothelium": ["Vascular Endothelial", "Capillary Aerocyte"],
            "Target_Parenchyma": ["AT1", "AT2", "Alveolar Epithelial"],
        },
    )

    # Matching provenance passes
    IntegrityFirewall.assert_matrix_provenance_and_ontology(prov, pack)

    # Formula-interpolated synthetic column MUST be blocked
    synthetic_prov = prov.model_copy(update={"derivation_method": "linear_weighted_interpolation"})
    with pytest.raises(FirewallViolationError, match="Firewall-2"):
        IntegrityFirewall.assert_matrix_provenance_and_ontology(synthetic_prov, pack)

    # Off-target excluded disease dataset (e.g., TBI dataset leaking into Sepsis pack) MUST be blocked
    offtarget_prov = prov.model_copy(update={"condition_summary": "Traumatic Brain Injury hippocampal cortex"})
    with pytest.raises(FirewallViolationError, match="Firewall-2"):
        IntegrityFirewall.assert_matrix_provenance_and_ontology(offtarget_prov, pack)


def test_firewall3_blocks_private_lab_paths_in_public_benchmark_mode(tmp_path: Path):
    private_path = Path(r"D:\AI工具\pr8_vgx_reanalysis\data\wetlab_il6.csv")
    neuro_archive_path = Path(r"D:\AI工具\NeuroImmune_Virtual_Biotech\flagship_paper\submission_package_v1.0\TableS1.csv")

    with pytest.raises(FirewallViolationError, match="Firewall-3"):
        IntegrityFirewall.assert_path_allowed(private_path, WorkspaceMode.PUBLIC_BENCHMARK)

    with pytest.raises(FirewallViolationError, match="Firewall-3"):
        IntegrityFirewall.assert_path_allowed(neuro_archive_path, WorkspaceMode.PUBLIC_BENCHMARK)

    # Clean public directory inside OmniNiche_Universal is allowed
    clean_public = tmp_path / "public_geo_cache" / "GSE12345.h5ad"
    IntegrityFirewall.assert_path_allowed(clean_public, WorkspaceMode.PUBLIC_BENCHMARK)

    # Private lab path is allowed ONLY in LAB_PRIVATE mode (except frozen nBME-26-6122 archive write check)
    IntegrityFirewall.assert_path_allowed(private_path, WorkspaceMode.LAB_PRIVATE)


def test_niche_grammar_and_cross_organ_transfer_schema_validation():
    grammar = NicheGrammar(
        grammar_id="fibrosis_core_v2",
        source_pack_id="fibrosis_ckm_ipf",
        organ_context="Kidney_Cortex",
        compartments=[
            CompartmentRule(name="Injured_Epithelium", baseline_production=0.05, decay_rate=0.35, noise_sigma=0.04, marker_genes=["HAVCR1", "VCAM1"]),
            CompartmentRule(name="Profibrotic_Macrophage", baseline_production=0.03, decay_rate=0.30, noise_sigma=0.04, marker_genes=["SPP1", "TREM2"]),
            CompartmentRule(name="Myofibroblast_ECM", baseline_production=0.02, decay_rate=0.18, noise_sigma=0.03, marker_genes=["POSTN", "COL1A1", "ACTA2"]),
        ],
        edges=[
            DirectedEdgeRule(
                source="Injured_Epithelium",
                target="Profibrotic_Macrophage",
                interaction_type="activation",
                weight_v=0.85,
                half_max_k=0.40,
                hill_n=2.0,
                ligand_gene="CSF1",
                receptor_gene="CSF1R",
                effector_tf="SPI1",
                partial_corr_gain=0.62,
            ),
            DirectedEdgeRule(
                source="Profibrotic_Macrophage",
                target="Myofibroblast_ECM",
                interaction_type="activation",
                weight_v=1.10,
                half_max_k=0.38,
                hill_n=3.0,
                ligand_gene="TGFB1",
                receptor_gene="TGFBR2",
                effector_tf="SMAD3",
                partial_corr_gain=0.71,
            ),
            DirectedEdgeRule(
                source="Myofibroblast_ECM",
                target="Injured_Epithelium",
                interaction_type="activation",
                weight_v=0.75,
                half_max_k=0.45,
                hill_n=3.0,
                ligand_gene="FN1",
                receptor_gene="ITGB6",
                effector_tf="YAP1",
                partial_corr_gain=0.58,
            ),
        ],
        remission_compartment="Myofibroblast_ECM",
        remission_threshold=0.25,
    )
    assert len(grammar.compartments) == 3
    assert len(grammar.edges) == 3

    tmap = CrossOrganTransferMap(
        source_organ="Kidney_Cortex",
        target_organ="Lung_Alveolar_Interstitium",
        compartment_alignment={
            "Injured_Epithelium": "Aberrant_Basaloid_AT2",
            "Profibrotic_Macrophage": "SPP1_Alveolar_Macrophage",
            "Myofibroblast_ECM": "CTHRC1_Myofibroblast",
        },
        gene_ortholog_or_functional_swap={
            "HAVCR1": "KRT17",
            "VCAM1": "GDF15",
        },
    )
    assert tmap.compartment_alignment["Injured_Epithelium"] == "Aberrant_Basaloid_AT2"
