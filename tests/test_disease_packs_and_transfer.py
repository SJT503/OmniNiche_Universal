from pathlib import Path
import pytest

from omniniche.core.schema import (
    CompartmentRule,
    DirectedEdgeRule,
    NicheGrammar,
    WorkspaceMode,
)
from omniniche.disease_packs.loader import (
    compute_barrier_access,
    list_builtin_disease_packs,
    load_disease_pack,
    transfer_grammar_cross_organ,
)
from omniniche.disease_packs.autopack_agent import AutoPackAgent


EXPECTED_FIVE_PACKS = {
    "sepsis_ards",
    "fibrosis_ckm_ipf",
    "solidtumor_time",
    "autoimmune_ibd_ra",
    "cardiometabolic_mash",
}


def test_five_orthogonal_builtin_packs_load_with_distinct_archetypes():
    available = set(list_builtin_disease_packs())
    assert EXPECTED_FIVE_PACKS.issubset(available)

    archetypes = set()
    for pid in EXPECTED_FIVE_PACKS:
        pack = load_disease_pack(pid)
        assert pack.pack_id == pid
        assert len(pack.clinical_query_spec.stages) >= 4
        assert len(pack.organ_purity_spec.ambient_contamination_genes) >= 3
        assert len(pack.cell_ontology_map) >= 4
        # Ensure strict zero CNS/TBI overlap in builtin benchmark packs
        joined_terms = " ".join(pack.clinical_query_spec.condition_terms).lower()
        assert "traumatic brain injury" not in joined_terms
        assert "epilepsy" not in joined_terms
        archetypes.add(pack.archetype)

    assert len(archetypes) == 5


def test_barrier_accessibility_reflects_organ_specific_biophysics():
    tumor_pack = load_disease_pack("solidtumor_time")
    # Early tumor vs late desmoplastic CAF stroma: mAb drops strongly, small_molecule / ADC retains higher access
    mab_early = compute_barrier_access(tumor_pack, modality="mAb", normalized_stage_progress=0.1)
    mab_late = compute_barrier_access(tumor_pack, modality="mAb", normalized_stage_progress=0.95)
    adc_late = compute_barrier_access(tumor_pack, modality="ADC", normalized_stage_progress=0.95)
    sm_late = compute_barrier_access(tumor_pack, modality="small_molecule", normalized_stage_progress=0.95)

    assert mab_early > mab_late
    assert adc_late > mab_late
    assert sm_late > mab_late


def test_cross_organ_zero_shot_fibrosis_grammar_transfer():
    fib_pack = load_disease_pack("fibrosis_ckm_ipf")
    assert fib_pack.cross_organ_transfer is not None

    kidney_grammar = NicheGrammar(
        grammar_id="kidney_aki_ckd_v2",
        source_pack_id="fibrosis_ckm_ipf",
        organ_context="Kidney_Cortex",
        compartments=[
            CompartmentRule(name="Injured_Epithelium", marker_genes=["HAVCR1", "VCAM1"]),
            CompartmentRule(name="Profibrotic_Macrophage", marker_genes=["SPP1", "TREM2"]),
            CompartmentRule(name="Myofibroblast_ECM", marker_genes=["POSTN", "COL1A1", "ACTA2"]),
            CompartmentRule(name="Microvascular_Endothelium", marker_genes=["PLVAP", "EMCN"]),
        ],
        edges=[
            DirectedEdgeRule(
                source="Injured_Epithelium",
                target="Profibrotic_Macrophage",
                interaction_type="activation",
                weight_v=0.90,
                half_max_k=0.38,
                hill_n=2.0,
                ligand_gene="CSF1",
                receptor_gene="CSF1R",
                effector_tf="SPI1",
            ),
            DirectedEdgeRule(
                source="Profibrotic_Macrophage",
                target="Myofibroblast_ECM",
                interaction_type="activation",
                weight_v=1.15,
                half_max_k=0.35,
                hill_n=3.0,
                ligand_gene="TGFB1",
                receptor_gene="TGFBR2",
                effector_tf="SMAD3",
            ),
            DirectedEdgeRule(
                source="Myofibroblast_ECM",
                target="Injured_Epithelium",
                interaction_type="activation",
                weight_v=0.80,
                half_max_k=0.42,
                hill_n=3.0,
                ligand_gene="FN1",
                receptor_gene="ITGB6",
                effector_tf="YAP1",
            ),
        ],
        remission_compartment="Myofibroblast_ECM",
        remission_threshold=0.28,
    )

    ipf_grammar = transfer_grammar_cross_organ(kidney_grammar, fib_pack.cross_organ_transfer)
    assert ipf_grammar.organ_context == "Lung_Alveolar_Interstitium"
    comp_names = [c.name for c in ipf_grammar.compartments]
    assert "Aberrant_Basaloid_AT2" in comp_names
    assert "CTHRC1_Myofibroblast" in comp_names
    assert ipf_grammar.remission_compartment == "CTHRC1_Myofibroblast"
    # Verify epithelial marker swap (HAVCR1 -> KRT17) while preserving dynamical Hill parameters
    ab_comp = next(c for c in ipf_grammar.compartments if c.name == "Aberrant_Basaloid_AT2")
    assert "KRT17" in ab_comp.marker_genes
    assert ipf_grammar.edges[1].weight_v == pytest.approx(1.15)
    assert ipf_grammar.edges[1].hill_n == pytest.approx(3.0)


def test_autopack_agent_generates_valid_new_disease_yaml(tmp_path: Path):
    out_yaml = tmp_path / "systemic_lupus_erythematosus.yaml"
    agent = AutoPackAgent()
    pack = agent.generate_disease_pack(
        disease_name="Systemic Lupus Erythematosus",
        organ_hint="Kidney_Blood_Immune",
        mode=WorkspaceMode.LAB_PRIVATE,
        output_yaml_path=out_yaml,
        use_network=False,
    )
    assert out_yaml.exists()
    assert pack.pack_id == "systemic_lupus_erythematosus"
    assert len(pack.organ_purity_spec.ambient_contamination_genes) >= 3
    reloaded = load_disease_pack(out_yaml)
    assert reloaded.display_name == pack.display_name
