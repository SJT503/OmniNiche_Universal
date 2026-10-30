from __future__ import annotations

from pathlib import Path
from typing import List
import numpy as np
import yaml

from omniniche.core.schema import (
    CompartmentRule,
    CrossOrganTransferMap,
    DirectedEdgeRule,
    DiseasePack,
    NicheGrammar,
)

PACKS_DIR = Path(__file__).resolve().parent


def list_builtin_disease_packs() -> List[str]:
    return sorted(p.stem for p in PACKS_DIR.glob("*.yaml"))


def load_disease_pack(pack_id_or_path: str | Path) -> DiseasePack:
    p = Path(pack_id_or_path)
    if p.suffix.lower() in (".yaml", ".yml") and p.exists():
        yaml_path = p
    else:
        yaml_path = PACKS_DIR / f"{pack_id_or_path}.yaml"
        if not yaml_path.exists():
            raise FileNotFoundError(
                f"DiseasePack {pack_id_or_path!r} not found. Available builtin packs: {list_builtin_disease_packs()}"
            )
    raw = yaml.safe_load(yaml_path.read_text(encoding="utf-8"))
    return DiseasePack.model_validate(raw)


def compute_barrier_access(
    pack: DiseasePack,
    modality: str,
    normalized_stage_progress: float,
) -> float:
    """Interpolate organ-specific pharmacokinetic/barrier delivery accessibility B(m, t) in [0, 1]."""
    curves = pack.barrier_delivery_spec.stage_permeability
    curve = curves.get(modality) or curves.get("small_molecule")
    if not curve:
        return 0.75
    prog = float(np.clip(normalized_stage_progress, 0.0, 1.0))
    if len(curve) == 1:
        return float(curve[0])
    xs = np.linspace(0.0, 1.0, len(curve))
    return float(np.interp(prog, xs, curve))


def transfer_grammar_cross_organ(
    source_grammar: NicheGrammar,
    transfer_map: CrossOrganTransferMap,
) -> NicheGrammar:
    """Zero-shot transfer of a compiled multicellular NicheGrammar from source_organ to target_organ
    (e.g., Kidney AKI-to-CKD fibrosis grammar -> Lung IPF alveolar-interstitial fibrosis grammar).
    """
    cmap = transfer_map.compartment_alignment
    gswap = transfer_map.gene_ortholog_or_functional_swap

    def _map_comp(name: str) -> str:
        return cmap.get(name, name)

    def _map_gene(gene: str) -> str:
        return gswap.get(gene, gene)

    new_compartments: List[CompartmentRule] = []
    for comp in source_grammar.compartments:
        new_compartments.append(
            comp.model_copy(
                update={
                    "name": _map_comp(comp.name),
                    "marker_genes": [_map_gene(g) for g in comp.marker_genes],
                }
            )
        )

    new_edges: List[DirectedEdgeRule] = []
    for edge in source_grammar.edges:
        new_edges.append(
            edge.model_copy(
                update={
                    "source": _map_comp(edge.source),
                    "target": _map_comp(edge.target),
                    "ligand_gene": _map_gene(edge.ligand_gene),
                    "receptor_gene": _map_gene(edge.receptor_gene),
                    "effector_tf": _map_gene(edge.effector_tf),
                }
            )
        )

    return NicheGrammar(
        grammar_id=f"{source_grammar.grammar_id}__transferred_to_{transfer_map.target_organ.lower()}",
        source_pack_id=source_grammar.source_pack_id,
        organ_context=transfer_map.target_organ,
        compartments=new_compartments,
        edges=new_edges,
        remission_compartment=_map_comp(source_grammar.remission_compartment),
        remission_threshold=source_grammar.remission_threshold,
        bic_score=source_grammar.bic_score,
    )
