from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Tuple

from omniniche.core.schema import (
    DiseasePack,
    MatrixProvenance,
    TrialRecord,
    WorkspaceMode,
)


class FirewallViolationError(RuntimeError):
    """Raised whenever any of OmniNiche's three data-authenticity or isolation firewalls is breached."""


class IntegrityFirewall:
    """Runtime hard assertions enforcing clinical-trial verbatim integrity, raw-matrix provenance,
    disease-ontology matching, and strict separation between public benchmarks and private lab projects.
    """

    ALLOWED_DERIVATION_METHODS: Tuple[str, ...] = (
        "raw_matrix_direct_computation",
        "pseudobulk_donor_aggregation",
        "backed_sparse_h5ad_direct",
    )

    FORBIDDEN_DERIVATION_KEYWORDS: Tuple[str, ...] = (
        "interpolat",
        "synthetic",
        "mock",
        "linear_weighted",
        "simulated_column",
    )

    PRIVATE_OR_FROZEN_PATH_TOKENS: Tuple[str, ...] = (
        "pr8_vgx_reanalysis",
        "neuroimmune_virtual_biotech",
        "nbme-26-6122",
        "submission_package_v1.0",
    )

    FROZEN_ARCHIVE_TOKENS: Tuple[str, ...] = (
        "neuroimmune_virtual_biotech\\flagship_paper\\submission_package_v1.0",
        "neuroimmune_virtual_biotech/flagship_paper/submission_package_v1.0",
        "neuroimmune_virtual_biotech_投稿包_v1.0",
    )

    @classmethod
    def assert_trial_verbatim_integrity(cls, record: TrialRecord, raw_ctg_study: dict) -> None:
        """Firewall-1: Guarantee 100% character-for-character fidelity against ClinicalTrials.gov API v2 JSON."""
        canonical = json.dumps(raw_ctg_study, sort_keys=True, ensure_ascii=False).encode("utf-8")
        expected_hash = hashlib.sha256(canonical).hexdigest()
        if record.raw_payload_sha256 != expected_hash:
            raise FirewallViolationError(
                f"Firewall-1 breach for {record.nct_id}: SHA-256 mismatch against raw ClinicalTrials.gov payload."
            )

        ref = TrialRecord.from_ctg_api_v2(raw_ctg_study)
        for field_name in ("nct_id", "official_title", "brief_title", "overall_status", "why_stopped", "enrollment"):
            actual_val = getattr(record, field_name)
            expected_val = getattr(ref, field_name)
            if actual_val != expected_val:
                raise FirewallViolationError(
                    f"Firewall-1 breach for {record.nct_id}: field '{field_name}' was altered "
                    f"(expected {expected_val!r}, got {actual_val!r}). Regex or manual rewriting is forbidden."
                )

    @classmethod
    def assert_matrix_provenance_and_ontology(cls, prov: MatrixProvenance, pack: DiseasePack) -> None:
        """Firewall-2: Block formula-interpolated columns, missing SHA-256, and off-target disease datasets."""
        if len(prov.source_matrix_sha256) != 64:
            raise FirewallViolationError(
                f"Firewall-2 breach for {prov.accession}: invalid 64-hex SHA-256 signature."
            )

        method_lower = prov.derivation_method.lower()
        if any(bad in method_lower for bad in cls.FORBIDDEN_DERIVATION_KEYWORDS) or (
            prov.derivation_method not in cls.ALLOWED_DERIVATION_METHODS
        ):
            raise FirewallViolationError(
                f"Firewall-2 breach for {prov.accession}: derivation_method={prov.derivation_method!r} is forbidden. "
                f"Only direct raw/pseudobulk matrix computations {cls.ALLOWED_DERIVATION_METHODS} are permitted."
            )

        summary_lower = f"{prov.tissue} {prov.condition_summary}".lower()
        for excl in pack.clinical_query_spec.exclusion_terms:
            if excl.lower() in summary_lower:
                raise FirewallViolationError(
                    f"Firewall-2 breach for {prov.accession}: dataset matches excluded ontology term "
                    f"{excl!r} in pack '{pack.pack_id}'."
                )

    @classmethod
    def assert_path_allowed(cls, path: Path | str, mode: WorkspaceMode, is_write: bool = False) -> None:
        """Firewall-3: Physically isolate public benchmarks from private lab data and frozen archives."""
        norm = str(path).replace("/", "\\").lower()

        # Never allow writing to frozen nBME-26-6122 submission archives in ANY mode
        if is_write and any(tok in norm for tok in cls.FROZEN_ARCHIVE_TOKENS):
            raise FirewallViolationError(
                f"Firewall-3 breach: write attempt to frozen nBME-26-6122 archive path {path!s} is strictly forbidden."
            )

        if mode == WorkspaceMode.PUBLIC_BENCHMARK:
            for tok in cls.PRIVATE_OR_FROZEN_PATH_TOKENS:
                if tok in norm:
                    raise FirewallViolationError(
                        f"Firewall-3 breach: path {path!s} contains restricted private/prior-project token "
                        f"{tok!r} while running in PUBLIC_BENCHMARK mode."
                    )
