"""Stage-11 expert registry with fail-closed scientific qualification gates."""
from __future__ import annotations
import csv
import json
from dataclasses import asdict, dataclass
from importlib.metadata import version
from typing import Any, Mapping
from .core import STAGE, file_sha256, stable_hash


@dataclass(frozen=True)
class AdapterSpec:
    algorithm: str; package: str; artifact_kind: str; ensemble_method: str; design_status: str; stage10_source: str; blocker: str = ""


_SPECS = (
    AdapterSpec("GIMME", "troppo", "context_specific_model", "common_optgp", "blocked_pending_full_qualification", "src/gimme_repair.py:fit_gimme_expert"),
    AdapterSpec("iMAT", "troppo", "context_specific_model", "common_optgp", "blocked_pending_full_qualification", "src/imat_canonical.py:fit_imat_expert", "Canonical Troppo iMAT requires current CT2A and GL261 model/smoke qualification."),
    AdapterSpec("INIT_tINIT", "troppo", "context_specific_model", "common_optgp", "primary_runnable", "src/contextualization.py:contextualize"),
    AdapterSpec("CORDA", "corda", "context_specific_model", "common_optgp", "blocked_pending_full_qualification", "src/corda_canonical.py:fit_corda_expert"),
    AdapterSpec("E-Flux", "cobra", "rna_bounded_model", "common_optgp", "blocked_pending_validation", "src/models.py:reaction_expression", "Full-parent E-Flux semantics/feasibility qualification is required."),
    AdapterSpec("RegrEx", "unavailable", "regrex_native_expert", "regrex_aos", "blocked_pending_validation", "new adapter", "No validated RegrEx-LAD/AOS implementation is available."),
    AdapterSpec("RIPTiDe", "riptide", "riptide_native_expert", "native_gapsplit", "blocked_pending_validation", "src/contextualization.py:contextualize", "Native persistence, set_bounds, and 50-to-250 reuse qualification is required."),
    AdapterSpec("FASTCORE_FASTCORMICS", "troppo", "context_specific_model", "common_optgp", "optional_secondary", "src/contextualization.py:_numpy_compatible_fastcore/contextualize"),
    AdapterSpec("mCADRE", "pymcadre", "none", "none", "retired_excluded", "src/adapters.py", "Incompatible non-RNA-only auxiliary-input/task requirements."),
    AdapterSpec("COMPASS", "compass", "none", "none", "retired_excluded", "src/adapters.py", "Does not naturally provide common parent-iMM1865 flux vectors."),
)
REGISTRY = {spec.algorithm: spec for spec in _SPECS}


_RIPTIDE_TUMORS = ("CT2A", "GL261")
_QUALIFICATION_PROTOCOL = "khalsa_all_mikolajewicz_validation"
_QUALIFICATION_FOLD = "all"
_QUALIFICATION_VARIANT = "rna_only_primary"


def _valid_marker(name: str, expected_fingerprint: str) -> bool:
    path = STAGE / "provenance" / name
    try: data = json.loads(path.read_text())
    except Exception: return False
    return bool(data.get("passed") is True and data.get("qualification_fingerprint") == expected_fingerprint and not data.get("skipped") and not data.get("xfailed"))


def expected_eflux_fingerprint() -> str:
    return stable_hash({"cobra_version": version("cobra"), "expert_fitting_sha256": file_sha256(STAGE / "src/expert_fitting.py"), "qualification_schema": 1})


def eflux_validated() -> bool:
    return _valid_marker("EFLUX_VALIDATED.json", expected_eflux_fingerprint())


def riptide_native_validated() -> bool:
    """Validate only the native API/source marker, independently of parent qualification."""
    try:
        from .riptide_native import qualification_fingerprint
        expected = qualification_fingerprint()
    except Exception:
        return False
    return _valid_marker("RIPTIDE_NATIVE_VALIDATED.json", expected)


def _current_manifest_reference(algorithm: str, tumor: str) -> dict[str, Any] | None:
    """Read the exact current administrative qualification reference, if available.

    This is used only for status queries outside manifest construction.  During
    manifest construction the in-memory current reference is passed explicitly,
    so an old on-disk manifest can never authorize a newly rebuilt manifest.
    """
    path = STAGE / "manifests/expert_references.tsv"
    try:
        with path.open(newline="") as handle:
            rows = [
                dict(row) for row in csv.DictReader(handle, delimiter="\t")
                if row.get("algorithm") == algorithm
                and row.get("tumor") == tumor
                and row.get("protocol") == _QUALIFICATION_PROTOCOL
                and str(row.get("fold")) == _QUALIFICATION_FOLD
                and row.get("expert_variant") == _QUALIFICATION_VARIANT
            ]
    except Exception:
        return None
    if len(rows) != 1:
        return None
    row = rows[0]
    if str(row.get("dmi_mouse_id", "")).strip().lower() in {"", "nan", "none"}:
        row["dmi_mouse_id"] = None
    return row


PRIMARY_ALGORITHMS = ("GIMME", "iMAT", "CORDA", "RIPTiDe")
PRIMARY_TUMORS = ("CT2A", "GL261")


def _qualification_ensemble(reference: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    from .cache_identity import ensemble_identity
    from .qualification import QUALIFICATION_VECTORS

    digest, seed, spec = ensemble_identity(
        str(reference["expert_hash"]), str(reference["algorithm"]), config, QUALIFICATION_VECTORS
    )
    return dict(reference) | {
        "execution_status": "planned",
        "ensemble_hash": digest,
        "ensemble_seed": seed,
        "ensemble_spec": json.dumps(spec, sort_keys=True, separators=(",", ":")),
        "requested_vectors": QUALIFICATION_VECTORS,
        "ensemble_replicate": 0,
        "sampling_thinning": int(config["sampling"]["thinning"]),
        "sampling_processes": int(config["sampling"]["processes"]),
    }


def qualification_evidence_matches(
    report: Mapping[str, Any],
    reference: Mapping[str, Any],
    metadata: Mapping[str, Any],
    persisted: Mapping[str, Any],
    *,
    fit_fingerprint: str,
    ensemble_fingerprint: str,
    sampling_fingerprint: str,
    controller_fingerprint: str,
    smoke_artifact_identity: str,
    smoke_ensemble_hash: str,
    qualification_identity: str,
) -> bool:
    """Pure identity/status check shared by all four production gates."""
    qualification = report.get("qualification_record")
    representative = report.get("representative_condition") or {}
    phases = report.get("phases") or {}
    exact_reference_keys = (
        "algorithm", "tumor", "condition_type", "context_id", "protocol", "fold",
        "expert_variant", "expert_hash", "generation_spec_identity", "rna_training_slice_hash",
    )
    return bool(
        isinstance(qualification, Mapping)
        and report.get("algorithm") == reference.get("algorithm")
        and report.get("tumor") == reference.get("tumor")
        and report.get("status") == "qualified_full_parent"
        and report.get("production") is False
        and report.get("all_workers_reaped") is True
        and report.get("vectors") == 2
        and all(phases.get(name, {}).get("status") == "passed" for name in ("fit", "ensemble", "projection", "reuse"))
        and report.get("algorithm_fit_fingerprint") == fit_fingerprint
        and report.get("ensemble_implementation_fingerprint") == ensemble_fingerprint
        and report.get("sampling_semantics_fingerprint") == sampling_fingerprint
        and report.get("qualification_implementation_fingerprint") == controller_fingerprint
        and report.get("qualification_smoke_artifact_identity") == smoke_artifact_identity
        and representative.get("tumor") == reference.get("tumor")
        and representative.get("expert_hash") == reference.get("expert_hash")
        and representative.get("ensemble_hash") == smoke_ensemble_hash
        and all(str(qualification.get(key)) == str(reference.get(key)) for key in exact_reference_keys)
        and _optional_reference(qualification.get("dmi_mouse_id")) == _optional_reference(reference.get("dmi_mouse_id"))
        and qualification.get("qualification_status") == "passed"
        and qualification.get("artifact_model_hash") == metadata.get("artifact_model_hash")
        and qualification.get("qualification_identity") == qualification_identity
        and persisted.get("qualification_identity") == qualification_identity
        and persisted.get("artifact_model_hash") == metadata.get("artifact_model_hash")
        and metadata.get("generation_spec_identity") == reference.get("generation_spec_identity")
    )


def _optional_reference(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return None if text.lower() in {"", "nan", "none"} else text


def primary_method_validated(
    algorithm: str,
    tumor: str,
    *,
    expected_record: Mapping[str, Any] | None = None,
) -> bool:
    """Validate one exact current model+smoke qualification, without production artifacts."""
    if algorithm not in PRIMARY_ALGORITHMS or tumor not in PRIMARY_TUMORS:
        return False
    if expected_record is None:
        expected_record = _current_manifest_reference(algorithm, tumor)
    if expected_record is None:
        return False
    if (str(expected_record.get("algorithm")) != algorithm
            or str(expected_record.get("tumor")) != tumor):
        return False
    try:
        from .cache_identity import (
            algorithm_fit_fingerprint, ensemble_implementation_fingerprint,
            qualification_identity, sampling_semantics_fingerprint,
        )
        from .core import load_configuration
        from .ensembles import valid_ensemble
        from .expert_artifacts import expert_directory, valid_expert
        from .projections import valid_projection
        from .qualification import (
            currently_qualified_expert,
            full_gem_qualification_fingerprint, projection_record,
            qualification_matches_reference, qualification_smoke_artifact_identity,
        )

        paths, config = load_configuration()
        report = json.loads(
            (STAGE / f"qc/full_gem_qualification/{algorithm}/{tumor}/qualification.json").read_text()
        )
        qualification = report.get("qualification_record")
        ensemble = _qualification_ensemble(expected_record, config)
        projection, _ = projection_record(ensemble, paths)
        metadata = json.loads((expert_directory(expected_record) / "metadata.json").read_text())
        persisted = json.loads(
            (expert_directory(expected_record) / "qc/model_qualification.json").read_text()
        )
        smoke_identity, _ = qualification_smoke_artifact_identity(
            expected_record, ensemble, paths, config
        )
        expected_qualification_identity, _ = qualification_identity(
            str(expected_record["generation_spec_identity"]),
            str(metadata["artifact_model_hash"]),
            config,
            qualification["expected_source_routes"],
            algorithm=algorithm,
        )
    except Exception:
        return False

    valid = bool(
        qualification_evidence_matches(
            report, expected_record, metadata, persisted,
            fit_fingerprint=algorithm_fit_fingerprint(algorithm),
            ensemble_fingerprint=ensemble_implementation_fingerprint(algorithm),
            sampling_fingerprint=sampling_semantics_fingerprint(algorithm, config),
            controller_fingerprint=full_gem_qualification_fingerprint(),
            smoke_artifact_identity=smoke_identity,
            smoke_ensemble_hash=str(ensemble["ensemble_hash"]),
            qualification_identity=expected_qualification_identity,
        )
        and qualification_matches_reference(expected_record, qualification)
        and currently_qualified_expert(qualification, qualification)
        and valid_expert(expected_record)
        and valid_ensemble(ensemble)
        and valid_projection(projection)
    )
    if algorithm == "RIPTiDe":
        try:
            from .riptide_native import qualification_fingerprint
            valid = bool(
                valid
                and riptide_native_validated()
                and report.get("native_qualification_fingerprint") == qualification_fingerprint()
            )
        except Exception:
            return False
    return valid


def gimme_validated(tumor: str = "CT2A", *, expected_record: Mapping[str, Any] | None = None) -> bool:
    return primary_method_validated("GIMME", tumor, expected_record=expected_record)


def imat_validated(tumor: str = "CT2A", *, expected_record: Mapping[str, Any] | None = None) -> bool:
    return primary_method_validated("iMAT", tumor, expected_record=expected_record)


def corda_validated(tumor: str = "CT2A", *, expected_record: Mapping[str, Any] | None = None) -> bool:
    return primary_method_validated("CORDA", tumor, expected_record=expected_record)


def riptide_validated(tumor: str = "CT2A", *, expected_record: Mapping[str, Any] | None = None) -> bool:
    return primary_method_validated("RIPTiDe", tumor, expected_record=expected_record)


def _all_current_references() -> dict[tuple[str, str], Mapping[str, Any]] | None:
    references = {
        (algorithm, tumor): _current_manifest_reference(algorithm, tumor)
        for algorithm in PRIMARY_ALGORITHMS for tumor in PRIMARY_TUMORS
    }
    return None if any(value is None for value in references.values()) else references


def four_method_barrier_validated(
    expected_records: Mapping[tuple[str, str], Mapping[str, Any] | None] | None = None,
) -> bool:
    """Require all eight exact current qualifications before production eligibility."""
    records = expected_records if expected_records is not None else _all_current_references()
    if records is None:
        return False
    for algorithm in PRIMARY_ALGORITHMS:
        for tumor in PRIMARY_TUMORS:
            reference = records.get((algorithm, tumor))
            if reference is None or not primary_method_validated(
                algorithm, tumor, expected_record=reference
            ):
                return False
    return True


def full_parent_rescue_validated(algorithm: str) -> bool:
    """Require a successful current rescue before failed secondary methods run."""
    try:
        report = json.loads((STAGE / f"qc/rescue_qualification/{algorithm}/qualification.json").read_text())
    except Exception:
        return False
    return bool(
        report.get("algorithm") == algorithm
        and report.get("status") == "qualified_full_parent"
        and report.get("all_workers_reaped") is True
        and report.get("expert_hash")
    )


def effective_status(
    algorithm: str,
    *,
    tumor: str | None = None,
    expected_record: Mapping[str, Any] | None = None,
    global_authorized: bool | None = None,
) -> str:
    spec = REGISTRY[algorithm]
    if algorithm in PRIMARY_ALGORITHMS:
        barrier = four_method_barrier_validated() if global_authorized is None else bool(global_authorized)
        case_valid = (
            primary_method_validated(algorithm, tumor, expected_record=expected_record)
            if tumor is not None else all(primary_method_validated(algorithm, item) for item in PRIMARY_TUMORS)
        )
        return "primary_runnable" if barrier and case_valid else "blocked_pending_full_qualification"
    if algorithm == "INIT_tINIT":
        return "primary_runnable" if full_parent_rescue_validated(algorithm) else "blocked_pending_validation"
    if algorithm == "E-Flux":
        return "blocked_pending_scientific_zero_policy"
    return spec.design_status


def adapter_records() -> list[dict[str, object]]:
    barrier = four_method_barrier_validated()
    return [
        asdict(spec) | {"effective_status": effective_status(spec.algorithm, global_authorized=barrier)}
        for spec in _SPECS
    ]
