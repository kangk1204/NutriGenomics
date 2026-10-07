"""Historical phase-one implementation guardrails, not verified final-contract targets."""
from __future__ import annotations

import math

PROTOCOL_ID = "phase1-final-20261003-v1"
# Each guard is an assertion backed by a checksum-frozen artifact. Missing is pending.
TARGETS = {
    "standardized_records": ("표준화 DB", ">=", 1_000_000, "deduplicated_records", ("loaded", "deduplicated", "units_audited")),
    "new_food_records": ("추가 식품 정보", ">=", 90_000, "unique_food_records", ("loaded", "baseline_delta", "deduplicated")),
    "nutriomics_database": ("NutriOmics DB", ">=", 1, "database", ("query_verified", "four_modalities")),
    "genomics": ("유전체 자료", ">=", 1, "dataset", ("actual_data", "modality_verified", "sample_mapping")),
    "epigenomics": ("후성유전체 자료", ">=", 1, "dataset", ("actual_data", "modality_verified", "sample_mapping")),
    "ffq": ("FFQ 자료", ">=", 1, "dataset", ("actual_data", "modality_verified", "sample_mapping")),
    "metabolomics": ("대사체 자료", ">=", 1, "dataset", ("actual_data", "modality_verified", "sample_mapping")),
    "omics_pipeline": ("오믹스 표준화 방법론", ">=", 1, "pipeline", ("raw_to_results", "independent_rerun")),
    "chronic_classifier": ("만성질환 판별 AI", ">", .8, "AUROC", ("held_out", "leakage_audited", "predictions_saved")),
    "evidence_recommender": ("문헌 근거 성분추천 AI", ">=", .8, "AUROC", ("sealed_gold", "document_overlap_audited", "predictions_saved", "endpoint_defined")),
    "measured_interactions": ("화합물–단백질 상호작용", ">=", 1_680_000, "measured_interactions", ("loaded", "deduplicated", "endpoint_separated", "units_audited")),
    "natural_products": ("천연물", ">=", 400_000, "unique_structures", ("loaded", "structure_verified", "deduplicated")),
    "food_compounds": ("단일 식품 성분", ">=", 2_000, "unique_food_compounds", ("loaded", "structure_verified", "food_occurrence_verified", "deduplicated")),
    "android": ("Android 앱", ">=", 1, "application", ("actual_build", "core_flow_verified", "version_saved")),
    "ios": ("iOS 앱", ">=", 1, "application", ("actual_build", "core_flow_verified", "version_saved")),
    "web": ("웹 서비스", ">=", 1, "application", ("actual_build", "core_flow_verified", "version_saved")),
    "nafld_candidates": ("지방간 후보", ">=", 17, "candidate", ("identity_verified", "evidence_cards", "uncertainty_recorded")),
    "immune_candidates": ("면역 관련 후보", ">=", 6, "candidate", ("identity_verified", "evidence_cards", "uncertainty_recorded")),
    "scie_paper": ("SCIE 논문", ">=", 1, "paper", ("bibliography_verified", "acknowledgement_verified", "counted_once")),
    "literature_precision": ("유전자 탐지 정밀도", ">=", .90, "precision", ("held_out", "predictions_saved")),
    "literature_recall": ("유전자 탐지 재현율", ">=", .70, "recall", ("held_out", "predictions_saved")),
    "relation_precision": ("관계 추출 정밀도", ">=", .85, "precision", ("held_out", "predictions_saved", "endpoint_compatible")),
    "relation_recall": ("관계 추출 재현율", ">=", .30, "recall", ("held_out", "predictions_saved", "endpoint_compatible")),
    "relation_precision_gain": ("관계 추출 정밀도 개선", ">=", .05, "absolute_gain", ("same_sealed_set", "endpoint_compatible")),
    "relation_recall_gain": ("관계 추출 재현율 비감소", ">=", 0., "absolute_gain", ("same_sealed_set", "endpoint_compatible")),
    "scope_abstract": ("초록 부정·추정 범위", ">=", .85, "token_F1", ("held_out", "predictions_saved")),
    "scope_fulltext": ("원문 부정·추정 범위", ">=", .85, "token_F1", ("held_out", "predictions_saved")),
    "provenance": ("문헌 출처 보존", ">=", 1., "fraction", ("audited",)),
    "operation_days": ("문헌 시스템 실제 운영", ">=", 7, "elapsed_days", ("continuous_logs", "gap_within_900_seconds", "resume_verified", "restore_verified", "rollback_verified")),
}


def evaluate(observations: dict, verified_artifacts: set[str]) -> dict:
    """Fail closed on missing evidence, wrong units, failed guards or train-set scores."""
    unknown = set(observations) - set(TARGETS)
    if unknown:
        raise ValueError(f"Unknown acceptance criteria: {sorted(unknown)}")
    items = []
    for key, (label, operator, threshold, unit, guards) in TARGETS.items():
        observation = observations.get(key, {})
        value = observation.get("value")
        numeric = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
        proof = observation.get("artifacts", [])
        proof_valid = isinstance(proof, list) and bool(proof) and all(isinstance(p, str) and p in verified_artifacts for p in proof)
        missing = [g for g in guards if observation.get("guards", {}).get(g) is not True]
        reasons = ([] if numeric else ["observation_missing_or_nonfinite"])
        if observation.get("unit") != unit:
            reasons.append("unit_missing_or_mismatch")
        if not proof_valid:
            reasons.append("verified_evidence_missing")
        if missing:
            reasons += ["guard:" + g for g in missing]
        complete = not reasons
        reached = numeric and (value > threshold if operator == ">" else value >= threshold)
        state = "achieved" if complete and reached else "not_achieved" if complete else "pending_verification"
        items.append({"id": key, "goal": label, "target": threshold, "operator": operator,
                      "unit": unit, "observed": value, "status": state, "reasons": reasons,
                      "artifacts": proof, "note": observation.get("note", "")})
    return {"schema_version": 3, "protocol_id": PROTOCOL_ID,
            "criteria_scope": "implementation_guardrails",
            "criterion_scopes": {"food_compounds": "additional_strict_identity_subset_not_overall_component_coverage"},
            "formal_contract_target_received": False,
            "target_basis": "historical_phase1_implementation_guardrails",
            "formal_acceptance_status": "pending_contract_target_verification",
            "formal_achievement_percentage": None,
            "all_metrics_satisfied": all(i["status"] == "achieved" for i in items),
            "achieved_criteria": sum(i["status"] == "achieved" for i in items),
            "required_criteria": len(items), "items": items,
            "interpretation": "Internal guardrail results only; the strict food subset does not determine overall component coverage, and 2026 goals require separate source-specific evaluation"}
