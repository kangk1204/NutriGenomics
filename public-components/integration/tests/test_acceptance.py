import copy
import hashlib
import json

import pytest

from nutriomics_final.acceptance import TARGETS, evaluate
from nutriomics_final.research_registry import ResearchRegistry


def observed(key, value):
    _, _, _, unit, guards = TARGETS[key]
    return {"value": value, "unit": unit, "artifacts": ["receipt"], "guards": dict.fromkeys(guards, True)}


def test_auc_boundary_and_no_training_score_acceptance():
    obs = {"chronic_classifier": observed("chronic_classifier", .8),
           "evidence_recommender": observed("evidence_recommender", .8)}
    result = {x['id']: x for x in evaluate(obs, {"receipt"})['items']}
    assert result['chronic_classifier']['status'] == 'not_achieved'
    assert result['evidence_recommender']['status'] == 'achieved'
    obs['chronic_classifier']['value'] = .99
    obs['chronic_classifier']['guards']['held_out'] = False
    assert evaluate(obs, {"receipt"})['all_metrics_satisfied'] is False
    assert next(x for x in evaluate(obs, {"receipt"})['items'] if x['id'] == 'chronic_classifier')['status'] == 'pending_verification'


def test_missing_proof_units_and_natural_food_distinction():
    natural = observed('natural_products', 500_000)
    assert next(x for x in evaluate({'natural_products': natural}, {'receipt'})['items'] if x['id']=='natural_products')['status']=='achieved'
    for artifacts, value in [(set(), 500_000), ({'receipt'}, float('nan')), ({'receipt'}, True)]:
        bad = copy.deepcopy(natural)
        bad['value'] = value
        assert not evaluate({'natural_products': bad}, artifacts)['all_metrics_satisfied']
    food = observed('food_compounds', 500_000)
    food['guards']['food_occurrence_verified'] = False
    assert next(x for x in evaluate({'food_compounds': food}, {'receipt'})['items'] if x['id']=='food_compounds')['status']=='pending_verification'
    with pytest.raises(ValueError, match='Unknown'):
        evaluate({'altered_threshold': {}}, set())


def test_status_rechecks_evidence_and_computed_result(tmp_path):
    observations = {'chronic_classifier': observed('chronic_classifier', .93)}
    values = {'receipt': {'prediction_sha256': 'a'*64}, 'observations': observations,
              'goals': evaluate(observations, {'receipt'})}
    manifest = {'schema_version':1,'snapshot_version':'phase1','generated_utc':'2026-10-03T00:00:00Z','artifacts':{}}
    for key, data in values.items():
        raw = json.dumps(data).encode()
        (tmp_path / (key+'.json')).write_bytes(raw)
        manifest['artifacts'][key] = {'path': key+'.json', 'sha256':hashlib.sha256(raw).hexdigest()}
    path = tmp_path/'manifest.json'
    path.write_text(json.dumps(manifest))
    registry = ResearchRegistry(tmp_path, path)
    assert registry.status()['goals']['formal_contract_target_received'] is False
    (tmp_path/'receipt.json').write_text('{}')
    with pytest.raises(ValueError, match='checksum'):
        registry.status()


def all_guardrails_met():
    return {
        key: observed(key, threshold + .01 if operator == ">" else threshold)
        for key, (_, operator, threshold, _, _) in TARGETS.items()
    }


def legacy_schema2(result):
    result = copy.deepcopy(result)
    for key in ("criteria_scope", "criterion_scopes", "formal_acceptance_status", "formal_achievement_percentage"):
        result.pop(key)
    result.update({
        "schema_version": 2,
        "formal_contract_target_received": True,
        "target_basis": "user_confirmed_phase1_final_table",
        "interpretation": "Research acceptance evidence; historical achievements require their own source documents",
    })
    return result


def install_fixture(tmp_path, observations, goals):
    values = {"receipt": {"fixture_only": True}, "observations": observations, "goals": goals}
    manifest = {"schema_version": 1, "snapshot_version": "test-only",
                "generated_utc": "2026-10-07T00:00:00Z", "artifacts": {}}
    for key, value in values.items():
        raw = json.dumps(value).encode()
        (tmp_path / (key + ".json")).write_bytes(raw)
        manifest["artifacts"][key] = {
            "path": key + ".json", "sha256": hashlib.sha256(raw).hexdigest()}
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    return ResearchRegistry(tmp_path, path)


def test_all_internal_passes_do_not_verify_final_contract():
    result = evaluate(all_guardrails_met(), {"receipt"})
    assert result["required_criteria"] == 29
    assert result["achieved_criteria"] == 29
    assert result["all_metrics_satisfied"] is True
    assert result["criteria_scope"] == "implementation_guardrails"
    assert result["criterion_scopes"]["food_compounds"] == "additional_strict_identity_subset_not_overall_component_coverage"
    assert result["formal_contract_target_received"] is False
    assert result["formal_achievement_percentage"] is None
    assert result["formal_acceptance_status"] == "pending_contract_target_verification"


def test_frozen_legacy_snapshot_is_checked_then_provenance_corrected(tmp_path):
    observations = {"chronic_classifier": observed("chronic_classifier", .93)}
    current = evaluate(observations, {"receipt"})
    registry = install_fixture(tmp_path, observations, legacy_schema2(current))
    frozen = (tmp_path / "goals.json").read_bytes()
    assert registry.status()["goals"] == current
    assert (tmp_path / "goals.json").read_bytes() == frozen
    (tmp_path / "receipt.json").write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        registry.status()


@pytest.mark.parametrize("field", ["observed", "target", "status", "unit"])
def test_legacy_compatibility_does_not_relax_scientific_comparison(tmp_path, field):
    observations = {"chronic_classifier": observed("chronic_classifier", .93)}
    legacy = legacy_schema2(evaluate(observations, {"receipt"}))
    item = next(x for x in legacy["items"] if x["id"] == "chronic_classifier")
    item[field] = {"observed": .99, "target": .7, "status": "not_achieved",
                   "unit": "accuracy"}[field]
    registry = install_fixture(tmp_path, observations, legacy)
    with pytest.raises(ValueError, match="differs"):
        registry.status()


def test_corrected_snapshot_provenance_cannot_be_tampered(tmp_path):
    observations = {"chronic_classifier": observed("chronic_classifier", .93)}
    goals = evaluate(observations, {"receipt"})
    goals["formal_contract_target_received"] = True
    registry = install_fixture(tmp_path, observations, goals)
    with pytest.raises(ValueError, match="differs"):
        registry.status()


def test_legacy_protocol_and_self_referential_proof_are_rejected(tmp_path):
    observations = {"chronic_classifier": observed("chronic_classifier", .93)}
    legacy = legacy_schema2(evaluate(observations, {"receipt"}))
    legacy["protocol_id"] = "unrecognized"
    registry = install_fixture(tmp_path, observations, legacy)
    with pytest.raises(ValueError, match="protocol"):
        registry.status()
    legacy["protocol_id"] = "phase1-final-20261003-v1"
    observations["chronic_classifier"]["artifacts"] = ["goals"]
    registry = install_fixture(tmp_path, observations, legacy)
    with pytest.raises(ValueError, match="self-referential"):
        registry.status()
