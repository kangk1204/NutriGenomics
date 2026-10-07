import json
import pytest
from nutriomics_evidence.quote_grounding import exact_occurrences, parse_quote_grounded


def example(text='Preface. Aspirin inhibits PTGS2. Tail.'):
    return {'pmid': '1', 'text': text, 'entities': [
        {'id': 'T1', 'type': 'CHEMICAL', 'start': text.index('Aspirin'), 'end': text.index('Aspirin') + 7, 'text': 'Aspirin'},
        {'id': 'T2', 'type': 'GENE', 'start': text.index('PTGS2'), 'end': text.index('PTGS2') + 5, 'text': 'PTGS2'}]}


def response(**changes):
    row = {'relation': 'INHIBITOR', 'arg1': 'T1', 'arg2': 'T2',
        'evidence_text': 'Aspirin inhibits PTGS2.', 'evidence_start': 0, 'evidence_end': 7}
    return {'relations': [{**row, **changes}]}


def reasons(value):
    return set(value['rejected'][0]['reasons'])


def test_original_json_deterministically_grounds_wrong_numeric_offsets():
    result = parse_quote_grounded(json.dumps(response()), example())
    row = result['accepted'][0]
    assert (row['evidence_start'], row['evidence_end']) == (9, 32)
    assert row['grounding_receipt']['original_offsets'] == {'start': 0, 'end': 7}
    assert row['grounding_receipt']['gold_relations_used'] is False
    assert row['binding_affinity'] is None and row['expert_approval'] is None
    assert row['semantic_entailment_validated'] is False and row['API_grounding_evaluated'] is False


def test_nullable_offsets_support_is_diagnostic_not_schema_or_api_validation():
    result = parse_quote_grounded(response(evidence_start=None, evidence_end=None), example())
    assert len(result['accepted']) == 1
    assert result['generation_schema_evaluated'] is False and result['API_grounding_evaluated'] is False


def test_valid_original_offsets_are_retained_without_change():
    row = parse_quote_grounded(response(evidence_start=9, evidence_end=32), example())['accepted'][0]
    assert row['grounding_receipt']['offsets_changed'] is False


def test_repeated_quote_is_ambiguous_even_when_supplied_offset_matches_first():
    source = example('Aspirin inhibits PTGS2. Aspirin inhibits PTGS2.')
    result = parse_quote_grounded(response(evidence_start=0, evidence_end=23), source)
    assert 'ambiguous_exact_quote' in reasons(result) and not result['accepted']


@pytest.mark.parametrize('quote', ['Aspirin  inhibits PTGS2.', 'aspirin inhibits PTGS2.', 'Aspirin inhibits PTGS2!', '  Aspirin inhibits PTGS2. '])
def test_no_whitespace_case_near_match_or_quote_shortening(quote):
    result = parse_quote_grounded(response(evidence_text=quote), example())
    assert 'quote_not_exact_source_substring' in reasons(result)


def test_quote_must_contain_native_mentions_not_merely_matching_names_elsewhere():
    source = example('Aspirin inhibits PTGS2. PTGS2 is elsewhere.')
    result = parse_quote_grounded(response(evidence_text='PTGS2 is elsewhere.'), source)
    assert 'evidence_does_not_contain_both_mentions' in reasons(result)


def test_negation_is_still_rejected_by_frozen_strict_gate():
    result = parse_quote_grounded(response(evidence_text='Aspirin does not inhibit PTGS2.'), example('Aspirin does not inhibit PTGS2.'))
    assert 'negation_requires_manual_review' in reasons(result)


@pytest.mark.parametrize('changes', [{'arg1': 'T2', 'arg2': 'T1'}, {'arg2': 'T9'}, {'arg2': ['T2']}, {'relation': 'THERAPEUTIC'}])
def test_wrong_entities_or_relation_are_not_repaired(changes):
    assert not parse_quote_grounded(response(**changes), example())['accepted']


@pytest.mark.parametrize('quote', [None, '', 4])
def test_missing_quote_abstains(quote):
    assert 'missing_or_invalid_exact_quote' in reasons(parse_quote_grounded(response(evidence_text=quote), example()))


@pytest.mark.parametrize('text', ['```json\n{"relations": []}\n```', '{"relations": []}```', '{"relations": [', '', '{"relations":[],"relations":[]}'])
def test_invalid_or_ambiguous_json_is_not_guessed(text):
    with pytest.raises(ValueError):
        parse_quote_grounded(text, example())


@pytest.mark.parametrize('offset', [True, '9', [9]])
def test_unsupported_offset_types_are_rejected(offset):
    assert 'unsupported_offset_type' in reasons(parse_quote_grounded(response(evidence_start=offset), example()))


def test_native_source_spans_are_verified():
    source = example();source['entities'][0]['start'] = 0
    with pytest.raises(ValueError, match='Native entity'):
        parse_quote_grounded(response(), source)


def test_overlapping_exact_occurrences_and_duplicate_predictions():
    assert exact_occurrences('aaaa', 'aaa') == [(0, 3), (1, 4)]
    payload = response();payload['relations'].append(dict(payload['relations'][0]))
    result = parse_quote_grounded(payload, example())
    assert len(result['accepted']) == 1 and 'duplicate_relation' in reasons(result)


def test_original_schema_missing_offsets_is_not_silently_changed():
    payload = response();del payload['relations'][0]['evidence_end']
    assert 'unsupported_or_missing_fields' in reasons(parse_quote_grounded(payload, example()))
