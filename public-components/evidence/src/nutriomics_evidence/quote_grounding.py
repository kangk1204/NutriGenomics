"""Exact-quotation grounding for a separately labelled post-hoc diagnostic.

This module neither reads reference relations nor repairs JSON, quotations,
entity IDs, labels or whitespace. A quotation must occur exactly once in the
unaltered source and contain both supplied native entity spans. Numeric offsets
are reconstructed from that occurrence, then the frozen strict parser is used
as the final gate. This is textual grounding, not semantic relation validation.
"""
from __future__ import annotations

import json
from .extraction import RELATIONS, parse_qwen

FIELDS = {'relation', 'arg1', 'arg2', 'evidence_text', 'evidence_start', 'evidence_end'}
POLICY_VERSION = 'exact-unique-quotation-native-entities-posthoc-v1'


def exact_occurrences(text: str, quote: str) -> list[tuple[int, int]]:
    """Find the complete literal quote, including overlapping occurrences."""
    if not isinstance(quote, str) or not quote:
        return []
    result, position = [], 0
    while True:
        start = text.find(quote, position)
        if start < 0:
            return result
        result.append((start, start + len(quote)))
        position = start + 1


def _unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError('Duplicate JSON object key; ambiguous response rejected')
        value[key] = item
    return value


def _entities(document):
    if not isinstance(document.get('text'), str):
        raise ValueError('Document requires original source text')
    lookup = {}
    for entity in document['entities']:
        start, end = entity['start'], entity['end']
        if (entity['id'] in lookup or not isinstance(entity['id'], str)
                or isinstance(start, bool) or isinstance(end, bool)
                or not isinstance(start, int) or not isinstance(end, int)
                or not 0 <= start < end <= len(document['text'])
                or document['text'][start:end] != entity['text']):
            raise ValueError('Native entity identity/span/text mismatch')
        lookup[entity['id']] = entity
    return lookup


def parse_quote_grounded(response: str | dict, document: dict) -> dict:
    """Ground one original six-field JSON response without reference labels.

    Supplied offsets may be integers or None; they are preserved in the receipt,
    not used to choose among multiple occurrences. Booleans/strings are rejected.
    The nullable-offset accommodation is diagnostic and is not a claim that a
    deployed generation schema or API has been evaluated.
    """
    try:
        payload = json.loads(response, object_pairs_hook=_unique_object) if isinstance(response, str) else response
    except json.JSONDecodeError as error:
        raise ValueError('Qwen output is not a single valid JSON value') from error
    if not isinstance(payload, dict) or set(payload) != {'relations'} or not isinstance(payload['relations'], list):
        raise ValueError('Qwen JSON must contain only a relations list')
    lookup = _entities(document)
    accepted, rejected, seen = [], [], set()
    for index, row in enumerate(payload['relations']):
        reasons = []
        if not isinstance(row, dict):
            rejected.append({'index': index, 'reasons': ['relation_not_object']})
            continue
        if set(row) != FIELDS:
            reasons.append('unsupported_or_missing_fields')
        relation = row.get('relation')
        if not isinstance(relation, str) or relation not in RELATIONS:
            reasons.append('unsupported_relation')
        arg1, arg2 = row.get('arg1'), row.get('arg2')
        a = lookup.get(arg1) if isinstance(arg1, str) else None
        b = lookup.get(arg2) if isinstance(arg2, str) else None
        if a is None or b is None or a['type'] != 'CHEMICAL' or not b['type'].startswith('GENE'):
            reasons.append('unknown_or_wrong_entity_types')
        original = {'start': row.get('evidence_start'), 'end': row.get('evidence_end')}
        if any(value is not None and (not isinstance(value, int) or isinstance(value, bool)) for value in original.values()):
            reasons.append('unsupported_offset_type')
        quote = row.get('evidence_text')
        occurrences = exact_occurrences(document['text'], quote)
        if not isinstance(quote, str) or not quote:
            reasons.append('missing_or_invalid_exact_quote')
        elif not occurrences:
            reasons.append('quote_not_exact_source_substring')
        elif len(occurrences) != 1:
            reasons.append('ambiguous_exact_quote')
        elif a and b and not all(occurrences[0][0] <= e['start'] and e['end'] <= occurrences[0][1] for e in (a, b)):
            reasons.append('evidence_does_not_contain_both_mentions')
        key = (relation, arg1, arg2)
        if all(isinstance(x, str) for x in key) and key in seen:
            reasons.append('duplicate_relation')
        if reasons:
            rejected.append({'index': index, 'reasons': reasons, 'exact_quote_occurrences': len(occurrences)})
            continue
        start, end = occurrences[0]
        grounded = {**row, 'evidence_start': start, 'evidence_end': end}
        # The unchanged parser rechecks the actual source slice, both native
        # mentions, allowed types, negation, and original six-field schema.
        checked = parse_qwen({'relations': [grounded]}, document)
        if checked['rejected']:
            rejected.append({'index': index, 'reasons': checked['rejected'][0]['reasons'], 'exact_quote_occurrences': 1})
            continue
        seen.add(key)
        accepted.append({**checked['accepted'][0],
            'method': 'qwen_exact_quote_posthoc_diagnostic',
            'grounding_receipt': {'policy_version': POLICY_VERSION, 'exact_quote_occurrences': 1,
                'original_offsets': original, 'resolved_offsets': {'start': start, 'end': end},
                'offsets_changed': original != {'start': start, 'end': end},
                'native_entity_ids': [arg1, arg2], 'gold_relations_used': False},
            'semantic_entailment_validated': False, 'binding_affinity': None,
            'relation_is_causal': False, 'expert_approval': None,
            'evaluation_status': 'post-hoc diagnostic on previously scored development responses',
            'generation_schema_evaluated': False, 'API_grounding_evaluated': False})
    return {'pmid': document['pmid'], 'accepted': accepted, 'rejected': rejected,
        'policy_version': POLICY_VERSION, 'gold_relations_used_for_grounding': False,
        'semantic_entailment_validated': False, 'generation_schema_evaluated': False,
        'API_grounding_evaluated': False}
