import pandas as pd
from nutriomics_dti.external_audit import audit, aliases


def row(**changes):
    return {"compound_id": "full-key-A", "target_id": "sequence-A", "doi": "10.1/a", "pmid": "123", "endpoint": "Kd", "relation": "=", "value_nm": 10, **changes}


def test_paper_alias_blocks_another_database_record():
    reference = pd.DataFrame([row()])
    external = pd.DataFrame([row(compound_id="new-key", target_id="new-sequence", doi="https://doi.org/10.1/A", pmid="")])
    retained, result = audit(reference, external, complete_reference=True)
    assert len(retained) == 0
    assert result['rejection_reason_counts']['paper_patent_assay_source_overlap'] == 1


def test_partial_ledger_and_missing_sources_never_certify_independence():
    reference = pd.DataFrame([row(doi="", pmid="")])
    external = pd.DataFrame([row(compound_id="new", doi="10.2/b", pmid="456")])
    _, result = audit(reference, external, complete_reference=True)
    assert not result['independent_evaluation_eligible']
    assert result['reference_without_source_alias'] == 1


def test_same_pair_and_censored_endpoint_rejected_even_with_new_paper():
    reference = pd.DataFrame([row()])
    external = pd.DataFrame([row(doi="10.2/b", pmid="456", relation="<")])
    _, result = audit(reference, external, complete_reference=True)
    assert set(result['rejection_reason_counts']) == {'unsupported_or_invalid_exact_positive_kd', 'compound_target_pair_seen_in_reference'}


def test_source_and_pair_nonoverlap_requires_complete_ledger():
    reference = pd.DataFrame([row()])
    external = pd.DataFrame([row(compound_id="new", target_id="new", doi="10.2/b", pmid="456")])
    _, result = audit(reference, external)
    assert result['source_and_pair_nonoverlapping_candidates'] == 1 and not result['independent_evaluation_eligible']
    _, result = audit(reference, external, complete_reference=True)
    assert result['independent_evaluation_eligible']
