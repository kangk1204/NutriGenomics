import gzip
import hashlib
import json
import sqlite3
import unicodedata

import pytest
from normalize_food_db import STRICT, add_term_seed, connect, import_projection, initialize, report

K1 = 'AAAAAAAAAAAAAA-BBBBBBBBBB-C'
K2 = 'AAAAAAAAAAAAAA-DDDDDDDDDD-C'

def row(identifier='r1', key=K1, name='Shared label', food_name='Plant raw leaf', part='leaf', status=STRICT):
    return {'id': identifier, 'full_inchikey': key, 'name': name, 'source_inchi': 'synthetic-test-inchi',
            'canonical_isomeric_smiles': 'test.salt' if key == K2 else 'test', 'standard_structure_verified': True,
            'single_component': key != K2, 'standard_unspecified_stereo': False, 'public_compound_id': key,
            'source_id': key, 'orig_source_name': name, 'food_id': 'food1', 'public_food_id': 'FDB-food1',
            'orig_food_common_name': food_name, 'orig_food_part': part, 'strict_status': status,
            'positive_value': 2.5, 'orig_unit': 'mg/100g', 'orig_method': None, 'citation': 'TEST-ONLY',
            'citation_type': 'FIXTURE', 'source_pmid': None, 'orig_source_id': identifier}

@pytest.fixture
def con():
    database = connect(':memory:'); initialize(database)
    yield database
    database.close()

def projection(tmp_path, audit, standards=None, footer=True):
    inputs = [{'remote_path': '/SYNTHETIC-TEST/strict_verified_food_compounds.parquet', 'sha256': '1'*64, 'bytes': 1},
              {'remote_path': '/SYNTHETIC-TEST/strict_occurrence_audit.parquet', 'sha256': '2'*64, 'bytes': 2}]
    standards = standards if standards is not None else [audit[0]]
    items = [{'kind': 'header', 'inputs': inputs}]
    items += [{'kind': 'row', 'table': table, 'record': r} for table, values in [('standards', standards), ('audit', audit)] for r in values]
    if footer: items += [{'kind': 'footer', 'counts': {'standards': len(standards), 'audit': len(audit)}, 'total_rows':len(standards)+len(audit)}]
    path = tmp_path / 'synthetic.ndjson.gz'
    with gzip.open(path, 'wt') as stream:
        for item in items: stream.write(json.dumps(item)+'\n')
    return path, hashlib.sha256(path.read_bytes()).hexdigest()

def test_full_key_not_name_or_connectivity_identity(con, tmp_path):
    path, sha = projection(tmp_path, [row(), row('r2', K2, status='form_ambiguous')])
    import_projection(con, path, sha)
    assert con.execute('SELECT count(*) FROM compound').fetchone()[0] == 2
    assert con.execute('SELECT count(*) FROM compound_alias WHERE original_name=\'Shared label\'').fetchone()[0] == 2
    assert con.execute('SELECT original_isomeric_smiles FROM structure_representation WHERE full_inchikey=?', (K2,)).fetchone()[0] == 'test.salt'

def test_idempotent_import_and_schema(con, tmp_path):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    changes = con.total_changes
    version = con.execute('PRAGMA schema_version').fetchone()[0]
    initialize(con)
    second = import_projection(con, path, sha)
    assert second['already_imported'] and second['new_raw_rows'] == 0
    assert con.total_changes == changes
    assert con.execute('PRAGMA schema_version').fetchone()[0] == version

def test_missing_footer_rolls_back_all_rows(con, tmp_path):
    path, sha = projection(tmp_path, [row()], footer=False)
    with pytest.raises(ValueError, match='footer'): import_projection(con, path, sha)
    assert con.execute('SELECT count(*) FROM raw_record').fetchone()[0] == 0
    assert con.execute('SELECT count(*) FROM source_snapshot').fetchone()[0] == 0

def test_duplicate_source_row_rolls_back(con, tmp_path):
    path, sha = projection(tmp_path, [row(), row()])
    with pytest.raises(sqlite3.IntegrityError): import_projection(con, path, sha)
    assert con.execute('SELECT count(*) FROM observation').fetchone()[0] == 0

def test_changed_projection_is_rejected_before_import(con, tmp_path):
    path, sha = projection(tmp_path, [row()])
    with pytest.raises(ValueError, match='SHA'): import_projection(con, path, '0'*64)
    assert con.execute('SELECT count(*) FROM source_snapshot').fetchone()[0] == 0

def test_held_analyte_is_not_promoted_as_synonym(con, tmp_path):
    held = row('r2', status='original_analyte_ambiguous'); held['orig_source_name'] = 'Unresolved class total'; held['orig_unit'] = 'unknown native unit'
    path, sha = projection(tmp_path, [row(), held]); import_projection(con, path, sha)
    assert con.execute('SELECT count(*) FROM compound_alias WHERE original_name=\'Unresolved class total\'').fetchone()[0] == 0
    observed = con.execute('SELECT quality_state,original_analyte_name,original_unit,normalized_unit_code FROM observation WHERE source_row_id=\'r2\'').fetchone()
    assert observed == ('retained_held', 'Unresolved class total', 'unknown native unit', None)
    assert 'Unresolved class total' in con.execute('SELECT raw_json FROM raw_record WHERE source_table=\'audit\' AND source_row_id=\'r2\'').fetchone()[0]

def test_food_parts_preparation_names_and_nulls_preserved(con, tmp_path):
    path, sha = projection(tmp_path, [row(), row('r2', food_name='Plant boiled root', part=None)])
    import_projection(con, path, sha)
    contexts = con.execute('SELECT original_common_name,original_part,preparation_state,taxon_state FROM food_context ORDER BY original_common_name').fetchall()
    assert contexts == [('Plant boiled root', None, 'not_provided', 'not_provided'), ('Plant raw leaf', 'leaf', 'not_provided', 'not_provided')]

def test_orphan_observation_is_rejected(con, tmp_path):
    path, sha = projection(tmp_path, [row()]); imported = import_projection(con, path, sha)
    with pytest.raises(sqlite3.IntegrityError):
        con.execute('UPDATE observation SET full_inchikey=?', (K2,))
    con.rollback()
    assert con.execute('PRAGMA foreign_key_check').fetchall() == []

def evidence(identifier='official', method='direct_primary_source', **extras):
    return {'evidence_id':identifier, 'authority':'SYNTHETIC-TEST-AUTHORITY', 'source_uri':'https://example.invalid/test-only',
            'source_version':'test-fixture', 'retrieved_utc':'2026-01-01T00:00:00Z', 'verification_method':method, **extras}

def seed(method='direct_primary_source', form=True, key=K1):
    return {'evidence':[evidence(method=method), evidence('identity', full_inchikey=key, form_review_passed=form)],
            'terms':[{'term_id':'term1','evidence_id':'official','original_ko':'시험명','original_en':'Test name','entity_scope':'compound_concept','scope_review_state':'reviewed'}]}

def binding(con, grade='exact_form', state='verified', key=K1):
    con.execute('INSERT INTO compound_term_binding VALUES(?,?,?,?,?,?,?)', (K1,'term1',grade,state,key,'identity','synthetic test'))

def test_language_presence_is_not_authority_verified(con, tmp_path):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    add_term_seed(con, seed(method='delegated_primary_source')); binding(con, state='candidate')
    coverage = report(con)['compound_language_coverage']
    assert coverage['source_en_present'] == 1 and coverage['source_ko_present'] == 0
    assert coverage['authority_ko_candidate_or_verified'] == 1 and coverage['authority_verified_exact_bilingual'] == 0
    strict = report(con)['strict_projection_language_coverage']
    assert strict['source_en_present'] == 1 and strict['source_ko_present'] == 0

@pytest.mark.parametrize('method,form,key,grade', [('delegated_primary_source',True,K1,'exact_form'),
                                                ('direct_primary_source',False,K1,'exact_form'),
                                                ('direct_primary_source',True,K2,'exact_form'),
                                                ('direct_primary_source',True,K1,'broader_concept')])
def test_bilingual_exact_proof_guards(con, tmp_path, method, form, key, grade):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    add_term_seed(con, seed(method,form,key))
    with pytest.raises(sqlite3.IntegrityError): binding(con, grade=grade)
    con.rollback()

def test_verified_pair_requires_matching_primary_structure(con, tmp_path):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    add_term_seed(con, seed()); binding(con)
    assert report(con)['compound_language_coverage']['authority_verified_exact_bilingual'] == 1

def test_authority_changes_cannot_silently_overwrite(con):
    data = seed(); add_term_seed(con, data); add_term_seed(con, data)
    data['evidence'][0]['source_version'] = 'different'
    with pytest.raises(ValueError, match='new evidence version'): add_term_seed(con, data)

@pytest.mark.parametrize('mutation', ["positive_value=9.0", "original_unit='g/100g'", "original_method='invented method'",
                                     "original_analyte_name='invented synonym'", "source_citation='invented citation'",
                                     "source_pmid='123'", "original_analyte_id='different-id'"])
def test_source_measurements_and_provenance_cannot_be_rewritten(con, tmp_path, mutation):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    with pytest.raises(sqlite3.IntegrityError): con.execute('UPDATE observation SET ' + mutation)
    con.rollback()
    assert con.execute('SELECT positive_value,original_unit,source_citation FROM observation').fetchone() == (2.5, 'mg/100g', 'TEST-ONLY')

def test_held_source_cannot_be_relabelled_strict(con, tmp_path):
    path, sha = projection(tmp_path, [row(status='ambiguous_analyte')]); import_projection(con, path, sha)
    with pytest.raises(sqlite3.IntegrityError):
        con.execute('UPDATE observation SET strict_status=?,quality_state=\'accepted_strict\'', (STRICT,))
    con.rollback()
    assert con.execute('SELECT quality_state FROM observation').fetchone()[0] == 'retained_held'

@pytest.mark.parametrize('sql', ['UPDATE raw_record SET raw_json=\'{}\'', 'DELETE FROM raw_record',
                               'UPDATE food_context SET original_part=\'seed\''])
def test_raw_source_and_food_context_immutable(con, tmp_path, sql):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    with pytest.raises(sqlite3.IntegrityError): con.execute(sql)
    con.rollback()

@pytest.mark.parametrize('sql', ["UPDATE authority_evidence SET verification_method='direct_primary_source'",
                               'DELETE FROM authority_evidence', "UPDATE official_term SET scope_review_state='reviewed'",
                               'DELETE FROM official_term'])
def test_authority_direct_sql_requires_new_version(con, sql):
    add_term_seed(con, seed(method='delegated_primary_source'))
    with pytest.raises(sqlite3.IntegrityError): con.execute(sql)
    con.rollback()

def test_candidate_cannot_be_promoted_by_update(con, tmp_path):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    add_term_seed(con, seed(method='delegated_primary_source')); binding(con, state='candidate'); con.commit()
    with pytest.raises(sqlite3.IntegrityError): con.execute("UPDATE compound_term_binding SET review_state='verified'")
    con.rollback()
    assert report(con)['compound_language_coverage']['authority_verified_exact_bilingual'] == 0

def test_blocked_term_is_not_counted_as_candidate_coverage(con, tmp_path):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    add_term_seed(con, seed()); binding(con, state='blocked')
    results = report(con)
    assert results['compound_language_coverage']['authority_ko_candidate_or_verified'] == 0
    assert results['strict_projection_language_coverage']['authority_ko_candidate_or_verified'] == 0

def test_decomposed_korean_detected_without_changing_original(con, tmp_path):
    original = unicodedata.normalize('NFD', '시험명')
    path, sha = projection(tmp_path, [row(name=original)]); import_projection(con, path, sha)
    assert con.execute('SELECT original_name,search_name,language FROM compound_alias').fetchone() == (original, '시험명', 'ko')

@pytest.mark.parametrize('scope,prep,context_match', [('organism', True, True), ('food_context', False, True), ('food_context', True, False)])
def test_nonblank_official_food_pair_needs_context_proof(con, tmp_path, scope, prep, context_match):
    path, sha = projection(tmp_path, [row()]); import_projection(con, path, sha)
    context = con.execute('SELECT context_id FROM food_context').fetchone()[0]
    data = {'evidence':[evidence(), evidence('food-identity',food_context_id=context if context_match else 'wrong-context',taxon_match=True,part_match=True,preparation_match=prep)],
            'terms':[{'term_id':'food-term','evidence_id':'official','original_ko':'정제 시험유','original_en':'Test plant','entity_scope':scope,'scope_review_state':'reviewed'}]}
    add_term_seed(con, data)
    con.execute('INSERT INTO food_term_binding VALUES(?,?,?,?,?,?,?,?,?)', (context,'food-term','exact_context','candidate',1,1,1,'food-identity','synthetic scope test'))
    con.commit()
    with pytest.raises(sqlite3.IntegrityError): con.execute("UPDATE food_term_binding SET review_state='verified'")
    con.rollback()
    assert con.execute("SELECT count(*) FROM food_term_binding WHERE review_state='verified'").fetchone()[0] == 0
